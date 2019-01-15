#ifndef GraphAlignerBitvectorBanded_h
#define GraphAlignerBitvectorBanded_h

#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include "AlignmentGraph.h"
#include "NodeSlice.h"
#include "CommonUtils.h"
#include "GraphAlignerWrapper.h"
#include "AlignmentCorrectnessEstimation.h"
#include "ThreadReadAssertion.h"
#include "WordSlice.h"
#include "GraphAlignerBitvectorCommon.h"
#include "GraphAlignerCommon.h"
#include "ArrayPriorityQueue.h"

#ifndef NDEBUG
thread_local int debugLastRowMinScore;
#endif

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerBitvectorBanded
{
private:
	using BV = GraphAlignerBitvectorCommon<LengthType, ScoreType, Word>;
	using Common = GraphAlignerCommon<LengthType, ScoreType, Word>;
	using AlignerGraphsizedState = typename Common::AlignerGraphsizedState;
	using Params = typename Common::Params;
	using MatrixPosition = typename Common::MatrixPosition;
	using Trace = typename Common::Trace;
	using OnewayTrace = typename Common::OnewayTrace;
	using WordSlice = typename BV::WordSlice;
	using EqVector = typename BV::EqVector;
	using EdgeWithPriority = typename Common::EdgeWithPriority;
	const Params& params;
	class DPSlice
	{
	public:
		DPSlice() :
		minScore(std::numeric_limits<ScoreType>::max()),
		minScoreNode(std::numeric_limits<LengthType>::max()),
		minScoreNodeOffset(std::numeric_limits<LengthType>::max()),
		scoresVectorMap(),
		scores(),
		correctness(),
		j(std::numeric_limits<LengthType>::max()),
		cellsProcessed(0),
		bandwidth(0),
		scoresNotValid(false)
#ifdef SLICEVERBOSE
		,nodesProcessed(0)
		,numCells(0)
#endif
		{}
		DPSlice(std::vector<typename NodeSlice<LengthType, ScoreType, Word, true>::MapItem>* vectorMap) :
		minScore(std::numeric_limits<ScoreType>::max()),
		minScoreNode(std::numeric_limits<LengthType>::max()),
		minScoreNodeOffset(std::numeric_limits<LengthType>::max()),
		scoresVectorMap(vectorMap),
		scores(),
		correctness(),
		j(std::numeric_limits<LengthType>::max()),
		cellsProcessed(0),
		bandwidth(0),
		scoresNotValid(false)
#ifdef SLICEVERBOSE
		,nodesProcessed(0)
		,numCells(0)
#endif
		{}
		ScoreType minScore;
		LengthType minScoreNode;
		LengthType minScoreNodeOffset;
		NodeSlice<LengthType, ScoreType, Word, true> scoresVectorMap;
		NodeSlice<LengthType, ScoreType, Word, false> scores;
		AlignmentCorrectnessEstimationState correctness;
		LengthType j;
		size_t cellsProcessed;
		size_t bandwidth;
		bool scoresNotValid;
#ifdef SLICEVERBOSE
		size_t nodesProcessed;
		size_t numCells;
#endif
		DPSlice getMapSlice() const
		{
			DPSlice result;
			result.minScore = minScore;
			result.minScoreNode = minScoreNode;
			result.minScoreNodeOffset = minScoreNodeOffset;
			assert(scores.size() != 0);
			result.scores = scores;
			result.correctness = correctness;
			result.j = j;
			result.cellsProcessed = cellsProcessed;
			result.bandwidth = bandwidth;
			result.scoresNotValid = scoresNotValid;
#ifdef SLICEVERBOSE
			result.nodesProcessed = nodesProcessed;
			result.numCells = numCells;
#endif
			return result;
		}
	};
	class DPTable
	{
	public:
		DPTable() :
		slices()
		{}
		std::vector<DPSlice> slices;
	};
public:

	GraphAlignerBitvectorBanded(const Params& params) :
	params(params)
	{
	}

	OnewayTrace getReverseTraceFromSeed(const std::string& sequence, int bigraphNodeId, size_t nodeOffset, bool forceGlobal, AlignerGraphsizedState& reusableState) const
	{
		size_t numSlices = (sequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		auto initialBandwidth = getInitialSliceExactPosition(bigraphNodeId, nodeOffset);
		auto slice = getSqrtSlices(sequence, initialBandwidth, numSlices, forceGlobal, reusableState);
		if (!forceGlobal) removeWronglyAlignedEnd(slice);
		if (slice.slices.size() <= 1)
		{
			return OnewayTrace::TraceFailed();
		}
		assert(sequence.size() <= std::numeric_limits<ScoreType>::max() - WordConfiguration<Word>::WordSize * 2);
		assert(slice.slices.back().minScore >= 0);
		assert(slice.slices.back().minScore <= (ScoreType)sequence.size() + (ScoreType)WordConfiguration<Word>::WordSize * 2);

		OnewayTrace result;

		result = getReverseTraceFromTable(sequence, slice, reusableState);

		return result;
	}

	OnewayTrace getBacktraceFullStart(std::string originalSequence, bool forceGlobal, AlignerGraphsizedState& reusableState) const
	{
		assert(originalSequence.size() > 1);
		DPSlice startSlice;
		startSlice.j = -WordConfiguration<Word>::WordSize;
		startSlice.scores.addEmptyNodeMap(params.graph.NodeSize());
		startSlice.bandwidth = 1;
		startSlice.minScore = 0;
		startSlice.minScoreNode = 0;
		startSlice.minScoreNodeOffset = 0;
		char firstChar = originalSequence[0];
		for (size_t i = 0; i < params.graph.NodeSize(); i++)
		{
			startSlice.scores.addNodeToMap(i);
			startSlice.scores.setMinScore(i, 0);
			auto& node = startSlice.scores.node(i);
			bool match = Common::characterMatch(firstChar, params.graph.NodeSequences(i, 0));
			node.startSlice = {0, 0, match ? 0 : 1};
			node.minScore = match ? 0 : 1;
			for (size_t j = 1; j < params.graph.NodeLength(i); j++)
			{
				bool oldMatch = match;
				match = Common::characterMatch(firstChar, params.graph.NodeSequences(i, j));
				if (oldMatch && !match)
				{
					node.HP[j / params.graph.SPLIT_NODE_SIZE] |= ((Word)1) << (j % params.graph.SPLIT_NODE_SIZE);
				}
				else if (match && !oldMatch)
				{
					node.HN[j / params.graph.SPLIT_NODE_SIZE] |= ((Word)1) << (j % params.graph.SPLIT_NODE_SIZE);
				}
				if (match) node.minScore = 0;
			}
			node.endSlice = {0, 0, match ? 0 : 1};
			node.exists = true;
		}
		std::string alignableSequence = originalSequence.substr(1);
		assert(alignableSequence.size() > 0);
		size_t numSlices = (alignableSequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		auto slice = getSqrtSlices(alignableSequence, startSlice, numSlices, forceGlobal, reusableState);
		if (!forceGlobal) removeWronglyAlignedEnd(slice);
		if (slice.slices.size() <= 1)
		{
			return OnewayTrace::TraceFailed();
		}

		auto result = getReverseTraceFromTable(alignableSequence, slice, reusableState);
		for (size_t i = 0; i < result.trace.size(); i++)
		{
			result.trace[i].first.seqPos += 1;
		}
		std::reverse(result.trace.begin(), result.trace.end());
		assert(result.trace[0].first.seqPos == 0);
		return result;
	}

private:

	OnewayTrace getReverseTraceFromTable(const std::string& sequence, const DPTable& slice, AlignerGraphsizedState& reusableState) const
	{
		assert(slice.slices.size() > 0);
		assert(slice.slices.back().minScoreNode != std::numeric_limits<LengthType>::max());
		assert(slice.slices.back().minScoreNodeOffset != std::numeric_limits<LengthType>::max());
		OnewayTrace result;
		result.score = slice.slices.back().minScore;
		result.trace.emplace_back(MatrixPosition {slice.slices.back().minScoreNode, slice.slices.back().minScoreNodeOffset, std::min(slice.slices.back().j + WordConfiguration<Word>::WordSize - 1, sequence.size()-1)}, false);
		LengthType currentNode = std::numeric_limits<LengthType>::max();
		size_t currentSlice = slice.slices.size();
		std::vector<WordSlice> nodeSlices;
		while (result.trace.back().first.seqPos != (size_t)-1)
		{
			size_t newSlice = result.trace.back().first.seqPos / WordConfiguration<Word>::WordSize + 1;
			assert(newSlice < slice.slices.size());
			assert(result.trace.back().first.seqPos >= slice.slices[newSlice].j);
			assert(result.trace.back().first.seqPos < slice.slices[newSlice].j + WordConfiguration<Word>::WordSize);
			LengthType newNode = result.trace.back().first.node;
			if (newSlice != currentSlice || newNode != currentNode)
			{
				currentSlice = newSlice;
				currentNode = newNode;
				assert(slice.slices[currentSlice].scores.hasNode(currentNode));
				assert(currentSlice > 0);
				typename NodeSlice<LengthType, ScoreType, Word, false>::NodeSliceMapItem previous;
				if (slice.slices[currentSlice-1].scores.hasNode(currentNode))
				{
					previous = slice.slices[currentSlice-1].scores.node(currentNode);
				}
				else
				{
					for (size_t i = 0; i < previous.NUM_CHUNKS; i++)
					{
						previous.HP[i] = WordConfiguration<Word>::AllOnes;
						previous.HN[i] = WordConfiguration<Word>::AllZeros;
					}
				}
				nodeSlices = recalcNodeWordslice(currentNode, slice.slices[currentSlice].scores.node(currentNode), previous, slice.slices[currentSlice].j, sequence);
#ifdef SLICEVERBOSE
				std::cerr << "j " << slice.slices[currentSlice].j << " firstbt-calc " << slice.slices[currentSlice].scores.node(currentNode).firstSlicesCalcedWhenCalced << " lastbt-calc " << slice.slices[currentSlice].scores.node(currentNode).slicesCalcedWhenCalced << std::endl;
#endif
			}
			assert(result.trace.back().first.node == currentNode);
			assert(result.trace.back().first.nodeOffset < params.graph.NodeLength(currentNode));
			assert(nodeSlices.size() == params.graph.NodeLength(currentNode));
			assert(result.trace.back().first.seqPos >= slice.slices[currentSlice].j);
			assert(result.trace.back().first.seqPos < slice.slices[currentSlice].j + WordConfiguration<Word>::WordSize);
			assert((ScoreType)slice.slices[currentSlice].bandwidth >= 0);
			assert((ScoreType)slice.slices[currentSlice].bandwidth < std::numeric_limits<ScoreType>::max());
			assert(slice.slices[currentSlice].minScore < std::numeric_limits<ScoreType>::max() - (ScoreType)slice.slices[currentSlice].bandwidth);
			if (result.trace.back().first.seqPos % WordConfiguration<Word>::WordSize == 0 && result.trace.back().first.nodeOffset == 0)
			{
				result.trace.emplace_back(pickBacktraceCorner(slice.slices[currentSlice].scores, slice.slices[currentSlice-1].scores, currentNode, slice.slices[currentSlice].j, sequence, slice.slices[currentSlice].minScore + slice.slices[currentSlice].bandwidth, slice.slices[currentSlice].scoresNotValid));
				checkBacktraceCircularity(result);
				continue;
			}
			if (result.trace.back().first.seqPos % WordConfiguration<Word>::WordSize == 0)
			{
				assert(currentSlice > 0);
				assert(result.trace.back().first.nodeOffset > 0);
				if (!slice.slices[currentSlice-1].scores.hasNode(currentNode))
				{
					result.trace.emplace_back(MatrixPosition {currentNode, 0, result.trace.back().first.seqPos}, false);
					continue;
				}
				auto crossing = pickBacktraceVerticalCrossing(slice.slices[currentSlice].scores, slice.slices[currentSlice-1].scores, nodeSlices, slice.slices[currentSlice].j, currentNode, result.trace.back().first, sequence, slice.slices[currentSlice].minScore + slice.slices[currentSlice].bandwidth, slice.slices[currentSlice].scoresNotValid);
				if (crossing.first.first != result.trace.back().first) result.trace.push_back(crossing.first);
				assert(crossing.second.first != result.trace.back().first);
				result.trace.push_back(crossing.second);
				continue;
			}
			if (result.trace.back().first.nodeOffset == 0)
			{
				assert(result.trace.back().first.seqPos % WordConfiguration<Word>::WordSize != 0);
				auto crossing = pickBacktraceHorizontalCrossing(slice.slices[currentSlice].scores, slice.slices[currentSlice-1].scores, slice.slices[currentSlice].j, currentNode, result.trace.back().first, sequence, slice.slices[currentSlice].minScore + slice.slices[currentSlice].bandwidth, slice.slices[currentSlice].scoresNotValid);
				if (crossing.first.first != result.trace.back().first) result.trace.push_back(crossing.first);
				assert(crossing.second.first != result.trace.back().first);
				result.trace.push_back(crossing.second);
				checkBacktraceCircularity(result);
				continue;
			}
			assert(result.trace.back().first.nodeOffset != 0);
			assert(result.trace.back().first.seqPos % WordConfiguration<Word>::WordSize != 0);
			result.trace.emplace_back(pickBacktraceInside(slice.slices[currentSlice].j, nodeSlices, result.trace.back().first, sequence), false);
		}
		do
		{
			assert(result.trace.back().first.seqPos == (size_t)-1);
			assert(slice.slices[0].scores.hasNode(result.trace.back().first.node));
			auto node = slice.slices[0].scores.node(result.trace.back().first.node);
			std::vector<ScoreType> beforeSliceScores;
			beforeSliceScores.resize(params.graph.NodeLength(result.trace.back().first.node));
			beforeSliceScores[0] = node.startSlice.scoreEnd;
			for (size_t i = 1; i < beforeSliceScores.size(); i++)
			{
				size_t chunk = i / params.graph.SPLIT_NODE_SIZE;
				size_t offset = i % params.graph.SPLIT_NODE_SIZE;
				Word mask = ((Word)1) << offset;
				beforeSliceScores[i] = beforeSliceScores[i-1] + ((node.HP[chunk] & mask) >> offset) - ((node.HN[chunk] & mask) >> offset);
			}
			assert(beforeSliceScores.back() == node.endSlice.scoreEnd);
			while (beforeSliceScores[result.trace.back().first.nodeOffset] != 0 && result.trace.back().first.nodeOffset > 0 && beforeSliceScores[result.trace.back().first.nodeOffset-1] == beforeSliceScores[result.trace.back().first.nodeOffset] - 1)
			{
				result.trace.emplace_back(MatrixPosition {result.trace.back().first.node, result.trace.back().first.nodeOffset-1, result.trace.back().first.seqPos}, false);
			}
			if (result.trace.back().first.nodeOffset == 0 && beforeSliceScores[result.trace.back().first.nodeOffset] != 0)
			{
				bool found = false;
				for (auto neighbor : params.graph.inNeighbors[result.trace.back().first.node])
				{
					if (slice.slices[0].scores.hasNode(neighbor) && slice.slices[0].scores.node(neighbor).endSlice.getScoreBeforeStart() == beforeSliceScores[result.trace.back().first.nodeOffset] - 1)
					{
						result.trace.emplace_back(MatrixPosition {neighbor, params.graph.NodeLength(neighbor)-1, result.trace.back().first.seqPos}, true);
						found = true;
						break;
					}
				}
				if (found) continue;
			}
		} while (false);
		return result;
	}

	void checkBacktraceCircularity(const OnewayTrace& result) const
	{
		for (size_t i = result.trace.size()-2; i < result.trace.size(); i--)
		{
			assert(result.trace[i].first != result.trace.back().first);
			if (result.trace[i].first == result.trace.back().first) std::abort();
			if (result.trace[i].first.seqPos != result.trace.back().first.seqPos) return;
		}
	}

	MatrixPosition pickBacktraceInside(LengthType verticalOffset, const std::vector<WordSlice>& nodeSlices, MatrixPosition pos, const std::string& sequence) const
	{
		assert(verticalOffset <= pos.seqPos);
		assert(verticalOffset + WordConfiguration<Word>::WordSize > pos.seqPos);
		assert((verticalOffset % WordConfiguration<Word>::WordSize) == 0);
		size_t hori = pos.nodeOffset;
		size_t vert = pos.seqPos - verticalOffset;
		assert(vert >= 0);
		assert(vert < WordConfiguration<Word>::WordSize);
		assert(hori >= 0);
		assert(hori < nodeSlices.size());
		while (hori > 0 && vert > 0)
		{
			ScoreType scoreHere = nodeSlices[hori].getValue(vert);
			ScoreType verticalScore = nodeSlices[hori].getValue(vert-1);
			ScoreType horizontalScore = nodeSlices[hori-1].getValue(vert);
			ScoreType diagonalScore = nodeSlices[hori-1].getValue(vert-1);
			bool eq = Common::characterMatch(sequence[vert + verticalOffset], params.graph.NodeSequences(pos.node, hori));
			assert(verticalScore >= scoreHere-1);
			assert(horizontalScore >= scoreHere-1);
			assert(diagonalScore >= scoreHere - (eq?0:1));
			if (verticalScore == scoreHere - 1)
			{
				vert--;
				continue;
			}
			if (diagonalScore == scoreHere - (eq?0:1))
			{
				hori--;
				vert--;
				continue;
			}
			assert(horizontalScore == scoreHere - 1);
			hori--;
			continue;
		}
		return MatrixPosition { pos.node, hori, vert + verticalOffset };
	}

	std::pair<std::pair<MatrixPosition, bool>, std::pair<MatrixPosition, bool>> pickBacktraceHorizontalCrossing(const NodeSlice<LengthType, ScoreType, Word, false>& current, const NodeSlice<LengthType, ScoreType, Word, false>& previous, size_t j, LengthType node, MatrixPosition pos, const std::string& sequence, ScoreType quitScore, bool scoresNotValid) const
	{
		assert(current.hasNode(node));
		auto startSlice = current.node(node).startSlice;
		while (pos.seqPos % WordConfiguration<Word>::WordSize != 0 && (startSlice.VP & ((Word)1 << (pos.seqPos % WordConfiguration<Word>::WordSize))))
		{
			pos.seqPos--;
		}
		size_t offset = pos.seqPos % WordConfiguration<Word>::WordSize;
		if (offset == 0)
		{
			return std::make_pair(std::make_pair(pos, false), pickBacktraceCorner(current, previous, node, j, sequence, quitScore, scoresNotValid));
		}
		bool eq = Common::characterMatch(sequence[pos.seqPos], params.graph.NodeSequences(pos.node, pos.nodeOffset));
		ScoreType scoreHere = startSlice.getValue(offset);
		if (scoresNotValid || scoreHere > quitScore)
		{
			//this location is out of the band so the usual horizontal and vertical score limits don't apply
			//just pick the smallest scoring in-neighbor
			assert(offset > 0);
			ScoreType smallestFound = startSlice.getValue(offset-1);
			MatrixPosition smallestPos { node, 0, pos.seqPos-1 };
			bool nodeChange = false;
			for (auto neighbor : params.graph.inNeighbors[node])
			{
				if (current.hasNode(neighbor))
				{
					auto neighborSlice = current.node(neighbor).endSlice;
					if (neighborSlice.getValue(offset-1) < smallestFound)
					{
						smallestFound = neighborSlice.getValue(offset-1);
						smallestPos = MatrixPosition { neighbor, params.graph.NodeLength(neighbor)-1, pos.seqPos-1 };
						nodeChange = true;
					}
					if (neighborSlice.getValue(offset) < smallestFound && neighbor != node)
					{
						smallestFound = neighborSlice.getValue(offset);
						smallestPos = MatrixPosition { neighbor, params.graph.NodeLength(neighbor)-1, pos.seqPos };
						nodeChange = true;
					}
				}
			}
			assert(smallestPos != pos);
			return std::make_pair(std::make_pair(pos, false), std::make_pair(smallestPos, nodeChange));
		}
		for (auto neighbor : params.graph.inNeighbors[node])
		{
			if (current.hasNode(neighbor))
			{
				auto neighborSlice = current.node(neighbor).endSlice;
				assert(neighborSlice.getValue(offset) >= scoreHere-1);
				assert(neighborSlice.getValue(offset-1) >= scoreHere - (eq?0:1));
				if (neighborSlice.getValue(offset) == scoreHere-1)
				{
					return std::make_pair(std::make_pair(pos, false), std::make_pair(MatrixPosition { neighbor, params.graph.NodeLength(neighbor)-1, pos.seqPos }, true));
				}
				if (neighborSlice.getValue(offset-1) == scoreHere - (eq?0:1))
				{
					return std::make_pair(std::make_pair(pos, false), std::make_pair(MatrixPosition { neighbor, params.graph.NodeLength(neighbor)-1, pos.seqPos-1 }, true));
				}
			}
		}
		assert(false);
		return std::make_pair(std::make_pair(MatrixPosition {0, 0, 0}, false), std::make_pair(MatrixPosition {0, 0, 0}, false));
	}

	std::pair<std::pair<MatrixPosition, bool>, std::pair<MatrixPosition, bool>> pickBacktraceVerticalCrossing(const NodeSlice<LengthType, ScoreType, Word, false>& current, const NodeSlice<LengthType, ScoreType, Word, false>& previous, const std::vector<WordSlice> nodeScores, size_t j, LengthType node, MatrixPosition pos, const std::string& sequence, ScoreType quitScore, bool scoresNotValid) const
	{
		assert(pos.nodeOffset > 0);
		assert(pos.nodeOffset < nodeScores.size());
		while (pos.nodeOffset > 0 && nodeScores[pos.nodeOffset-1].getValue(0) == nodeScores[pos.nodeOffset].getValue(0) - 1)
		{
			pos.nodeOffset--;
		}
		if (pos.nodeOffset == 0)
		{
			return std::make_pair(std::make_pair(pos, false), pickBacktraceCorner(current, previous, node, j, sequence, quitScore, scoresNotValid));
		}
		assert(previous.hasNode(node));
		bool eq = Common::characterMatch(sequence[pos.seqPos], params.graph.NodeSequences(pos.node, pos.nodeOffset));
		auto previousNode = previous.node(node);
		ScoreType scoreHere = nodeScores[pos.nodeOffset].getValue(0);
		ScoreType scoreDiagonal = previousNode.startSlice.scoreEnd;
		for (size_t i = 1; i <= pos.nodeOffset - 1; i++)
		{
			scoreDiagonal += (previousNode.HP[i / WordConfiguration<Word>::WordSize] >> (i % WordConfiguration<Word>::WordSize)) & 1;
			scoreDiagonal -= (previousNode.HN[i / WordConfiguration<Word>::WordSize] >> (i % WordConfiguration<Word>::WordSize)) & 1;
		}
		ScoreType scoreUp = scoreDiagonal;
		scoreUp += (previousNode.HP[(pos.nodeOffset) / WordConfiguration<Word>::WordSize] >> ((pos.nodeOffset) % WordConfiguration<Word>::WordSize)) & 1;
		scoreUp -= (previousNode.HN[(pos.nodeOffset) / WordConfiguration<Word>::WordSize] >> ((pos.nodeOffset) % WordConfiguration<Word>::WordSize)) & 1;
		if (scoresNotValid || scoreHere > quitScore)
		{
			//this location is out of the band so the usual horizontal and vertical score limits don't apply
			//just pick the smallest scoring in-neighbor
			if (scoreDiagonal < scoreUp)
			{
				return std::make_pair(std::make_pair(pos, false), std::make_pair(MatrixPosition{pos.node, pos.nodeOffset - 1, pos.seqPos-1}, false));
			}
			else
			{
				return std::make_pair(std::make_pair(pos, false), std::make_pair(MatrixPosition{pos.node, pos.nodeOffset, pos.seqPos-1}, false));
			}
		}
		assert(scoreUp >= scoreHere - 1);
		assert(scoreDiagonal >= scoreHere - (eq?0:1));
		if (scoreUp == scoreHere - 1) return std::make_pair(std::make_pair(pos, false), std::make_pair(MatrixPosition{pos.node, pos.nodeOffset, pos.seqPos-1}, false));
		assert(scoreDiagonal == scoreHere - (eq?0:1));
		return std::make_pair(std::make_pair(pos, false), std::make_pair(MatrixPosition{pos.node, pos.nodeOffset - 1, pos.seqPos-1}, false));
	}

	std::pair<MatrixPosition, bool> pickBacktraceCorner(const NodeSlice<LengthType, ScoreType, Word, false>& current, const NodeSlice<LengthType, ScoreType, Word, false>& previous, LengthType node, size_t j, const std::string& sequence, ScoreType quitScore, bool scoresNotValid) const
	{
		ScoreType scoreHere = current.node(node).startSlice.getValue(0);
		if (scoresNotValid || scoreHere > quitScore)
		{
			//this location is out of the band so the usual horizontal and vertical score limits don't apply
			//just pick the smallest scoring in-neighbor
			ScoreType smallestFound = scoreHere+1;
			MatrixPosition smallestPos { 0, 0, 0 };
			bool nodeChange = false;
			if (previous.hasNode(node))
			{
				auto slice = previous.node(node).startSlice;
				smallestFound = slice.scoreEnd;
				smallestPos = MatrixPosition { node, 0, j-1 };
				nodeChange = false;
			}
			for (auto neighbor : params.graph.inNeighbors[node])
			{
				if (previous.hasNode(neighbor))
				{
					auto neighborSlice = previous.node(neighbor).endSlice;
					if (neighborSlice.scoreEnd < smallestFound)
					{
						smallestFound = neighborSlice.scoreEnd;
						smallestPos = MatrixPosition { neighbor, params.graph.NodeLength(neighbor)-1, j-1 };
						nodeChange = true;
					}
				}
				if (current.hasNode(neighbor) && neighbor != node)
				{
					auto neighborSlice = current.node(neighbor).endSlice;
					if (neighborSlice.getValue(0) < smallestFound)
					{
						smallestFound = neighborSlice.getValue(0);
						smallestPos = MatrixPosition { neighbor, params.graph.NodeLength(neighbor)-1, j };
						nodeChange = true;
					}
				}
			}
			return std::make_pair(smallestPos, nodeChange);
		}
		assert(scoreHere <= quitScore);
		bool eq = Common::characterMatch(sequence[j], params.graph.NodeSequences(node, 0));
		if (previous.hasNode(node))
		{
			assert(previous.node(node).startSlice.scoreEnd >= scoreHere-1);
			if (previous.node(node).startSlice.scoreEnd == scoreHere-1)
			{
				return std::make_pair(MatrixPosition {node, 0, j-1 }, false);
			}
		}
		for (auto neighbor : params.graph.inNeighbors[node])
		{
			if (current.hasNode(neighbor))
			{
				assert(current.node(neighbor).endSlice.getValue(0) >= scoreHere-1);
				if (current.node(neighbor).endSlice.getValue(0) == scoreHere-1)
				{
					return std::make_pair(MatrixPosition {neighbor, params.graph.NodeLength(neighbor)-1, j}, true);
				}
			}
			if (previous.hasNode(neighbor))
			{
				assert(previous.node(neighbor).endSlice.scoreEnd >= scoreHere-(eq?0:1));
				if (previous.node(neighbor).endSlice.scoreEnd == scoreHere-(eq?0:1))
				{
					return std::make_pair(MatrixPosition {neighbor, params.graph.NodeLength(neighbor)-1, j-1 }, true);
				}
			}
		}
		assert(false);
		return std::make_pair(MatrixPosition {0, 0, 0}, false);
	}

	WordSlice getSourceSliceFromScore(ScoreType previousScore) const
	{
		WordSlice result { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, previousScore+WordConfiguration<Word>::WordSize };
		return result;
	}

	class NodeCalculationResult
	{
	public:
		ScoreType minScore;
		LengthType minScoreNode;
		LengthType minScoreNodeOffset;
		size_t cellsProcessed;
#ifdef SLICEVERBOSE
		size_t nodesProcessed;
#endif
	};

	void assertSliceCorrectness(WordSlice oldSlice, WordSlice newSlice, Word Eq, int hin) const
	{
		ScoreType foundMinScore = newSlice.getScoreBeforeStart() + 1;
		foundMinScore = std::min(foundMinScore, oldSlice.getValue(0) + 1);
		foundMinScore = std::min(foundMinScore, oldSlice.getScoreBeforeStart() + ((Eq & 1) ? 0 : 1));
		assert(newSlice.getScoreBeforeStart() == oldSlice.getScoreBeforeStart() + hin);
		assert(newSlice.getValue(0) == foundMinScore);
		for (size_t i = 1; i < WordConfiguration<Word>::WordSize; i++)
		{
			foundMinScore = newSlice.getValue(i-1)+1;
			foundMinScore = std::min(foundMinScore, oldSlice.getValue(i)+1);
			foundMinScore = std::min(foundMinScore, oldSlice.getValue(i-1) + ((Eq & ((Word)1 << i)) ? 0 : 1));
			assert(newSlice.getValue(i) == foundMinScore);
		}
	}

	std::vector<WordSlice> recalcNodeWordslice(LengthType node, const typename NodeSlice<LengthType, ScoreType, Word, false>::NodeSliceMapItem& slice, const typename NodeSlice<LengthType, ScoreType, Word, false>::NodeSliceMapItem& previousSlice, LengthType j, const std::string& sequence) const
	{
		EqVector EqV = BV::getEqVector(sequence, j);
		std::vector<WordSlice> result;
		WordSlice ws = slice.startSlice;
		result.push_back(ws);

		ScoreType scoreBefore = ws.getScoreBeforeStart();
		ScoreType scoreComparison = 0;
		Word forceVP = WordConfiguration<Word>::AllOnes;
		Word forceVN = WordConfiguration<Word>::AllZeros;
		Word forceEq = WordConfiguration<Word>::AllOnes;
		if (!previousSlice.exists)
		{
			forceVP ^= 1;
			forceVN = 1;
			forceEq ^= 1;
		}
		else
		{
			scoreComparison = previousSlice.startSlice.scoreEnd;
		}

		size_t chunk = 0;
		size_t offset = 1;
		Word hinN, hinP, Eq;
		size_t pos;
		WordSlice newWs;
		size_t nodeLength = params.graph.NodeLength(node);
		for (; chunk < slice.NUM_CHUNKS; chunk++)
		{
			Word HP = previousSlice.HP[chunk];
			Word HN = previousSlice.HN[chunk];
			HP >>= offset;
			HN >>= offset;
			for (; offset < WordConfiguration<Word>::WordSize; offset++)
			{
				pos = chunk * WordConfiguration<Word>::WordSize + offset;
				if (pos >= nodeLength) break;
				Eq = EqV.getEqC(params.graph.NodeSequences(node, pos));
				Eq &= forceEq;
				if ((HN & 1) && (scoreBefore == scoreComparison - 1))
				{
					assert(previousSlice.exists);
					assert(ws.getScoreBeforeStart() == scoreBefore);
					std::tie(newWs, hinP, hinN) = BV::getNextSlice(Eq, ws, 0, 0);
#ifdef EXTRACORRECTNESSASSERTIONS
					assertSliceCorrectness(ws, newWs, Eq, 0);
#endif
				}
				else if (scoreBefore < scoreComparison)
				{
					assert(previousSlice.exists);
					assert(ws.getScoreBeforeStart() == scoreBefore);
					std::tie(newWs, hinP, hinN) = BV::getNextSlice(Eq, ws, 1, 0);
					newWs.VP &= ~(Word)1;
					newWs.VN |= 1;
#ifdef EXTRACORRECTNESSASSERTIONS
					assertSliceCorrectness(ws, newWs, Eq, 1);
#endif
				}
				else
				{
					assert(!previousSlice.exists || ws.getScoreBeforeStart() == scoreComparison);
					std::tie(newWs, hinP, hinN) = BV::getNextSlice(Eq, ws, HP & 1, HN & 1);
					newWs.VP &= forceVP;
					newWs.VN |= forceVN;
#ifdef EXTRACORRECTNESSASSERTIONS
					assertSliceCorrectness(ws, newWs, Eq, (HP & 1) - (HN & 1));
#endif
				}
				assert(forceVN == 0 || scoreBefore < scoreComparison || newWs.getScoreBeforeStart() == ws.getScoreBeforeStart() + 1);
				ws = newWs;
				result.push_back(ws);
				scoreComparison += HP & 1;
				scoreComparison -= HN & 1;
				HP >>= 1;
				HN >>= 1;
				scoreBefore++;
			}
			offset = 0;
		}
		assert(result.back().VP == slice.endSlice.VP);
		assert(result.back().VN == slice.endSlice.VN);
		assert(result.back().scoreEnd == slice.endSlice.scoreEnd);
		return result;
	}

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	template <typename NodeChunkType>
	NodeCalculationResult calculateNode(size_t i, typename NodeSlice<LengthType, ScoreType, Word, true>::NodeSliceMapItem& slice, const EqVector& EqV, typename NodeSlice<LengthType, ScoreType, Word, true>::NodeSliceMapItem previousSlice, const std::vector<EdgeWithPriority>& incoming, const std::vector<bool>& previousBand, NodeChunkType nodeChunks) const
	{
		assert(incoming.size() > 0);
		WordSlice newWs;
		WordSlice ws;
		bool hasWs = false;
		NodeCalculationResult result;
		result.minScore = std::numeric_limits<ScoreType>::max();
		result.minScoreNode = std::numeric_limits<LengthType>::max();
		result.minScoreNodeOffset = std::numeric_limits<LengthType>::max();
		result.cellsProcessed = 0;
		auto nodeLength = params.graph.NodeLength(i);

		Word Eq = EqV.getEqI(nodeChunks[0] & 3);
		bool hasSkipless = false;

		for (auto inc : incoming)
		{
			result.cellsProcessed++;
			if (inc.skipFirst)
			{
				if (!hasWs)
				{
					ws = inc.incoming;
					hasWs = true;
				}
				else
				{
					ws = ws.mergeWith(inc.incoming);
				}
				continue;
			}
			hasSkipless = true;
			Word hinP;
			Word hinN;
			if (previousSlice.exists)
			{
				ScoreType incomingScoreBeforeStart = inc.incoming.getScoreBeforeStart();
				if (previousSlice.startSlice.scoreEnd < incomingScoreBeforeStart)
				{
					hinP = 0;
					hinN = 1;
				}
				else if (previousSlice.startSlice.scoreEnd > incomingScoreBeforeStart)
				{
					hinP = 1;
					hinN = 0;
				}
				else
				{
					hinP = 0;
					hinN = 0;
				}
			}
			else
			{
				hinP = 1;
				hinN = 0;
			}

			WordSlice newWs;
			std::tie(newWs, hinP, hinN) = BV::getNextSlice(Eq, inc.incoming, hinP, hinN);
			if (!previousSlice.exists || newWs.getScoreBeforeStart() < previousSlice.startSlice.scoreEnd)
			{
				newWs.VP &= WordConfiguration<Word>::AllOnes ^ 1;
				newWs.VN |= 1;
			}
			assert(newWs.getScoreBeforeStart() >= debugLastRowMinScore);
			if (!hasWs)
			{
				ws = newWs;
				hasWs = true;
			}
			else
			{
				ws = ws.mergeWith(newWs);
			}
		}

		assert(hasWs);

		result.minScore = ws.scoreEnd;
		result.minScoreNode = i;
		result.minScoreNodeOffset = 0;

		if (slice.exists)
		{
			if (hasSkipless && params.graph.inNeighbors[i].size() == 1 && previousBand[params.graph.inNeighbors[i][0]])
			{
				if (ws.scoreEnd > slice.startSlice.scoreEnd)
				{
#ifdef EXTRACORRECTNESSASSERTIONS
					auto debugTest = ws.mergeWith(slice.startSlice);
					assert(debugTest.VP == slice.startSlice.VP);
					assert(debugTest.VN == slice.startSlice.VN);
					assert(debugTest.scoreEnd == slice.startSlice.scoreEnd);
#endif
					return result;
				}
				else if (ws.scoreEnd < slice.startSlice.scoreEnd)
				{
#ifdef EXTRACORRECTNESSASSERTIONS
					auto debugTest = ws.mergeWith(slice.startSlice);
					assert(debugTest.VP == ws.VP);
					assert(debugTest.VN == ws.VN);
					assert(debugTest.scoreEnd == ws.scoreEnd);
#endif
				}
				else
				{
					Word newBigger = (ws.VP & ~slice.startSlice.VP) | (slice.startSlice.VN & ~ws.VN);
					Word oldBigger = (slice.startSlice.VP & ~ws.VP) | (ws.VN & ~slice.startSlice.VN);
					if (newBigger > oldBigger)
					{
#ifdef EXTRACORRECTNESSASSERTIONS
						auto debugTest = ws.mergeWith(slice.startSlice);
						assert(debugTest.VP == ws.VP);
						assert(debugTest.VN == ws.VN);
						assert(debugTest.scoreEnd == ws.scoreEnd);
#endif
					}
					else if (oldBigger > newBigger)
					{
#ifdef EXTRACORRECTNESSASSERTIONS
						auto debugTest = ws.mergeWith(slice.startSlice);
						assert(debugTest.VP == slice.startSlice.VP);
						assert(debugTest.VN == slice.startSlice.VN);
						assert(debugTest.scoreEnd == slice.startSlice.scoreEnd);
#endif
						return result;
					}
					else if (newBigger == 0 && oldBigger == 0)
					{
						assert(ws.VP == slice.startSlice.VP);
						assert(ws.VN == slice.startSlice.VN);
						assert(ws.scoreEnd == slice.startSlice.scoreEnd);
						return result;
					}
					else
					{
						WordSlice test = ws.mergeWith(slice.startSlice);
						if (test.scoreEnd == slice.startSlice.scoreEnd && test.VP == slice.startSlice.VP && test.VN == slice.startSlice.VN)
						{
							return result;
						}
						ws = test;
					}
				}
			}
			else
			{
				WordSlice test = ws.mergeWith(slice.startSlice);
				if (test.scoreEnd == slice.startSlice.scoreEnd && test.VP == slice.startSlice.VP && test.VP == slice.startSlice.VN)
				{
					return result;
				}
				ws = test;
			}
		}

		if (hasSkipless && previousSlice.exists)
		{
			if (ws.getScoreBeforeStart() > previousSlice.startSlice.scoreEnd)
			{
				//todo vertical merge
				ws = ws.mergeWith(getSourceSliceFromScore(previousSlice.startSlice.scoreEnd));
			}
		}

		for (size_t i = 0; i < slice.NUM_CHUNKS; i++)
		{
			slice.HP[i] = WordConfiguration<Word>::AllZeros;
			slice.HN[i] = WordConfiguration<Word>::AllZeros;
		}

		LengthType pos = 1;
		size_t forceUntil = 0;
		if (previousSlice.exists)
		{
			ScoreType scoreBefore = ws.getScoreBeforeStart();
			ScoreType scoreComparison = previousSlice.startSlice.scoreEnd;
			assert(scoreBefore <= scoreComparison);
			if (scoreBefore < scoreComparison)
			{
				size_t fixoffset = 1;
				for (size_t fixchunk = 0; fixchunk < slice.NUM_CHUNKS; fixchunk++)
				{
					for (; fixoffset < WordConfiguration<Word>::WordSize; fixoffset++)
					{
						ScoreType newScoreComparison = scoreComparison;
						newScoreComparison += (previousSlice.HP[fixchunk] >> fixoffset) & 1;
						newScoreComparison -= (previousSlice.HN[fixchunk] >> fixoffset) & 1;
						Word mask = ((Word)1) << fixoffset;
						assert(scoreBefore <= newScoreComparison);
						if (scoreBefore < newScoreComparison)
						{
							previousSlice.HP[fixchunk] |= mask;
							previousSlice.HN[fixchunk] &= ~mask;
							forceUntil = fixchunk * WordConfiguration<Word>::WordSize + fixoffset;
						}
						if (scoreBefore == newScoreComparison)
						{
							previousSlice.HP[fixchunk] &= ~mask;
							previousSlice.HN[fixchunk] &= ~mask;
						}
						scoreBefore++;
						scoreComparison = newScoreComparison;
						if (scoreBefore >= scoreComparison) break;
					}
					if (scoreBefore >= scoreComparison) break;
					fixoffset = 0;
				}
			}
		}
		else
		{
			forceUntil = nodeLength;
		}
		slice.startSlice = ws;
		slice.exists = true;
		Word forceEq = WordConfiguration<Word>::AllOnes;
		Word hinP, hinN;
		if (!previousSlice.exists) forceEq ^= 1;
		size_t smallChunk = 0;
		size_t offset = 1;
		pos = smallChunk * (WordConfiguration<Word>::WordSize / 2) + offset;
		for (; smallChunk < params.graph.CHUNKS_IN_NODE; smallChunk++)
		{
			size_t bigChunk = smallChunk / 2;
			size_t bigChunkOffset = (smallChunk % 2) * (WordConfiguration<Word>::WordSize / 2);
			Word HP = previousSlice.HP[bigChunk] >> bigChunkOffset;
			Word HN = previousSlice.HN[bigChunk] >> bigChunkOffset;
			auto charChunk = nodeChunks[smallChunk];
			HP >>= offset;
			HN >>= offset;
			charChunk >>= offset * 2;
			for (; offset < WordConfiguration<Word>::WordSize / 2 && pos < nodeLength; offset++)
			{
				Eq = EqV.getEqI(charChunk & 3);
				Eq &= forceEq;
				std::tie(newWs, hinP, hinN) = BV::getNextSlice(Eq, ws, HP & 1, HN & 1);
				if (forceUntil >= pos)
				{
					newWs.VP &= WordConfiguration<Word>::AllOnes ^ 1;
					newWs.VN |= 1;
				}
#ifdef EXTRACORRECTNESSASSERTIONS
				assertSliceCorrectness(ws, newWs, Eq, (HP & 1) - (HN & 1));
#endif
				assert(newWs.getScoreBeforeStart() >= debugLastRowMinScore);
				ws = newWs;
				if (ws.scoreEnd < result.minScore)
				{
					result.minScore = ws.scoreEnd;
					result.minScoreNodeOffset = pos;
				}
				charChunk >>= 2;
				HP >>= 1;
				HN >>= 1;
				pos++;
				slice.HP[bigChunk] |= hinP << (offset + bigChunkOffset);
				slice.HN[bigChunk] |= hinN << (offset + bigChunkOffset);
			}
			offset = 0;
		}
		result.cellsProcessed = pos;
		slice.endSlice = ws;
#ifndef NDEBUG
		if (previousSlice.exists && forceUntil < nodeLength - 1) assert(slice.endSlice.getScoreBeforeStart() == previousSlice.endSlice.scoreEnd);
#endif
		return result;
	}

#ifdef EXTRACORRECTNESSASSERTIONS
	void assertBitvectorConfirmedAreConsistent(WordSlice newslice, WordSlice oldslice, ScoreType quitScore) const
	{
		assert(newslice.scoreBeforeStart <= oldslice.scoreBeforeStart);
		for (int i = 0; i < 64; i++)
		{
			auto newScore = newslice.getValue(i);
			auto oldScore = oldslice.getValue(i);
			if (oldScore <= quitScore) assert(newslice.getValue(i) <= oldslice.getValue(i));
		}
	}

	template <bool HasVectorMap, bool PreviousHasVectorMap>
	void checkNodeBoundaryCorrectness(const NodeSlice<LengthType, ScoreType, Word, HasVectorMap>& currentSlice, const NodeSlice<LengthType, ScoreType, Word, PreviousHasVectorMap>& previousSlice, const std::string& sequence, size_t j, ScoreType maxScore, ScoreType previousMaxScore) const
	{
		for (auto pair : currentSlice)
		{
			auto node = pair.first;
			if (previousSlice.hasNode(node) && previousSlice.node(node).exists)
			{
				assert(pair.second.startSlice.getScoreBeforeStart() <= previousSlice.node(node).startSlice.scoreEnd);
			}
			bool eq = Common::characterMatch(sequence[j], params.graph.NodeSequences(node, 0));
			ScoreType foundMinScore = std::numeric_limits<ScoreType>::max();
			foundMinScore = std::min(foundMinScore, currentSlice.node(node).startSlice.getScoreBeforeStart()+1);
			if (previousSlice.hasNode(node) && previousSlice.node(node).exists && previousSlice.node(node).minScore <= previousMaxScore)
			{
				foundMinScore = std::min(foundMinScore, previousSlice.node(node).startSlice.scoreEnd+1);
			}
			for (auto neighbor : params.graph.inNeighbors[node])
			{
				if (currentSlice.hasNode(neighbor) && currentSlice.node(neighbor).exists)
				{
					foundMinScore = std::min(foundMinScore, currentSlice.node(neighbor).endSlice.getValue(0)+1);
					foundMinScore = std::min(foundMinScore, currentSlice.node(neighbor).endSlice.getScoreBeforeStart() + (eq?0:1));
				}
				if (previousSlice.hasNode(neighbor) && previousSlice.node(neighbor).exists && previousSlice.node(neighbor).minScore <= previousMaxScore)
				{
					foundMinScore = std::min(foundMinScore, previousSlice.node(neighbor).endSlice.scoreEnd + (eq ? 0 : 1));
				}
			}
			if (pair.second.startSlice.getValue(0) <= maxScore || foundMinScore <= maxScore)
			{
				assert(foundMinScore != std::numeric_limits<ScoreType>::max());
				assert(pair.second.startSlice.getValue(0) == foundMinScore);
			}
			bool hasMatch = Common::characterMatch(sequence[0], params.graph.NodeSequences(node, 0));
			for (size_t i = 1; i < WordConfiguration<Word>::WordSize; i++)
			{
				if (j+i >= sequence.size()) break;
				eq = Common::characterMatch(sequence[j+i], params.graph.NodeSequences(node, 0));
				if (eq) hasMatch = true;
				ScoreType foundMinScore = pair.second.startSlice.getValue(i-1)+1;
				for (auto neighbor : params.graph.inNeighbors[node])
				{
					if (!currentSlice.hasNode(neighbor)) continue;
					if (!currentSlice.node(neighbor).exists) continue;
					foundMinScore = std::min(foundMinScore, currentSlice.node(neighbor).endSlice.getValue(i)+1);
					foundMinScore = std::min(foundMinScore, currentSlice.node(neighbor).endSlice.getValue(i-1) + (eq ? 0 : 1));
				}
				if (pair.second.startSlice.getValue(i) <= maxScore || foundMinScore <= maxScore)
				{
					assert(pair.second.startSlice.getValue(i) == foundMinScore);
				}
			}
		}
	}
#endif

	template <bool HasVectorMap, bool PreviousHasVectorMap, typename PriorityQueue>
	NodeCalculationResult calculateSlice(const std::string& sequence, size_t j, NodeSlice<LengthType, ScoreType, Word, HasVectorMap>& currentSlice, const NodeSlice<LengthType, ScoreType, Word, PreviousHasVectorMap>& previousSlice, std::vector<bool>& currentBand, const std::vector<bool>& previousBand, PriorityQueue& calculableQueue, ScoreType previousQuitScore, int bandwidth, ScoreType previousMinScore) const
	{
		ScoreType currentMinimumScore = std::numeric_limits<ScoreType>::max() - bandwidth - 1;
		LengthType currentMinimumNode = -1;
		LengthType currentMinimumNodeOffset = -1;
		size_t cellsProcessed = 0;
#ifdef SLICEVERBOSE
		size_t nodesProcessed = 0;
#endif

		EqVector EqV = BV::getEqVector(sequence, j);

		assert(previousSlice.size() > 0);
		ScoreType zeroScore = previousMinScore - j / 2 - 32;
		if (j == 0)
		{
			for (auto node : previousSlice)
			{
				assert(node.second.minScore <= previousQuitScore);
				WordSlice startSlice = getSourceSliceFromScore(node.second.startSlice.scoreEnd);
				if (std::is_same<decltype(calculableQueue), ComponentPriorityQueue<EdgeWithPriority>&>::value)
				{
					calculableQueue.insert(params.graph.componentNumber[node.first], node.second.minScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
				}
				else
				{
					calculableQueue.insert(node.second.minScore - j/2 - zeroScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
				}
			}
		}
		else
		{
			for (auto node : previousSlice)
			{
				assert(node.second.exists);
				if (node.second.minScore > previousQuitScore) continue;
				if (params.graph.inNeighbors[node.first].size() == 1)
				{
					auto neighbor = params.graph.inNeighbors[node.first][0];
				 	if (previousBand[params.graph.inNeighbors[node.first][0]] && previousSlice.node(neighbor).endSlice.scoreEnd < previousQuitScore)
				 	{
				 		//linear area, no need to add the later node into the queue 
				 		//because calculating the earlier node will guarantee that the later node will get added
				 		continue;
				 	}
				}
				WordSlice startSlice = getSourceSliceFromScore(node.second.startSlice.scoreEnd);
				if (std::is_same<decltype(calculableQueue), ComponentPriorityQueue<EdgeWithPriority>&>::value)
				{
					calculableQueue.insert(params.graph.componentNumber[node.first], node.second.minScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
				}
				else
				{
					calculableQueue.insert(node.second.minScore - j/2 - zeroScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
				}
			}
		}
		assert(calculableQueue.size() > 0);
		
		ScoreType currentMinScoreAtEndRow = currentMinimumScore;
		while (calculableQueue.size() > 0)
		{
			auto pair = calculableQueue.top();
			if (!std::is_same<decltype(calculableQueue), ComponentPriorityQueue<EdgeWithPriority>&>::value)
			{
				if (pair.priority > currentMinScoreAtEndRow + bandwidth) break;
			}
			if (calculableQueue.extraSize(pair.target) == 0)
			{
				calculableQueue.pop();
				continue;
			}
			auto i = pair.target;
			if (!currentBand[i])
			{
				assert(!currentSlice.hasNode(i));
				currentSlice.addNode(i);
				currentBand[i] = true;
			}
			assert(currentBand[i]);
			const std::vector<EdgeWithPriority>* extras;
			extras = &calculableQueue.getExtras(i);
			auto& thisNode = currentSlice.node(i);
			auto oldEnd = thisNode.endSlice;
			if (!thisNode.exists) oldEnd = { 0, 0, std::numeric_limits<ScoreType>::max() };
			typename NodeSlice<LengthType, ScoreType, Word, HasVectorMap>::NodeSliceMapItem previousThisNode;

			if (previousBand[i])
			{
				previousThisNode = previousSlice.node(i);
				assert(previousThisNode.exists);
			}
			else
			{
				for (size_t chunk = 0; chunk < previousThisNode.NUM_CHUNKS; chunk++)
				{
					previousThisNode.HP[chunk] = WordConfiguration<Word>::AllOnes;
					previousThisNode.HN[chunk] = WordConfiguration<Word>::AllZeros;
				}
				previousThisNode.exists = false;
			}
			NodeCalculationResult nodeCalc;
			if (i < params.graph.firstAmbiguous)
			{
				nodeCalc = calculateNode(i, thisNode, EqV, previousThisNode, *extras, previousBand, params.graph.NodeChunks(i));
			}
			else
			{
				nodeCalc = calculateNode(i, thisNode, EqV, previousThisNode, *extras, previousBand, params.graph.AmbiguousNodeChunks(i));
			}
			calculableQueue.pop();
			if (!std::is_same<decltype(calculableQueue), ComponentPriorityQueue<EdgeWithPriority>&>::value)
			{
				calculableQueue.removeExtras(i);
			}
			assert(nodeCalc.minScore <= previousQuitScore + bandwidth + params.graph.SPLIT_NODE_SIZE + WordConfiguration<Word>::WordSize);
			currentMinScoreAtEndRow = std::min(currentMinScoreAtEndRow, nodeCalc.minScore);
			currentSlice.setMinScoreIfSmaller(i, nodeCalc.minScore);
#ifdef SLICEVERBOSE
			volatile size_t firstslices = currentSlice.node(i).firstSlicesCalcedWhenCalced;
			volatile size_t calcedslices = currentSlice.node(i).slicesCalcedWhenCalced;
			if (currentSlice.node(i).firstSlicesCalcedWhenCalced == std::numeric_limits<size_t>::max()) currentSlice.node(i).firstSlicesCalcedWhenCalced = cellsProcessed;
			if (currentSlice.node(i).slicesCalcedWhenCalced != std::numeric_limits<size_t>::max()) assert(currentSlice.node(i).slicesCalcedWhenCalced < cellsProcessed);
			currentSlice.node(i).slicesCalcedWhenCalced = cellsProcessed;
			assert(currentSlice.node(i).firstSlicesCalcedWhenCalced <= currentSlice.node(i).slicesCalcedWhenCalced);
#endif
			auto newEnd = thisNode.endSlice;

			if (newEnd.scoreEnd != oldEnd.scoreEnd || newEnd.VP != oldEnd.VP || newEnd.VN != oldEnd.VN)
			{
				ScoreType newEndMinScore = newEnd.changedMinScore(oldEnd);
				assert(newEndMinScore >= previousMinScore);
				assert(newEndMinScore != std::numeric_limits<ScoreType>::max());
				if (newEndMinScore <= currentMinScoreAtEndRow + bandwidth)
				{
					for (auto neighbor : params.graph.outNeighbors[i])
					{
						if (std::is_same<decltype(calculableQueue), ComponentPriorityQueue<EdgeWithPriority>&>::value)
						{
							calculableQueue.insert(params.graph.componentNumber[neighbor], newEndMinScore, EdgeWithPriority { neighbor, newEndMinScore - previousMinScore, newEnd, false });
						}
						else
						{
							ScoreType newEndPriorityScore = newEnd.getPriorityScore(j);
							assert(newEndPriorityScore >= zeroScore);
							calculableQueue.insert(newEndPriorityScore - zeroScore, EdgeWithPriority { neighbor, newEndMinScore - previousMinScore, newEnd, false });
						}
					}
				}
			}
			if (nodeCalc.minScore < currentMinimumScore)
			{
				currentMinimumScore = nodeCalc.minScore;
				currentMinimumNode = nodeCalc.minScoreNode;
				currentMinimumNodeOffset = nodeCalc.minScoreNodeOffset;
			}
			assert(currentMinimumScore == currentMinScoreAtEndRow);
			cellsProcessed += nodeCalc.cellsProcessed;
			assert(nodeCalc.cellsProcessed > 0);
#ifdef SLICEVERBOSE
			nodesProcessed++;
#endif
			if (cellsProcessed > params.maxCellsPerSlice) break;
		}

#ifdef EXTRACORRECTNESSASSERTIONS
		checkNodeBoundaryCorrectness<HasVectorMap, PreviousHasVectorMap>(currentSlice, previousSlice, sequence, j, currentMinScoreAtEndRow + bandwidth, previousQuitScore);
#endif

		assert(currentMinimumNode != std::numeric_limits<LengthType>::max());
		NodeCalculationResult result;
		result.minScore = currentMinimumScore;
		result.minScoreNode = currentMinimumNode;
		result.minScoreNodeOffset = currentMinimumNodeOffset;
		result.cellsProcessed = cellsProcessed;
#ifdef SLICEVERBOSE
		result.nodesProcessed = nodesProcessed;
#endif

		if (j + WordConfiguration<Word>::WordSize > sequence.size())
		{
			flattenLastSliceEnd<HasVectorMap, PreviousHasVectorMap>(currentSlice, previousSlice, result, j, sequence);
		}

#ifdef SLICEVERBOSE
		std::cerr << "prefilternodes " << currentSlice.size() << " ";
#endif

		//this sometimes removes nodes which must be used for out-of-band backtrace
		//comment until we figure out a solution
		// finalizeSlice(currentSlice, currentBand, currentMinScoreAtEndRow + bandwidth);

		calculableQueue.clear();

		return result;
	}

	template <bool HasVectorMap>
	void finalizeSlice(NodeSlice<LengthType, ScoreType, Word, HasVectorMap>& slice, std::vector<bool>& currentBand, ScoreType maxScore) const
	{
		for (auto node : slice)
		{
			if (node.second.minScore > maxScore && node.second.endSlice.getMinScore() > maxScore)
			{
				currentBand[node.first] = false;
				slice.node(node.first).exists = false;
			}
		}
		slice.removeNonExistant();
	}

	template <bool HasVectorMap, bool PreviousHasVectorMap>
	void flattenLastSliceEnd(NodeSlice<LengthType, ScoreType, Word, HasVectorMap>& slice, const NodeSlice<LengthType, ScoreType, Word, PreviousHasVectorMap>& previousSlice, NodeCalculationResult& sliceCalc, LengthType j, const std::string& sequence) const
	{
		assert(j < sequence.size());
		assert(sequence.size() - j < WordConfiguration<Word>::WordSize);
		sliceCalc.minScore = std::numeric_limits<ScoreType>::max();
		sliceCalc.minScoreNode = std::numeric_limits<LengthType>::max();
		sliceCalc.minScoreNodeOffset = std::numeric_limits<LengthType>::max();
		auto offset = sequence.size() - j;
		assert(offset >= 0);
		assert(offset < WordConfiguration<Word>::WordSize);
		for (auto node : slice)
		{
			auto current = node.second;
			decltype(current) old;
			if (previousSlice.hasNode(node.first))
			{
				old = previousSlice.node(node.first);
			}
			else
			{
				old.exists = false;
				for (size_t i = 0; i < old.NUM_CHUNKS; i++)
				{
					old.HP[i] = WordConfiguration<Word>::AllOnes;
					old.HN[i] = WordConfiguration<Word>::AllZeros;
				}
			}
			auto nodeSlices = recalcNodeWordslice(node.first, current, old, j, sequence);
			assert(nodeSlices[0].VP == node.second.startSlice.VP);
			assert(nodeSlices[0].VN == node.second.startSlice.VN);
			assert(nodeSlices[0].scoreEnd == node.second.startSlice.scoreEnd);
			assert(nodeSlices.back().VP == node.second.endSlice.VP);
			assert(nodeSlices.back().VN == node.second.endSlice.VN);
			assert(nodeSlices.back().scoreEnd == node.second.endSlice.scoreEnd);
			for (size_t i = 0; i < nodeSlices.size(); i++)
			{
				auto wordSliceResult = BV::flattenWordSlice(nodeSlices[i], offset);
				if (wordSliceResult.scoreEnd < sliceCalc.minScore)
				{
					sliceCalc.minScore = wordSliceResult.scoreEnd;
					sliceCalc.minScoreNode = node.first;
					sliceCalc.minScoreNodeOffset = i;
				}
			}
		}
		assert(sliceCalc.minScore != std::numeric_limits<ScoreType>::max());
		assert(sliceCalc.minScoreNode < params.graph.NodeSize());
		assert(sliceCalc.minScoreNodeOffset < params.graph.NodeLength(sliceCalc.minScoreNode));
	}

	template <typename PriorityQueue>
	void fillDPSlice(const std::string& sequence, DPSlice& slice, const DPSlice& previousSlice, const std::vector<bool>& previousBand, std::vector<bool>& currentBand, PriorityQueue& calculableQueue, int bandwidth) const
	{
		NodeCalculationResult sliceResult;
		assert((ScoreType)previousSlice.bandwidth < std::numeric_limits<ScoreType>::max());
		assert((ScoreType)previousSlice.bandwidth >= 0);
		assert(previousSlice.minScore < std::numeric_limits<ScoreType>::max() - (ScoreType)previousSlice.bandwidth);
		if (slice.scoresVectorMap.hasVectorMapCurrently())
		{
			if (previousSlice.scoresVectorMap.hasVectorMapCurrently())
			{
				sliceResult = calculateSlice<true, true>(sequence, slice.j, slice.scoresVectorMap, previousSlice.scoresVectorMap, currentBand, previousBand, calculableQueue, previousSlice.minScore + previousSlice.bandwidth, bandwidth, previousSlice.minScore);
			}
			else
			{
				sliceResult = calculateSlice<true, false>(sequence, slice.j, slice.scoresVectorMap, previousSlice.scores, currentBand, previousBand, calculableQueue, previousSlice.minScore + previousSlice.bandwidth, bandwidth, previousSlice.minScore);
			}
			slice.scores = slice.scoresVectorMap.getMapSlice();
		}
		else
		{
			assert(!previousSlice.scoresVectorMap.hasVectorMapCurrently());
			sliceResult = calculateSlice<false, false>(sequence, slice.j, slice.scores, previousSlice.scores, currentBand, previousBand, calculableQueue, previousSlice.minScore + previousSlice.bandwidth, bandwidth, previousSlice.minScore);
		}
		slice.cellsProcessed = sliceResult.cellsProcessed;
		slice.minScoreNode = sliceResult.minScoreNode;
		slice.minScoreNodeOffset = sliceResult.minScoreNodeOffset;
		slice.minScore = sliceResult.minScore;
		assert(slice.minScore >= previousSlice.minScore);
		slice.correctness = slice.correctness.NextState(slice.minScore - previousSlice.minScore, WordConfiguration<Word>::WordSize);
		slice.bandwidth = bandwidth;
#ifdef SLICEVERBOSE
		slice.nodesProcessed = sliceResult.nodesProcessed;
		for (auto node : slice.scores)
		{
			if (currentBand[node.first])
			{
				slice.numCells += params.graph.NodeLength(node.first);
			}
		}
#endif
	}

	template <typename PriorityQueue>
	DPSlice pickMethodAndExtendFill(const std::string& sequence, const DPSlice& previous, const std::vector<bool>& previousBand, std::vector<bool>& currentBand, std::vector<typename NodeSlice<LengthType, ScoreType, Word, true>::MapItem>& nodesliceMap, PriorityQueue& calculableQueue, int bandwidth) const
	{
		if (!params.lowMemory)
		{
			DPSlice bandTest { &nodesliceMap };
			bandTest.j = previous.j + WordConfiguration<Word>::WordSize;
			bandTest.correctness = previous.correctness;
			fillDPSlice(sequence, bandTest, previous, previousBand, currentBand, calculableQueue, bandwidth);
			return bandTest;
		}
		else
		{
			DPSlice bandTest;
			bandTest.scores.addEmptyNodeMap(previous.scores.size());
			bandTest.j = previous.j + WordConfiguration<Word>::WordSize;
			bandTest.correctness = previous.correctness;
			fillDPSlice(sequence, bandTest, previous, previousBand, currentBand, calculableQueue, bandwidth);
			return bandTest;
		}
	}

	void removeWronglyAlignedEnd(DPTable& table) const
	{
		if (table.slices.size() == 0) return;
		bool currentlyCorrect = table.slices.back().correctness.CurrentlyCorrect();
		while (!currentlyCorrect)
		{
			currentlyCorrect = table.slices.back().correctness.FalseFromCorrect();
			table.slices.pop_back();
			if (table.slices.size() == 0) break;
		}
	}

	DPTable getSqrtSlices(const std::string& sequence, const DPSlice& initialSlice, size_t numSlices, bool forceGlobal, AlignerGraphsizedState& reusableState) const
	{
		assert(initialSlice.j == (size_t)-WordConfiguration<Word>::WordSize);
		assert(initialSlice.j + numSlices * WordConfiguration<Word>::WordSize <= sequence.size() + WordConfiguration<Word>::WordSize);
		DPTable result;
		result.slices.reserve(numSlices + 1);
		size_t cellsProcessed = 0;
		std::vector<size_t> partOfComponent;
		{
			for (auto node : initialSlice.scores)
			{
				reusableState.previousBand[node.first] = true;
			}
		}
#ifndef NDEBUG
		debugLastRowMinScore = 0;
#endif
		DPSlice lastSlice = initialSlice;
		result.slices.push_back(initialSlice);
		assert(lastSlice.correctness.CurrentlyCorrect());
		DPSlice rampSlice = lastSlice;
		size_t rampRedoIndex = -1;
		size_t rampUntil = 0;
#ifndef NDEBUG
		volatile size_t debugLastProcessedSlice;
		// we want to keep this variable for debugging purposes
		// useless self-assignment to prevent unused variable compilation warning
		debugLastProcessedSlice = debugLastProcessedSlice;
#endif
		for (size_t slice = 0; slice < numSlices; slice++)
		{
			int bandwidth = (params.rampBandwidth > params.initialBandwidth && rampUntil >= slice) ? params.rampBandwidth : params.initialBandwidth;
#ifndef NDEBUG
			debugLastProcessedSlice = slice;
			debugLastRowMinScore = lastSlice.minScore;
#endif
#ifdef SLICEVERBOSE
			auto timeStart = std::chrono::system_clock::now();
#endif
			DPSlice newSlice;
			if (reusableState.componentQueue.valid())
			{
				newSlice = pickMethodAndExtendFill(sequence, lastSlice, reusableState.previousBand, reusableState.currentBand, (slice % 2 == 0) ? reusableState.evenNodesliceMap : reusableState.oddNodesliceMap, reusableState.componentQueue, bandwidth);
			}
			else
			{
				newSlice = pickMethodAndExtendFill(sequence, lastSlice, reusableState.previousBand, reusableState.currentBand, (slice % 2 == 0) ? reusableState.evenNodesliceMap : reusableState.oddNodesliceMap, reusableState.calculableQueue, bandwidth);
			}
#ifdef SLICEVERBOSE
			auto timeEnd = std::chrono::system_clock::now();
			auto time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
			std::cerr << "slice " << slice << " bandwidth " << bandwidth << " minscore " << newSlice.minScore << " diff " << (newSlice.minScore - lastSlice.minScore) << " time " << time << " nodes " << newSlice.scores.size() << " slices " << newSlice.numCells << " nodesprocessed " << newSlice.nodesProcessed << " cellsprocessed " << newSlice.cellsProcessed << " overhead " << (100 * (int)(newSlice.cellsProcessed - newSlice.numCells) / (int)(newSlice.numCells)) << "%";
#endif
			assert(newSlice.minScore != std::numeric_limits<ScoreType>::max());
			assert(newSlice.minScoreNode != std::numeric_limits<LengthType>::max());
			assert(newSlice.minScoreNodeOffset != std::numeric_limits<LengthType>::max());
			assert(newSlice.scores.hasNode(newSlice.minScoreNode));
			assert(newSlice.minScoreNodeOffset < params.graph.NodeLength(newSlice.minScoreNode));

			if ((rampUntil == slice-1 || (rampUntil < slice && newSlice.correctness.CurrentlyCorrect() && newSlice.correctness.FalseFromCorrect())))
			{
				rampSlice = lastSlice.getMapSlice();
				rampRedoIndex = slice-1;
			}
			assert(newSlice.j == lastSlice.j + WordConfiguration<Word>::WordSize);

			cellsProcessed += newSlice.cellsProcessed;

			if (newSlice.cellsProcessed > params.maxCellsPerSlice)
			{
				newSlice.scoresNotValid = true;
			}

			// if forcing global alignment, the whole read must be aligned -> don't allow breaking even if the alignment looks wrong
			if (!forceGlobal)
			{
				if (!newSlice.correctness.CorrectFromCorrect())
				{
#ifndef NDEBUG
					debugLastProcessedSlice = slice-1;
#endif
					for (auto node : lastSlice.scores)
					{
						assert(reusableState.previousBand[node.first]);
						reusableState.previousBand[node.first] = false;
					}
					for (auto node : newSlice.scores)
					{
						assert(reusableState.currentBand[node.first]);
						reusableState.currentBand[node.first] = false;
					}
					lastSlice.scoresVectorMap.removeVectorArray();
					newSlice.scoresVectorMap.removeVectorArray();
					break;
				}
				if (!newSlice.correctness.CurrentlyCorrect() && rampUntil < slice && params.rampBandwidth > params.initialBandwidth)
				{
					for (auto node : newSlice.scores)
					{
						assert(reusableState.currentBand[node.first]);
						reusableState.currentBand[node.first] = false;
					}
					for (auto node : lastSlice.scores)
					{
						assert(reusableState.previousBand[node.first]);
						reusableState.previousBand[node.first] = false;
					}
					lastSlice.scoresVectorMap.removeVectorArray();
					newSlice.scoresVectorMap.removeVectorArray();
					rampUntil = slice;
					std::swap(slice, rampRedoIndex);
					std::swap(lastSlice, rampSlice);
					for (auto node : lastSlice.scores)
					{
						assert(!reusableState.previousBand[node.first]);
						reusableState.previousBand[node.first] = true;
					}
					if (slice == (size_t)-1)
					{
						result.slices.clear();
					}
					while (result.slices.size() > 1 && result.slices.back().j > slice * WordConfiguration<Word>::WordSize) result.slices.pop_back();
					assert(slice == (size_t)-1 || result.slices.size() == slice+2);
					assert(result.slices.back().j == lastSlice.j);
#ifdef SLICEVERBOSE
					std::cerr << " ramp to " << slice;
					std::cerr << std::endl;
					if (result.slices.size() > 0) std::cerr << " slices.back().j " << result.slices.back().j; else std::cerr << " slices.size() 0";
					std::cerr << std::endl;
#endif
					continue;
				}
			}

#ifdef SLICEVERBOSE
			std::cerr << std::endl;
#endif

			result.slices.push_back(newSlice.getMapSlice());
			for (auto node : lastSlice.scores)
			{
				assert(reusableState.previousBand[node.first]);
				reusableState.previousBand[node.first] = false;
			}
			assert(newSlice.minScore != std::numeric_limits<ScoreType>::max());
			assert(newSlice.minScore >= lastSlice.minScore);
			if (slice == numSlices - 1)
			{
				for (auto node : newSlice.scores)
				{
					assert(reusableState.currentBand[node.first]);
					reusableState.currentBand[node.first] = false;
				}
			}
			else
			{
				std::swap(reusableState.previousBand, reusableState.currentBand);
			}
			lastSlice.scoresVectorMap.removeVectorArray();
			lastSlice = std::move(newSlice);
		}
		lastSlice.scoresVectorMap.removeVectorArray();

		assert(result.slices.size() <= numSlices + 1);

#ifdef EXTRACORRECTNESSASSERTIONS
		assert(reusableState.calculableQueue.size() == 0);
		for (size_t i = 0; i < reusableState.currentBand.size(); i++)
		{
			assert(!reusableState.currentBand[i]);
			assert(!reusableState.previousBand[i]);
		}
#endif

#ifndef NDEBUG
		if (result.slices.size() > 0)
		{
			for (size_t i = 1; i < result.slices.size(); i++)
			{
				assert(result.slices[i].j == result.slices[i-1].j + WordConfiguration<Word>::WordSize);
			}
			for (size_t i = 1; i < result.slices.size(); i++)
			{
				assert(result.slices[i].minScore >= result.slices[i-1].minScore);
			}
		}
#endif
		return result;
	}

	DPSlice getInitialSliceExactPosition(LengthType bigraphNodeId, size_t offset) const
	{
		DPSlice result;
		result.j = -WordConfiguration<Word>::WordSize;
		result.bandwidth = 1;
		result.minScore = 0;
		result.scores.addEmptyNodeMap(1);
		auto nodes = params.graph.nodeLookup.at(bigraphNodeId);
		assert(offset < params.graph.originalNodeSize.at(bigraphNodeId));
		size_t nodeIndex = params.graph.GetUnitigNode(bigraphNodeId, offset);
		assert(params.graph.nodeOffset[nodeIndex] <= offset);
		assert(params.graph.nodeOffset[nodeIndex] + params.graph.NodeLength(nodeIndex) > offset);
		size_t offsetInNode = offset - params.graph.nodeOffset[nodeIndex];
		assert(offsetInNode < params.graph.NodeLength(nodeIndex));
		result.scores.addNodeToMap(nodeIndex);
		result.minScoreNode = nodeIndex;
		result.minScoreNodeOffset = offsetInNode;
		auto& node = result.scores.node(nodeIndex);
		node.startSlice = {0, 0, (int)offsetInNode};
		node.endSlice = {0, 0, (int)params.graph.NodeLength(nodeIndex) - 1 - (int)offsetInNode};
		node.minScore = 0;
		node.exists = true;
		for (size_t i = 1; i <= offsetInNode; i++)
		{
			size_t chunkIndex = i / (sizeof(Word) * 8);
			size_t chunkOffset = i % (sizeof(Word) * 8);
			node.HN[chunkIndex] |= ((Word)1) << chunkOffset;
		}
		for (size_t i = offsetInNode+1; i < params.graph.NodeLength(nodeIndex); i++)
		{
			size_t chunkIndex = i / (sizeof(Word) * 8);
			size_t chunkOffset = i % (sizeof(Word) * 8);
			node.HP[chunkIndex] |= ((Word)1) << chunkOffset;
		}
		return result;
	}

	DPSlice getInitialSliceOneNodeGroup(const std::vector<LengthType>& nodeIndices) const
	{
		DPSlice result;
		result.j = -WordConfiguration<Word>::WordSize;
		result.bandwidth = 1;
		result.minScore = 0;
		result.scores.addEmptyNodeMap(nodeIndices.size());
		for (auto nodeIndex : nodeIndices)
		{
			result.scores.addNodeToMap(nodeIndex);
			result.minScoreNode = nodeIndex;
			result.minScoreNodeOffset = params.graph.NodeLength(nodeIndex) - 1;
			auto& node = result.scores.node(nodeIndex);
			node.startSlice = {0, 0, 0};
			node.endSlice = {0, 0, 0};
			node.minScore = 0;
			node.exists = true;
		}
		return result;
	}

};

#endif
