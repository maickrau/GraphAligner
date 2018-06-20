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
	using SeedHit = typename Common::SeedHit;
	using WordSlice = typename BV::WordSlice;
	using EqVector = typename BV::EqVector;
	using EdgeWithPriority = typename Common::EdgeWithPriority;
	const Params& params;
	class DPSlice
	{
	public:
		DPSlice() :
		minScore(std::numeric_limits<ScoreType>::min()),
		minScoreIndex(-1),
		scores(),
		nodes(),
		correctness(),
		j(std::numeric_limits<LengthType>::max()),
		bandwidth(0),
		cellsProcessed(0),
		numCells(0)
		{}
		DPSlice(std::vector<typename NodeSlice<LengthType, ScoreType, Word>::MapItem>* vectorMap) :
		minScore(std::numeric_limits<ScoreType>::min()),
		minScoreIndex(-1),
		scores(vectorMap),
		nodes(),
		correctness(),
		j(std::numeric_limits<LengthType>::max()),
		bandwidth(0),
		cellsProcessed(0),
		numCells(0)
		{}
		ScoreType minScore;
		LengthType minScoreIndex;
		NodeSlice<LengthType, ScoreType, Word> scores;
		std::vector<size_t> nodes;
		AlignmentCorrectnessEstimationState correctness;
		LengthType j;
		int bandwidth;
		size_t cellsProcessed;
		size_t numCells;
	};
	class DPTable
	{
	public:
		DPTable() :
		slices(),
		bandwidthPerSlice(),
		correctness()
		{}
		std::vector<DPSlice> slices;
		std::vector<ScoreType> bandwidthPerSlice;
		std::vector<AlignmentCorrectnessEstimationState> correctness;
	};
	class TwoDirectionalSplitAlignment
	{
	public:
		size_t EstimatedCorrectlyAligned() const
		{
			return (forward.bandwidthPerSlice.size() + backward.bandwidthPerSlice.size()) * WordConfiguration<Word>::WordSize;
		}
		size_t sequenceSplitIndex;
		DPTable forward;
		DPTable backward;
	};
public:

	GraphAlignerBitvectorBanded(const Params& params) :
	params(params)
	{
	}

	OnewayTrace getTraceFromSeed(const std::string& sequence, int bigraphNodeId, AlignerGraphsizedState& reusableState) const
	{
		std::vector<size_t> nodes;
		nodes = params.graph.nodeLookup.at(bigraphNodeId);
		assert(sequence.size() >= params.graph.DBGOverlap);
		size_t numSlices = (sequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		auto initialBandwidth = getInitialSliceOneNodeGroup(nodes);
		auto slice = getSqrtSlices(sequence, initialBandwidth, numSlices, reusableState);
		removeWronglyAlignedEnd(slice);
		if (slice.slices.size() == 0)
		{
			return OnewayTrace::TraceFailed();
		}
		assert(slice.slices.back().minScore <= sequence.size() + WordConfiguration<Word>::WordSize * 2);

		OnewayTrace result;

		result = getTraceFromTable(sequence, slice, reusableState);
		size_t backtraceableSize = sequence.size() - params.graph.DBGOverlap;
		while (result.trace.size() > 0 && result.trace.back().second >= backtraceableSize)
		{
			result.trace.pop_back();
		}

		return result;
	}

private:

	OnewayTrace getTraceFromTable(const std::string& sequence, const DPTable& slice, AlignerGraphsizedState& reusableState) const
	{
		assert(slice.slices.size() > 0);
		assert(slice.slices.back().minScoreIndex < params.graph.NodeSequencesSize());
		OnewayTrace result;
		result.score = slice.slices.back().minScore;
		result.trace.emplace_back(slice.slices.back().minScoreIndex, std::min(slice.slices.back().j + WordConfiguration<Word>::WordSize - 1, sequence.size()-1));
		LengthType currentNode = -1;
		size_t currentSlice = slice.slices.size();
		std::vector<WordSlice> nodeSlices;
		while (result.trace.back().second != 0)
		{
			size_t newSlice = result.trace.back().second / WordConfiguration<Word>::WordSize + 1;
			assert(newSlice < slice.slices.size());
			assert(result.trace.back().second >= slice.slices[newSlice].j);
			assert(result.trace.back().second < slice.slices[newSlice].j + WordConfiguration<Word>::WordSize);
			LengthType newNode = params.graph.IndexToNode(result.trace.back().first);
			if (newSlice != currentSlice || newNode != currentNode)
			{
				currentSlice = newSlice;
				currentNode = newNode;
				assert(slice.slices[currentSlice].scores.hasNode(currentNode));
				assert(currentSlice > 0);
				typename NodeSlice<LengthType, ScoreType, Word>::NodeSliceMapItem previous;
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
			}
			assert(result.trace.back().first >= params.graph.NodeStart(currentNode));
			assert(result.trace.back().first < params.graph.NodeEnd(currentNode));
			if (result.trace.back().second % WordConfiguration<Word>::WordSize == 0 && result.trace.back().first == params.graph.NodeStart(currentNode))
			{
				result.trace.emplace_back(pickBacktraceCorner(slice.slices[currentSlice].scores, slice.slices[currentSlice-1].scores, currentNode, slice.slices[currentSlice].j, sequence));
				continue;
			}
			if (result.trace.back().second % WordConfiguration<Word>::WordSize == 0)
			{
				assert(currentSlice > 0);
				assert(result.trace.back().first != params.graph.NodeStart(currentNode));
				if (!slice.slices[currentSlice-1].scores.hasNode(currentNode))
				{
					result.trace.emplace_back(params.graph.NodeStart(currentNode), result.trace.back().second);
					continue;
				}
				auto crossing = pickBacktraceVerticalCrossing(slice.slices[currentSlice].scores, slice.slices[currentSlice-1].scores, nodeSlices, slice.slices[currentSlice].j, currentNode, result.trace.back(), sequence);
				if (crossing.first != result.trace.back()) result.trace.push_back(crossing.first);
				assert(crossing.second != result.trace.back());
				result.trace.push_back(crossing.second);
				continue;
			}
			if (result.trace.back().first == params.graph.NodeStart(currentNode))
			{
				assert(result.trace.back().second % WordConfiguration<Word>::WordSize != 0);
				auto crossing = pickBacktraceHorizontalCrossing(slice.slices[currentSlice].scores, slice.slices[currentSlice-1].scores, slice.slices[currentSlice].j, currentNode, result.trace.back(), sequence);
				if (crossing.first != result.trace.back()) result.trace.push_back(crossing.first);
				assert(crossing.second != result.trace.back());
				result.trace.push_back(crossing.second);
				continue;
			}
			assert(result.trace.back().first != params.graph.NodeStart(currentNode));
			assert(result.trace.back().second % WordConfiguration<Word>::WordSize != 0);
			result.trace.push_back(pickBacktraceInside(params.graph.NodeStart(currentNode), slice.slices[currentSlice].j, nodeSlices, result.trace.back(), sequence));
		}
		std::reverse(result.trace.begin(), result.trace.end());
		return result;
	}

	MatrixPosition pickBacktraceInside(LengthType horizontalOffset, LengthType verticalOffset, const std::vector<WordSlice>& nodeSlices, MatrixPosition pos, const std::string& sequence) const
	{
		assert(verticalOffset <= pos.second);
		assert(verticalOffset + WordConfiguration<Word>::WordSize > pos.second);
		assert((verticalOffset % WordConfiguration<Word>::WordSize) == 0);
		assert(horizontalOffset == params.graph.NodeStart(params.graph.IndexToNode(pos.first)));
		pos.first -= horizontalOffset;
		pos.second -= verticalOffset;
		assert(pos.second >= 0);
		assert(pos.second < WordConfiguration<Word>::WordSize);
		assert(pos.first >= 0);
		assert(pos.first < nodeSlices.size());
		while (pos.first > 0 && pos.second > 0)
		{
			ScoreType scoreHere = nodeSlices[pos.first].getValue(pos.second);
			ScoreType verticalScore = nodeSlices[pos.first].getValue(pos.second-1);
			ScoreType horizontalScore = nodeSlices[pos.first-1].getValue(pos.second);
			ScoreType diagonalScore = nodeSlices[pos.first-1].getValue(pos.second-1);
			bool eq = Common::characterMatch(sequence[pos.second + verticalOffset], params.graph.NodeSequences(pos.first + horizontalOffset));
			assert(verticalScore >= scoreHere-1);
			assert(horizontalScore >= scoreHere-1);
			assert(diagonalScore >= scoreHere - (eq?0:1));
			if (diagonalScore == scoreHere - (eq?0:1))
			{
				pos.first--;
				pos.second--;
				continue;
			}
			if (verticalScore == scoreHere - 1)
			{
				pos.second--;
				continue;
			}
			assert(horizontalScore == scoreHere - 1);
			pos.first--;
			continue;
		}
		return std::make_pair(pos.first + horizontalOffset, pos.second + verticalOffset);
	}

	std::pair<MatrixPosition, MatrixPosition> pickBacktraceHorizontalCrossing(const NodeSlice<LengthType, ScoreType, Word>& current, const NodeSlice<LengthType, ScoreType, Word>& previous, size_t j, LengthType node, MatrixPosition pos, const std::string& sequence) const
	{
		assert(current.hasNode(node));
		auto startSlice = current.node(node).startSlice;
		while (pos.second % WordConfiguration<Word>::WordSize != 0 && (startSlice.VP & ((Word)1 << (pos.second % WordConfiguration<Word>::WordSize))))
		{
			pos.second--;
		}
		size_t offset = pos.second % WordConfiguration<Word>::WordSize;
		if (offset == 0)
		{
			return std::make_pair(pos, pickBacktraceCorner(current, previous, node, j, sequence));
		}
		bool eq = Common::characterMatch(sequence[pos.second], params.graph.NodeSequences(pos.first));
		ScoreType scoreHere = startSlice.getValue(offset);
		for (auto neighbor : params.graph.inNeighbors[node])
		{
			if (current.hasNode(neighbor))
			{
				auto neighborSlice = current.node(neighbor).endSlice;
				assert(neighborSlice.getValue(offset) >= scoreHere-1);
				assert(neighborSlice.getValue(offset-1) >= scoreHere - (eq?0:1));
				if (neighborSlice.getValue(offset) == scoreHere-1)
				{
					return std::make_pair(pos, std::make_pair(params.graph.NodeEnd(neighbor)-1, pos.second));
				}
				if (neighborSlice.getValue(offset-1) == scoreHere - (eq?0:1))
				{
					return std::make_pair(pos, std::make_pair(params.graph.NodeEnd(neighbor)-1, pos.second-1));
				}
			}
		}
		assert(false);
		return std::make_pair(std::make_pair(0, 0), std::make_pair(0, 0));
	}

	std::pair<MatrixPosition, MatrixPosition> pickBacktraceVerticalCrossing(const NodeSlice<LengthType, ScoreType, Word>& current, const NodeSlice<LengthType, ScoreType, Word>& previous, const std::vector<WordSlice> nodeScores, size_t j, LengthType node, MatrixPosition pos, const std::string& sequence) const
	{
		LengthType nodeStart = params.graph.NodeStart(node);
		assert(pos.first > nodeStart);
		assert(pos.first - nodeStart < nodeScores.size());
		while (pos.first > nodeStart && nodeScores[pos.first-nodeStart-1].getValue(0) == nodeScores[pos.first-nodeStart].getValue(0) - 1)
		{
			pos.first--;
		}
		if (pos.first == nodeStart)
		{
			return std::make_pair(pos, pickBacktraceCorner(current, previous, node, j, sequence));
		}
		assert(previous.hasNode(node));
		bool eq = Common::characterMatch(sequence[pos.second], params.graph.NodeSequences(pos.first));
		auto previousNode = previous.node(node);
		ScoreType scoreHere = nodeScores[pos.first-nodeStart].getValue(0);
		ScoreType scoreDiagonal = previousNode.startSlice.scoreEnd;
		for (size_t i = 1; i <= pos.first - nodeStart - 1; i++)
		{
			scoreDiagonal += (previousNode.HP[i / WordConfiguration<Word>::WordSize] >> (i % WordConfiguration<Word>::WordSize)) & 1;
			scoreDiagonal -= (previousNode.HN[i / WordConfiguration<Word>::WordSize] >> (i % WordConfiguration<Word>::WordSize)) & 1;
		}
		ScoreType scoreUp = scoreDiagonal;
		scoreUp += (previousNode.HP[(pos.first - nodeStart) / WordConfiguration<Word>::WordSize] >> ((pos.first - nodeStart) % WordConfiguration<Word>::WordSize)) & 1;
		scoreUp -= (previousNode.HN[(pos.first - nodeStart) / WordConfiguration<Word>::WordSize] >> ((pos.first - nodeStart) % WordConfiguration<Word>::WordSize)) & 1;
		assert(scoreUp >= scoreHere - 1);
		assert(scoreDiagonal >= scoreHere - (eq?0:1));
		if (scoreDiagonal == scoreHere - (eq?0:1)) return std::make_pair(pos, std::make_pair(pos.first-1, pos.second-1));
		assert(scoreUp == scoreHere - 1);
		return std::make_pair(pos, std::make_pair(pos.first, pos.second-1));
	}

	MatrixPosition pickBacktraceCorner(const NodeSlice<LengthType, ScoreType, Word>& current, const NodeSlice<LengthType, ScoreType, Word>& previous, LengthType node, size_t j, const std::string& sequence) const
	{
		ScoreType scoreHere = current.node(node).startSlice.getValue(0);
		bool eq = Common::characterMatch(sequence[j], params.graph.NodeSequences(params.graph.NodeStart(node)));
		if (previous.hasNode(node))
		{
			assert(previous.node(node).startSlice.scoreEnd >= scoreHere-1);
			if (previous.node(node).startSlice.scoreEnd == scoreHere-1)
			{
				return std::make_pair(params.graph.NodeStart(node), j-1);
			}
		}
		for (auto neighbor : params.graph.inNeighbors[node])
		{
			if (current.hasNode(neighbor))
			{
				assert(current.node(neighbor).endSlice.getValue(0) >= scoreHere-1);
				if (current.node(neighbor).endSlice.getValue(0) == scoreHere-1)
				{
					return std::make_pair(params.graph.NodeEnd(neighbor)-1, j);
				}
			}
			if (previous.hasNode(neighbor))
			{
				assert(previous.node(neighbor).endSlice.scoreEnd >= scoreHere-(eq?0:1));
				if (previous.node(neighbor).endSlice.scoreEnd == scoreHere-(eq?0:1))
				{
					return std::make_pair(params.graph.NodeEnd(neighbor)-1, j-1);
				}
			}
		}
		assert(false);
		return std::make_pair(0, 0);
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
		LengthType minScoreIndex;
		size_t cellsProcessed;
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

	std::vector<WordSlice> recalcNodeWordslice(LengthType node, const typename NodeSlice<LengthType, ScoreType, Word>::NodeSliceMapItem& slice, const typename NodeSlice<LengthType, ScoreType, Word>::NodeSliceMapItem& previousSlice, LengthType j, const std::string& sequence) const
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
		char graphChar;
		size_t pos;
		WordSlice newWs;
		size_t nodeLength = params.graph.NodeLength(node);
		size_t nodeStart = params.graph.NodeStart(node);
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
				graphChar = params.graph.NodeSequences(nodeStart+pos);
				Eq = EqV.getEq(graphChar);
				Eq &= forceEq;
				if ((HN & 1) && (scoreBefore == scoreComparison - 1))
				{
					assert(previousSlice.exists);
					assert(ws.getScoreBeforeStart() == scoreBefore);
					std::tie(newWs, hinP, hinN) = BV::getNextSlice(Eq, ws, 0, 0);
// #ifdef EXTRACORRECTNESSASSERTIONS
					assertSliceCorrectness(ws, newWs, Eq, 0);
// #endif
				}
				else if (scoreBefore < scoreComparison)
				{
					assert(previousSlice.exists);
					assert(ws.getScoreBeforeStart() == scoreBefore);
					std::tie(newWs, hinP, hinN) = BV::getNextSlice(Eq, ws, 1, 0);
					newWs.VP &= ~(Word)1;
					newWs.VN |= 1;
// #ifdef EXTRACORRECTNESSASSERTIONS
					assertSliceCorrectness(ws, newWs, Eq, 1);
// #endif
				}
				else
				{
					assert(!previousSlice.exists || ws.getScoreBeforeStart() == scoreComparison);
					std::tie(newWs, hinP, hinN) = BV::getNextSlice(Eq, ws, HP & 1, HN & 1);
					newWs.VP &= forceVP;
					newWs.VN |= forceVN;
// #ifdef EXTRACORRECTNESSASSERTIONS
					assertSliceCorrectness(ws, newWs, Eq, (HP & 1) - (HN & 1));
// #endif
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
		return result;
	}

	NodeCalculationResult calculateNode(size_t i, typename NodeSlice<LengthType, ScoreType, Word>::NodeSliceMapItem& slice, const EqVector& EqV, typename NodeSlice<LengthType, ScoreType, Word>::NodeSliceMapItem previousSlice, WordSlice ws, bool skipFirst) const
	{
		NodeCalculationResult result;
		result.minScore = std::numeric_limits<ScoreType>::max();
		result.minScoreIndex = -1;
		result.cellsProcessed = 0;
		auto nodeStart = params.graph.NodeStart(i);
		auto nodeLength = params.graph.NodeLength(i);

		char graphChar = params.graph.NodeSequences(nodeStart);
		Word Eq = EqV.getEq(graphChar);

		Word hinP;
		Word hinN;
		if (skipFirst)
		{
			hinP = 0;
			hinN = 0;
		}
		else
		{
			if (previousSlice.exists)
			{
				ScoreType incomingScoreBeforeStart = ws.getScoreBeforeStart();
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
			std::tie(newWs, hinP, hinN) = BV::getNextSlice(Eq, ws, hinP, hinN);
			if (!previousSlice.exists || newWs.getScoreBeforeStart() < previousSlice.startSlice.scoreEnd)
			{
				newWs.VP &= WordConfiguration<Word>::AllOnes ^ 1;
				newWs.VN |= 1;
			}
			assert(newWs.getScoreBeforeStart() >= debugLastRowMinScore);
			ws = newWs;
		}
		result.cellsProcessed++;
		result.minScore = ws.scoreEnd;
		result.minScoreIndex = nodeStart;

		if (slice.exists)
		{
			//todo fix
			if (true)
			// if (params.graph.inNeighbors[i].size() > 1)
			{
				WordSlice test = ws.mergeWith(slice.startSlice);
				if (test.scoreEnd == slice.startSlice.scoreEnd && test.VP == slice.startSlice.VP && test.VP == slice.startSlice.VN)
				{
					return result;
				}
				ws = test;
			}
			else
			{
#ifndef NDEBUG
				WordSlice test = ws.mergeWith(slice.startSlice);
#endif
				if (slice.startSlice.scoreEnd < ws.scoreEnd)
				{
					assert(test.VP == slice.startSlice.VP);
					assert(test.VN == slice.startSlice.VN);
					assert(test.scoreEnd == slice.startSlice.scoreEnd);
					return result;
				}
				else
				{
					assert(test.VP == ws.VP);
					assert(test.VN == ws.VN);
					assert(test.scoreEnd == ws.scoreEnd);
				}
			}
		}
		else if (!skipFirst && previousSlice.exists)
		{
			if (ws.getValue(0) > previousSlice.startSlice.scoreEnd + 1)
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
		size_t chunk = 0;
		size_t offset = 1;
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
							previousSlice.HP[chunk] |= mask;
							previousSlice.HN[chunk] &= ~mask;
							forceUntil = fixchunk * WordConfiguration<Word>::WordSize + fixoffset;
						}
						if (scoreBefore == newScoreComparison)
						{
							previousSlice.HP[chunk] &= ~mask;
							previousSlice.HN[chunk] &= ~mask;
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
		WordSlice newWs;
		Word forceEq = WordConfiguration<Word>::AllOnes;
		if (!previousSlice.exists) Eq ^= 1;
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
				graphChar = params.graph.NodeSequences(nodeStart+pos);
				Eq = EqV.getEq(graphChar);
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
				HP >>= 1;
				HN >>= 1;
				slice.HP[chunk] |= hinP << offset;
				slice.HN[chunk] |= hinN << offset;
				result.cellsProcessed++;
				if (ws.scoreEnd < result.minScore)
				{
					result.minScore = ws.scoreEnd;
					result.minScoreIndex = nodeStart + pos;
				}
			}
			offset = 0;
		}
		slice.endSlice = ws;
#ifndef NDEBUG
		if (previousSlice.exists && forceUntil < nodeLength - 1) assert(slice.endSlice.getScoreBeforeStart() == previousSlice.endSlice.scoreEnd);
#endif
		return result;
	}

	static bool cellExists(const Params& params, const DPSlice& slice, LengthType row, LengthType w)
	{
		auto node = params.graph.IndexToNode(w);
		auto offset = w - params.graph.NodeStart(node);
		if (!slice.scores.hasNode(node)) return false;
		auto wordslice = slice.scores.node(node)[offset];
		return wordslice.sliceExists;
	}

	bool cellExists(const DPSlice& slice, LengthType row, LengthType w) const
	{
		return cellExists(params, slice, row, w);
	}

	static ScoreType getValueIfExists(const Params& params, const DPSlice& slice, int row, LengthType cell, ScoreType defaultValue)
	{
		auto nodeIndex = params.graph.IndexToNode(cell);
		if (!slice.scores.hasNode(nodeIndex)) return defaultValue;
		auto wordslice = slice.scores.node(nodeIndex)[cell - params.graph.NodeStart(nodeIndex)];
		if (cellExists(wordslice, row, slice.minScore + slice.bandwidth)) return wordslice.getValue(row);
		return defaultValue;
	}

	ScoreType getValueIfExists(const DPSlice& slice, int row, LengthType cell, ScoreType defaultValue) const
	{
		return getValueIfExists(params, slice, row, cell, defaultValue);
	}

	static bool cellExists(WordSlice slice, int row, int maxScore)
	{
		if (!slice.sliceExists) return false;
		auto score = slice.getValue(row);
		return score <= maxScore;
	}

	static ScoreType getValueOrMax(const Params& params, const DPTable& band, LengthType j, LengthType w, ScoreType max)
	{
		//todo fix
		// auto node = params.graph.IndexToNode(w);
		// auto slice = (j - band.startj) / WordConfiguration<Word>::WordSize / band.samplingFrequency;
		// assert(band.slices[slice].j == (j / WordConfiguration<Word>::WordSize) * WordConfiguration<Word>::WordSize);
		// if (!band.slices[slice].scores.hasNode(node)) return max;
		// auto word = band.slices[slice].scores.node(node)[w - params.graph.NodeStart(node)];
		// auto off = j % WordConfiguration<Word>::WordSize;
		// return word.getValue(off);
	}

	static ScoreType getValueOrMax(const Params& params, const DPSlice& slice, LengthType j, LengthType w, ScoreType max)
	{
		//todo fix
		// assert(j >= 0);
		// assert(j < WordConfiguration<Word>::WordSize);
		// auto node = params.graph.IndexToNode(w);
		// if (!slice.scores.hasNode(node)) return max;
		// auto word = slice.scores.node(node)[w - params.graph.NodeStart(node)];
		// auto off = j % WordConfiguration<Word>::WordSize;
		// return word.getValue(off);
	}

	static ScoreType getValue(const Params& params, const DPSlice& slice, LengthType j, LengthType w)
	{
		//todo fix
		// assert(j >= 0);
		// assert(j < WordConfiguration<Word>::WordSize);
		// auto node = params.graph.IndexToNode(w);
		// auto word = slice.scores.node(node)[w - params.graph.NodeStart(node)];
		// auto off = j % WordConfiguration<Word>::WordSize;
		// return word.getValue(off);
	}

	static ScoreType getValue(const Params& params, const DPTable& band, LengthType j, LengthType w)
	{
		// auto node = params.graph.IndexToNode(w);
		// auto slice = (j - band.startj) / WordConfiguration<Word>::WordSize / band.samplingFrequency;
		// assert(band.slices[slice].j == (j / WordConfiguration<Word>::WordSize) * WordConfiguration<Word>::WordSize);
		// auto word = band.slices[slice].scores.node(node)[w - params.graph.NodeStart(node)];
		// auto off = j % WordConfiguration<Word>::WordSize;
		// return word.getValue(off);
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
#endif

#ifdef EXTRACORRECTNESSASSERTIONS
	void checkNodeBoundaryCorrectness(const NodeSlice<LengthType, ScoreType, Word>& currentSlice, const NodeSlice<LengthType, ScoreType, Word>& previousSlice, const std::string& sequence, size_t j, ScoreType maxScore, ScoreType previousMaxScore) const
	{
		for (auto pair : currentSlice)
		{
			auto node = pair.first;
			if (previousSlice.hasNode(node) && previousSlice.node(node).exists)
			{
				assert(pair.second.startSlice.getScoreBeforeStart() <= previousSlice.node(node).startSlice.scoreEnd);
			}
			bool eq = Common::characterMatch(sequence[j], params.graph.NodeSequences(params.graph.NodeStart(node)));
			if (j == 0 && previousSlice.hasNode(node))
			{
				assert(pair.second.startSlice.getValue(0) == (eq ? 0 : 1));
			}
			else
			{
				ScoreType foundMinScore = std::numeric_limits<ScoreType>::max();
				if (previousSlice.hasNode(node) && previousSlice.node(node).exists)
				{
					foundMinScore = std::min(foundMinScore, previousSlice.node(node).startSlice.scoreEnd+1);
				}
				for (auto neighbor : params.graph.inNeighbors[node])
				{
					if (currentSlice.hasNode(neighbor) && currentSlice.node(neighbor).exists)
					{
						foundMinScore = std::min(foundMinScore, currentSlice.node(neighbor).endSlice.getValue(0)+1);
					}
					if (previousSlice.hasNode(neighbor) && previousSlice.node(neighbor).exists)
					{
						foundMinScore = std::min(foundMinScore, previousSlice.node(neighbor).endSlice.scoreEnd + (eq ? 0 : 1));
					}
				}
				if (pair.second.startSlice.getValue(0) <= maxScore || foundMinScore <= maxScore)
				{
					assert(foundMinScore != std::numeric_limits<ScoreType>::max());
					assert(pair.second.startSlice.getValue(0) == foundMinScore);
				}
			}
			for (size_t i = 1; i < WordConfiguration<Word>::WordSize; i++)
			{
				eq = Common::characterMatch(sequence[j+i], params.graph.NodeSequences(params.graph.NodeStart(node)));
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

	NodeCalculationResult calculateSlice(const std::string& sequence, size_t j, NodeSlice<LengthType, ScoreType, Word>& currentSlice, const NodeSlice<LengthType, ScoreType, Word>& previousSlice, const std::vector<LengthType>& previousNodes, std::vector<bool>& currentBand, const std::vector<bool>& previousBand, ArrayPriorityQueue<EdgeWithPriority>& calculableQueue, ScoreType previousQuitScore, int bandwidth, ScoreType previousMinScore) const
	{
		ScoreType currentMinimumScore = std::numeric_limits<ScoreType>::max() - bandwidth - 1;
		LengthType currentMinimumIndex;
		size_t cellsProcessed = 0;

		EqVector EqV = BV::getEqVector(sequence, j);

		assert(previousNodes.size() > 0);
		if (j == 0)
		{
			for (auto node : previousSlice)
			{
				assert(node.second.minScore <= previousQuitScore);
				WordSlice startSlice = getSourceSliceFromScore(node.second.startSlice.scoreEnd);
				calculableQueue.insert(node.second.minScore - previousMinScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, false });
			}
		}
		else
		{
			for (auto node : previousSlice)
			{
				assert(node.second.exists);
				// if (node.second.minScore <= previousQuitScore)
				// {
					WordSlice startSlice = getSourceSliceFromScore(node.second.startSlice.scoreEnd);
					calculableQueue.insert(node.second.minScore - previousMinScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
				// }
			}
		}
		assert(calculableQueue.size() > 0);
		
		ScoreType currentMinScoreAtEndRow = currentMinimumScore;
		while (calculableQueue.size() > 0)
		{
			auto pair = calculableQueue.top();
			if (pair.priority > currentMinScoreAtEndRow + bandwidth) break;
			auto i = pair.target;
			WordSlice incoming = pair.incoming;
			bool skipFirst = pair.skipFirst;
			if (!currentBand[i])
			{
				assert(!currentSlice.hasNode(i));
				currentSlice.addNode(i);
				currentBand[i] = true;
			}
			assert(currentBand[i]);
			calculableQueue.pop();
			auto& thisNode = currentSlice.node(i);
			auto oldEnd = thisNode.endSlice;
			typename NodeSlice<LengthType, ScoreType, Word>::NodeSliceMapItem previousThisNode;

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
			auto nodeCalc = calculateNode(i, thisNode, EqV, previousThisNode, incoming, skipFirst);
			// assert(nodeCalc.minScore <= previousQuitScore + WordConfiguration<Word>::WordSize);
			currentMinScoreAtEndRow = std::min(currentMinScoreAtEndRow, nodeCalc.minScore);
			currentSlice.setMinScoreIfSmaller(i, nodeCalc.minScore);
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
						calculableQueue.insert(newEndMinScore - previousMinScore, EdgeWithPriority { neighbor, newEndMinScore - previousMinScore, newEnd, false });
					}
				}
			}
			if (nodeCalc.minScore < currentMinimumScore)
			{
				currentMinimumScore = nodeCalc.minScore;
				currentMinimumIndex = nodeCalc.minScoreIndex;
			}
			assert(currentMinimumScore == currentMinScoreAtEndRow);
			cellsProcessed += nodeCalc.cellsProcessed;
			if (cellsProcessed > params.maxCellsPerSlice) break;
		}

#ifdef EXTRACORRECTNESSASSERTIONS
		checkNodeBoundaryCorrectness(currentSlice, previousSlice, sequence, j, currentMinScoreAtEndRow + bandwidth, previousQuitScore);
#endif

		NodeCalculationResult result;
		result.minScore = currentMinimumScore;
		result.minScoreIndex = currentMinimumIndex;
		result.cellsProcessed = cellsProcessed;

		if (j + WordConfiguration<Word>::WordSize > sequence.size())
		{
			flattenLastSliceEnd(currentSlice, result, j, sequence.size());
		}

		finalizeSlice(currentSlice, currentBand, currentMinScoreAtEndRow + bandwidth);

		calculableQueue.clear();

		return result;
	}

	void finalizeSlice(NodeSlice<LengthType, ScoreType, Word>& slice, std::vector<bool>& currentBand, ScoreType maxScore) const
	{
		for (auto node : slice)
		{
			if (node.second.minScore > maxScore && node.second.endSlice.getMinScore() > maxScore)
			{
				currentBand[node.first] = false;
				slice.node(node.first).exists = false;
			}
		}
	}

	void flattenLastSliceEnd(NodeSlice<LengthType, ScoreType, Word>& slice, NodeCalculationResult& sliceCalc, LengthType j, size_t sequenceSize) const
	{
		//todo fix
		// assert(j < sequenceSize);
		// assert(sequenceSize - j < WordConfiguration<Word>::WordSize);
		// sliceCalc.minScore = std::numeric_limits<ScoreType>::max();
		// sliceCalc.minScoreIndex = -1;
		// for (auto node : slice)
		// {
		// 	for (size_t i = 0; i < node.second.size(); i++)
		// 	{
		// 		auto wordSliceResult = BV::flattenWordSlice(node.second[i], sequenceSize - j);
		// 		node.second[i] = wordSliceResult;
		// 		if (wordSliceResult.scoreEnd < sliceCalc.minScore)
		// 		{
		// 			sliceCalc.minScore = wordSliceResult.scoreEnd;
		// 			sliceCalc.minScoreIndex = params.graph.NodeStart(node.first)+i;
		// 		}
		// 	}
		// }
	}

	void fillDPSlice(const std::string& sequence, DPSlice& slice, const DPSlice& previousSlice, const std::vector<bool>& previousBand, std::vector<bool>& currentBand, ArrayPriorityQueue<EdgeWithPriority>& calculableQueue, int bandwidth) const
	{
		auto sliceResult = calculateSlice(sequence, slice.j, slice.scores, previousSlice.scores, previousSlice.nodes, currentBand, previousBand, calculableQueue, previousSlice.minScore + previousSlice.bandwidth, bandwidth, previousSlice.minScore);
		slice.bandwidth = bandwidth;
		slice.cellsProcessed = sliceResult.cellsProcessed;
		slice.minScoreIndex = sliceResult.minScoreIndex;
		slice.minScore = sliceResult.minScore;
		assert(slice.minScore >= previousSlice.minScore);
		for (auto node : slice.scores)
		{
			if (currentBand[node.first])
			{
				slice.nodes.push_back(node.first);
				slice.numCells += params.graph.NodeLength(node.first);
			}
		}
		slice.correctness = slice.correctness.NextState(slice.minScore - previousSlice.minScore, WordConfiguration<Word>::WordSize);
	}

	DPSlice pickMethodAndExtendFill(const std::string& sequence, const DPSlice& previous, const std::vector<bool>& previousBand, std::vector<bool>& currentBand, std::vector<typename NodeSlice<LengthType, ScoreType, Word>::MapItem>& nodesliceMap, ArrayPriorityQueue<EdgeWithPriority>& calculableQueue, int bandwidth) const
	{
		DPSlice bandTest { &nodesliceMap };
		bandTest.j = previous.j + WordConfiguration<Word>::WordSize;
		bandTest.correctness = previous.correctness;

		fillDPSlice(sequence, bandTest, previous, previousBand, currentBand, calculableQueue, bandwidth);

		bandTest.scores.convertVectorArrayToMap();

		return bandTest;
	}

	void removeWronglyAlignedEnd(DPTable& table) const
	{
		if (table.slices.size() == 0) return;
		bool currentlyCorrect = table.correctness.back().CurrentlyCorrect();
		while (!currentlyCorrect)
		{
			currentlyCorrect = table.correctness.back().FalseFromCorrect();
			table.correctness.pop_back();
			table.bandwidthPerSlice.pop_back();
			if (table.correctness.size() == 0) break;
		}
		if (table.correctness.size() == 0)
		{
			table.slices.clear();
		}
		while (table.slices.size() > 1 && table.slices.back().j >= table.correctness.size() * WordConfiguration<Word>::WordSize) table.slices.pop_back();
	}

#ifndef NDEBUG
	void printPathExtensions(LengthType startpos, std::string prefix) const
	{
		auto node = params.graph.IndexToNode(startpos);
		auto end = params.graph.NodeEnd(node);
		for (size_t i = startpos; i < end && prefix.size() < 64; i++)
		{
			prefix += params.graph.NodeSequences(i);
		}
		if (prefix.size() == 64)
		{
			std::cerr << prefix << " " << node << std::endl;
		}
		else
		{
			if (params.graph.outNeighbors[node].size() == 0)
			{
				std::cerr << prefix << " TIP! " << node << std::endl;
			}
			else
			{
				for (auto neighbor : params.graph.outNeighbors[node])
				{
					printPathExtensions(params.graph.NodeStart(neighbor), prefix);
				}
			}
		}
	}

	void __attribute__ ((noinline)) printPathExtensions(LengthType startpos) const
	{
		printPathExtensions(startpos, "");
		asm ("");
	}
#endif

	DPTable getSqrtSlices(const std::string& sequence, const DPSlice& initialSlice, size_t numSlices, AlignerGraphsizedState& reusableState) const
	{
		assert(initialSlice.j == -WordConfiguration<Word>::WordSize);
		assert(initialSlice.j + numSlices * WordConfiguration<Word>::WordSize <= sequence.size() + WordConfiguration<Word>::WordSize);
		DPTable result;
		size_t realCells = 0;
		size_t cellsProcessed = 0;
		std::vector<size_t> partOfComponent;
		{
			auto initialOrder = initialSlice.nodes;
			for (auto node : initialOrder)
			{
				reusableState.previousBand[node] = true;
			}
		}
#ifndef NDEBUG
		debugLastRowMinScore = 0;
#endif
		DPSlice lastSlice = initialSlice;
		result.slices.push_back(initialSlice);
		result.correctness.push_back(initialSlice.correctness);
		result.bandwidthPerSlice.push_back(1);
		assert(lastSlice.correctness.CurrentlyCorrect());
		DPSlice rampSlice = lastSlice;
		size_t rampRedoIndex = -1;
		size_t rampUntil = 0;
#ifndef NDEBUG
		volatile size_t debugLastProcessedSlice;
#endif
		for (size_t slice = 0; slice < numSlices; slice++)
		{
			int bandwidth = (params.rampBandwidth > params.initialBandwidth && rampUntil >= slice) ? params.rampBandwidth : params.initialBandwidth;
#ifndef NDEBUG
			debugLastProcessedSlice = slice;
			debugLastRowMinScore = lastSlice.minScore;
#endif
			auto timeStart = std::chrono::system_clock::now();
			auto newSlice = pickMethodAndExtendFill(sequence, lastSlice, reusableState.previousBand, reusableState.currentBand, reusableState.nodesliceMap, reusableState.calculableQueue, bandwidth);
			auto timeEnd = std::chrono::system_clock::now();
			auto time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
#ifdef SLICEVERBOSE
			std::cerr << "slice " << slice << " bandwidth " << bandwidth << " minscore " << newSlice.minScore << " diff " << (newSlice.minScore - lastSlice.minScore) << " time " << time << " slices " << newSlice.numCells << " cellsprocessed " << newSlice.cellsProcessed << " overhead " << (100 * (int)(newSlice.cellsProcessed - newSlice.numCells) / (int)(newSlice.numCells)) << "%";
#endif

			if (rampUntil == slice && newSlice.cellsProcessed >= params.BacktraceOverrideCutoff)
			{
				rampUntil++;
			}
			if ((rampUntil == slice-1 || (rampUntil < slice && newSlice.correctness.CurrentlyCorrect() && newSlice.correctness.FalseFromCorrect())) && newSlice.cellsProcessed < params.BacktraceOverrideCutoff)
			{
				rampSlice = lastSlice;
				rampRedoIndex = slice-1;
			}
			assert(newSlice.j == lastSlice.j + WordConfiguration<Word>::WordSize);

			realCells += newSlice.numCells;
			cellsProcessed += newSlice.cellsProcessed;

			if (newSlice.cellsProcessed > params.maxCellsPerSlice)
			{
#ifndef NDEBUG
				debugLastProcessedSlice = slice-1;
#endif
				for (auto node : lastSlice.nodes)
				{
					assert(reusableState.previousBand[node]);
					reusableState.previousBand[node] = false;
				}
				for (auto node : newSlice.nodes)
				{
					assert(reusableState.currentBand[node]);
					reusableState.currentBand[node] = false;
				}
				break;
			}

			if (!newSlice.correctness.CorrectFromCorrect())
			{
#ifndef NDEBUG
				debugLastProcessedSlice = slice-1;
#endif
				for (auto node : lastSlice.nodes)
				{
					assert(reusableState.previousBand[node]);
					reusableState.previousBand[node] = false;
				}
				for (auto node : newSlice.nodes)
				{
					assert(reusableState.currentBand[node]);
					reusableState.currentBand[node] = false;
				}
				break;
			}
			if (!newSlice.correctness.CurrentlyCorrect() && rampUntil < slice && params.rampBandwidth > params.initialBandwidth)
			{
				for (auto node : newSlice.nodes)
				{
					assert(reusableState.currentBand[node]);
					reusableState.currentBand[node] = false;
				}
				for (auto node : lastSlice.nodes)
				{
					assert(reusableState.previousBand[node]);
					reusableState.previousBand[node] = false;
				}
				rampUntil = slice;
				std::swap(slice, rampRedoIndex);
				std::swap(lastSlice, rampSlice);
				for (auto node : lastSlice.nodes)
				{
					assert(!reusableState.previousBand[node]);
					reusableState.previousBand[node] = true;
				}
				if (slice == -1)
				{
					result.bandwidthPerSlice.clear();
					result.correctness.clear();
					result.slices.clear();
				}
				while (result.bandwidthPerSlice.size() > slice+2) result.bandwidthPerSlice.pop_back();
				while (result.correctness.size() > slice+2) result.correctness.pop_back();
				while (result.slices.size() > 1 && result.slices.back().j > slice * WordConfiguration<Word>::WordSize) result.slices.pop_back();
				assert(result.slices.size() == result.bandwidthPerSlice.size());
				assert(result.correctness.size() == result.slices.size());
				assert(result.slices.back().j == lastSlice.j);
#ifdef SLICEVERBOSE
				std::cerr << " ramp to " << slice;
				std::cerr << std::endl;
				std::cerr << "bandwidthPerSlice.size() " << result.bandwidthPerSlice.size();
				if (result.slices.size() > 0) std::cerr << " slices.back().j " << result.slices.back().j; else std::cerr << " slices.size() 0";
				std::cerr << std::endl;
#endif
				continue;
			}

#ifdef SLICEVERBOSE
			std::cerr << std::endl;
#endif

			assert(result.bandwidthPerSlice.size() == slice+1);
			result.bandwidthPerSlice.push_back(bandwidth);
			result.correctness.push_back(newSlice.correctness);
			result.slices.push_back(newSlice);
			for (auto node : lastSlice.nodes)
			{
				assert(reusableState.previousBand[node]);
				reusableState.previousBand[node] = false;
			}
			assert(newSlice.minScore != std::numeric_limits<LengthType>::max());
			assert(newSlice.minScore >= lastSlice.minScore);
			if (slice == numSlices - 1)
			{
				for (auto node : newSlice.nodes)
				{
					assert(reusableState.currentBand[node]);
					reusableState.currentBand[node] = false;
				}
			}
			else
			{
				std::swap(reusableState.previousBand, reusableState.currentBand);
			}
			lastSlice = newSlice;
		}

#ifdef EXTRACORRECTNESSASSERTIONS
		assert(reusableState.calculableQueue.size() == 0);
		for (size_t i = 0; i < reusableState.currentBand.size(); i++)
		{
			assert(!reusableState.currentBand[i]);
			assert(!reusableState.previousBand[i]);
			assert(!reusableState.nodesliceMap[i].exists);
		}
#endif

#ifndef NDEBUG
		assert(result.slices.size() == result.bandwidthPerSlice.size());
		assert(result.correctness.size() == result.slices.size());
		if (result.slices.size() > 0)
		{
			volatile size_t lastExisting = 0;
			assert(result.bandwidthPerSlice.size() == debugLastProcessedSlice + 2);
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

	std::vector<DPSlice> getSlicesFromTable(const std::string& sequence, LengthType overrideLastJ, const DPTable& table, size_t startIndex, AlignerGraphsizedState& reusableState) const
	{
		assert(startIndex < table.slices.size());
		size_t startSlice = (table.slices[startIndex].j + WordConfiguration<Word>::WordSize) / WordConfiguration<Word>::WordSize;
		assert(overrideLastJ > startSlice * WordConfiguration<Word>::WordSize);
		size_t endSlice;
		if (startIndex == table.slices.size()-1) endSlice = table.bandwidthPerSlice.size(); else endSlice = (table.slices[startIndex+1].j + WordConfiguration<Word>::WordSize) / WordConfiguration<Word>::WordSize;
		if (endSlice * WordConfiguration<Word>::WordSize >= overrideLastJ) endSlice = (overrideLastJ / WordConfiguration<Word>::WordSize);
		assert(endSlice > startSlice);
		assert(endSlice <= table.bandwidthPerSlice.size());
		assert(startIndex < table.slices.size());
		const auto& initialSlice = table.slices[startIndex];
		std::vector<DPSlice> result;
		size_t realCells = 0;
		size_t cellsProcessed = 0;
		{
			auto initialOrder = initialSlice.nodes;
			for (auto node : initialOrder)
			{
				reusableState.previousBand[node] = true;
			}
		}
#ifndef NDEBUG
		debugLastRowMinScore = 0;
#endif
		DPSlice lastSlice = initialSlice;
		// assert(lastSlice.correctness.CurrentlyCorrect());
		DPSlice rampSlice = lastSlice;
		for (size_t slice = startSlice; slice < endSlice; slice++)
		{
			int bandwidth = table.bandwidthPerSlice[slice];
#ifndef NDEBUG
			debugLastRowMinScore = lastSlice.minScore;
#endif
			auto newSlice = pickMethodAndExtendFill(sequence, lastSlice, reusableState.previousBand, reusableState.currentBand, reusableState.nodesliceMap, reusableState.calculableQueue, bandwidth);
			assert(result.size() == 0 || newSlice.j == result.back().j + WordConfiguration<Word>::WordSize);

			size_t sliceCells = 0;
			for (auto node : newSlice.nodes)
			{
				sliceCells += params.graph.NodeEnd(node) - params.graph.NodeStart(node);
			}
			realCells += sliceCells;
			cellsProcessed += newSlice.cellsProcessed;

			// assert(slice == endSlice-1 || newSlice.correctness.CurrentlyCorrect());
			result.push_back(newSlice);
			for (auto node : lastSlice.nodes)
			{
				assert(reusableState.previousBand[node]);
				reusableState.previousBand[node] = false;
			}
			assert(newSlice.minScore != std::numeric_limits<LengthType>::max());
			assert(newSlice.minScore >= lastSlice.minScore);
#ifndef NDEBUG
			auto debugMinimumNode = params.graph.IndexToNode(newSlice.minScoreIndex);
			assert(newSlice.scores.hasNode(debugMinimumNode));
			auto debugslice = newSlice.scores.node(debugMinimumNode);
			assert(newSlice.minScoreIndex >= params.graph.NodeStart(debugMinimumNode));
			assert(newSlice.minScoreIndex < params.graph.NodeEnd(debugMinimumNode));
			assert(debugslice[newSlice.minScoreIndex - params.graph.NodeStart(debugMinimumNode)].scoreEnd == newSlice.minScore);
#endif
			if (slice == endSlice-1)
			{
				for (auto node : newSlice.nodes)
				{
					assert(reusableState.currentBand[node]);
					reusableState.currentBand[node] = false;
				}
			}
			else
			{
				std::swap(reusableState.previousBand, reusableState.currentBand);
			}
			lastSlice = std::move(newSlice);
		}
#ifndef NDEBUG
		for (size_t i = 1; i < result.size(); i++)
		{
			assert(result[i].minScore >= result[i-1].minScore);
		}
#endif

#ifdef EXTRACORRECTNESSASSERTIONS
		assert(reusableState.calculableQueue.size() == 0);
		for (size_t i = 0; i < reusableState.currentBand.size(); i++)
		{
			assert(!reusableState.currentBand[i]);
			assert(!reusableState.previousBand[i]);
			assert(reusableState.nodesliceMap[i].start == 0);
			assert(reusableState.nodesliceMap[i].end == 0);
			assert(!reusableState.nodesliceMap[i].exists);
		}
#endif

		return result;
	}

	DPSlice getInitialSliceOneNodeGroup(const std::vector<LengthType>& nodeIndices) const
	{
		DPSlice result;
		result.j = -WordConfiguration<Word>::WordSize;
		result.bandwidth = 1;
		result.minScore = 0;
		for (auto nodeIndex : nodeIndices)
		{
			result.scores.addNodeToMap(nodeIndex);
			result.minScoreIndex = params.graph.NodeEnd(nodeIndex) - 1;
			auto& node = result.scores.node(nodeIndex);
			node.startSlice = {0, 0, 0};
			node.endSlice = {0, 0, 0};
			node.minScore = 0;
			node.exists = true;
			result.nodes.push_back(nodeIndex);
		}
		return result;
	}

	OnewayTrace getBacktraceFullStart(std::string sequence, AlignerGraphsizedState& reusableState) const
	{
		int padding = (WordConfiguration<Word>::WordSize - (sequence.size() % WordConfiguration<Word>::WordSize)) % WordConfiguration<Word>::WordSize;
		for (int i = 0; i < padding; i++)
		{
			sequence += 'N';
		}
		DPSlice startSlice;
		for (size_t i = 0; i < params.graph.nodeStart.size(); i++)
		{
			startSlice.scores.addNodeToMap(i);
			startSlice.scores.setMinScore(i, 0);
			startSlice.j = -WordConfiguration<Word>::WordSize;
			startSlice.nodes.push_back(i);
		}
		auto slice = getSqrtSlices(sequence, startSlice, sequence.size() / WordConfiguration<Word>::WordSize, reusableState);
		removeWronglyAlignedEnd(slice);
		if (slice.slices.size() == 0)
		{
			return OnewayTrace::TraceFailed();
		}
		// std::cerr << "score: " << slice.slices.back().minScore << std::endl;

		auto result = getTraceFromTable(sequence, slice, reusableState);
		while (result.trace.back().second >= sequence.size() - padding)
		{
			result.trace.pop_back();
		}
		assert(result.trace[0].second == 0);
		return result;
	}

};

#endif
