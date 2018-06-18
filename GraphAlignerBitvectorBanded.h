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

	static MatrixPosition pickBacktracePredecessor(const Params& params, const std::string& sequence, const DPSlice& slice, const MatrixPosition pos, const DPSlice& previousSlice)
	{
		const ScoreType defaultValue = sequence.size() + slice.bandwidth + 1;
		assert(pos.second >= slice.j);
		assert(pos.second < slice.j + WordConfiguration<Word>::WordSize);
		assert(cellExists(params, slice, pos.second - slice.j, pos.first));
		auto nodeIndex = params.graph.IndexToNode(pos.first);
		assert(slice.scores.hasNode(nodeIndex));
		auto scoreHere = getValue(params, slice, pos.second - slice.j, pos.first);
		if (pos.second == 0 && previousSlice.scores.hasNode(nodeIndex) && (scoreHere == 0 || scoreHere == 1)) return { pos.first, pos.second - 1 };
		if (pos.first == params.graph.NodeStart(nodeIndex))
		{
			for (auto neighbor : params.graph.inNeighbors[nodeIndex])
			{
				LengthType u = params.graph.NodeEnd(neighbor)-1;
				auto horizontalScore = getValueIfExists(params, slice, pos.second - slice.j, u, defaultValue);
				assert(horizontalScore >= scoreHere-1);
				if (horizontalScore == scoreHere-1)
				{
					assert(cellExists(params, slice, pos.second - slice.j, u));
					return { u, pos.second };
				}
				ScoreType diagonalScore;
				if (pos.second == slice.j)
				{
					diagonalScore = getValueIfExists(params, previousSlice, WordConfiguration<Word>::WordSize - 1, u, defaultValue);
				}
				else
				{
					diagonalScore = getValueIfExists(params, slice, pos.second - 1 - slice.j, u, defaultValue);
				}
				if (Common::characterMatch(sequence[pos.second], params.graph.NodeSequences(pos.first)))
				{
					assert(diagonalScore >= scoreHere);
					if (diagonalScore == scoreHere)
					{
						assert(pos.second == slice.j || cellExists(params, slice, pos.second - slice.j - 1, u));
						assert(pos.second != slice.j || cellExists(params, previousSlice, pos.second - previousSlice.j - 1, u));
						return { u, pos.second - 1 };
					}
				}
				else
				{
					assert(diagonalScore >= scoreHere-1);
					if (diagonalScore == scoreHere-1)
					{
						assert(pos.second == slice.j || cellExists(params, slice, pos.second - slice.j - 1, u));
						assert(pos.second != slice.j || cellExists(params, previousSlice, pos.second - previousSlice.j - 1, u));
						return { u, pos.second - 1 };
					}
				}
			}
		}
		else
		{
			auto horizontalScore = getValueIfExists(params, slice, pos.second - slice.j, pos.first-1, defaultValue);
			assert(horizontalScore >= scoreHere-1);
			if (horizontalScore == scoreHere-1)
			{
				assert(cellExists(params, slice, pos.second - slice.j, pos.first - 1));
				return { pos.first - 1, pos.second };
			}
			ScoreType diagonalScore;
			if (pos.second == slice.j)
			{
				diagonalScore = getValueIfExists(params, previousSlice, WordConfiguration<Word>::WordSize - 1, pos.first-1, defaultValue);
			}
			else
			{
				diagonalScore = getValueIfExists(params, slice, pos.second - 1 - slice.j, pos.first-1, defaultValue);
			}
			if (Common::characterMatch(sequence[pos.second], params.graph.NodeSequences(pos.first)))
			{
				assert(diagonalScore >= scoreHere);
				if (diagonalScore == scoreHere)
				{
					assert(pos.second == slice.j || cellExists(params, slice, pos.second - slice.j - 1, pos.first - 1));
					assert(pos.second != slice.j || cellExists(params, previousSlice, pos.second - previousSlice.j - 1, pos.first - 1));
					return { pos.first - 1, pos.second - 1 };
				}
			}
			else
			{
				assert(diagonalScore >= scoreHere-1);
				if (diagonalScore == scoreHere-1)
				{
					assert(pos.second == slice.j || cellExists(params, slice, pos.second - slice.j - 1, pos.first - 1));
					assert(pos.second != slice.j || cellExists(params, previousSlice, pos.second - previousSlice.j - 1, pos.first - 1));
					return { pos.first - 1, pos.second - 1 };
				}
			}
		}
		ScoreType scoreUp;
		if (pos.second == slice.j)
		{
			assert(previousSlice.j + WordConfiguration<Word>::WordSize == slice.j);
			scoreUp = getValueIfExists(params, previousSlice, WordConfiguration<Word>::WordSize - 1, pos.first, defaultValue);
		}
		else
		{
			scoreUp = getValueIfExists(params, slice, pos.second - 1 - slice.j, pos.first, defaultValue);
		}
		assert(scoreUp >= scoreHere-1);
		if (scoreUp == scoreHere - 1)
		{
			assert(pos.second == slice.j || cellExists(params, slice, pos.second - slice.j - 1, pos.first));
			assert(pos.second != slice.j || cellExists(params, previousSlice, pos.second - previousSlice.j - 1, pos.first));
			return { pos.first, pos.second - 1 };
		}
		assert(false);
		std::abort();
		return pos;
	}
private:

	OnewayTrace getTraceFromTable(const std::string& sequence, const DPTable& slice, AlignerGraphsizedState& reusableState) const
	{
		//todo fix
		// assert(slice.bandwidthPerSlice.size() == slice.correctness.size());
		// if (slice.slices.size() == 0)
		// {
		// 	return OnewayTrace::TraceFailed();
		// }
		// if (slice.bandwidthPerSlice.size() == 0)
		// {
		// 	return OnewayTrace::TraceFailed();
		// }
		// assert(slice.samplingFrequency > 1);
		// OnewayTrace result;
		// result.score = 0;
		// size_t backtraceOverrideIndex = -1;
		// LengthType lastBacktraceOverrideStartJ = -1;
		// LengthType nextBacktraceOverrideEndJ = -1;
		// if (slice.backtraceOverrides.size() > 0)
		// {
		// 	backtraceOverrideIndex = slice.backtraceOverrides.size()-1;
		// 	nextBacktraceOverrideEndJ = slice.backtraceOverrides.back().endj;
		// }
		// for (size_t i = slice.slices.size()-1; i < slice.slices.size(); i--)
		// {
		// 	if ((slice.slices[i].j + WordConfiguration<Word>::WordSize) / WordConfiguration<Word>::WordSize == slice.bandwidthPerSlice.size())
		// 	{
		// 		assert(i == slice.slices.size() - 1);
		// 		result.score = slice.slices.back().minScore;
		// 		result.trace.emplace_back(slice.slices.back().minScoreIndex, std::min( slice.slices.back().j + WordConfiguration<Word>::WordSize - 1, sequence.size()-1));
		// 		continue;
		// 	}
		// 	if (lastBacktraceOverrideStartJ == slice.slices[i].j + WordConfiguration<Word>::WordSize) continue;
		// 	auto partTable = getSlicesFromTable(sequence, lastBacktraceOverrideStartJ, slice, i, reusableState);
		// 	assert(partTable.size() > 0);
		// 	if (i == slice.slices.size() - 1)
		// 	{
		// 		result.score = partTable.back().minScore;
		// 		assert(partTable.back().minScoreIndex != -1);
		// 		result.trace.emplace_back(partTable.back().minScoreIndex, std::min(partTable.back().j + WordConfiguration<Word>::WordSize - 1, sequence.size()-1));
		// 		if (result.trace.back().second == partTable.back().j)
		// 		{
		// 			if (partTable.size() == 1)
		// 			{
		// 				auto boundaryTrace = getSliceBoundaryTrace(sequence, partTable.back(), slice.slices[i], result.trace.back().first);
		// 				result.trace.insert(result.trace.end(), boundaryTrace.begin(), boundaryTrace.end());
		// 				continue;
		// 			}
		// 			else
		// 			{
		// 				auto boundaryTrace = getSliceBoundaryTrace(sequence, partTable.back(), partTable[partTable.size()-2], result.trace.back().first);
		// 				result.trace.insert(result.trace.end(), boundaryTrace.begin(), boundaryTrace.end());
		// 				partTable.pop_back();
		// 			}
		// 		}
		// 	}
		// 	auto partTrace = getTraceFromTableInner(sequence, partTable, result.trace.back());
		// 	assert(partTrace.size() >= 1);
		// 	//begin()+1 because the starting position was already inserted earlier
		// 	result.trace.insert(result.trace.end(), partTrace.begin()+1, partTrace.end());
		// 	auto boundaryTrace = getSliceBoundaryTrace(sequence, partTable[0], slice.slices[i], result.trace.back().first);
		// 	result.trace.insert(result.trace.end(), boundaryTrace.begin(), boundaryTrace.end());
		// 	assert(boundaryTrace.size() > 0);
		// 	if (slice.slices[i].j == nextBacktraceOverrideEndJ)
		// 	{
		// 		auto trace = slice.backtraceOverrides[backtraceOverrideIndex].GetBacktrace(result.trace.back());
		// 		result.trace.insert(result.trace.end(), trace.begin()+1, trace.end());
		// 		lastBacktraceOverrideStartJ = slice.backtraceOverrides[backtraceOverrideIndex].startj;
		// 		backtraceOverrideIndex--;
		// 		if (backtraceOverrideIndex != -1) nextBacktraceOverrideEndJ = slice.backtraceOverrides[backtraceOverrideIndex].endj;
		// 		if (lastBacktraceOverrideStartJ == 0) break;
		// 	}
		// }
		// assert(result.trace.back().second == -1);
		// result.trace.pop_back();
		// assert(result.trace.back().second == 0);
		// std::reverse(result.trace.begin(), result.trace.end());
		// return result;
	}

	//returns the trace backwards, aka result[0] is at the bottom of the slice and result.back() at the top
	std::vector<MatrixPosition> getTraceFromSlice(const std::string& sequence, const DPSlice& slice, MatrixPosition pos) const
	{
		assert(pos.second >= slice.j);
		assert(pos.second < slice.j + WordConfiguration<Word>::WordSize);
		// auto distance = params.graph.MinDistance(startColumn, slice.minScoreIndex);
		// std::cerr << "distance from min: " << distance << std::endl;
		// auto score = getValue(slice, WordConfiguration<Word>::WordSize-1, startColumn);
		// std::cerr << "score from min: " << (score - slice.minScore) << std::endl;
		std::vector<MatrixPosition> result;
		while (pos.second != slice.j)
		{
			assert(slice.scores.hasNode(params.graph.IndexToNode(pos.first)));
			pos = pickBacktracePredecessor(params, sequence, slice, pos, slice);
			result.push_back(pos);
		}
		assert(slice.scores.hasNode(params.graph.IndexToNode(pos.first)));
		return result;
	}

	//returns the trace backwards, aka result[0] is after the boundary (later slice) and result.back() over it (earlier slice)
	std::vector<MatrixPosition> getSliceBoundaryTrace(const std::string& sequence, const DPSlice& after, const DPSlice& before, LengthType afterColumn) const
	{
		MatrixPosition pos { afterColumn, after.j };
		assert(after.j == before.j + WordConfiguration<Word>::WordSize);
		std::vector<MatrixPosition> result;
		while (pos.second == after.j)
		{
			assert(after.scores.hasNode(params.graph.IndexToNode(pos.first)));
			pos = pickBacktracePredecessor(params, sequence, after, pos, before);
			result.push_back(pos);
		}
		assert(before.scores.hasNode(params.graph.IndexToNode(pos.first)));
		return result;
	}

	//returns the trace backwards, aka result[0] is at the bottom of the table and result.back() at the top
	std::vector<MatrixPosition> getTraceFromTableInner(const std::string& sequence, const std::vector<DPSlice>& table, MatrixPosition pos) const
	{
		assert(table.size() > 0);
		assert(pos.second >= table.back().j);
		assert(pos.second < table.back().j + WordConfiguration<Word>::WordSize);
		std::vector<MatrixPosition> result;
		result.push_back(pos);
		for (size_t slice = table.size()-1; slice < table.size(); slice--)
		{
#ifdef SLICEVERBOSE
			auto node = params.graph.IndexToNode(result.back().first);
			auto offset = result.back().first - params.graph.NodeStart(node);
			assert(table[slice].scores.hasNode(node));
			std::cerr << "1: j " << result.back().second << " score " << table[slice].scores.node(node)[offset].scoreEnd << std::endl;
#endif
			assert(table[slice].j <= result.back().second);
			assert(table[slice].j + WordConfiguration<Word>::WordSize > result.back().second);
			auto partialTrace = getTraceFromSlice(sequence, table[slice], result.back());
			assert(partialTrace.size() >= result.back().second - table[slice].j);
			result.insert(result.end(), partialTrace.begin(), partialTrace.end());
			assert(result.back().second == table[slice].j);
			if (slice > 0)
			{
				auto boundaryTrace = getSliceBoundaryTrace(sequence, table[slice], table[slice-1], result.back().first);
				result.insert(result.end(), boundaryTrace.begin(), boundaryTrace.end());
				assert(result.back().second == table[slice-1].j + WordConfiguration<Word>::WordSize - 1);
			}
		}
		assert(result.back().second == table[0].j);
		assert(table[0].scores.hasNode(params.graph.IndexToNode(result.back().first)));
		return result;
	}

#ifdef EXTRABITVECTORASSERTIONS

	WordSlice getWordSliceCellByCell(size_t j, size_t w, const std::string& sequence, const NodeSlice<LengthType, ScoreType, Word>& currentSlice, const NodeSlice<LengthType, ScoreType, Word>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand) const
	{
		const auto lastBitMask = ((Word)1) << (WordConfiguration<Word>::WordSize-1);
		WordSlice result;
		auto nodeIndex = params.graph.IndexToNode(w);
		assert(currentBand[nodeIndex]);
		const std::vector<WordSlice>& oldNode = previousBand[nodeIndex] ? previousSlice.node(nodeIndex) : currentSlice.node(nodeIndex);
		assert(currentBand[nodeIndex]);
		ScoreType current[66];
		current[0] = j+1;
		current[1] = j;
		if (j > 0 && previousBand[nodeIndex]) current[1] = std::min(current[1], oldNode[w-params.graph.NodeStart(nodeIndex)].scoreEnd);
		if (j > 0 && previousBand[nodeIndex]) current[0] = std::min(current[0], oldNode[w-params.graph.NodeStart(nodeIndex)].scoreEnd - ((oldNode[w-params.graph.NodeStart(nodeIndex)].VP & lastBitMask) ? 1 : 0) + ((oldNode[w-params.graph.NodeStart(nodeIndex)].VN & lastBitMask) ? 1 : 0));
		for (int i = 1; i < 65; i++)
		{
			current[i+1] = current[i]+1;
		}
		if (w == params.graph.NodeStart(nodeIndex))
		{
			for (auto neighbor : params.graph.inNeighbors[nodeIndex])
			{
				if (!previousBand[neighbor] && !currentBand[neighbor]) continue;
				const std::vector<WordSlice>& neighborSlice = currentBand[neighbor] ? currentSlice.node(neighbor) : previousSlice.node(neighbor);
				const std::vector<WordSlice>& oldNeighborSlice = previousBand[neighbor] ? previousSlice.node(neighbor) : currentSlice.node(neighbor);
				auto u = params.graph.NodeEnd(neighbor)-1;
				ScoreType previous[66];
				previous[0] = j+1;
				previous[1] = j;
				if (j > 0 && previousBand[neighbor]) previous[1] = std::min(previous[1], oldNeighborSlice.back().scoreEnd);
				if (j > 0 && previousBand[neighbor]) previous[0] = std::min(previous[0], oldNeighborSlice.back().scoreEnd - ((oldNeighborSlice.back().VP & lastBitMask) ? 1 : 0) + ((oldNeighborSlice.back().VN & lastBitMask) ? 1 : 0));
				if (currentBand[neighbor]) previous[1] = std::min(previous[1], neighborSlice.back().scoreBeforeStart);
				for (int i = 1; i < 65; i++)
				{
					if (currentBand[neighbor])
					{
						previous[i+1] = previous[i];
						previous[i+1] += (neighborSlice.back().VP & (((Word)1) << (i-1)) ? 1 : 0);
						previous[i+1] -= (neighborSlice.back().VN & (((Word)1) << (i-1)) ? 1 : 0);
					}
					else
					{
						previous[i+1] = previous[i]+1;
					}
				}
				current[0] = std::min(current[0], previous[0]+1);
				for (int i = 0; i < 65; i++)
				{
					current[i+1] = std::min(current[i+1], previous[i+1]+1);
					current[i+1] = std::min(current[i+1], current[i]+1);
					if (j+i > 0 && (sequence[j+i-1] == params.graph.NodeSequences(w) || sequence[j+i-1] == 'N'))
					{
						current[i+1] = std::min(current[i+1], previous[i]);
					}
					else
					{
						current[i+1] = std::min(current[i+1], previous[i]+1);
					}
				}
			}
		}
		else
		{
			const std::vector<WordSlice>& slice = currentSlice.node(nodeIndex);
			const std::vector<WordSlice>& oldSlice = previousBand[nodeIndex] ? previousSlice.node(nodeIndex) : slice;
			auto u = w-1;
			ScoreType previous[66];
			previous[0] = slice[u-params.graph.NodeStart(nodeIndex)].scoreBeforeStart+1;
			previous[1] = slice[u-params.graph.NodeStart(nodeIndex)].scoreBeforeStart;
			if (previousBand[nodeIndex]) previous[0] = std::min(previous[0], oldSlice[u-params.graph.NodeStart(nodeIndex)].scoreEnd - ((oldSlice[u-params.graph.NodeStart(nodeIndex)].VP & lastBitMask) ? 1 : 0) + ((oldSlice[u-params.graph.NodeStart(nodeIndex)].VN & lastBitMask) ? 1 : 0));
			if (previousBand[nodeIndex]) previous[1] = std::min(previous[1], oldSlice[u-params.graph.NodeStart(nodeIndex)].scoreEnd);
			for (int i = 1; i < 65; i++)
			{
				previous[i+1] = previous[i];
				previous[i+1] += (slice[u-params.graph.NodeStart(nodeIndex)].VP & (((Word)1) << (i-1)) ? 1 : 0);
				previous[i+1] -= (slice[u-params.graph.NodeStart(nodeIndex)].VN & (((Word)1) << (i-1)) ? 1 : 0);
			}
			current[0] = std::min(current[0], previous[0]+1);
			for (int i = 0; i < 65; i++)
			{
				current[i+1] = std::min(current[i+1], current[i]+1);
				current[i+1] = std::min(current[i+1], previous[i+1]+1);
				if (j+i > 0 && (sequence[j+i-1] == params.graph.NodeSequences(w) || sequence[j+i-1] == 'N'))
				{
					current[i+1] = std::min(current[i+1], previous[i]);
				}
				else
				{
					current[i+1] = std::min(current[i+1], previous[i]+1);
				}
			}
		}
		for (int i = 1; i < 65; i++)
		{
			assert(current[i+1] >= debugLastRowMinScore);
			assert(current[i+1] >= current[i]-1);
			assert(current[i+1] <= current[i]+1);
			if (current[i+1] == current[i]+1) result.VP |= ((Word)1) << (i-1);
			if (current[i+1] == current[i]-1) result.VN |= ((Word)1) << (i-1);
		}
		result.scoreBeforeStart = current[1];
		result.scoreEnd = current[65];
		assert(result.scoreEnd == result.scoreBeforeStart + WordConfiguration<Word>::popcount(result.VP) - WordConfiguration<Word>::popcount(result.VN));
		return result;
	}

#endif

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

	NodeCalculationResult calculateNode(size_t i, typename NodeSlice<LengthType, ScoreType, Word>::NodeSliceMapItem& slice, const EqVector& EqV, typename NodeSlice<LengthType, ScoreType, Word>::NodeSliceMapItem& previousSlice, WordSlice ws, Word forceVN) const
	{
		NodeCalculationResult result;
		result.minScore = std::numeric_limits<ScoreType>::max();
		result.minScoreIndex = -1;
		result.cellsProcessed = 0;
		auto nodeStart = params.graph.NodeStart(i);
		auto nodeLength = params.graph.NodeLength(i);

		LengthType pos = 0;
		char graphChar = params.graph.NodeSequences(nodeStart+pos);
		Word Eq = EqV.getEq(graphChar);
		size_t chunk = pos / WordConfiguration<Word>::WordSize;
		size_t offset = pos % WordConfiguration<Word>::WordSize;

		//todo fix
		Word hinP = 0;
		Word hinN = 0;
		if (slice.exists)
		{
		}

		std::tie(ws, hinP, hinP) = BV::getNextSlice(Eq, ws, hinP, hinN);
		assert((ws.VP & forceVN) == WordConfiguration<Word>::AllZeros);
		ws.VN |= forceVN;
		result.cellsProcessed++;
		result.minScore = ws.scoreEnd;
		result.minScoreIndex = nodeStart+pos;

		if (slice.exists)
		{
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

		for (size_t i = 0; i < slice.NUM_CHUNKS; i++)
		{
			slice.HP[i] = WordConfiguration<Word>::AllZeros;
			slice.HN[i] = WordConfiguration<Word>::AllZeros;
		}
		slice.HP[chunk] |= hinP << offset;
		slice.HN[chunk] |= hinN << offset;

		slice.startSlice = ws;
		slice.exists = true;
		offset++;
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
				std::tie(ws, hinP, hinN) = BV::getNextSlice(Eq, ws, HN & 1, HN & 1);
				assert((ws.VP & forceVN) == WordConfiguration<Word>::AllZeros);
				ws.VN |= forceVN;
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

#ifdef EXTRACORRECTNESSASSERTIONS

	template <typename T>
	T volmin(volatile T& a, T b) const
	{
		return std::min((T)a, b);
	}

	void verifySliceBitvector(const std::string& sequence, const DPSlice& current, const DPSlice& previous) const
	{
		const ScoreType uninitScore = sequence.size() + 10000;
		const auto lastrow = WordConfiguration<Word>::WordSize - 1;
		for (auto pair : current.scores)
		{
			auto start = params.graph.NodeStart(pair.first);
			volatile ScoreType foundMinScore = uninitScore;
			volatile bool match = Common::characterMatch(sequence[current.j], params.graph.NodeSequences(start));
			if (current.j == 0 && previous.scores.hasNode(pair.first))
			{
				foundMinScore = volmin(foundMinScore, match ? 0 : 1);
			}
			if (previous.scores.hasNode(pair.first))
			{
				foundMinScore = volmin(foundMinScore, getValueIfExists(params, previous, lastrow, start, uninitScore)+1);
			}
			for (auto neighbor : params.graph.inNeighbors[pair.first])
			{
				if (current.scores.hasNode(neighbor))
				{
					foundMinScore = volmin(foundMinScore, getValueIfExists(params, current, 0, params.graph.NodeEnd(neighbor)-1, uninitScore)+1);
				}
				if (previous.scores.hasNode(neighbor))
				{
					foundMinScore = volmin(foundMinScore, getValueIfExists(params, previous, lastrow, params.graph.NodeEnd(neighbor)-1, uninitScore) + (match ? 0 : 1));
				}
			}
			if (cellExists(pair.second[0], 0, current.minScore + current.bandwidth)) assert(getValueIfExists(params, current, 0, start, uninitScore) == foundMinScore);
			for (int j = 1; j < WordConfiguration<Word>::WordSize; j++)
			{
				foundMinScore = uninitScore;
				match = Common::characterMatch(sequence[current.j+j], params.graph.NodeSequences(start));
				foundMinScore = volmin(foundMinScore, getValueIfExists(params, current, j-1, start, uninitScore)+1);
				for (auto neighbor : params.graph.inNeighbors[pair.first])
				{
					if (!current.scores.hasNode(neighbor)) continue;
					foundMinScore = volmin(foundMinScore, getValueIfExists(params, current, j, params.graph.NodeEnd(neighbor)-1, uninitScore)+1);
					foundMinScore = volmin(foundMinScore, getValueIfExists(params, current, j-1, params.graph.NodeEnd(neighbor)-1, uninitScore)+(match ? 0 : 1));
				}
				if (cellExists(pair.second[0], j, current.minScore + current.bandwidth)) assert(getValueIfExists(params, current, j, start, uninitScore) == foundMinScore);
			}
			for (size_t i = 1; i < pair.second.size(); i++)
			{
				volatile bool match = Common::characterMatch(sequence[current.j], params.graph.NodeSequences(start+i));
				volatile ScoreType foundMinScore = uninitScore;
				foundMinScore = volmin(foundMinScore, getValueIfExists(params, current, 0, start+i-1, uninitScore)+1);
				if (previous.scores.hasNode(pair.first))
				{
					foundMinScore = volmin(foundMinScore, getValueIfExists(params, previous, lastrow, start+i, uninitScore)+1);
					foundMinScore = volmin(foundMinScore, getValueIfExists(params, previous, lastrow, start+i-1, uninitScore) + (match ? 0 : 1));
				}
				if (cellExists(pair.second[i], 0, current.minScore + current.bandwidth)) assert(getValueIfExists(params, current, 0, start+i, uninitScore) == foundMinScore);
				for (int j = 1; j < WordConfiguration<Word>::WordSize; j++)
				{
					match = Common::characterMatch(sequence[current.j+j], params.graph.NodeSequences(start+i));
					foundMinScore = uninitScore;
					foundMinScore = volmin(foundMinScore, getValueIfExists(params, current, j-1, start+i, uninitScore)+1);
					foundMinScore = volmin(foundMinScore, getValueIfExists(params, current, j, start+i-1, uninitScore)+1);
					foundMinScore = volmin(foundMinScore, getValueIfExists(params, current, j-1, start+i-1, uninitScore)+(match ? 0 : 1));
					if (cellExists(pair.second[i], j, current.minScore + current.bandwidth)) assert(getValueIfExists(params, current, j, start+i, uninitScore) == foundMinScore);
				}
			}
		}
	}

#endif

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

	NodeCalculationResult calculateSlice(const std::string& sequence, size_t j, NodeSlice<LengthType, ScoreType, Word>& currentSlice, const NodeSlice<LengthType, ScoreType, Word>& previousSlice, const std::vector<LengthType>& previousNodes, std::vector<bool>& currentBand, const std::vector<bool>& previousBand, ArrayPriorityQueue<EdgeWithPriority>& calculableQueue, ScoreType previousQuitScore, int bandwidth, ScoreType previousMinScore) const
	{
		ScoreType currentMinimumScore = std::numeric_limits<ScoreType>::max() - bandwidth - 1;
		LengthType currentMinimumIndex;
		size_t cellsProcessed = 0;

		//preprocessed bitvectors for character equality
		Word BA = WordConfiguration<Word>::AllZeros;
		Word BT = WordConfiguration<Word>::AllZeros;
		Word BC = WordConfiguration<Word>::AllZeros;
		Word BG = WordConfiguration<Word>::AllZeros;
		for (int i = 0; i < WordConfiguration<Word>::WordSize && j+i < sequence.size(); i++)
		{
			Word mask = ((Word)1) << i;
			if (Common::characterMatch(sequence[j+i], 'A')) BA |= mask;
			if (Common::characterMatch(sequence[j+i], 'C')) BC |= mask;
			if (Common::characterMatch(sequence[j+i], 'T')) BT |= mask;
			if (Common::characterMatch(sequence[j+i], 'G')) BG |= mask;
		}
		if (j + WordConfiguration<Word>::WordSize > sequence.size())
		{
			Word mask = WordConfiguration<Word>::AllOnes << (WordConfiguration<Word>::WordSize - j + sequence.size());
			assert((BA & mask) == 0);
			assert((BT & mask) == 0);
			assert((BC & mask) == 0);
			assert((BG & mask) == 0);
			BA |= mask;
			BT |= mask;
			BC |= mask;
			BG |= mask;
		}
		assert((BA | BC | BT | BG) == WordConfiguration<Word>::AllOnes);
		EqVector EqV {BA, BT, BC, BG};

		assert(previousNodes.size() > 0);
		for (auto node : previousSlice)
		{
			if (node.second.minScore <= previousQuitScore)
			{
				WordSlice startSlice = getSourceSliceFromScore(node.second.startSlice.scoreEnd);
				calculableQueue.insert(node.second.minScore, EdgeWithPriority { node.first, node.second.minScore, startSlice });
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
			Word forceVN = 0;

			if (previousBand[i])
			{
				previousThisNode = previousSlice.node(i);
			}
			else
			{
				forceVN = 1;
				for (size_t chunk = 0; chunk < previousThisNode.NUM_CHUNKS; chunk++)
				{
					previousThisNode.HP[chunk] = WordConfiguration<Word>::AllOnes;
					previousThisNode.HN[chunk] = WordConfiguration<Word>::AllZeros;
				}
			}
			auto nodeCalc = calculateNode(i, thisNode, EqV, previousThisNode, incoming, forceVN);
			assert(nodeCalc.minScore <= previousQuitScore + WordConfiguration<Word>::WordSize);
			currentMinScoreAtEndRow = std::min(currentMinScoreAtEndRow, nodeCalc.minScore);
			currentSlice.setMinScoreIfSmaller(i, nodeCalc.minScore);
			auto newEnd = thisNode.endSlice;

			if (newEnd.scoreEnd != oldEnd.scoreEnd || newEnd.VP != oldEnd.VP || newEnd.VN != oldEnd.VN)
			{
				ScoreType newEndMinScore = newEnd.changedMinScore(oldEnd);
				assert(newEndMinScore != std::numeric_limits<ScoreType>::max());
				if (newEndMinScore <= currentMinScoreAtEndRow + bandwidth)
				{
					Word hinP = 1;
					Word hinN = 0;
					if (previousSlice.hasNode(i))
					{
						auto prevNode = previousSlice.node(i);
						auto size = params.graph.NodeLength(i);
						size_t chunk = size / WordConfiguration<Word>::WordSize;
						size_t offset = size % WordConfiguration<Word>::WordSize;
						hinP = (prevNode.HP[chunk] >> offset) & 1;
						hinN = (prevNode.HN[chunk] >> offset) & 1;
					}
					for (auto neighbor : params.graph.outNeighbors[i])
					{
						calculableQueue.insert(newEndMinScore - previousMinScore, EdgeWithPriority { neighbor, newEndMinScore - previousMinScore, newEnd });
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

		NodeCalculationResult result;
		result.minScore = currentMinimumScore;
		result.minScoreIndex = currentMinimumIndex;
		result.cellsProcessed = cellsProcessed;

		if (j + WordConfiguration<Word>::WordSize > sequence.size())
		{
			flattenLastSliceEnd(currentSlice, result, j, sequence.size());
		}

		calculableQueue.clear();

		return result;
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
			slice.nodes.push_back(node.first);
			slice.numCells += params.graph.NodeLength(node.first);
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

#ifdef EXTRACORRECTNESSASSERTIONS
		if (bandTest.cellsProcessed <= params.maxCellsPerSlice && sequence.size() >= bandTest.j + WordConfiguration<Word>::WordSize) verifySliceBitvector(sequence, bandTest, previous);
#endif
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
			size_t debugSmallCells = 0;
			for (auto node : newSlice.scores)
			{
				for (auto cell : node.second)
				{
					if (cell.scoreEnd <= newSlice.minScore + bandwidth) debugSmallCells++;
				}
			}
			std::cerr << " small endcells " << debugSmallCells;
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
				while (result.bandwidthPerSlice.size() > slice+1) result.bandwidthPerSlice.pop_back();
				while (result.correctness.size() > slice+1) result.correctness.pop_back();
				while (result.slices.size() > 1 && result.slices.back().j > slice * WordConfiguration<Word>::WordSize) result.slices.pop_back();
#ifdef SLICEVERBOSE
				std::cerr << " ramp to " << slice;
#endif
#ifdef SLICEVERBOSE
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

			assert(result.bandwidthPerSlice.size() == slice);
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
			assert(reusableState.nodesliceMap[i].start == 0);
			assert(reusableState.nodesliceMap[i].end == 0);
		}
#endif

#ifndef NDEBUG
		assert(result.slices.size() == result.bandwidthPerSlice.size());
		assert(result.correctness.size() == result.slices.size());
		if (result.slices.size() > 0)
		{
			volatile size_t lastExisting = 0;
			assert(result.bandwidthPerSlice.size() == debugLastProcessedSlice + 1);
			for (size_t i = 0; i < result.slices.size(); i++)
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
