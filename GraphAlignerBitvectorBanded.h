#ifndef GraphAlignerBitvectorBanded_h
#define GraphAlignerBitvectorBanded_h
/*
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
	using Params = typename Common::Params;
	using MatrixPosition = typename Common::MatrixPosition;
	using Trace = typename Common::Trace;
	using OnewayTrace = typename Common::OnewayTrace;
	using SeedHit = typename Common::SeedHit;
	using WordSlice = typename WordContainer<LengthType, ScoreType, Word>::Slice;
	using EqVector = typename BV::EqVector;
	using NodeWithPriority = typename BV::NodeWithPriority;
	mutable ArrayPriorityQueue<NodeWithPriority> calculableQueue;
	mutable std::vector<typename NodeSlice<WordSlice>::MapItem> nodesliceMap;
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
		DPSlice(std::vector<typename NodeSlice<WordSlice>::MapItem>* vectorMap) :
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
		NodeSlice<WordSlice> scores;
		std::vector<size_t> nodes;
		AlignmentCorrectnessEstimationState correctness;
		LengthType j;
		int bandwidth;
		size_t cellsProcessed;
		size_t numCells;
		size_t EstimatedMemoryUsage() const
		{
			return numCells * sizeof(typename WordContainer<LengthType, ScoreType, Word>::TinySlice) + scores.size() * (sizeof(size_t) * 3 + sizeof(int));
		}
		DPSlice getFrozenSqrtEndScores() const
		{
			DPSlice result;
			result.scores = scores.getFrozenSqrtEndScores();
			result.minScore = minScore;
			result.minScoreIndex = minScoreIndex;
			result.nodes = nodes;
			result.correctness = correctness;
			result.j = j;
			result.bandwidth = bandwidth;
			result.cellsProcessed = cellsProcessed;
			result.numCells = numCells;
			return result;
		}
		DPSlice getFrozenScores() const
		{
			DPSlice result;
			result.scores = scores.getFrozenScores();
			result.minScore = minScore;
			result.minScoreIndex = minScoreIndex;
			result.nodes = nodes;
			result.correctness = correctness;
			result.j = j;
			result.bandwidth = bandwidth;
			result.cellsProcessed = cellsProcessed;
			result.numCells = numCells;
			return result;
		}
	};
	class BacktraceOverride
	{
	public:
		class BacktraceItem
		{
		public:
			BacktraceItem(bool previousInSameRow, size_t previousIndex, MatrixPosition pos) :
			end(false),
			previousInSameRow(previousInSameRow),
			previousIndex(previousIndex),
			pos(pos)
#ifdef SLICEVERBOSE
			,score(0)
#endif
			{}
			bool end;
			bool previousInSameRow;
			size_t previousIndex;
			MatrixPosition pos;
#ifdef SLICEVERBOSE
			ScoreType score;
#endif
		};
		BacktraceOverride()
		{
		}
		BacktraceOverride(const Params& params, const std::string& sequence, const DPSlice& previous, const std::vector<DPSlice>& slices)
		{
			assert(slices.size() > 0);
			startj = slices[0].j;
			endj = slices.back().j;
			assert(endj == startj + (slices.size()-1) * WordConfiguration<Word>::WordSize);
			if (startj + WordConfiguration<Word>::WordSize * slices.size() <= sequence.size())
			{
				items.resize(WordConfiguration<Word>::WordSize * slices.size());
			}
			else
			{
				items.resize(sequence.size() - startj);
			}
			makeTrace(params, sequence, previous, slices);
		};
		//returns the trace backwards, aka result[0] is at the bottom of the slice and result.back() at the top
		std::vector<MatrixPosition> GetBacktrace(MatrixPosition start) const
		{
			assert(items.size() > 0);
			assert(items.back().size() > 0);
			assert(items.back()[0].pos.second == start.second);
			size_t currentIndex = -1;
			size_t currentRow = items.size()-1;
			std::vector<MatrixPosition> result;
			for (size_t i = 0; i < items.back().size(); i++)
			{
				if (items.back()[i].pos == start)
				{
					currentIndex = i;
					break;
				}
			}
			assert(currentIndex != -1);
#ifdef SLICEVERBOSE
			std::cerr << "2: j " << items[currentRow][currentIndex].pos.second << " score " << items[currentRow][currentIndex].score << std::endl;
#endif
			while (true)
			{
				auto current = items[currentRow][currentIndex];
				assert(!current.end);
				result.push_back(current.pos);
#ifdef SLICEVERBOSE
				if (current.pos.second % WordConfiguration<Word>::WordSize == WordConfiguration<Word>::WordSize - 1)
				{
					std::cerr << "2: j " << items[currentRow][currentIndex].pos.second << " score " << items[currentRow][currentIndex].score << std::endl;
				}
#endif
				size_t nextIndex = current.previousIndex;
				size_t nextRow = current.previousInSameRow ? currentRow : currentRow - 1;
				if (nextRow == -1)
				{
					result.emplace_back(nextIndex, current.pos.second-1);
					break;
				}
				currentIndex = nextIndex;
				currentRow = nextRow;
			}
			return result;
		}
		LengthType startj;
		LengthType endj;
	private:
		void addReachableRec(const Params& params, MatrixPosition pos, size_t row, const std::string& sequence, const DPSlice& previous, const std::vector<DPSlice>& slices, std::vector<std::unordered_map<LengthType, size_t>>& indices)
		{
		//manual tail call
		//because gcc doesn't optimize it for some reason and recursion runs into a stack overflow
		start:
			assert(slices.size() > 0);
			assert(row < indices.size());
			if (indices[row].count(pos.first) == 1) return;
			size_t size = indices[row].size();
			indices[row][pos.first] = size;
			if (row > 0 && row % WordConfiguration<Word>::WordSize == WordConfiguration<Word>::WordSize - 1)
			{
				size_t sliceIndex = row / WordConfiguration<Word>::WordSize;
				assert(sliceIndex < slices.size());
				auto nodeIndex = params.graph.IndexToNode(pos.first);
				assert(slices[sliceIndex].scores.hasNode(nodeIndex));
				auto nodeStart = params.graph.NodeStart(nodeIndex);
				auto offset = pos.first - nodeStart;
				assert(slices[sliceIndex].scores.node(nodeIndex).size() > offset);
				if (slices[sliceIndex].scores.node(nodeIndex)[offset].scoreEnd > slices[sliceIndex].minScore + slices[sliceIndex].bandwidth) return;
			}
			assert(row == pos.second - slices[0].j);
			size_t sliceIndex = row / WordConfiguration<Word>::WordSize;
			assert(sliceIndex < slices.size());
			MatrixPosition predecessor;
			if (sliceIndex > 0)
			{
				predecessor = pickBacktracePredecessor(params, sequence, slices[sliceIndex], pos, slices[sliceIndex-1]);
			}
			else
			{
				predecessor = pickBacktracePredecessor(params, sequence, slices[0], pos, previous);
			}
			assert(predecessor.second == pos.second || predecessor.second == pos.second-1);
			if (predecessor.second >= slices[0].j && predecessor.second != -1)
			{
				//manual tail call
				//because gcc doesn't optimize it for some reason and recursion runs into a stack overflow
				pos = predecessor;
				row = predecessor.second - slices[0].j;
				goto start;
			}
		}
		void makeTrace(const Params& params, const std::string& sequence, const DPSlice& previous, const std::vector<DPSlice>& slices)
		{
			assert(slices.size() > 0);
			std::vector<std::unordered_map<LengthType, size_t>> indexOfPos;
			indexOfPos.resize(items.size());
			size_t endrow = items.size()-1;
			{
				std::vector<MatrixPosition> endPositions;
				for (const auto& pair : slices.back().scores)
				{
					LengthType nodeStart = params.graph.NodeStart(pair.first);
					LengthType endj = slices[0].j + endrow;
					for (size_t i = 0; i < pair.second.size(); i++)
					{
						if (pair.second[i].scoreEnd <= slices.back().minScore + slices.back().bandwidth)
						{
							endPositions.emplace_back(nodeStart+i, endj);
						}
					}
				}
				for (size_t i = 0; i < indexOfPos.size(); i++)
				{
					indexOfPos[i].reserve(endPositions.size());
				}
				for (auto pos : endPositions)
				{
					addReachableRec(params, pos, endrow, sequence, previous, slices, indexOfPos);
				}
#ifdef SLICEVERBOSE
				std::cerr << " endcells " << endPositions.size();
#endif
			}

			for (size_t row = items.size()-1; row < items.size(); row--)
			{
				items[row].resize(indexOfPos[row].size(), BacktraceItem {false, 0, MatrixPosition {0, 0}});
				for (auto pair : indexOfPos[row])
				{
					auto w = pair.first;
					auto index = pair.second;
					MatrixPosition pos { w, slices[0].j + row };
					items[row][index].pos = pos;
#ifdef SLICEVERBOSE
					if (pos.second % WordConfiguration<Word>::WordSize == WordConfiguration<Word>::WordSize - 1)
					{
						auto sliceIndex = row / WordConfiguration<Word>::WordSize;
						auto node = params.graph.IndexToNode(pos.first);
						auto offset = pos.first - params.graph.NodeStart(node);
						assert(sliceIndex < slices.size());
						assert(slices[sliceIndex].scores.hasNode(node));
						items[row][index].score = slices[sliceIndex].scores.node(node)[offset].scoreEnd;
					}
#endif
					MatrixPosition predecessor;
					size_t sliceIndex = row / WordConfiguration<Word>::WordSize;
					if (row % WordConfiguration<Word>::WordSize == WordConfiguration<Word>::WordSize - 1)
					{
						auto nodeIndex = params.graph.IndexToNode(w);
						auto offset = w - params.graph.NodeStart(nodeIndex);
						assert(slices[sliceIndex].scores.hasNode(nodeIndex));
						if (slices[sliceIndex].scores.node(nodeIndex)[offset].scoreEnd > slices[sliceIndex].minScore + slices[sliceIndex].bandwidth)
						{
							items[row][index].end = true;
							continue;
						}
					}
					if (sliceIndex > 0)
					{
						predecessor = pickBacktracePredecessor(params, sequence, slices[sliceIndex], pos, slices[sliceIndex-1]);
					}
					else
					{
						predecessor = pickBacktracePredecessor(params, sequence, slices[0], pos, previous);
					}
					if (predecessor.second == pos.second)
					{
						items[row][index].previousInSameRow = true;
						items[row][index].previousIndex = indexOfPos[row].at(predecessor.first);
					}
					else
					{
						items[row][index].previousInSameRow = false;
						if (row != 0)
						{
							items[row][index].previousIndex = indexOfPos[row-1].at(predecessor.first);
						}
						else
						{
							items[row][index].previousIndex = predecessor.first;
						}
					}
				}
#ifndef NDEBUG
				for (size_t i = 0; i < items[row].size(); i++)
				{
					assert(items[row][i].end || items[row][i].pos.first != 0);
				}
#endif
			}
		}
		std::vector<std::vector<BacktraceItem>> items;
	};
	class DPTable
	{
	public:
		DPTable() :
		slices(),
		samplingFrequency(0)
		{}
		std::vector<DPSlice> slices;
		size_t samplingFrequency;
		std::vector<LengthType> bandwidthPerSlice;
		std::vector<AlignmentCorrectnessEstimationState> correctness;
		std::vector<BacktraceOverride> backtraceOverrides;
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
	calculableQueue(WordConfiguration<Word>::WordSize + std::max(params.initialBandwidth, params.rampBandwidth) + 1),
	nodesliceMap(),
	params(params)
	{
		nodesliceMap.resize(params.graph.NodeSize(), {0, 0, 0, 0, 0});
	}

	OnewayTrace getTraceFromSeed(const std::string& sequence, int bigraphNodeId) const
	{
		std::vector<size_t> nodes;
		nodes = params.graph.nodeLookup.at(bigraphNodeId);
		assert(sequence.size() >= params.graph.DBGOverlap);
		size_t samplingFrequency = getSamplingFrequency(sequence.size());
		size_t numSlices = (sequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		auto initialBandwidth = getInitialSliceOneNodeGroup(nodes);
		auto slice = getSqrtSlices(sequence, initialBandwidth, numSlices, samplingFrequency);
		removeWronglyAlignedEnd(slice);
		if (slice.slices.size() == 0)
		{
			assert(slice.bandwidthPerSlice.size() == 0);
			return OnewayTrace::TraceFailed();
		}
		assert(slice.slices.back().minScore <= sequence.size() + WordConfiguration<Word>::WordSize * 2);

		OnewayTrace result;

		result = getTraceFromTable(sequence, slice);
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

	OnewayTrace getTraceFromTable(const std::string& sequence, const DPTable& slice) const
	{
		assert(slice.bandwidthPerSlice.size() == slice.correctness.size());
		if (slice.slices.size() == 0)
		{
			return OnewayTrace::TraceFailed();
		}
		if (slice.bandwidthPerSlice.size() == 0)
		{
			return OnewayTrace::TraceFailed();
		}
		assert(slice.samplingFrequency > 1);
		OnewayTrace result;
		result.score = 0;
		size_t backtraceOverrideIndex = -1;
		LengthType lastBacktraceOverrideStartJ = -1;
		LengthType nextBacktraceOverrideEndJ = -1;
		if (slice.backtraceOverrides.size() > 0)
		{
			backtraceOverrideIndex = slice.backtraceOverrides.size()-1;
			nextBacktraceOverrideEndJ = slice.backtraceOverrides.back().endj;
		}
		for (size_t i = slice.slices.size()-1; i < slice.slices.size(); i--)
		{
			if ((slice.slices[i].j + WordConfiguration<Word>::WordSize) / WordConfiguration<Word>::WordSize == slice.bandwidthPerSlice.size())
			{
				assert(i == slice.slices.size() - 1);
				result.score = slice.slices.back().minScore;
				result.trace.emplace_back(slice.slices.back().minScoreIndex, std::min( slice.slices.back().j + WordConfiguration<Word>::WordSize - 1, sequence.size()-1));
				continue;
			}
			if (lastBacktraceOverrideStartJ == slice.slices[i].j + WordConfiguration<Word>::WordSize) continue;
			auto partTable = getSlicesFromTable(sequence, lastBacktraceOverrideStartJ, slice, i);
			assert(partTable.size() > 0);
			if (i == slice.slices.size() - 1)
			{
				result.score = partTable.back().minScore;
				assert(partTable.back().minScoreIndex != -1);
				result.trace.emplace_back(partTable.back().minScoreIndex, std::min(partTable.back().j + WordConfiguration<Word>::WordSize - 1, sequence.size()-1));
				if (result.trace.back().second == partTable.back().j)
				{
					if (partTable.size() == 1)
					{
						auto boundaryTrace = getSliceBoundaryTrace(sequence, partTable.back(), slice.slices[i], result.trace.back().first);
						result.trace.insert(result.trace.end(), boundaryTrace.begin(), boundaryTrace.end());
						continue;
					}
					else
					{
						auto boundaryTrace = getSliceBoundaryTrace(sequence, partTable.back(), partTable[partTable.size()-2], result.trace.back().first);
						result.trace.insert(result.trace.end(), boundaryTrace.begin(), boundaryTrace.end());
						partTable.pop_back();
					}
				}
			}
			auto partTrace = getTraceFromTableInner(sequence, partTable, result.trace.back());
			assert(partTrace.size() >= 1);
			//begin()+1 because the starting position was already inserted earlier
			result.trace.insert(result.trace.end(), partTrace.begin()+1, partTrace.end());
			auto boundaryTrace = getSliceBoundaryTrace(sequence, partTable[0], slice.slices[i], result.trace.back().first);
			result.trace.insert(result.trace.end(), boundaryTrace.begin(), boundaryTrace.end());
			assert(boundaryTrace.size() > 0);
			if (slice.slices[i].j == nextBacktraceOverrideEndJ)
			{
				auto trace = slice.backtraceOverrides[backtraceOverrideIndex].GetBacktrace(result.trace.back());
				result.trace.insert(result.trace.end(), trace.begin()+1, trace.end());
				lastBacktraceOverrideStartJ = slice.backtraceOverrides[backtraceOverrideIndex].startj;
				backtraceOverrideIndex--;
				if (backtraceOverrideIndex != -1) nextBacktraceOverrideEndJ = slice.backtraceOverrides[backtraceOverrideIndex].endj;
				if (lastBacktraceOverrideStartJ == 0) break;
			}
		}
		assert(result.trace.back().second == -1);
		result.trace.pop_back();
		assert(result.trace.back().second == 0);
		std::reverse(result.trace.begin(), result.trace.end());
		return result;
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

	WordSlice getWordSliceCellByCell(size_t j, size_t w, const std::string& sequence, const NodeSlice<WordSlice>& currentSlice, const NodeSlice<WordSlice>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand) const
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

	WordSlice getNodeStartSlice(const Word Eq, const size_t nodeIndex, const NodeSlice<WordSlice>& previousSlice, const NodeSlice<WordSlice>& currentSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand, const bool previousEq, ScoreType quitScore, ScoreType previousSliceQuitScore) const
	{
		WordSlice result;
		WordSlice up;
		if (previousBand[nodeIndex]) up = previousSlice.node(nodeIndex)[0];
		if (!up.sliceExists) up = WordSlice {};
		bool foundOne = false;
		for (auto neighbor : params.graph.inNeighbors[nodeIndex])
		{
			if (!currentBand[neighbor] && !previousBand[neighbor]) continue;
			Word EqHere = Eq;
			WordSlice previous;
			WordSlice previousUp;
			bool foundOneUp = false;
			bool hasRealNeighbor = false;
			bool foundSomething = false;
#ifndef NDEBUG
			if (currentBand[neighbor] && previousBand[neighbor])
			{
				WordSlice neighborSlice = currentSlice.node(neighbor).back();
				WordSlice upNeighborSlice = previousSlice.node(neighbor).back();
				if (neighborSlice.sliceExists && neighborSlice.minScore <= quitScore && upNeighborSlice.sliceExists && upNeighborSlice.scoreEnd <= previousSliceQuitScore)
				{
					BV::assertSliceCorrectness(neighborSlice, upNeighborSlice, previousBand[neighbor]);
				}
			}
#endif
			if (currentBand[neighbor])
			{
				auto possiblePrevious = currentSlice.node(neighbor).back();
				if (possiblePrevious.sliceExists && possiblePrevious.minScore <= quitScore)
				{
					previous = possiblePrevious;
					hasRealNeighbor = true;
					foundSomething = true;
				}
			}
			if (previousBand[neighbor])
			{
				auto possiblePreviousUp = previousSlice.node(neighbor).back();
				if (possiblePreviousUp.sliceExists && possiblePreviousUp.scoreEnd <= previousSliceQuitScore)
				{
					previousUp = possiblePreviousUp;
					foundOneUp = true;
					foundSomething = true;
					auto possiblePrevious = getSourceSliceFromScore(possiblePreviousUp.scoreEnd);
					if (!previous.sliceExists)
					{
						previous = possiblePrevious;
						previous.scoreBeforeExists = true;
						previous.sliceExists = true;
					}
					else if (possiblePrevious.scoreEnd < previous.scoreBeforeStart)
					{
						previous = previous.mergeWithVertical(possiblePrevious);
						previous.scoreBeforeExists = true;
						previous.sliceExists = true;
					}
					assert(previousBand[neighbor]);
					foundSomething = true;
				}
			}
			if (!foundSomething) continue;
			BV::assertSliceCorrectness(previous, previousUp, foundOneUp);
			if (!hasRealNeighbor) EqHere &= 1;
			auto resultHere = BV::getNextSlice(EqHere, previous, up.sliceExists, up.sliceExists && foundOneUp, foundOneUp, previousEq, previousUp, up);
			if (!foundOne)
			{
				result = resultHere;
				foundOne = true;
			}
			else
			{
				result = result.mergeWith(resultHere);
			}
		}
		assert(foundOne);
		return result;
	}

	WordSlice getSourceSliceFromScore(ScoreType previousScore) const
	{
		WordSlice result { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, previousScore+WordConfiguration<Word>::WordSize, previousScore, false };
		result.sliceExists = true;
		return result;
	}

	WordSlice getSourceSliceFromStartMatch(char sequenceChar, char graphChar, ScoreType previousScore) const
	{
		Word firstVP = Common::characterMatch(sequenceChar, graphChar) ? 0 : 1;
		WordSlice result { WordConfiguration<Word>::AllOnes & ~(Word)1 | firstVP, WordConfiguration<Word>::AllZeros, previousScore+WordConfiguration<Word>::WordSize - 1 + firstVP, previousScore, true };
		result.sliceExists = true;
		return result;
	}

	WordSlice getSourceSliceFromBefore(size_t nodeIndex, const NodeSlice<WordSlice>& previousSlice) const
	{
		auto previousWordSlice = previousSlice.node(nodeIndex)[0];
		WordSlice result { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, previousWordSlice.scoreEnd+WordConfiguration<Word>::WordSize, previousWordSlice.scoreEnd, true };
		result.sliceExists = true;
		return result;
	}

	bool isSource(size_t nodeIndex, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand, const NodeSlice<WordSlice>& currentSlice, const NodeSlice<WordSlice>& previousSlice, ScoreType quitScore, ScoreType previousSliceQuitScore) const
	{
		for (auto neighbor : params.graph.inNeighbors[nodeIndex])
		{
			if (currentBand[neighbor])
			{
				auto slice = currentSlice.node(neighbor).back();
				if (slice.sliceExists && slice.minScore <= quitScore) return false;
			}
			if (previousBand[neighbor])
			{
				auto slice = previousSlice.node(neighbor).back();
				if (slice.sliceExists && slice.scoreEnd <= previousSliceQuitScore) return false;
			}
		}
		return true;
	}

	class NodeCalculationResult
	{
	public:
		ScoreType minScore;
		LengthType minScoreIndex;
		size_t cellsProcessed;
	};

	NodeCalculationResult calculateNode(size_t i, size_t j, size_t startIndex, size_t endIndex, const std::string& sequence, const EqVector& EqV, NodeSlice<WordSlice>& currentSlice, const NodeSlice<WordSlice>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand, ScoreType previousSliceQuitScore, ScoreType quitScore, int bandwidth) const
	{
		NodeCalculationResult result;
		result.minScore = std::numeric_limits<ScoreType>::max();
		result.minScoreIndex = -1;
		result.cellsProcessed = 0;
		auto slice = currentSlice.node(i);
		const auto oldSlice = previousBand[i] ? previousSlice.node(i) : slice;
		assert(slice.size() == params.graph.NodeEnd(i) - params.graph.NodeStart(i));
		assert(startIndex < slice.size());
		auto nodeStart = params.graph.NodeStart(i);

		WordSlice oldWordSlice;
		WordSlice currentWordSlice;
		WordSlice upWordSlice;

		oldWordSlice = slice[startIndex];
		upWordSlice = oldSlice[startIndex];

		bool upExists = previousBand[i];
		bool upWordSliceExists = upExists && upWordSlice.sliceExists && upWordSlice.scoreEnd <= previousSliceQuitScore;

		if (startIndex == 0)
		{
			if (isSource(i, currentBand, previousBand, currentSlice, previousSlice, quitScore, previousSliceQuitScore))
			{
				assert(upWordSliceExists);
				if (j == 0)
				{
					currentWordSlice = getSourceSliceFromStartMatch(sequence[0], params.graph.NodeSequences(nodeStart), previousSlice.node(i)[0].scoreEnd);
				}
				else
				{
					currentWordSlice = getSourceSliceFromBefore(i, previousSlice);
				}
				if (currentWordSlice.scoreEnd < result.minScore)
				{
					result.minScore = currentWordSlice.scoreEnd;
					result.minScoreIndex = nodeStart;
					quitScore = std::min(quitScore, result.minScore + bandwidth);
				}
				BV::assertSliceCorrectness(currentWordSlice, upWordSlice, upWordSliceExists);
			}
			else
			{
				Word Eq = EqV.getEq(params.graph.NodeSequences(nodeStart));
				currentWordSlice = getNodeStartSlice(Eq, i, previousSlice, currentSlice, currentBand, previousBand, (j == 0 && previousBand[i]) || (j > 0 && params.graph.NodeSequences(params.graph.NodeStart(i)) == sequence[j-1]), quitScore, previousSliceQuitScore);
				if (upWordSliceExists && currentWordSlice.scoreBeforeStart > upWordSlice.scoreEnd)
				{
					auto mergable = getSourceSliceFromScore(upWordSlice.scoreEnd);
					mergable.scoreBeforeExists = true;
					currentWordSlice = currentWordSlice.mergeWithVertical(mergable);
				}
				if (currentWordSlice.scoreEnd < result.minScore)
				{
					result.minScore = currentWordSlice.scoreEnd;
					result.minScoreIndex = nodeStart;
					quitScore = std::min(quitScore, result.minScore + bandwidth);
				}
				BV::assertSliceCorrectness(currentWordSlice, upWordSlice, upWordSliceExists);
			}
		}
		else
		{
			assert(upWordSliceExists);
			if (j == 0)
			{
				currentWordSlice = getSourceSliceFromStartMatch(sequence[0], params.graph.NodeSequences(nodeStart), upWordSlice.scoreEnd);
			}
			else
			{
				currentWordSlice = getSourceSliceFromScore(upWordSlice.scoreEnd);
				currentWordSlice.scoreBeforeExists = true;
			}
			if (currentWordSlice.scoreEnd < result.minScore)
			{
				result.minScore = currentWordSlice.scoreEnd;
				result.minScoreIndex = nodeStart + startIndex;
				quitScore = std::min(quitScore, result.minScore + bandwidth);
			}
			BV::assertSliceCorrectness(currentWordSlice, upWordSlice, upWordSliceExists);
		}

		if (oldWordSlice.sliceExists)
		{
			currentWordSlice = currentWordSlice.mergeWith(oldWordSlice);
		}
		currentWordSlice.calcMinScore();
		slice[startIndex] = currentWordSlice;
#ifdef SLICEVERBOSE
		slice[startIndex].debugCalcCount = oldWordSlice.debugCalcCount+1;
#endif

		int timeUntilNextScoreCheck = quitScore - currentWordSlice.minScore;
		if (endIndex <= startIndex && timeUntilNextScoreCheck < 0)
		{
			result.cellsProcessed = 1;
			return result;
		}
		if (endIndex <= startIndex && oldWordSlice.sliceExists && currentWordSlice.VP == oldWordSlice.VP && currentWordSlice.VN == oldWordSlice.VN && currentWordSlice.scoreBeforeStart == oldWordSlice.scoreBeforeStart)
		{
			result.cellsProcessed = 1;
			return result;
		}

		WordSlice previousWordSlice = currentWordSlice;
		WordSlice upPreviousWordSlice = upWordSlice;
		size_t nodeSize = params.graph.NodeEnd(i) - params.graph.NodeStart(i);
		for (LengthType w = startIndex+1; w < nodeSize; w++)
		{
			char graphChar = params.graph.NodeSequences(nodeStart+w);
			upWordSlice = oldSlice[w];
			upWordSliceExists = upExists && upWordSlice.sliceExists && upWordSlice.scoreEnd <= previousSliceQuitScore;
			bool leftUpWordSliceExists = upExists && upPreviousWordSlice.sliceExists && upPreviousWordSlice.scoreEnd <= previousSliceQuitScore;
			Word Eq = EqV.getEq(graphChar);

			oldWordSlice = slice[w];

			currentWordSlice = BV::getNextSlice(Eq, previousWordSlice, upWordSliceExists, leftUpWordSliceExists && upWordSliceExists, leftUpWordSliceExists, (j == 0 && previousBand[i]) || (j > 0 && graphChar == sequence[j-1]), upPreviousWordSlice, upWordSlice);
			if (upWordSliceExists && currentWordSlice.scoreBeforeStart > upWordSlice.scoreEnd)
			{
				auto mergable = getSourceSliceFromScore(upWordSlice.scoreEnd);
				mergable.scoreBeforeExists = true;
				currentWordSlice = currentWordSlice.mergeWithVertical(mergable);
			}

			assert(previousBand[i] || currentWordSlice.scoreBeforeStart == j || currentWordSlice.scoreBeforeStart == previousWordSlice.scoreBeforeStart + 1);
			BV::assertSliceCorrectness(currentWordSlice, upWordSlice, upWordSliceExists);

			if (currentWordSlice.scoreEnd < result.minScore)
			{
				result.minScore = currentWordSlice.scoreEnd;
				result.minScoreIndex = nodeStart + w;
				quitScore = std::min(quitScore, result.minScore + bandwidth);
			}

#ifdef EXTRACORRECTNESSASSERTIONS
			if (oldWordSlice.sliceExists)
			{
				auto debugTestSlice = currentWordSlice.mergeWith(oldWordSlice);
				assert(debugTestSlice.scoreBeforeStart == currentWordSlice.scoreBeforeStart);
				assert(debugTestSlice.VP == currentWordSlice.VP);
				assert(debugTestSlice.VN == currentWordSlice.VN);
			}
#endif

			timeUntilNextScoreCheck--;
			if (timeUntilNextScoreCheck < 0 || w == nodeSize-1)
			{
				currentWordSlice.calcMinScore();
				timeUntilNextScoreCheck = quitScore - currentWordSlice.minScore;
			}
			slice[w] = currentWordSlice;
#ifdef SLICEVERBOSE
			slice[w].debugCalcCount = oldWordSlice.debugCalcCount+1;
#endif

			if (endIndex <= w && timeUntilNextScoreCheck < 0)
			{
				result.cellsProcessed = w - startIndex + 1;
				return result;
			}
			if (endIndex <= w && oldWordSlice.sliceExists && currentWordSlice.VP == oldWordSlice.VP && currentWordSlice.VN == oldWordSlice.VN && currentWordSlice.scoreBeforeStart == oldWordSlice.scoreBeforeStart)
			{
				result.cellsProcessed = w - startIndex + 1;
				return result;
			}

#ifdef EXTRABITVECTORASSERTIONS
			auto correctslice = getWordSliceCellByCell(j, nodeStart+w, sequence, currentSlice, previousSlice, currentBand, previousBand);
			assert(currentWordSlice.scoreBeforeStart == correctslice.scoreBeforeStart);
			assert(currentWordSlice.scoreEnd == correctslice.scoreEnd);
			assert(currentWordSlice.VP == correctslice.VP);
			assert(currentWordSlice.VN == correctslice.VN);
#endif
			previousWordSlice = currentWordSlice;
			upPreviousWordSlice = upWordSlice;
		}
		result.cellsProcessed = nodeSize - startIndex;
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
		}
	}

#endif

	static ScoreType getValueOrMax(const Params& params, const DPTable& band, LengthType j, LengthType w, ScoreType max)
	{
		auto node = params.graph.IndexToNode(w);
		auto slice = (j - band.startj) / WordConfiguration<Word>::WordSize / band.samplingFrequency;
		assert(band.slices[slice].j == (j / WordConfiguration<Word>::WordSize) * WordConfiguration<Word>::WordSize);
		if (!band.slices[slice].scores.hasNode(node)) return max;
		auto word = band.slices[slice].scores.node(node)[w - params.graph.NodeStart(node)];
		auto off = j % WordConfiguration<Word>::WordSize;
		return word.getValue(off);
	}

	static ScoreType getValueOrMax(const Params& params, const DPSlice& slice, LengthType j, LengthType w, ScoreType max)
	{
		assert(j >= 0);
		assert(j < WordConfiguration<Word>::WordSize);
		auto node = params.graph.IndexToNode(w);
		if (!slice.scores.hasNode(node)) return max;
		auto word = slice.scores.node(node)[w - params.graph.NodeStart(node)];
		auto off = j % WordConfiguration<Word>::WordSize;
		return word.getValue(off);
	}

	static ScoreType getValue(const Params& params, const DPSlice& slice, LengthType j, LengthType w)
	{
		assert(j >= 0);
		assert(j < WordConfiguration<Word>::WordSize);
		auto node = params.graph.IndexToNode(w);
		auto word = slice.scores.node(node)[w - params.graph.NodeStart(node)];
		auto off = j % WordConfiguration<Word>::WordSize;
		return word.getValue(off);
	}

	static ScoreType getValue(const Params& params, const DPTable& band, LengthType j, LengthType w)
	{
		auto node = params.graph.IndexToNode(w);
		auto slice = (j - band.startj) / WordConfiguration<Word>::WordSize / band.samplingFrequency;
		assert(band.slices[slice].j == (j / WordConfiguration<Word>::WordSize) * WordConfiguration<Word>::WordSize);
		auto word = band.slices[slice].scores.node(node)[w - params.graph.NodeStart(node)];
		auto off = j % WordConfiguration<Word>::WordSize;
		return word.getValue(off);
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

	NodeCalculationResult calculateSlice(const std::string& sequence, size_t j, NodeSlice<WordSlice>& currentSlice, const NodeSlice<WordSlice>& previousSlice, const std::vector<LengthType>& previousNodes, std::vector<bool>& currentBand, const std::vector<bool>& previousBand, std::vector<size_t>& partOfComponent, ScoreType previousQuitScore, int bandwidth, ScoreType previousMinScore) const
	{
		ScoreType currentMinimumScore = std::numeric_limits<ScoreType>::max();
		LengthType currentMinimumIndex;
		size_t cellsProcessed = 0;

		EqVector EqV = BV::getEqVector(sequence, j);

		for (auto node : previousNodes)
		{
			if (previousSlice.minScore(node) <= previousQuitScore)
			{
				calculableQueue.insert(previousSlice.minScore(node) - previousMinScore, NodeWithPriority { node, previousSlice.startIndex(node), previousSlice.endIndex(node), previousSlice.minScore(node) });
			}
		}
		assert(calculableQueue.size() != 0);
		
		ScoreType currentMinScoreAtEndRow = std::numeric_limits<ScoreType>::max() - bandwidth - 1;
		while (calculableQueue.size() > 0)
		{
			auto pair = calculableQueue.top();
			if (pair.priority > currentMinScoreAtEndRow + bandwidth) break;
			auto i = pair.node;
			size_t offset = pair.offset;
			size_t endOffset = pair.endOffset;
			if (!currentSlice.hasNode(i))
			{
				assert(!currentBand[i]);
				currentSlice.addNode(i, params.graph.NodeLength(i), WordSlice {0, 0, std::numeric_limits<ScoreType>::max(), std::numeric_limits<ScoreType>::max(), 0 });
				currentSlice.setMinScore(i, std::numeric_limits<ScoreType>::max());
				currentSlice.setStartIndex(i, offset);
				currentBand[i] = true;
			}
			assert(currentBand[i]);
			calculableQueue.pop();
			auto oldEnd = currentSlice.node(i).back();
#ifdef EXTRACORRECTNESSASSERTIONS
			std::vector<WordSlice> debugOldNode;
			auto debugNode = currentSlice.node(i);
			for (size_t ii = 0; ii < debugNode.size(); ii++)
			{
				debugOldNode.push_back(debugNode[ii]);
			}
#endif
			auto nodeCalc = calculateNode(i, j, offset, endOffset, sequence, EqV, currentSlice, previousSlice, currentBand, previousBand, previousQuitScore, currentMinScoreAtEndRow + bandwidth, bandwidth);
			currentMinScoreAtEndRow = std::min(currentMinScoreAtEndRow, nodeCalc.minScore);
			currentSlice.setMinScoreIfSmaller(i, nodeCalc.minScore);
			auto newEnd = currentSlice.node(i).back();
#ifdef EXTRACORRECTNESSASSERTIONS
			auto debugNewNode = currentSlice.node(i);
			for (size_t debugi = 0; debugi < debugOldNode.size(); debugi++)
			{
				assertBitvectorConfirmedAreConsistent(debugNewNode[debugi], debugOldNode[debugi], currentMinScoreAtEndRow + bandwidth);
			}
#endif
			if (newEnd.scoreBeforeStart != oldEnd.scoreBeforeStart || newEnd.VP != oldEnd.VP || newEnd.VN != oldEnd.VN)
			{
				ScoreType newEndMinScore = newEnd.changedMinScore(oldEnd);
				if (newEndMinScore <= currentMinScoreAtEndRow + bandwidth)
				{
					for (auto neighbor : params.graph.outNeighbors[i])
					{
						calculableQueue.insert(newEndMinScore - previousMinScore, NodeWithPriority { neighbor, 0, 0, newEndMinScore });
					}
				}
			}
#ifndef NDEBUG
			auto debugslice = currentSlice.node(i);
			if (nodeCalc.minScore != std::numeric_limits<ScoreType>::max() && nodeCalc.minScore <= currentMinScoreAtEndRow + bandwidth)
			{
				assert(nodeCalc.minScoreIndex >= params.graph.NodeStart(i));
				assert(nodeCalc.minScoreIndex < params.graph.NodeEnd(i));
				assert(debugslice[nodeCalc.minScoreIndex - params.graph.NodeStart(i)].scoreEnd == nodeCalc.minScore);
			}
#endif
			if (nodeCalc.minScore < currentMinimumScore)
			{
				currentMinimumScore = nodeCalc.minScore;
				currentMinimumIndex = nodeCalc.minScoreIndex;
			}
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

		finalizeSlice(currentSlice, result, bandwidth);

		calculableQueue.clear();

		return result;
	}

	void flattenLastSliceEnd(NodeSlice<WordSlice>& slice, NodeCalculationResult& sliceCalc, LengthType j, size_t sequenceSize) const
	{
		assert(j < sequenceSize);
		assert(sequenceSize - j < WordConfiguration<Word>::WordSize);
		sliceCalc.minScore = std::numeric_limits<ScoreType>::max();
		sliceCalc.minScoreIndex = -1;
		for (auto node : slice)
		{
			for (size_t i = 0; i < node.second.size(); i++)
			{
				auto wordSliceResult = BV::flattenWordSlice(node.second[i], sequenceSize - j);
				node.second[i] = wordSliceResult;
				if (wordSliceResult.scoreEnd < sliceCalc.minScore)
				{
					sliceCalc.minScore = wordSliceResult.scoreEnd;
					sliceCalc.minScoreIndex = params.graph.NodeStart(node.first)+i;
				}
			}
		}
	}

	void finalizeSlice(NodeSlice<WordSlice>& slice, NodeCalculationResult sliceCalc, int bandwidth) const
	{
#ifdef SLICEVERBOSE
		size_t uselessCells = 0;
		size_t doubleCalcs = 0;
#endif
		ScoreType uninitScore = std::numeric_limits<ScoreType>::max();
		ScoreType minScore = sliceCalc.minScore;
		assert(minScore < uninitScore);
		assert(minScore < uninitScore - 2 * bandwidth - 1);

		for (auto node : slice)
		{
			size_t startOffset = -1;
			size_t endOffset = -1;
#ifndef NDEBUG
			ScoreType debugNodeMinScore = std::numeric_limits<ScoreType>::max();
#endif
			ScoreType cellMinScore = minScore+bandwidth+1;
			for (size_t i = 0; i < node.second.size(); i++)
			{
				WordSlice& cell = node.second[i];
#ifdef SLICEVERBOSE
				if (cell.debugCalcCount > 1)
				{
					doubleCalcs += cell.debugCalcCount - 1;
				}
#endif
#ifndef NDEBUG
				debugNodeMinScore = std::min(debugNodeMinScore, cell.scoreEnd);
#endif
				if (cell.scoreBeforeStart != uninitScore)
				{
					if (cellMinScore < minScore+bandwidth)
					{
						cellMinScore++;
					}
					else
					{
						cell.calcMinScore();
						cellMinScore = cell.minScore;
					}
				}
				if (cell.scoreBeforeStart == uninitScore || cellMinScore > minScore+bandwidth)
				{
#ifdef SLICEVERBOSE
					if (cell.scoreEnd != uninitScore) uselessCells++;
#endif
					cell.scoreEnd = minScore + 2 * bandwidth + 1;
					cell.scoreBeforeStart = minScore + 2 * bandwidth + 1;
					cell.VP = 0;
					cell.VN = 0;
					cell.sliceExists = false;
				}
				else if (cell.scoreEnd <= minScore + bandwidth)
				{
					if (startOffset == -1) startOffset = i;
					endOffset = i;
				}
			}
			assert(startOffset != -1 || debugNodeMinScore > minScore + bandwidth);
			slice.setStartIndex(node.first, startOffset);
			slice.setEndIndex(node.first, endOffset);
		}
#ifdef SLICEVERBOSE
		std::cerr << "useless cells " << uselessCells << " ";
		std::cerr << "doublecounts " << doubleCalcs << " ";
#endif
	}

	void fillDPSlice(const std::string& sequence, DPSlice& slice, const DPSlice& previousSlice, const std::vector<bool>& previousBand, std::vector<size_t>& partOfComponent, std::vector<bool>& currentBand, int bandwidth) const
	{
		auto sliceResult = calculateSlice(sequence, slice.j, slice.scores, previousSlice.scores, previousSlice.nodes, currentBand, previousBand, partOfComponent, previousSlice.minScore + previousSlice.bandwidth, bandwidth, previousSlice.minScore);
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

	DPSlice pickMethodAndExtendFill(const std::string& sequence, const DPSlice& previous, const std::vector<bool>& previousBand, std::vector<bool>& currentBand, std::vector<size_t>& partOfComponent, std::vector<bool>& processed, int bandwidth) const
	{
		DPSlice bandTest { &nodesliceMap };
		bandTest.j = previous.j + WordConfiguration<Word>::WordSize;
		bandTest.correctness = previous.correctness;
		bandTest.scores.reserve(previous.numCells);

		fillDPSlice(sequence, bandTest, previous, previousBand, partOfComponent, currentBand, bandwidth);

#ifdef EXTRACORRECTNESSASSERTIONS
		if (bandTest.cellsProcessed <= params.maxCellsPerSlice && sequence.size() >= bandTest.j + WordConfiguration<Word>::WordSize) verifySliceBitvector(sequence, bandTest, previous);
#endif
		return bandTest;
	}

	void removeWronglyAlignedEnd(DPTable& table) const
	{
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

	DPTable getSqrtSlices(const std::string& sequence, const DPSlice& initialSlice, size_t numSlices, size_t samplingFrequency) const
	{
		assert(initialSlice.j == -WordConfiguration<Word>::WordSize);
		assert(initialSlice.j + numSlices * WordConfiguration<Word>::WordSize <= sequence.size() + WordConfiguration<Word>::WordSize);
		DPTable result;
		size_t realCells = 0;
		size_t cellsProcessed = 0;
		result.samplingFrequency = samplingFrequency;
		std::vector<bool> previousBand;
		std::vector<bool> currentBand;
		std::vector<size_t> partOfComponent;
		previousBand.resize(params.graph.NodeSize(), false);
		currentBand.resize(params.graph.NodeSize(), false);
		partOfComponent.resize(params.graph.NodeSize(), std::numeric_limits<size_t>::max());
		{
			auto initialOrder = initialSlice.nodes;
			for (auto node : initialOrder)
			{
				previousBand[node] = true;
			}
		}
#ifndef NDEBUG
		debugLastRowMinScore = 0;
#endif
		DPSlice lastSlice = initialSlice.getFrozenSqrtEndScores();
		DPSlice storeSlice = lastSlice;
		assert(lastSlice.correctness.CurrentlyCorrect());
		DPSlice rampSlice = lastSlice;
		std::vector<bool> processed;
		processed.resize(params.graph.SizeInBp(), false);
		size_t rampRedoIndex = -1;
		size_t rampUntil = 0;
		DPSlice backtraceOverridePreslice = lastSlice;
		std::vector<DPSlice> backtraceOverrideTemps;
		bool backtraceOverriding = false;
#ifndef NDEBUG
		volatile size_t debugLastProcessedSlice;
#endif
		for (size_t slice = 0; slice < numSlices; slice++)
		{
			int bandwidth = (params.rampBandwidth > params.initialBandwidth && rampUntil >= slice) ? params.rampBandwidth : params.initialBandwidth;
			size_t storeSliceIndex = slice / samplingFrequency + 1;
#ifndef NDEBUG
			debugLastProcessedSlice = slice;
			debugLastRowMinScore = lastSlice.minScore;
#endif
			auto timeStart = std::chrono::system_clock::now();
			auto newSlice = pickMethodAndExtendFill(sequence, lastSlice, previousBand, currentBand, partOfComponent, processed, bandwidth);
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

			if (!newSlice.correctness.CorrectFromCorrect())
			{
				newSlice.scores.clearVectorMap();
#ifndef NDEBUG
				debugLastProcessedSlice = slice-1;
#endif
				break;
			}
			if (!newSlice.correctness.CurrentlyCorrect() && rampUntil < slice && params.rampBandwidth > params.initialBandwidth)
			{
				for (auto node : newSlice.nodes)
				{
					assert(currentBand[node]);
					currentBand[node] = false;
				}
				for (auto node : lastSlice.nodes)
				{
					assert(previousBand[node]);
					previousBand[node] = false;
				}
				newSlice.scores.clearVectorMap();
				rampUntil = slice;
				std::swap(slice, rampRedoIndex);
				std::swap(lastSlice, rampSlice);
				for (auto node : lastSlice.nodes)
				{
					assert(!previousBand[node]);
					previousBand[node] = true;
				}
				while (result.bandwidthPerSlice.size() > slice+1) result.bandwidthPerSlice.pop_back();
				while (result.correctness.size() > slice+1) result.correctness.pop_back();
				while (result.slices.size() > 1 && result.slices.back().j > slice * WordConfiguration<Word>::WordSize) result.slices.pop_back();
				storeSlice = lastSlice;
#ifdef SLICEVERBOSE
				std::cerr << " ramp to " << slice;
#endif
				if (backtraceOverriding)
				{
#ifdef SLICEVERBOSE
					std::cerr << " preslicej " << backtraceOverridePreslice.j << " lastslicej " << lastSlice.j;
#endif
					if (backtraceOverridePreslice.j > lastSlice.j)
					{
#ifdef SLICEVERBOSE
						std::cerr << " empty backtrace override";
#endif
						backtraceOverriding = false;
						//empty memory
						{
							decltype(backtraceOverrideTemps) tmp;
							std::swap(backtraceOverrideTemps, tmp);
						}
					}
					else
					{
#ifdef SLICEVERBOSE
						std::cerr << " shorten backtrace override";
#endif
						while (backtraceOverrideTemps.size() > 0 && backtraceOverrideTemps.back().j > lastSlice.j)
						{
							DPSlice tmp;
							std::swap(backtraceOverrideTemps.back(), tmp);
							backtraceOverrideTemps.pop_back();
						}
#ifdef SLICEVERBOSE
						std::cerr << " to " << backtraceOverrideTemps.size() << " temps";
#endif
					}
				}
				while (result.backtraceOverrides.size() > 0 && result.backtraceOverrides.back().endj > lastSlice.j)
				{
					result.backtraceOverrides.pop_back();
				}
#ifdef SLICEVERBOSE
				std::cerr << std::endl;
				std::cerr << "bandwidthPerSlice.size() " << result.bandwidthPerSlice.size();
				if (result.slices.size() > 0) std::cerr << " slices.back().j " << result.slices.back().j; else std::cerr << " slices.size() 0";
				if (result.backtraceOverrides.size() > 0) std::cerr << " backtraceOverrides.back().endj " << result.backtraceOverrides.back().endj; else std::cerr << " backtraceOverrides.size() 0";
				std::cerr << std::endl;
#endif
				continue;
			}

			if (!backtraceOverriding && newSlice.cellsProcessed >= params.BacktraceOverrideCutoff)
			{
#ifdef SLICEVERBOSE
				std::cerr << " start backtrace override";
#endif
				assert(!lastSlice.cellsProcessed < params.BacktraceOverrideCutoff);
				backtraceOverridePreslice = lastSlice;
				backtraceOverriding = true;
				backtraceOverrideTemps.push_back(newSlice.getFrozenScores());
			}
			else if (backtraceOverriding)
			{
				if (newSlice.cellsProcessed < params.BacktraceOverrideCutoff)
				{
#ifdef SLICEVERBOSE
					std::cerr << " end backtrace override";
#endif
					assert(lastSlice.j == backtraceOverrideTemps.back().j);
					assert(backtraceOverrideTemps.size() > 0);
					result.backtraceOverrides.emplace_back(params, sequence, backtraceOverridePreslice, backtraceOverrideTemps);
					backtraceOverriding = false;
					while (result.slices.size() > 0 && result.slices.back().j >= result.backtraceOverrides.back().startj && result.slices.back().j <= result.backtraceOverrides.back().endj)
					{
						result.slices.pop_back();
					}
					result.slices.push_back(lastSlice);
#ifdef SLICEVERBOSE
					std::cerr << " push slice j " << lastSlice.j;
#endif
					storeSlice = newSlice.getFrozenSqrtEndScores();
					//empty memory
					{
						decltype(backtraceOverrideTemps) tmp;
						std::swap(backtraceOverrideTemps, tmp);
					}
				}
				else
				{
#ifdef SLICEVERBOSE
					std::cerr << " continue backtrace override";
#endif
					backtraceOverrideTemps.push_back(newSlice.getFrozenScores());
				}
			}
#ifdef SLICEVERBOSE
			std::cerr << std::endl;
#endif

			assert(result.bandwidthPerSlice.size() == slice);
			result.bandwidthPerSlice.push_back(bandwidth);
			result.correctness.push_back(newSlice.correctness);
			if (slice % samplingFrequency == 0)
			{
				if (result.slices.size() == 0 || storeSlice.j != result.slices.back().j)
				{
					assert(result.slices.size() == 0 || result.slices.back().j == -WordConfiguration<Word>::WordSize || storeSlice.j > result.slices.back().j);
					result.slices.push_back(storeSlice);
#ifdef SLICEVERBOSE
					std::cerr << " push slice j " << storeSlice.j;
#endif
					storeSlice = newSlice.getFrozenSqrtEndScores();
				}
			}
			if (newSlice.EstimatedMemoryUsage() < storeSlice.EstimatedMemoryUsage())
			{
				storeSlice = newSlice.getFrozenSqrtEndScores();
			}
			for (auto node : lastSlice.nodes)
			{
				assert(previousBand[node]);
				previousBand[node] = false;
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
			lastSlice = newSlice.getFrozenSqrtEndScores();
			newSlice.scores.clearVectorMap();
			std::swap(previousBand, currentBand);
		}

		if (backtraceOverriding)
		{
			assert(backtraceOverrideTemps.size() > 0);
			assert(lastSlice.j == backtraceOverrideTemps.back().j);
			result.backtraceOverrides.emplace_back(params, sequence, backtraceOverridePreslice, backtraceOverrideTemps);
			backtraceOverriding = false;
			//empty memory
			{
				decltype(backtraceOverrideTemps) tmp;
				std::swap(backtraceOverrideTemps, tmp);
			}
			while (result.slices.size() > 0 && result.slices.back().j >= result.backtraceOverrides.back().startj && result.slices.back().j <= result.backtraceOverrides.back().endj)
			{
				result.slices.pop_back();
			}
		}
#ifndef NDEBUF
		volatile
#endif
		size_t lastExisting = 0;
		// assert(debugLastProcessedSlice == -1 || lastExisting == debugLastProcessedSlice / samplingFrequency || lastExisting == debugLastProcessedSlice / samplingFrequency + 1);
		// storeSlices.erase(storeSlices.begin() + lastExisting + 1, storeSlices.end());
		// result.slices = storeSlices;
		assert(result.bandwidthPerSlice.size() == debugLastProcessedSlice + 1);
#ifndef NDEBUG
		assert(result.slices.size() > 0);
		for (size_t i = 0; i < result.slices.size(); i++)
		{
			// assert(i == 0 || result.slices[i].j / WordConfiguration<Word>::WordSize / samplingFrequency == i-1);
			assert(i <= 1 || result.slices[i].j > result.slices[i-1].j);
			// assert(i != 1 || result.slices[i].j >= 0);
		}
		for (size_t i = 1; i < result.slices.size(); i++)
		{
			assert(result.slices[i].minScore >= result.slices[i-1].minScore);
		}
		for (size_t i = 0; i < result.backtraceOverrides.size(); i++)
		{
			assert(result.backtraceOverrides[i].endj >= result.backtraceOverrides[i].startj);
		}
		for (size_t i = 1; i < result.backtraceOverrides.size(); i++)
		{
			assert(result.backtraceOverrides[i].startj > result.backtraceOverrides[i-1].endj);
		}
#endif
		return result;
	}

	std::vector<DPSlice> getSlicesFromTable(const std::string& sequence, LengthType overrideLastJ, const DPTable& table, size_t startIndex) const
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
		std::vector<bool> previousBand;
		std::vector<bool> currentBand;
		std::vector<size_t> partOfComponent;
		previousBand.resize(params.graph.NodeSize(), false);
		currentBand.resize(params.graph.NodeSize(), false);
		partOfComponent.resize(params.graph.NodeSize(), std::numeric_limits<size_t>::max());
		{
			auto initialOrder = initialSlice.nodes;
			for (auto node : initialOrder)
			{
				previousBand[node] = true;
			}
		}
#ifndef NDEBUG
		debugLastRowMinScore = 0;
#endif
		DPSlice lastSlice = initialSlice.getFrozenSqrtEndScores();
		// assert(lastSlice.correctness.CurrentlyCorrect());
		DPSlice rampSlice = lastSlice;
		std::vector<bool> processed;
		processed.resize(params.graph.SizeInBp(), false);
		for (size_t slice = startSlice; slice < endSlice; slice++)
		{
			int bandwidth = table.bandwidthPerSlice[slice];
#ifndef NDEBUG
			debugLastRowMinScore = lastSlice.minScore;
#endif
			auto newSlice = pickMethodAndExtendFill(sequence, lastSlice, previousBand, currentBand, partOfComponent, processed, bandwidth);
			assert(result.size() == 0 || newSlice.j == result.back().j + WordConfiguration<Word>::WordSize);

			size_t sliceCells = 0;
			for (auto node : newSlice.nodes)
			{
				sliceCells += params.graph.NodeEnd(node) - params.graph.NodeStart(node);
			}
			realCells += sliceCells;
			cellsProcessed += newSlice.cellsProcessed;

			// assert(slice == endSlice-1 || newSlice.correctness.CurrentlyCorrect());
			result.push_back(newSlice.getFrozenScores());
			for (auto node : lastSlice.nodes)
			{
				assert(previousBand[node]);
				previousBand[node] = false;
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
			lastSlice = newSlice.getFrozenSqrtEndScores();
			newSlice.scores.clearVectorMap();
			std::swap(previousBand, currentBand);
		}
#ifndef NDEBUG
		for (size_t i = 1; i < result.size(); i++)
		{
			assert(result[i].minScore >= result[i-1].minScore);
		}
#endif
		return result;
	}

	DPSlice getInitialSliceOnlyOneNode(LengthType nodeIndex) const
	{
		DPSlice result;
		result.j = -WordConfiguration<Word>::WordSize;
		result.bandwidth = 1;
		result.minScore = 0;
		result.scores.addNode(nodeIndex, params.graph.NodeLength(nodeIndex));
		result.scores.setMinScore(nodeIndex, 0);
		result.scores.setStartIndex(nodeIndex, 0);
		result.scores.setEndIndex(nodeIndex, params.graph.NodeLength(nodeIndex));
		result.minScoreIndex = params.graph.NodeEnd(nodeIndex) - 1;
		result.nodes.push_back(nodeIndex);
		auto slice = result.scores.node(nodeIndex);
		for (size_t i = 0; i < slice.size(); i++)
		{
			slice[i] = {0, 0, 0, 0, false};
			slice[i].sliceExists = true;
		}
		return result.getFrozenSqrtEndScores();
	}

	DPSlice getInitialSliceOneNodeGroup(const std::vector<LengthType>& nodeIndices) const
	{
		DPSlice result;
		result.j = -WordConfiguration<Word>::WordSize;
		result.bandwidth = 1;
		result.minScore = 0;
		for (auto nodeIndex : nodeIndices)
		{
			result.scores.addNode(nodeIndex, params.graph.NodeLength(nodeIndex));
			result.scores.setMinScore(nodeIndex, 0);
			result.scores.setStartIndex(nodeIndex, 0);
			result.scores.setEndIndex(nodeIndex, params.graph.NodeLength(nodeIndex));
			result.minScoreIndex = params.graph.NodeEnd(nodeIndex) - 1;
			result.nodes.push_back(nodeIndex);
			auto slice = result.scores.node(nodeIndex);
			for (size_t i = 0; i < slice.size(); i++)
			{
				slice[i] = {0, 0, 0, 0, false};
				slice[i].sliceExists = true;
			}
		}
		return result.getFrozenSqrtEndScores();
	}

	int getSamplingFrequency(size_t sequenceLen) const
	{
		size_t samplingFrequency = 1;
		samplingFrequency = (int)(sqrt(sequenceLen / WordConfiguration<Word>::WordSize));
		if (samplingFrequency <= 1) samplingFrequency = 2;
		return samplingFrequency;
	}

	OnewayTrace getBacktraceFullStart(std::string sequence) const
	{
		int padding = (WordConfiguration<Word>::WordSize - (sequence.size() % WordConfiguration<Word>::WordSize)) % WordConfiguration<Word>::WordSize;
		for (int i = 0; i < padding; i++)
		{
			sequence += 'N';
		}
		DPSlice startSlice;
		for (size_t i = 0; i < params.graph.nodeStart.size(); i++)
		{
			startSlice.scores.addNode(i, params.graph.NodeEnd(i) - params.graph.NodeStart(i), WordSlice {0, 0, 0, 0, WordConfiguration<Word>::WordSize, false});
			startSlice.scores.setMinScore(i, 0);
			startSlice.j = -WordConfiguration<Word>::WordSize;
			startSlice.nodes.push_back(i);
		}
		size_t samplingFrequency = getSamplingFrequency(sequence.size());
		auto slice = getSqrtSlices(sequence, startSlice, sequence.size() / WordConfiguration<Word>::WordSize, samplingFrequency);
		removeWronglyAlignedEnd(slice);
		// std::cerr << "score: " << slice.slices.back().minScore << std::endl;

		auto result = getTraceFromTable(sequence, slice);
		while (result.trace.back().second >= sequence.size() - padding)
		{
			result.trace.pop_back();
		}
		assert(result.trace[0].second == 0);
		return result;
	}

};
*/
#endif
