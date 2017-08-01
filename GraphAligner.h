#ifndef GraphAligner_H
#define GraphAligner_H

#include <chrono>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <boost/config.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/container/flat_set.hpp>
#include <unordered_set>
#include <queue>
#include "SliceRow.h"
#include "SparseBoolMatrix.h"
#include "AlignmentGraph.h"
#include "vg.pb.h"
#include "ThreadReadAssertion.h"

using namespace boost;

void printtime(const char* msg)
{
	static auto time = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
	auto newtime = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
	std::cout << msg << " " << newtime << " (" << (newtime - time) << ")" << std::endl;
	time = newtime;
}

template <typename Word>
class WordConfiguration
{
};

template <>
class WordConfiguration<uint64_t>
{
public:
	static constexpr int WordSize = 64;
	//number of bits per chunk
	//prefix sum differences are calculated in chunks of log w bits
	static constexpr int ChunkBits = 8;
	static constexpr uint64_t AllZeros = 0x0000000000000000;
	static constexpr uint64_t AllOnes = 0xFFFFFFFFFFFFFFFF;
	//positions of the sign bits for each chunk
	static constexpr uint64_t SignMask = 0x8080808080808080;
	//constant for multiplying the chunk popcounts into prefix sums
	//this should be 1 at the start of each chunk
	static constexpr uint64_t PrefixSumMultiplierConstant = 0x0101010101010101;
	//positions of the least significant bits for each chunk
	static constexpr uint64_t LSBMask = 0x0101010101010101;

	static int popcount(uint64_t x)
	{
		//https://en.wikipedia.org/wiki/Hamming_weight
		x -= (x >> 1) & 0x5555555555555555;
		x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
		x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f;
		return (x * 0x0101010101010101) >> 56;
	}

	static uint64_t ChunkPopcounts(uint64_t value)
	{
		uint64_t x = value;
		x -= (x >> 1) & 0x5555555555555555;
		x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
		x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f;
		return x;
	}
};

template <typename LengthType, typename ScoreType, typename Word>
class GraphAligner
{
private:
	typedef std::pair<LengthType, LengthType> MatrixPosition;
	class MatrixSlice
	{
	public:
		std::vector<ScoreType> minScorePerWordSlice;
		LengthType finalMinScoreColumn;
		ScoreType finalMinScore;
		size_t cellsProcessed;
	};
	class WordSlice
	{
	public:
		WordSlice() :
		VP(WordConfiguration<Word>::AllZeros),
		VN(WordConfiguration<Word>::AllZeros),
		scoreEnd(0),
		scoreBeforeStart(0)
		{}
		WordSlice(Word VP, Word VN, ScoreType scoreEnd, ScoreType scoreBeforeStart) :
		VP(VP),
		VN(VN),
		scoreEnd(scoreEnd),
		scoreBeforeStart(scoreBeforeStart)
		{}
		Word VP;
		Word VN;
		ScoreType scoreEnd;
		ScoreType scoreBeforeStart;
	};
public:
	class AlignmentResult
	{
	public:
		AlignmentResult()
		{
		}
		AlignmentResult(vg::Alignment alignment, bool alignmentFailed, size_t cellsProcessed) :
		alignment(alignment),
		alignmentFailed(alignmentFailed),
		cellsProcessed(cellsProcessed)
		{
		}
		vg::Alignment alignment;
		bool alignmentFailed;
		size_t cellsProcessed;
	};

	GraphAligner(const AlignmentGraph& graph) :
	graph(graph)
	{
	}
	
	AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence, int dynamicWidth, LengthType dynamicRowStart) const
	{
		assert(graph.finalized);
		auto band = getFullBand(sequence.size(), dynamicRowStart);
		auto trace = getBacktrace(sequence, dynamicWidth, dynamicRowStart, band);
		//failed alignment, don't output
		if (std::get<0>(trace) == std::numeric_limits<ScoreType>::min()) return emptyAlignment();
		auto result = traceToAlignment(seq_id, sequence, std::get<0>(trace), std::get<1>(trace), std::get<2>(trace));
		return result;
	}

	AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence, int dynamicWidth, LengthType dynamicRowStart, const std::vector<AlignmentGraph::SeedHit>& seedHits, int startBandwidth) const
	{
		assert(graph.finalized);
		auto band = getSeededStartBand(seedHits, dynamicRowStart, startBandwidth, sequence);
		auto trace = getBacktrace(sequence, dynamicWidth, dynamicRowStart, band);
		//failed alignment, don't output
		if (std::get<0>(trace) == std::numeric_limits<ScoreType>::min()) return emptyAlignment();
		auto result = traceToAlignment(seq_id, sequence, std::get<0>(trace), std::get<1>(trace), std::get<2>(trace));
		return result;
	}

private:

	class ExpandoCell
	{
	public:
		ExpandoCell(LengthType w, LengthType j, size_t bt) :
		position(w, j),
		backtraceIndex(bt)
		{}
		MatrixPosition position;
		size_t backtraceIndex;
	};

	std::vector<MatrixPosition> backtrace(MatrixPosition endPosition, const std::string& sequence, const std::vector<ScoreType>& minScorePerWordSlice) const
	{
		//first try the faster way
		auto result = backtraceInner(endPosition, sequence, minScorePerWordSlice, false);
		//if that breaks then do the full backtrace
		if (result.size() == 0) result = backtraceInner(endPosition, sequence, minScorePerWordSlice, true);
		assert(result.size() != 0);
		return result;
	}

	std::vector<MatrixPosition> backtraceInner(MatrixPosition endPosition, const std::string& sequence, const std::vector<ScoreType>& minScorePerWordSlice, bool fullBacktrace) const
	{
		ScoreType scoreAtEnd = minScorePerWordSlice.back();
		ScoreType currentDistance = 0;
		std::vector<ExpandoCell> visitedExpandos;
		std::vector<ExpandoCell> currentDistanceQueue;
		std::vector<ExpandoCell> currentDistancePlusOneQueue;
		currentDistanceQueue.emplace_back(endPosition.first, endPosition.second, 0);
		SparseBoolMatrix<SliceRow<LengthType>> visitedCells {graph.nodeSequences.size(), sequence.size()+1};
		while (true)
		{
			if (currentDistanceQueue.size() == 0)
			{
				if (currentDistancePlusOneQueue.size() == 0)
				{
					return std::vector<MatrixPosition> {};
				}
				assert(currentDistancePlusOneQueue.size() > 0);
				std::swap(currentDistanceQueue, currentDistancePlusOneQueue);
				currentDistance++;
				if (currentDistance > scoreAtEnd)
				{
					return std::vector<MatrixPosition> {};
				}
			}
			auto current = currentDistanceQueue.back();
			currentDistanceQueue.pop_back();
			auto w = current.position.first;
			auto j = current.position.second;
			if (j == 0)
			{
				visitedExpandos.push_back(current);
				break;
			}
			auto sliceIndex = (j-1) / WordConfiguration<Word>::WordSize;
			assert(sliceIndex < minScorePerWordSlice.size());
			ScoreType maxDistanceHere;
			if (fullBacktrace)
			{
				maxDistanceHere = scoreAtEnd;
			}
			else
			{
				maxDistanceHere = scoreAtEnd - minScorePerWordSlice[sliceIndex];
			}
			if (currentDistance > maxDistanceHere) continue;
			if (visitedCells.get(w, j)) continue;
			visitedCells.set(w, j);
			visitedExpandos.push_back(current);
			auto nodeIndex = graph.indexToNode[w];
			auto backtraceIndexToCurrent = visitedExpandos.size()-1;
			currentDistancePlusOneQueue.emplace_back(w, j-1, backtraceIndexToCurrent);
			if (w == graph.nodeStart[nodeIndex])
			{
				for (size_t i = 0; i < graph.inNeighbors[nodeIndex].size(); i++)
				{
					auto u = graph.nodeEnd[graph.inNeighbors[nodeIndex][i]]-1;
					currentDistancePlusOneQueue.emplace_back(u, j, backtraceIndexToCurrent);
					if (sequence[j-1] == 'N' || graph.nodeSequences[w] == sequence[j-1])
					{
						currentDistanceQueue.emplace_back(u, j-1, backtraceIndexToCurrent);
					}
					else
					{
						currentDistancePlusOneQueue.emplace_back(u, j-1, backtraceIndexToCurrent);
					}
				}
			}
			else
			{
				auto u = w-1;
				currentDistancePlusOneQueue.emplace_back(u, j, backtraceIndexToCurrent);
				if (sequence[j-1] == 'N' || graph.nodeSequences[w] == sequence[j-1])
				{
					currentDistanceQueue.emplace_back(u, j-1, backtraceIndexToCurrent);
				}
				else
				{
					currentDistancePlusOneQueue.emplace_back(u, j-1, backtraceIndexToCurrent);
				}
			}
		}
		std::cerr << "backtrace visited " << visitedCells.totalOnes() << " cells" << std::endl;
		assert(currentDistance <= scoreAtEnd);
		auto index = visitedExpandos.size()-1;
		std::vector<MatrixPosition> result;
		while (index > 0)
		{
			result.push_back(visitedExpandos[index].position);
			assert(visitedExpandos[index].backtraceIndex < index);
			index = visitedExpandos[index].backtraceIndex;
		}
		return result;
	}

	std::vector<std::vector<bool>> getFullBand(size_t sequenceSize, LengthType dynamicRowStart) const
	{
		std::vector<std::vector<bool>> result;
		result.resize(dynamicRowStart/WordConfiguration<Word>::WordSize);
		for (size_t i = 0; i < dynamicRowStart/WordConfiguration<Word>::WordSize; i++)
		{
			result[i].resize(graph.nodeStart.size(), true);
		}
		return result;
	}

	std::vector<std::vector<bool>> getSeededStartBand(const std::vector<AlignmentGraph::SeedHit>& originalSeedHits, int dynamicRowStart, int startBandwidth, const std::string& sequence) const
	{
		auto seedHitsInMatrix = graph.GetSeedHitPositionsInMatrix(sequence, originalSeedHits);
		auto bandLocations = getBandLocations(sequence.size(), seedHitsInMatrix, dynamicRowStart);
		std::vector<std::vector<bool>> result;
		result.resize(dynamicRowStart/WordConfiguration<Word>::WordSize);
		for (LengthType j = 0; j < dynamicRowStart && j < sequence.size()+1; j += WordConfiguration<Word>::WordSize)
		{
			auto index = j/WordConfiguration<Word>::WordSize;
			result[index].resize(graph.nodeStart.size(), false);
			for (size_t i = 0; i < bandLocations[j].size(); i++)
			{
				expandBandDynamically(result[index], bandLocations[j][i], startBandwidth);
			}
		}
		return result;
	}

	std::vector<std::vector<LengthType>> getBandLocations(int sequenceLength, const std::vector<MatrixPosition>& seedHits, LengthType maxRow) const
	{
		if (maxRow > sequenceLength+1) maxRow = sequenceLength+1;
		std::vector<std::vector<LengthType>> forwardResult;
		std::vector<std::vector<LengthType>> backwardResult;
		backwardResult.resize(sequenceLength+1);
		forwardResult.resize(sequenceLength+1);
		backwardResult[0].emplace_back(0);
		forwardResult[0].emplace_back(0);
		for (auto hit : seedHits)
		{
			if (hit.second < maxRow) expandBandForwards(forwardResult, hit.first, hit.second, sequenceLength, maxRow);
			expandBandBackwards(backwardResult, hit.first, hit.second, sequenceLength);
		}
		std::vector<std::vector<LengthType>> result;
		result.resize(maxRow+1);
		for (size_t j = 0; j < maxRow; j++)
		{
			std::set<LengthType> rowResult;
			rowResult.insert(forwardResult[j].begin(), forwardResult[j].end());
			rowResult.insert(backwardResult[j].begin(), backwardResult[j].end());
			result[j].insert(result[j].end(), rowResult.begin(), rowResult.end());
		}
		return result;
	}

	AlignmentResult emptyAlignment() const
	{
		vg::Alignment result;
		result.set_score(std::numeric_limits<decltype(result.score())>::min());
		return AlignmentResult { result, true, 0 };
	}

	AlignmentResult traceToAlignment(const std::string& seq_id, const std::string& sequence, ScoreType score, const std::vector<MatrixPosition>& trace, size_t cellsProcessed) const
	{
		vg::Alignment result;
		result.set_name(seq_id);
		auto path = new vg::Path;
		result.set_allocated_path(path);
		size_t pos = 0;
		size_t oldNode = graph.indexToNode[trace[0].first];
		while (oldNode == graph.dummyNodeStart)
		{
			pos++;
			if (pos == trace.size()) return emptyAlignment();
			assert(pos < trace.size());
			oldNode = graph.indexToNode[trace[pos].first];
			assert(oldNode < graph.nodeIDs.size());
		}
		if (oldNode == graph.dummyNodeEnd) return emptyAlignment();
		int rank = 0;
		auto vgmapping = path->add_mapping();
		auto position = new vg::Position;
		vgmapping->set_allocated_position(position);
		vgmapping->set_rank(rank);
		position->set_node_id(graph.nodeIDs[oldNode]);
		position->set_is_reverse(graph.reverse[oldNode]);
		position->set_offset(trace[pos].first - graph.nodeStart[oldNode]);
		MatrixPosition btNodeStart = trace[pos];
		MatrixPosition btNodeEnd = trace[pos];
		for (; pos < trace.size(); pos++)
		{
			if (graph.indexToNode[trace[pos].first] == graph.dummyNodeEnd) break;
			if (graph.indexToNode[trace[pos].first] == oldNode)
			{
				btNodeEnd = trace[pos];
				continue;
			}
			assert(graph.indexToNode[btNodeEnd.first] == graph.indexToNode[btNodeStart.first]);
			assert(btNodeEnd.second >= btNodeStart.second);
			assert(btNodeEnd.first >= btNodeStart.first);
			auto edit = vgmapping->add_edit();
			edit->set_from_length(btNodeEnd.first - btNodeStart.first + 1);
			edit->set_to_length(btNodeEnd.second - btNodeStart.second + 1);
			edit->set_sequence(sequence.substr(btNodeStart.second, btNodeEnd.second - btNodeStart.second + 1));
			oldNode = graph.indexToNode[trace[pos].first];
			btNodeStart = trace[pos];
			btNodeEnd = trace[pos];
			rank++;
			vgmapping = path->add_mapping();
			position = new vg::Position;
			vgmapping->set_allocated_position(position);
			vgmapping->set_rank(rank);
			position->set_node_id(graph.nodeIDs[oldNode]);
			position->set_is_reverse(graph.reverse[oldNode]);
		}
		auto edit = vgmapping->add_edit();
		edit->set_from_length(btNodeEnd.first - btNodeStart.first);
		edit->set_to_length(btNodeEnd.second - btNodeStart.second);
		edit->set_sequence(sequence.substr(btNodeStart.second, btNodeEnd.second - btNodeStart.second));
		result.set_score(score);
		result.set_sequence(sequence);
		return AlignmentResult { result, false, cellsProcessed };
	}

	void expandBandForwards(std::vector<std::vector<LengthType>>& result, LengthType w, LengthType j, size_t sequenceLength, LengthType maxRow) const
	{
		if (std::find(result[j].begin(), result[j].end(), w) != result[j].end()) return;
		auto nodeIndex = graph.indexToNode[w];
		auto end = graph.nodeEnd[nodeIndex];
		while (w != end && j < sequenceLength+1 && j < maxRow)
		{
			result[j].emplace_back(w);
			w++;
			j++;
		}
		if (w == end && j < sequenceLength+1 && j < maxRow)
		{
			for (size_t i = 0; i < graph.outNeighbors[nodeIndex].size(); i++)
			{
				expandBandForwards(result, graph.nodeStart[graph.outNeighbors[nodeIndex][i]], j, sequenceLength, maxRow);
			}
		}
	}

	void expandBandBackwards(std::vector<std::vector<LengthType>>& result, LengthType w, LengthType j, size_t sequenceLength) const
	{
		if (std::find(result[j].begin(), result[j].end(), w) != result[j].end()) return;
		auto nodeIndex = graph.indexToNode[w];
		auto start = graph.nodeStart[nodeIndex];
		while (w != start && j > 0)
		{
			result[j].emplace_back(w);
			w--;
			j--;
		}
		result[j].emplace_back(w);
		if (w == start && j > 0)
		{
			for (size_t i = 0; i < graph.inNeighbors[nodeIndex].size(); i++)
			{
				expandBandBackwards(result, graph.nodeEnd[graph.inNeighbors[nodeIndex][i]] - 1, j-1, sequenceLength);
			}
		}
	}

	class IndexWithScore
	{
	public:
		IndexWithScore(LengthType index, ScoreType score) :
		index(index),
		score(score)
		{}
		LengthType index;
		ScoreType score;
		bool operator>(const IndexWithScore& other) const
		{
			return score > other.score;
		}
	};

	void expandBandFromPrevious(std::vector<bool>& currentBand, const std::vector<bool>& previousBand, const std::vector<WordSlice>& previousSlice, LengthType dynamicWidth, const std::vector<ScoreType>& nodeMinScores) const
	{
		//todo optimization: 18% inclusive 17% exclusive. can this be improved?
		assert(dynamicWidth < graph.nodeSequences.size());
		assert(currentBand.size() == previousBand.size());
		assert(currentBand.size() == graph.nodeStart.size());
		std::priority_queue<IndexWithScore, std::vector<IndexWithScore>, std::greater<IndexWithScore>> nodeQueue;
		std::priority_queue<IndexWithScore, std::vector<IndexWithScore>, std::greater<IndexWithScore>> endQueue;
		for (size_t i = 0; i < graph.nodeStart.size(); i++)
		{
			if (!previousBand[i]) continue;
			nodeQueue.emplace(i, nodeMinScores[i]);
			endQueue.emplace(i, previousSlice[graph.nodeEnd[i]-1].scoreEnd);
		}
		LengthType currentWidth = 0;
		while (currentWidth < dynamicWidth)
		{
			assert(nodeQueue.size() > 0 || endQueue.size() > 0);
			if (nodeQueue.size() == 0)
			{
				auto nextNode = endQueue.top();
				endQueue.pop();
				for (size_t i = 0; i < graph.outNeighbors[nextNode.index].size(); i++)
				{
					auto neighbor = graph.outNeighbors[nextNode.index][i];
					if (!currentBand[neighbor])
					{
						currentBand[neighbor] = true;
						currentWidth += graph.nodeEnd[neighbor] - graph.nodeStart[neighbor];
						endQueue.emplace(neighbor, nextNode.score + graph.nodeEnd[neighbor] - graph.nodeStart[neighbor]);
					}
				}
				continue;
			}
			if (endQueue.size() == 0)
			{
				auto nextNode = nodeQueue.top();
				nodeQueue.pop();
				if (!currentBand[nextNode.index])
				{
					currentBand[nextNode.index] = true;
					currentWidth += graph.nodeEnd[nextNode.index] - graph.nodeStart[nextNode.index];
				}
				continue;
			}
			auto nodeBest = nodeQueue.top();
			auto endBest = endQueue.top();
			if (nodeBest.score <= endBest.score)
			{
				nodeQueue.pop();
				assert(previousBand[nodeBest.index]);
				if (!currentBand[nodeBest.index])
				{
					currentBand[nodeBest.index] = true;
					currentWidth += graph.nodeEnd[nodeBest.index] - graph.nodeStart[nodeBest.index];
				}
			}
			else
			{
				endQueue.pop();
				assert(currentBand[endBest.index]);
				for (size_t i = 0; i < graph.outNeighbors[endBest.index].size(); i++)
				{
					auto neighbor = graph.outNeighbors[endBest.index][i];
					if (!currentBand[neighbor])
					{
						currentBand[neighbor] = true;
						currentWidth += graph.nodeEnd[neighbor] - graph.nodeStart[neighbor];
						endQueue.emplace(neighbor, endBest.score + graph.nodeEnd[neighbor] - graph.nodeStart[neighbor]);
					}
				}
			}
		}
	}

	void expandBandDynamically(std::vector<bool>& band, LengthType previousMinimumIndex, LengthType dynamicWidth) const
	{
		assert(previousMinimumIndex < graph.nodeSequences.size());
		auto nodeIndex = graph.indexToNode[previousMinimumIndex];
		band[nodeIndex] = true;
		LengthType start = graph.nodeStart[nodeIndex];
		LengthType end = graph.nodeEnd[nodeIndex];
		if (dynamicWidth > previousMinimumIndex - start)
		{
			for (size_t i = 0; i < graph.inNeighbors[nodeIndex].size(); i++)
			{
				expandDynamicBandBackward(band, graph.inNeighbors[nodeIndex][i], dynamicWidth - (previousMinimumIndex - start));
			}
		}
		if (dynamicWidth > end - previousMinimumIndex)
		{
			for (size_t i = 0; i < graph.outNeighbors[nodeIndex].size(); i++)
			{
				expandDynamicBandForward(band, graph.outNeighbors[nodeIndex][i], dynamicWidth - (end - previousMinimumIndex));
			}
		}
	}

	void expandDynamicBandBackward(std::vector<bool>& band, LengthType nodeIndex, LengthType dynamicWidth) const
	{
		if (dynamicWidth == 0) return;
		assert(nodeIndex < graph.nodeStart.size());
		//todo fix: currently only the first path that reaches the node is considered
		//this means that it might not reach some nodes with distance < dynamicwidth
		if (band[nodeIndex]) return;
		band[nodeIndex] = true;
		for (size_t i = 0; i < graph.outNeighbors[nodeIndex].size(); i++)
		{
			expandDynamicBandForward(band, graph.outNeighbors[nodeIndex][i], dynamicWidth - 1);
		}
		auto nodeSize = graph.nodeEnd[nodeIndex] - graph.nodeStart[nodeIndex];
		if (dynamicWidth > nodeSize)
		{
			for (size_t i = 0; i < graph.inNeighbors[nodeIndex].size(); i++)
			{
				expandDynamicBandBackward(band, graph.inNeighbors[nodeIndex][i], dynamicWidth - nodeSize);
			}
		}
	}

	template <typename MatrixType>
	void expandDynamicBandForward(MatrixType& band, LengthType nodeIndex, LengthType dynamicWidth) const
	{
		if (dynamicWidth == 0) return;
		assert(nodeIndex < graph.nodeStart.size());
		//todo fix: currently only the first path that reaches the node is considered
		//this means that it might not reach some nodes with distance < dynamicwidth if there's multiple paths and it arbitrarily picks the longest one first
		if (band[nodeIndex]) return;
		band[nodeIndex] = true;
		for (size_t i = 0; i < graph.inNeighbors[nodeIndex].size(); i++)
		{
			expandDynamicBandBackward(band, graph.inNeighbors[nodeIndex][i], dynamicWidth - 1);
		}
		auto nodeSize = graph.nodeEnd[nodeIndex] - graph.nodeStart[nodeIndex];
		if (dynamicWidth > nodeSize)
		{
			for (size_t i = 0; i < graph.outNeighbors[nodeIndex].size(); i++)
			{
				expandDynamicBandForward(band, graph.outNeighbors[nodeIndex][i], dynamicWidth - nodeSize);
			}
		}
	}

	uint64_t bytePrefixSums(uint64_t value, int addition) const
	{
		value <<= WordConfiguration<Word>::ChunkBits;
		assert(addition >= 0);
		value += addition;
		return value * WordConfiguration<Word>::PrefixSumMultiplierConstant;
	}

	uint64_t byteVPVNSum(uint64_t prefixSumVP, uint64_t prefixSumVN) const
	{
		uint64_t result = WordConfiguration<Word>::SignMask;
		assert((prefixSumVP & result) == 0);
		assert((prefixSumVN & result) == 0);
		result += prefixSumVP;
		result -= prefixSumVN;
		result ^= WordConfiguration<Word>::SignMask;
		return result;
	}

#ifdef EXTRAASSERTIONS

	WordSlice getWordSliceCellByCell(size_t j, size_t w, const std::string& sequence, const std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand) const
	{
		const auto lastBitMask = ((Word)1) << (WordConfiguration<Word>::WordSize-1);
		WordSlice result;
		auto nodeIndex = graph.indexToNode[w];
		assert(currentBand[nodeIndex]);
		ScoreType current[66];
		current[0] = j+1;
		current[1] = j;
		if (j > 0 && previousBand[nodeIndex]) current[1] = std::min(current[1], previousSlice[w].scoreEnd);
		if (j > 0 && previousBand[nodeIndex]) current[0] = std::min(current[0], previousSlice[w].scoreEnd - ((previousSlice[w].VP & lastBitMask) ? 1 : 0) + ((previousSlice[w].VN & lastBitMask) ? 1 : 0));
		for (int i = 1; i < 65; i++)
		{
			current[i+1] = current[i]+1;
		}
		if (w == graph.nodeStart[nodeIndex])
		{
			for (size_t neighbori = 0; neighbori < graph.inNeighbors[nodeIndex].size(); neighbori++)
			{
				auto neighbor = graph.inNeighbors[nodeIndex][neighbori];
				if (!previousBand[neighbor] && !currentBand[neighbor]) continue;
				auto u = graph.nodeEnd[neighbor]-1;
				ScoreType previous[66];
				previous[0] = j+1;
				previous[1] = j;
				if (j > 0 && previousBand[neighbor]) previous[1] = std::min(previous[1], previousSlice[u].scoreEnd);
				if (j > 0 && previousBand[neighbor]) previous[0] = std::min(previous[0], previousSlice[u].scoreEnd - ((previousSlice[u].VP & lastBitMask) ? 1 : 0) + ((previousSlice[u].VN & lastBitMask) ? 1 : 0));
				if (currentBand[neighbor]) previous[1] = std::min(previous[1], currentSlice[u].scoreBeforeStart);
				for (int i = 1; i < 65; i++)
				{
					if (currentBand[neighbor])
					{
						previous[i+1] = previous[i];
						previous[i+1] += (currentSlice[u].VP & (((Word)1) << (i-1)) ? 1 : 0);
						previous[i+1] -= (currentSlice[u].VN & (((Word)1) << (i-1)) ? 1 : 0);
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
					if (j+i > 0 && (sequence[j+i-1] == graph.nodeSequences[w] || sequence[j+i-1] == 'N'))
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
			auto u = w-1;
			ScoreType previous[66];
			previous[0] = currentSlice[u].scoreBeforeStart+1;
			previous[1] = currentSlice[u].scoreBeforeStart;
			if (previousBand[nodeIndex]) previous[0] = std::min(previous[0], previousSlice[u].scoreEnd - ((previousSlice[u].VP & lastBitMask) ? 1 : 0) + ((previousSlice[u].VN & lastBitMask) ? 1 : 0));
			if (previousBand[nodeIndex]) previous[1] = std::min(previous[1], previousSlice[u].scoreEnd);
			for (int i = 1; i < 65; i++)
			{
				previous[i+1] = previous[i];
				previous[i+1] += (currentSlice[u].VP & (((Word)1) << (i-1)) ? 1 : 0);
				previous[i+1] -= (currentSlice[u].VN & (((Word)1) << (i-1)) ? 1 : 0);
			}
			current[0] = std::min(current[0], previous[0]+1);
			for (int i = 0; i < 65; i++)
			{
				current[i+1] = std::min(current[i+1], current[i]+1);
				current[i+1] = std::min(current[i+1], previous[i+1]+1);
				if (j+i > 0 && (sequence[j+i-1] == graph.nodeSequences[w] || sequence[j+i-1] == 'N'))
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

#ifdef EXTRAASSERTIONS
	std::pair<uint64_t, uint64_t> differenceMasksCellByCell(uint64_t leftVP, uint64_t leftVN, uint64_t rightVP, uint64_t rightVN, int scoreDifference) const
	{
		int leftscore = 0;
		int rightscore = scoreDifference;
		uint64_t leftSmaller = 0;
		uint64_t rightSmaller = 0;
		for (int i = 0; i < 64; i++)
		{
			leftscore += leftVP & 1;
			leftscore -= leftVN & 1;
			rightscore += rightVP & 1;
			rightscore -= rightVN & 1;
			leftVP >>= 1;
			leftVN >>= 1;
			rightVP >>= 1;
			rightVN >>= 1;
			if (leftscore < rightscore) leftSmaller |= ((Word)1) << i;
			if (rightscore < leftscore) rightSmaller |= ((Word)1) << i;
		}
		return std::make_pair(leftSmaller, rightSmaller);
	}
#endif

	std::pair<uint64_t, uint64_t> differenceMasks(uint64_t leftVP, uint64_t leftVN, uint64_t rightVP, uint64_t rightVN, int scoreDifference) const
	{
#ifdef EXTRAASSERTIONS
		auto correctValue = differenceMasksCellByCell(leftVP, leftVN, rightVP, rightVN, scoreDifference);
#endif
		const uint64_t signmask = WordConfiguration<Word>::SignMask;
		const uint64_t lsbmask = WordConfiguration<Word>::LSBMask;
		const int chunksize = WordConfiguration<Word>::ChunkBits;
		const uint64_t allones = WordConfiguration<Word>::AllOnes;
		const uint64_t allzeros = WordConfiguration<Word>::AllZeros;
		uint64_t VPcommon = ~(leftVP & rightVP);
		uint64_t VNcommon = ~(leftVN & rightVN);
		leftVP &= VPcommon;
		leftVN &= VNcommon;
		rightVP &= VPcommon;
		rightVN &= VNcommon;
		//left is lower everywhere
		if (scoreDifference > WordConfiguration<Word>::popcount(rightVN) + WordConfiguration<Word>::popcount(leftVP))
		{
			return std::make_pair(allones, allzeros);
		}
		assert(scoreDifference >= 0);
		uint64_t byteVPVNSumLeft = byteVPVNSum(bytePrefixSums(WordConfiguration<Word>::ChunkPopcounts(leftVP), 0), bytePrefixSums(WordConfiguration<Word>::ChunkPopcounts(leftVN), 0));
		uint64_t byteVPVNSumRight = byteVPVNSum(bytePrefixSums(WordConfiguration<Word>::ChunkPopcounts(rightVP), scoreDifference), bytePrefixSums(WordConfiguration<Word>::ChunkPopcounts(rightVN), 0));
		uint64_t difference = byteVPVNSumLeft;
		{
			//take the bytvpvnsumright and split it from positive/negative values into two vectors with positive values, one which needs to be added and the other deducted
			//smearmask is 1 where the number needs to be deducted, and 0 where it needs to be added
			//except sign bits which are all 0
			uint64_t smearmask = ((byteVPVNSumRight & signmask) >> (chunksize-1)) * ((((Word)1) << (chunksize-1))-1);
			assert((smearmask & signmask) == 0);
			uint64_t deductions = ~smearmask & byteVPVNSumRight & ~signmask;
			//byteVPVNSumRight is in one's complement so take the not-value + 1
			uint64_t additions = (smearmask & ~byteVPVNSumRight) + (smearmask & lsbmask);
			assert((deductions & signmask) == 0);
			uint64_t signsBefore = difference & signmask;
			//unset the sign bits so additions don't interfere with other chunks
			difference &= ~signmask;
			difference += additions;
			//the sign bit is 1 if the value went from <0 to >=0
			//so in that case we need to flip it
			difference ^= signsBefore;
			signsBefore = difference & signmask;
			//set the sign bits so that deductions don't interfere with other chunks
			difference |= signmask;
			difference -= deductions;
			//sign bit is 0 if the value went from >=0 to <0
			//so flip them to the correct values
			signsBefore ^= signmask & ~difference;
			difference &= ~signmask;
			difference |= signsBefore;
		}
		//difference now contains the prefix sum difference (left-right) at each chunk
		uint64_t resultLeftSmallerThanRight = 0;
		uint64_t resultRightSmallerThanLeft = 0;
		for (int bit = 0; bit < chunksize; bit++)
		{
			uint64_t signsBefore = difference & signmask;
			//unset the sign bits so additions don't interfere with other chunks
			difference &= ~signmask;
			difference += leftVP & lsbmask;
			difference += rightVN & lsbmask;
			//the sign bit is 1 if the value went from <0 to >=0
			//so in that case we need to flip it
			difference ^= signsBefore;
			signsBefore = difference & signmask;
			//set the sign bits so that deductions don't interfere with other chunks
			difference |= signmask;
			difference -= leftVN & lsbmask;
			difference -= rightVP & lsbmask;
			//sign bit is 0 if the value went from >=0 to <0
			//so flip them to the correct values
			signsBefore ^= signmask & ~difference;
			difference &= ~signmask;
			difference |= signsBefore;
			leftVN >>= 1;
			leftVP >>= 1;
			rightVN >>= 1;
			rightVP >>= 1;
			//difference now contains the prefix sums difference (left-right) at each byte at (bit)'th bit
			//left < right when the prefix sum difference is negative (sign bit is set)
			uint64_t negative = (difference & signmask);
			resultLeftSmallerThanRight |= negative >> (WordConfiguration<Word>::ChunkBits - 1 - bit);
			//Test equality to zero. If it's zero, substracting one will make the sign bit 0, otherwise 1
			uint64_t notEqualToZero = ((difference | signmask) - lsbmask) & signmask;
			//right > left when the prefix sum difference is positive (not zero and not negative)
			resultRightSmallerThanLeft |= (notEqualToZero & ~negative) >> (WordConfiguration<Word>::ChunkBits - 1 - bit);
		}
#ifdef EXTRAASSERTIONS
		assert(resultLeftSmallerThanRight == correctValue.first);
		assert(resultRightSmallerThanLeft == correctValue.second);
#endif
		return std::make_pair(resultLeftSmallerThanRight, resultRightSmallerThanLeft);
	}

	WordSlice mergeTwoSlices(WordSlice left, WordSlice right) const
	{
		//optimization: 11% time inclusive 9% exclusive. can this be improved?
		//O(log w), because prefix sums need log w chunks of log w bits
		static_assert(std::is_same<Word, uint64_t>::value);
#ifdef EXTRAASSERTIONS
		auto correctValue = mergeTwoSlicesCellByCell(left, right);
#endif
		if (left.scoreBeforeStart > right.scoreBeforeStart) std::swap(left, right);
		WordSlice result;
		assert((left.VP & left.VN) == WordConfiguration<Word>::AllZeros);
		assert((right.VP & right.VN) == WordConfiguration<Word>::AllZeros);
		auto masks = differenceMasks(left.VP, left.VN, right.VP, right.VN, right.scoreBeforeStart - left.scoreBeforeStart);
		auto leftSmaller = masks.first;
		auto rightSmaller = masks.second;
		auto mask = (rightSmaller | ((leftSmaller | rightSmaller) - (rightSmaller << 1))) & ~leftSmaller;
		uint64_t leftReduction = leftSmaller & (rightSmaller << 1);
		uint64_t rightReduction = rightSmaller & (leftSmaller << 1);
		if ((rightSmaller & 1) && left.scoreBeforeStart < right.scoreBeforeStart)
		{
			rightReduction |= 1;
		}
		assert((leftReduction & right.VP) == leftReduction);
		assert((rightReduction & left.VP) == rightReduction);
		assert((leftReduction & left.VN) == leftReduction);
		assert((rightReduction & right.VN) == rightReduction);
		left.VN &= ~leftReduction;
		right.VN &= ~rightReduction;
		result.VN = (left.VN & ~mask) | (right.VN & mask);
		result.VP = (left.VP & ~mask) | (right.VP & mask);
		assert((result.VP & result.VN) == 0);
		result.scoreBeforeStart = std::min(left.scoreBeforeStart, right.scoreBeforeStart);
		result.scoreEnd = std::min(left.scoreEnd, right.scoreEnd);
		assert(result.scoreEnd == result.scoreBeforeStart + WordConfiguration<Word>::popcount(result.VP) - WordConfiguration<Word>::popcount(result.VN));
#ifdef EXTRAASSERTIONS
		assert(result.VP == correctValue.VP);
		assert(result.VN == correctValue.VN);
		assert(result.scoreBeforeStart == correctValue.scoreBeforeStart);
		assert(result.scoreEnd == correctValue.scoreEnd);
#endif
		return result;
	}

#ifdef EXTRAASSERTIONS
	WordSlice mergeTwoSlicesCellByCell(WordSlice left, WordSlice right) const
	{
		assert((left.VP & left.VN) == WordConfiguration<Word>::AllZeros);
		assert((right.VP & right.VN) == WordConfiguration<Word>::AllZeros);
		ScoreType leftScore = left.scoreBeforeStart;
		WordSlice merged;
		merged.scoreBeforeStart = std::min(left.scoreBeforeStart, right.scoreBeforeStart);
		merged.VP = WordConfiguration<Word>::AllZeros;
		merged.VN = WordConfiguration<Word>::AllZeros;
		ScoreType rightScore = right.scoreBeforeStart;
		ScoreType previousScore = merged.scoreBeforeStart;
		for (size_t j = 0; j < WordConfiguration<Word>::WordSize; j++)
		{
			Word mask = ((Word)1) << j;
			if (left.VP & mask) leftScore++;
			else if (left.VN & mask) leftScore--;
			if (right.VN & mask) rightScore--;
			else if (right.VP & mask) rightScore++;
			ScoreType betterScore = std::min(leftScore, rightScore);
			if (betterScore == previousScore+1) merged.VP |= mask;
			else if (betterScore == previousScore-1) merged.VN |= mask;
			assert((merged.VP & merged.VN) == WordConfiguration<Word>::AllZeros);
			assert(betterScore >= previousScore-1);
			assert(betterScore <= previousScore+1);
			previousScore = betterScore;
		}
		merged.scoreEnd = previousScore;
		assert((merged.VP & merged.VN) == WordConfiguration<Word>::AllZeros);
		assert(merged.scoreEnd <= left.scoreEnd);
		assert(merged.scoreEnd <= right.scoreEnd);
		assert(merged.scoreBeforeStart <= left.scoreBeforeStart);
		assert(merged.scoreBeforeStart <= right.scoreBeforeStart);
		return merged;
	}
#endif

	WordSlice getNodeStartSlice(Word Eq, size_t nodeIndex, const std::vector<WordSlice>& previousSlice, const std::vector<WordSlice>& currentSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand, bool previousEq) const
	{
		//todo optimization: 10% time inclusive 4% exclusive. can this be improved?
		WordSlice previous;
		WordSlice previousUp;
		bool foundOne = false;
		bool foundOneUp = false;
		for (size_t i = 0; i < graph.inNeighbors[nodeIndex].size(); i++)
		{
			auto neighbor = graph.inNeighbors[nodeIndex][i];
			if (previousBand[neighbor])
			{
				if (!foundOneUp)
				{
					previousUp = previousSlice[graph.nodeEnd[neighbor]-1];
					foundOneUp = true;
				}
				else
				{
					auto competitor = previousSlice[graph.nodeEnd[neighbor]-1];
					previousUp = mergeTwoSlices(previousUp, competitor);
				}
			}
			if (previousBand[neighbor] && !currentBand[neighbor])
			{
				if (!foundOne)
				{
					previous = getSourceSliceFromColumn(graph.nodeEnd[neighbor]-1, previousSlice);
					foundOne = true;
				}
				else
				{
					auto competitor = getSourceSliceFromColumn(graph.nodeEnd[neighbor]-1, previousSlice);
					previous = mergeTwoSlices(previous, competitor);
				}
			}
			if (!currentBand[neighbor]) continue;
			if (!foundOne)
			{
				previous = currentSlice[graph.nodeEnd[neighbor]-1];
				foundOne = true;
			}
			else
			{
				auto competitor = currentSlice[graph.nodeEnd[neighbor]-1];
				previous = mergeTwoSlices(previous, competitor);
			}
		}
		assert(foundOne);
		auto result = getNextSlice(Eq, previous, foundOneUp, previousEq, previousUp);
		return result;
	}

	WordSlice getSourceSliceWithoutBefore(size_t row) const
	{
		return { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, row+WordConfiguration<Word>::WordSize, row };
	}

	WordSlice getSourceSliceFromColumn(size_t column, const std::vector<WordSlice>& previousSlice) const
	{
		return { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, previousSlice[column].scoreEnd+WordConfiguration<Word>::WordSize, previousSlice[column].scoreEnd };
	}

	WordSlice getSourceSlice(size_t nodeIndex, const std::vector<WordSlice>& previousSlice) const
	{
		auto start = graph.nodeStart[nodeIndex];
		return getSourceSliceFromColumn(start, previousSlice);
	}

	bool isSource(size_t nodeIndex, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand) const
	{
		for (size_t i = 0; i < graph.inNeighbors[nodeIndex].size(); i++)
		{
			if (currentBand[graph.inNeighbors[nodeIndex][i]]) return false;
			if (previousBand[graph.inNeighbors[nodeIndex][i]]) return false;
		}
		return true;
	}

	Word getEq(Word BA, Word BT, Word BC, Word BG, LengthType w) const
	{
		switch(graph.nodeSequences[w])
		{
			case 'A':
			return BA;
			break;
			case 'T':
			return BT;
			break;
			case 'C':
			return BC;
			break;
			case 'G':
			return BG;
			break;
			case '-':
			assert(false);
			break;
			default:
			assert(false);
			break;
		}
		assert(false);
		return 0;
	}

	WordSlice getNextSlice(Word Eq, WordSlice slice, bool previousInsideBand, bool previousEq, WordSlice previous) const
	{
		//optimization: 13% of time. probably can't be improved easily.
		//http://www.gersteinlab.org/courses/452/09-spring/pdf/Myers.pdf
		//pages 405 and 408

		auto oldValue = slice.scoreBeforeStart;
		if (!previousInsideBand)
		{
			slice.scoreBeforeStart += 1;
		}
		else
		{
			const auto lastBitMask = ((Word)1) << (WordConfiguration<Word>::WordSize - 1);
			assert(slice.scoreBeforeStart <= previous.scoreEnd);
			slice.scoreBeforeStart = std::min(slice.scoreBeforeStart + 1, previous.scoreEnd - ((previous.VP & lastBitMask) ? 1 : 0) + ((previous.VN & lastBitMask) ? 1 : 0) + (previousEq ? 0 : 1));
		}
		auto hin = slice.scoreBeforeStart - oldValue;

		Word Xv = Eq | slice.VN;
		//between 7 and 8
		if (hin < 0) Eq |= 1;
		Word Xh = (((Eq & slice.VP) + slice.VP) ^ slice.VP) | Eq;
		Word Ph = slice.VN | ~(Xh | slice.VP);
		Word Mh = slice.VP & Xh;
		Word lastBitMask = (((Word)1) << (WordConfiguration<Word>::WordSize - 1));
		if (Ph & lastBitMask)
		{
			slice.scoreEnd += 1;
		}
		else if (Mh & lastBitMask)
		{
			slice.scoreEnd -= 1;
		}
		Ph <<= 1;
		Mh <<= 1;
		//between 16 and 17
		if (hin < 0) Mh |= 1; else if (hin > 0) Ph |= 1;
		slice.VP = Mh | ~(Xv | Ph);
		slice.VN = Ph & Xv;

#ifndef NDEBUG
		auto wcvp = WordConfiguration<Word>::popcount(slice.VP);
		auto wcvn = WordConfiguration<Word>::popcount(slice.VN);
		assert(slice.scoreEnd == slice.scoreBeforeStart + wcvp - wcvn);
#endif

		return slice;
	}

	class NodeCalculationResult
	{
	public:
		ScoreType minScore;
		LengthType minScoreIndex;
		size_t cellsProcessed;
	};

	bool firstZeroForced(const std::vector<bool>& previousBand, const std::vector<bool>& currentBand, LengthType neighborNodeIndex, LengthType neighborColumn, Word currentEq, const std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice) const
	{
		if (previousBand[neighborNodeIndex] && currentBand[neighborNodeIndex])
		{
			if (currentSlice[neighborColumn].VN & 1)
			{
				return true;
			}
			if (!(currentSlice[neighborColumn].VP & 1) && !(currentSlice[neighborColumn].VN & 1) && !(currentEq & 1))
			{
				return true;
			}
			return false;
		}
		else if (previousBand[neighborNodeIndex] && !currentBand[neighborNodeIndex])
		{
			return false;
		}
		else
		{
			return true;
		}
		assert(false);
	}

	NodeCalculationResult calculateNode(size_t i, size_t j, const std::string& sequence, Word BA, Word BT, Word BC, Word BG, std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand, bool forceSource) const
	{
		//todo optimization: 42% inclusive 15% exclusive. can this be improved?
		NodeCalculationResult result;
		result.minScore = std::numeric_limits<ScoreType>::max();
		result.minScoreIndex = 0;
		result.cellsProcessed = 0;

		LengthType start = graph.nodeStart[i];
		if (forceSource || isSource(i, currentBand, previousBand))
		{
			if (previousBand[i])
			{
				currentSlice[start] = getSourceSlice(i, previousSlice);
			}
			else
			{
				currentSlice[start] = getSourceSliceWithoutBefore(j);
			}
			if (currentSlice[start].scoreEnd < result.minScore)
			{
				result.minScore = currentSlice[start].scoreEnd;
				result.minScoreIndex = start;
			}
			start++;
		}
		else
		{
			Word Eq = getEq(BA, BT, BC, BG, graph.nodeStart[i]);
			currentSlice[start] = getNodeStartSlice(Eq, i, previousSlice, currentSlice, currentBand, previousBand, j == 0 || graph.nodeSequences[graph.nodeStart[i]] == sequence[j-1]);
			if (previousBand[i] && currentSlice[start].scoreBeforeStart > previousSlice[start].scoreEnd)
			{
				currentSlice[start] = mergeTwoSlices(getSourceSliceFromColumn(start, previousSlice), currentSlice[start]);
			}
			if (currentSlice[start].scoreBeforeStart > j)
			{
				currentSlice[start] = mergeTwoSlices(getSourceSliceWithoutBefore(j), currentSlice[start]);
			}
			if (currentSlice[start].scoreEnd < result.minScore)
			{
				result.minScore = currentSlice[start].scoreEnd;
				result.minScoreIndex = start;
			}
			//note: currentSlice[start].score - optimalInNeighborEndScore IS NOT within {-1, 0, 1} always because of the band
			start++;
		}

#ifdef EXTRAASSERTIONS
		if (!forceSource)
		{
			auto correctstart = getWordSliceCellByCell(j, graph.nodeStart[i], sequence, currentSlice, previousSlice, currentBand, previousBand);
			assert(currentSlice[graph.nodeStart[i]].scoreBeforeStart == correctstart.scoreBeforeStart);
			assert(currentSlice[graph.nodeStart[i]].scoreEnd == correctstart.scoreEnd);
			assert(currentSlice[graph.nodeStart[i]].VP == correctstart.VP);
			assert(currentSlice[graph.nodeStart[i]].VN == correctstart.VN);
		}
#endif

		assert(start == graph.nodeStart[i]+1);

		for (LengthType w = start; w < graph.nodeEnd[i]; w++)
		{
			Word Eq = getEq(BA, BT, BC, BG, w);
			bool forceFirstHorizontalPositive = firstZeroForced(previousBand, currentBand, i, w-1, Eq, currentSlice, previousSlice);

			currentSlice[w] = getNextSlice(Eq, currentSlice[w-1], previousBand[i], j == 0 || graph.nodeSequences[w] == sequence[j-1], previousSlice[w-1]);

			if (previousBand[i] && currentSlice[w].scoreBeforeStart > previousSlice[w].scoreEnd)
			{
				currentSlice[w] = mergeTwoSlices(getSourceSliceFromColumn(w, previousSlice), currentSlice[w]);
			}
			if (currentSlice[w].scoreBeforeStart > j)
			{
				currentSlice[w] = mergeTwoSlices(getSourceSliceWithoutBefore(j), currentSlice[w]);
			}

#ifndef NDEBUG
			auto wcvp = WordConfiguration<Word>::popcount(currentSlice[w].VP);
			auto wcvn = WordConfiguration<Word>::popcount(currentSlice[w].VN);
			assert(currentSlice[w].scoreEnd == currentSlice[w].scoreBeforeStart + wcvp - wcvn);
#endif
			assert(previousBand[i] || currentSlice[w].scoreBeforeStart == j || currentSlice[w].scoreBeforeStart == currentSlice[w-1].scoreBeforeStart + 1);
			assert(currentSlice[w].scoreBeforeStart >= 0);
			assert(currentSlice[w].scoreEnd >= 0);
			assert(currentSlice[w].scoreBeforeStart <= currentSlice[w].scoreEnd + WordConfiguration<Word>::WordSize);
			assert(currentSlice[w].scoreEnd <= currentSlice[w].scoreBeforeStart + WordConfiguration<Word>::WordSize);
			assert((currentSlice[w].VP & currentSlice[w].VN) == WordConfiguration<Word>::AllZeros);

			assert(!previousBand[i] || currentSlice[w].scoreBeforeStart <= previousSlice[w].scoreEnd);
			assert(currentSlice[w].scoreBeforeStart >= 0);
			if (currentSlice[w].scoreEnd < result.minScore)
			{
				result.minScore = currentSlice[w].scoreEnd;
				result.minScoreIndex = w;
			}

#ifdef EXTRAASSERTIONS
			if (!forceSource)
			{
				auto correctslice = getWordSliceCellByCell(j, w, sequence, currentSlice, previousSlice, currentBand, previousBand);
				assert(currentSlice[w].scoreBeforeStart == correctslice.scoreBeforeStart);
				assert(currentSlice[w].scoreEnd == correctslice.scoreEnd);
				assert(currentSlice[w].VP == correctslice.VP);
				assert(currentSlice[w].VN == correctslice.VN);
			}
#endif
		}
		result.cellsProcessed = (graph.nodeEnd[i] - graph.nodeStart[i]) * WordConfiguration<Word>::WordSize;
		return result;
	}

	void calculateDoublesliceBackwardsRec(size_t i, size_t j, const std::string sequence, Word BA, Word BT, Word BC, Word BG, LengthType sizeLeft, std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand) const
	{
		assert(sizeLeft > 0);
		bool forceSource = true;
		if (graph.nodeEnd[i]-graph.nodeStart[i] < sizeLeft)
		{
			sizeLeft -= graph.nodeEnd[i]-graph.nodeStart[i];
			for (size_t neighbori = 0; neighbori < graph.inNeighbors[i].size(); neighbori++)
			{
				auto neighbor = graph.inNeighbors[i][neighbori];
				if (!currentBand[neighbor]) continue;
				calculateDoublesliceBackwardsRec(neighbor, j, sequence, BA, BT, BC, BG, sizeLeft, currentSlice, previousSlice, currentBand, previousBand);
				forceSource = false;
			}
		}
		calculateNode(i, j, sequence, BA, BT, BC, BG, currentSlice, previousSlice, currentBand, previousBand, forceSource);
	}

	void calculateDoublesliceBackwards(size_t j, const std::string& sequence, Word BA, Word BT, Word BC, Word BG, std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand) const
	{
		if (graph.firstInOrder == 0) return;
		//if there are cycles within 2*w of eachothers, calculating a latter slice may overwrite the earlier slice's value
		//store the correct values here and then merge them at the end
		std::vector<WordSlice> correctEndValues;
		correctEndValues.resize(graph.firstInOrder);
		for (size_t i = 1; i < graph.firstInOrder; i++)
		{
			assert(graph.notInOrder[i]);
			if (!currentBand[i]) continue;
			calculateDoublesliceBackwardsRec(i, j, sequence, BA, BT, BC, BG, WordConfiguration<Word>::WordSize*2, currentSlice, previousSlice, currentBand, previousBand);
			correctEndValues[i] = currentSlice[graph.nodeEnd[i]-1];
		}
		for (size_t i = 1; i < graph.firstInOrder; i++)
		{
			if (!currentBand[i]) continue;
			currentSlice[graph.nodeEnd[i]-1] = correctEndValues[i];
		}
	}

	MatrixSlice getBitvectorSliceScoresAndFinalPosition(const std::string& sequence, int dynamicWidth, std::vector<std::vector<bool>>& startBand, LengthType dynamicRowStart) const
	{
		//todo optimization: 82% inclusive 17% exclusive. can this be improved?
		MatrixSlice result;
		result.cellsProcessed = 0;
		result.finalMinScore = 0;
		result.finalMinScoreColumn = 0;
		result.minScorePerWordSlice.emplace_back(0);

		std::vector<WordSlice> slice1;
		std::vector<WordSlice> slice2;
		slice1.resize(graph.nodeSequences.size());
		slice2.resize(graph.nodeSequences.size());

		std::vector<WordSlice>& currentSlice = slice1;
		std::vector<WordSlice>& previousSlice = slice2;

		std::vector<ScoreType> nodeMinScores;
		nodeMinScores.resize(graph.nodeStart.size(), std::numeric_limits<ScoreType>::max());

		LengthType previousMinimumIndex;
		std::vector<bool> currentBand;
		std::vector<bool> previousBand;
		assert(startBand.size() > 0);
		assert(startBand[0].size() == graph.nodeStart.size());
		currentBand.resize(graph.nodeStart.size());
		previousBand = startBand[0];

		for (size_t j = 0; j < sequence.size(); j += WordConfiguration<Word>::WordSize)
		{
			ScoreType currentMinimumScore = std::numeric_limits<ScoreType>::max();
			LengthType currentMinimumIndex = 0;
			//preprocessed bitvectors for character equality
			Word BA = WordConfiguration<Word>::AllZeros;
			Word BT = WordConfiguration<Word>::AllZeros;
			Word BC = WordConfiguration<Word>::AllZeros;
			Word BG = WordConfiguration<Word>::AllZeros;
			for (int i = 0; i < WordConfiguration<Word>::WordSize && j+i < sequence.size(); i++)
			{
				Word mask = ((Word)1) << i;
				switch(sequence[j+i])
				{
					case 'A':
					case 'a':
					BA |= mask;
					break;
					case 'T':
					case 't':
					BT |= mask;
					break;
					case 'C':
					case 'c':
					BC |= mask;
					break;
					case 'G':
					case 'g':
					BG |= mask;
					break;
					case 'N':
					case 'n':
					BA |= mask;
					BC |= mask;
					BT |= mask;
					BG |= mask;
					break;
				}
			}
			size_t slice = j / WordConfiguration<Word>::WordSize;
			if (startBand.size() > slice)
			{
				if (slice > 0) previousBand = currentBand;
				currentBand = startBand[slice];
			}
			else
			{
				std::swap(currentBand, previousBand);
				currentBand.assign(currentBand.size(), false);
				expandBandFromPrevious(currentBand, previousBand, previousSlice, dynamicWidth, nodeMinScores);
			}
			calculateDoublesliceBackwards(j, sequence, BA, BT, BC, BG, currentSlice, previousSlice, currentBand, previousBand);
			for (size_t i = graph.firstInOrder; i < graph.nodeStart.size(); i++)
			{
				nodeMinScores[i] = std::numeric_limits<ScoreType>::max();
				if (!currentBand[i]) continue;
				auto nodeCalc = calculateNode(i, j, sequence, BA, BT, BC, BG, currentSlice, previousSlice, currentBand, previousBand, false);
				nodeMinScores[i] = nodeCalc.minScore;
				if (nodeCalc.minScore < currentMinimumScore)
				{
					currentMinimumScore = nodeCalc.minScore;
					currentMinimumIndex = nodeCalc.minScoreIndex;
				}
				result.cellsProcessed += nodeCalc.cellsProcessed;
			}
			for (size_t i = 0; i < graph.firstInOrder; i++)
			{
				nodeMinScores[i] = std::numeric_limits<ScoreType>::max();
				if (!currentBand[i]) continue;
				auto nodeCalc = calculateNode(i, j, sequence, BA, BT, BC, BG, currentSlice, previousSlice, currentBand, previousBand, false);
				nodeMinScores[i] = nodeCalc.minScore;
				if (nodeCalc.minScore < currentMinimumScore)
				{
					currentMinimumScore = nodeCalc.minScore;
					currentMinimumIndex = nodeCalc.minScoreIndex;
				}
				result.cellsProcessed += nodeCalc.cellsProcessed;
			}
			std::swap(currentSlice, previousSlice);
			previousMinimumIndex = currentMinimumIndex;
			result.minScorePerWordSlice.emplace_back(currentMinimumScore);
		}
		result.finalMinScoreColumn = previousMinimumIndex;
		result.finalMinScore = result.minScorePerWordSlice.back();
		return result;
	}

	std::tuple<ScoreType, std::vector<MatrixPosition>, size_t> getBacktrace(std::string sequence, int dynamicWidth, LengthType dynamicRowStart, std::vector<std::vector<bool>>& startBand) const
	{
		int padding = (WordConfiguration<Word>::WordSize - (sequence.size() % WordConfiguration<Word>::WordSize)) % WordConfiguration<Word>::WordSize;
		for (int i = 0; i < padding; i++)
		{
			sequence += 'N';
		}
		auto slice = getBitvectorSliceScoresAndFinalPosition(sequence, dynamicWidth, startBand, dynamicRowStart);
		std::cerr << "score: " << slice.minScorePerWordSlice.back() << std::endl;
		auto backtraceresult = backtrace(std::make_pair(slice.finalMinScoreColumn, sequence.size()), sequence, slice.minScorePerWordSlice);
		for (int i = 0; i < padding; i++)
		{
			assert(backtraceresult.back().second >= sequence.size() - padding);
			backtraceresult.pop_back();
		}
		//if there's a mismatch at the last base, the backtrace might be spending one more jump in the padding
		if (backtraceresult.back().second == sequence.size() - padding && graph.nodeSequences[backtraceresult.back().first] != sequence[backtraceresult.back().second])
		{
			backtraceresult.pop_back();
		}
		assert(backtraceresult.back().second == sequence.size() - padding - 1);
		return std::make_tuple(slice.finalMinScore, backtraceresult, slice.cellsProcessed);
	}


	const AlignmentGraph& graph;
};

#endif