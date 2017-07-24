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
	static constexpr uint64_t AllZeros = 0x0000000000000000;
	static constexpr uint64_t AllOnes = 0xFFFFFFFFFFFFFFFF;
	static constexpr int WordSize = 64;
	static int popcount(uint64_t x)
	{
		//https://en.wikipedia.org/wiki/Hamming_weight
		x -= (x >> 1) & 0x5555555555555555;
		x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
		x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f;
		return (x * 0x0101010101010101) >> 56;
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
		scoreStart(0),
		scoreEnd(0),
		scoreBeforeStart(0),
		hout(0)
		{}
		WordSlice(Word VP, Word VN, ScoreType scoreStart, ScoreType scoreEnd, ScoreType scoreBeforeStart) :
		VP(VP),
		VN(VN),
		scoreStart(scoreStart),
		scoreEnd(scoreEnd),
		scoreBeforeStart(scoreBeforeStart),
		hout(0)
		{}
		Word VP;
		Word VN;
		ScoreType scoreStart;
		ScoreType scoreEnd;
		ScoreType scoreBeforeStart;
		int hout;
	};
public:
	class AlignmentResult
	{
	public:
		AlignmentResult()
		{
		}
		AlignmentResult(vg::Alignment alignment, int maxDistanceFromBand, bool alignmentFailed, size_t cellsProcessed) :
		alignment(alignment),
		maxDistanceFromBand(maxDistanceFromBand),
		alignmentFailed(alignmentFailed),
		cellsProcessed(cellsProcessed)
		{
		}
		vg::Alignment alignment;
		int maxDistanceFromBand;
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
		auto result = traceToAlignment(seq_id, sequence, std::get<0>(trace), std::get<2>(trace), std::get<1>(trace), std::get<3>(trace));
		return result;
	}

	AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence, int dynamicWidth, LengthType dynamicRowStart, const std::vector<AlignmentGraph::SeedHit>& seedHits, int startBandwidth) const
	{
		assert(graph.finalized);
		auto band = getSeededStartBand(seedHits, dynamicRowStart, startBandwidth, sequence);
		auto trace = getBacktrace(sequence, dynamicWidth, dynamicRowStart, band);
		//failed alignment, don't output
		if (std::get<0>(trace) == std::numeric_limits<ScoreType>::min()) return emptyAlignment();
		auto result = traceToAlignment(seq_id, sequence, std::get<0>(trace), std::get<2>(trace), std::get<1>(trace), std::get<3>(trace));
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
		ScoreType scoreAtEnd = minScorePerWordSlice.back()+1;
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
		return AlignmentResult { result, 0, true, 0 };
	}

	AlignmentResult traceToAlignment(const std::string& seq_id, const std::string& sequence, ScoreType score, const std::vector<MatrixPosition>& trace, int maxDistanceFromBand, size_t cellsProcessed) const
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
		return AlignmentResult { result, maxDistanceFromBand, false, cellsProcessed };
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

	uint64_t bytePopcounts(uint64_t value) const
	{
		uint64_t x = value;
		x -= (x >> 1) & 0x5555555555555555;
		x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
		x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f;
		return x;
	}

	uint64_t bytePrefixSums(uint64_t value, int addition) const
	{
		value <<= 8;
		assert(addition >= 0);
		value += addition;
		return value * 0x0101010101010101;
	}

	uint64_t byteVPVNSum(uint64_t prefixSumVP, uint64_t prefixSumVN) const
	{
		uint64_t result = 0x8080808080808080;
		assert((prefixSumVP & result) == 0);
		assert((prefixSumVN & result) == 0);
		result += prefixSumVP;
		result -= prefixSumVN;
		return result;
	}

	std::string wordToStr(uint64_t value) const
	{
		std::string result;
		for (int i = 63; i >= 0; i--)
		{
			if (value & (((uint64_t)1) << i))
			{
				result += "1";
			}
			else
			{
				result += "0";
			}
		}
		return result;
	}

	std::pair<uint64_t, uint64_t> differenceMasks(uint64_t leftVP, uint64_t leftVN, uint64_t rightVP, uint64_t rightVN, int scoreDifference) const
	{
		//do the addition in chunks to prevent calculations from interfering with the other bytes
		const uint64_t chunkmask1 = 0xFF00FF00FF00FF00;
		const uint64_t chunkmask2 = 0x00FF00FF00FF00FF;
		//the fences make sure that when moving positive<->negative the other bytes are not effected
		const uint64_t fencemask1 = 0x0055005500550055;
		const uint64_t fencemask2 = 0x5500550055005500;
		const uint64_t signmask = 0x8080808080808080;
		uint64_t VPcommon = ~(leftVP & rightVP);
		uint64_t VNcommon = ~(leftVN & rightVN);
		leftVP &= VPcommon;
		leftVN &= VNcommon;
		rightVP &= VPcommon;
		rightVN &= VNcommon;
		//left is lower everywhere
		if (scoreDifference > WordConfiguration<Word>::popcount(rightVN) + WordConfiguration<Word>::popcount(leftVP))
		{
			return std::make_pair(0xFFFFFFFFFFFFFFFF, 0x0000000000000000);
		}
		assert(scoreDifference >= 0);
		uint64_t byteVPVNSumLeft = byteVPVNSum(bytePrefixSums(bytePopcounts(leftVP), 0), bytePrefixSums(bytePopcounts(leftVN), 0));
		uint64_t byteVPVNSumRight = byteVPVNSum(bytePrefixSums(bytePopcounts(rightVP), scoreDifference), bytePrefixSums(bytePopcounts(rightVN), 0));
		uint64_t left = (byteVPVNSumLeft & chunkmask1) | fencemask1;
		uint64_t right = (byteVPVNSumRight & chunkmask1);
		uint64_t difference = (left - right) & chunkmask1;
		left = (byteVPVNSumLeft & chunkmask2) | fencemask2;
		right = (byteVPVNSumRight & chunkmask2);
		difference |= (left - right) & chunkmask2;
		//difference now contains the prefix sum difference (left-right) at each byte
		uint64_t resultLeftSmallerThanRight = 0;
		uint64_t resultRightSmallerThanLeft = 0;
		const uint64_t LSBmask1 = 0x0101010101010101 & chunkmask1;
		const uint64_t LSBmask2 = 0x0101010101010101 & chunkmask2;
		for (int bit = 0; bit < 8; bit++)
		{
			uint64_t difference1 = (difference & chunkmask1) | fencemask1;
			uint64_t difference2 = (difference & chunkmask2) | fencemask2;
			difference1 += (leftVP & LSBmask1);
			difference1 -= (leftVN & LSBmask1);
			difference1 -= (rightVP & LSBmask1);
			difference1 += (rightVN & LSBmask1);
			difference2 += (leftVP & LSBmask2);
			difference2 -= (leftVN & LSBmask2);
			difference2 -= (rightVP & LSBmask2);
			difference2 += (rightVN & LSBmask2);
			leftVN >>= 1;
			leftVP >>= 1;
			rightVN >>= 1;
			rightVP >>= 1;
			difference = (difference1 & chunkmask1) | (difference2 & chunkmask2);
			//difference now contains the prefix sums difference (left-right) at each byte at (bit)'th bit
			//left < right when the prefix sum difference is negative (sign bit is set)
			uint64_t negative = (difference & signmask);
			resultLeftSmallerThanRight |= negative >> (7 - bit);
			//Test equality to zero. If it's zero, substracting one will make the sign bit 0, otherwise 1
			uint64_t notEqualToZero = ((difference | signmask) - 0x0101010101010101) & signmask;
			//right > left when the prefix sum difference is positive (not zero and not negative)
			resultRightSmallerThanLeft |= (notEqualToZero & ~negative) >> (7 - bit);
		}
		return std::make_pair(resultLeftSmallerThanRight, resultRightSmallerThanLeft);
	}

	WordSlice mergeTwoSlices(WordSlice left, WordSlice right) const
	{
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
		result.scoreStart = std::min(left.scoreStart, right.scoreStart);
		result.scoreEnd = std::min(left.scoreEnd, right.scoreEnd);
		assert(result.scoreEnd == result.scoreBeforeStart + WordConfiguration<Word>::popcount(result.VP) - WordConfiguration<Word>::popcount(result.VN));
		assert(result.scoreStart == result.scoreBeforeStart + (result.VP & 1) - (result.VN & 1));
#ifdef EXTRAASSERTIONS
		assert(result.VP == correctValue.VP);
		assert(result.VN == correctValue.VN);
		assert(result.scoreBeforeStart == correctValue.scoreBeforeStart);
		assert(result.scoreStart == correctValue.scoreStart);
		assert(result.scoreEnd == correctValue.scoreEnd);
#endif
		return result;
	}

	WordSlice getNodeStartSlice(Word Eq, size_t nodeIndex, const std::vector<WordSlice>& previousSlice, const std::vector<WordSlice>& currentSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand) const
	{
		WordSlice previous;
		bool foundOne = false;
		bool forceFirstHorizontalPositive = false;
		int hin = 0;
		if (previousBand[nodeIndex]) hin = previousSlice[graph.nodeStart[nodeIndex]].hout;
		for (size_t i = 0; i < graph.inNeighbors[nodeIndex].size(); i++)
		{
			auto neighbor = graph.inNeighbors[nodeIndex][i];
			if (!currentBand[neighbor]) continue;
			if (!foundOne)
			{
				previous = currentSlice[graph.nodeEnd[neighbor]-1];
				forceFirstHorizontalPositive = firstZeroForced(previousBand, neighbor, graph.nodeEnd[neighbor]-1, Eq, currentSlice, previousSlice);
				foundOne = true;
			}
			else
			{
				auto competitor = currentSlice[graph.nodeEnd[neighbor]-1];
				if (competitor.scoreBeforeStart < previous.scoreBeforeStart)
				{
					forceFirstHorizontalPositive = firstZeroForced(previousBand, neighbor, graph.nodeEnd[neighbor]-1, Eq, currentSlice, previousSlice);
				}
				else if (competitor.scoreBeforeStart == previous.scoreBeforeStart)
				{
					forceFirstHorizontalPositive = forceFirstHorizontalPositive & firstZeroForced(previousBand, neighbor, graph.nodeEnd[neighbor]-1, Eq, currentSlice, previousSlice);
				}
				previous = mergeTwoSlices(previous, competitor);
			}
		}
		assert(foundOne);
		auto result = getNextSlice(Eq, previous, hin, forceFirstHorizontalPositive);
		return result;
	}

	WordSlice getSourceSliceWithoutBefore(size_t row) const
	{
		return { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, row+1, row+WordConfiguration<Word>::WordSize, row };
	}

	WordSlice getSourceSliceFromColumn(size_t column, const std::vector<WordSlice>& previousSlice) const
	{
		return { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, previousSlice[column].scoreEnd+1, previousSlice[column].scoreEnd+WordConfiguration<Word>::WordSize, previousSlice[column].scoreEnd };
	}

	WordSlice getSourceSlice(size_t nodeIndex, const std::vector<WordSlice>& previousSlice) const
	{
		auto start = graph.nodeStart[nodeIndex];
		return getSourceSliceFromColumn(start, previousSlice);
	}

	bool isSource(size_t nodeIndex, const std::vector<bool>& band) const
	{
		for (size_t i = 0; i < graph.inNeighbors[nodeIndex].size(); i++)
		{
			if (band[graph.inNeighbors[nodeIndex][i]]) return false;
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

	WordSlice getNextSlice(Word Eq, WordSlice slice, int hin, bool forceFirstHorizontalPositive) const
	{
		//http://www.gersteinlab.org/courses/452/09-spring/pdf/Myers.pdf
		//pages 405 and 408
		Word Xv = Eq | slice.VN;
		//between 7 and 8
		if (hin < 0) Eq |= 1;
		if (forceFirstHorizontalPositive)
		{
			Eq &= ~((Word)1);
		}
		Word Xh = (((Eq & slice.VP) + slice.VP) ^ slice.VP) | Eq;
		Word Ph = slice.VN | ~(Xh | slice.VP);
		Word Mh = slice.VP & Xh;
		if (forceFirstHorizontalPositive)
		{
			Ph |= 1;
			Mh &= ~((Word)1);
		}
		if (Ph & ((Word)1))
		{
			slice.scoreStart += 1;
		}
		else if (Mh & ((Word)1))
		{
			slice.scoreStart -= 1;
		}
		Word lastBitMask = (((Word)1) << (WordConfiguration<Word>::WordSize - 1));
		if (Ph & lastBitMask)
		{
			slice.scoreEnd += 1;
			slice.hout = 1;
		}
		else if (Mh & lastBitMask)
		{
			slice.scoreEnd -= 1;
			slice.hout = -1;
		}
		else
		{
			slice.hout = 0;
		}
		Ph <<= 1;
		Mh <<= 1;
		//between 16 and 17
		if (hin < 0) Mh |= 1; else if (hin > 0) Ph |= 1;
		slice.VP = Mh | ~(Xv | Ph);
		slice.VN = Ph & Xv;

#ifndef NDEBUG
		auto wcvpExceptFirst = WordConfiguration<Word>::popcount(slice.VP & ~((Word)1));
		auto wcvnExceptFirst = WordConfiguration<Word>::popcount(slice.VN & ~((Word)1));
		assert(slice.scoreEnd == slice.scoreStart + wcvpExceptFirst - wcvnExceptFirst);
#endif
		slice.scoreBeforeStart = slice.scoreStart;
		if (slice.VP & 1) slice.scoreBeforeStart -= 1;
		if (slice.VN & 1) slice.scoreBeforeStart += 1;

		return slice;
	}

	class NodeCalculationResult
	{
	public:
		ScoreType minScore;
		LengthType minScoreIndex;
		size_t cellsProcessed;
	};

	bool firstZeroForced(const std::vector<bool>& previousBand, LengthType neighborNodeIndex, LengthType neighborColumn, Word currentEq, const std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice) const
	{
		if (previousBand[neighborNodeIndex])
		{
			if (currentSlice[neighborColumn].scoreStart < previousSlice[neighborColumn].scoreEnd)
			{
				return true;
			}
			if (currentSlice[neighborColumn].scoreStart == previousSlice[neighborColumn].scoreEnd && !(currentEq & 1))
			{
				return true;
			}
			return false;
		}
		else
		{
			return true;
		}
		assert(false);
	}

	NodeCalculationResult calculateNode(size_t i, size_t j, Word BA, Word BT, Word BC, Word BG, std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand, bool forceSource) const
	{
		NodeCalculationResult result;
		result.minScore = std::numeric_limits<ScoreType>::max();
		result.minScoreIndex = 0;
		result.cellsProcessed = 0;

		LengthType start = graph.nodeStart[i];
		if (forceSource || isSource(i, currentBand))
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
			currentSlice[start] = getNodeStartSlice(Eq, i, previousSlice, currentSlice, currentBand, previousBand);
			if (previousBand[i] && currentSlice[start].scoreBeforeStart > previousSlice[start].scoreEnd)
			{
				currentSlice[start] = mergeTwoSlices(getSourceSliceFromColumn(start, previousSlice), currentSlice[start]);
			}
			if (currentSlice[start].scoreEnd < result.minScore)
			{
				result.minScore = currentSlice[start].scoreEnd;
				result.minScoreIndex = start;
			}
			start++;
		}

		assert(start == graph.nodeStart[i]+1);

		for (LengthType w = start; w < graph.nodeEnd[i]; w++)
		{
			Word Eq = getEq(BA, BT, BC, BG, w);
			int hin = 0;
			if (previousBand[i]) hin = previousSlice[w].hout;
			bool forceFirstHorizontalPositive = firstZeroForced(previousBand, i, w-1, Eq, currentSlice, previousSlice);

			currentSlice[w] = getNextSlice(Eq, currentSlice[w-1], hin, forceFirstHorizontalPositive);

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
			assert(previousBand[i] || currentSlice[w].scoreBeforeStart == j || currentSlice[w].scoreStart == currentSlice[w-1].scoreStart + 1);
			assert(currentSlice[w].scoreStart >= 0);
			assert(currentSlice[w].scoreEnd >= 0);
			assert(currentSlice[w].scoreStart <= currentSlice[w].scoreEnd + WordConfiguration<Word>::WordSize);
			assert(currentSlice[w].scoreEnd <= currentSlice[w].scoreStart + WordConfiguration<Word>::WordSize);
			assert((currentSlice[w].VP & currentSlice[w].VN) == WordConfiguration<Word>::AllZeros);

			assert(!previousBand[i] || currentSlice[w].scoreBeforeStart <= previousSlice[w].scoreEnd);
			assert(currentSlice[w].scoreBeforeStart >= 0);
			if (currentSlice[w].scoreEnd < result.minScore)
			{
				result.minScore = currentSlice[w].scoreEnd;
				result.minScoreIndex = w;
			}
		}
		result.cellsProcessed = (graph.nodeEnd[i] - graph.nodeStart[i]) * WordConfiguration<Word>::WordSize;
		return result;
	}

	void calculateDoublesliceBackwardsRec(size_t i, size_t j, size_t start, Word BA, Word BT, Word BC, Word BG, LengthType sizeLeft, std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand) const
	{
		assert(sizeLeft > 0);
		bool forceSource = true;
		if (graph.nodeEnd[i]-graph.nodeStart[i] < sizeLeft)
		{
			sizeLeft -= graph.nodeEnd[i]-graph.nodeStart[i];
			for (size_t neighbori = 0; neighbori < graph.inNeighbors[i].size(); neighbori++)
			{
				//todo what happens when two cycles are with 2*w of eachothers?
				//the slice calculation may overwrite a previous correct value with a new wrong value
				//let's assume the cycles are far enough from each others
				//cycling over itself is fine
				assert(graph.inNeighbors[i][neighbori] == start || !graph.notInOrder[graph.inNeighbors[i][neighbori]]);
				auto neighbor = graph.inNeighbors[i][neighbori];
				if (!currentBand[neighbor]) continue;
				calculateDoublesliceBackwardsRec(neighbor, j, start, BA, BT, BC, BG, sizeLeft, currentSlice, previousSlice, currentBand, previousBand);
				forceSource = false;
			}
		}
		calculateNode(i, j, BA, BT, BC, BG, currentSlice, previousSlice, currentBand, previousBand, forceSource);
	}

	void calculateDoublesliceBackwards(size_t j, Word BA, Word BT, Word BC, Word BG, std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand) const
	{
		if (graph.firstInOrder == 0) return;
		for (size_t i = graph.firstInOrder-1; i > 0; i--)
		{
			std::set<size_t> currentCalculables;
			std::vector<size_t> currentCalculableOrder;
			assert(graph.notInOrder[i]);
			if (!currentBand[i]) continue;
			calculateDoublesliceBackwardsRec(i, j, i, BA, BT, BC, BG, WordConfiguration<Word>::WordSize*2, currentSlice, previousSlice, currentBand, previousBand);
		}
	}

	MatrixSlice getBitvectorSliceScoresAndFinalPosition(const std::string& sequence, int dynamicWidth, std::vector<std::vector<bool>>& startBand, LengthType dynamicRowStart) const
	{
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
			calculateDoublesliceBackwards(j, BA, BT, BC, BG, currentSlice, previousSlice, currentBand, previousBand);
			for (size_t i = graph.firstInOrder; i < graph.nodeStart.size(); i++)
			{
				nodeMinScores[i] = std::numeric_limits<ScoreType>::max();
				if (!currentBand[i]) continue;
				auto nodeCalc = calculateNode(i, j, BA, BT, BC, BG, currentSlice, previousSlice, currentBand, previousBand, false);
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
				auto nodeCalc = calculateNode(i, j, BA, BT, BC, BG, currentSlice, previousSlice, currentBand, previousBand, false);
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

	std::tuple<ScoreType, int, std::vector<MatrixPosition>, size_t> getBacktrace(std::string sequence, int dynamicWidth, LengthType dynamicRowStart, std::vector<std::vector<bool>>& startBand) const
	{
		int padding = (WordConfiguration<Word>::WordSize - (sequence.size() % WordConfiguration<Word>::WordSize)) % WordConfiguration<Word>::WordSize;
		for (int i = 0; i < padding; i++)
		{
			sequence += 'N';
		}
		auto slice = getBitvectorSliceScoresAndFinalPosition(sequence, dynamicWidth, startBand, dynamicRowStart);
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
		return std::make_tuple(slice.finalMinScore, 0, backtraceresult, slice.cellsProcessed);
	}


	const AlignmentGraph& graph;
};

#endif