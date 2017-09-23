#ifndef GraphAligner_H
#define GraphAligner_H

#include <chrono>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <unordered_set>
#include <queue>
#include <boost/config.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/rational.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include "SliceRow.h"
#include "SparseBoolMatrix.h"
#include "AlignmentGraph.h"
#include "vg.pb.h"
#include "ThreadReadAssertion.h"
#include "NodeSlice.h"
#include "CommonUtils.h"

void printtime(const char* msg)
{
	static auto time = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
	auto newtime = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
	std::cout << msg << " " << newtime << " (" << (newtime - time) << ")" << std::endl;
	time = newtime;
}

size_t getset(const std::set<size_t>& set, size_t pos)
{
	auto iter = set.begin();
	for (auto i = 0; i < pos; i++)
	{
		++iter;;
	}
	return *iter;
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

#ifndef NDEBUG
thread_local int debugLastRowMinScore;
#endif

template <typename LengthType, typename ScoreType, typename Word>
class GraphAligner
{
private:
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
	typedef std::pair<LengthType, LengthType> MatrixPosition;
	class MatrixSlice
	{
	public:
		std::vector<ScoreType> minScorePerWordSlice;
		std::vector<LengthType> minScoreIndexPerWordSlice;
		std::vector<NodeSlice<WordSlice>> scoreSlices;
		size_t cellsProcessed;
	};
	class TwoDirectionalSplitAlignment
	{
	public:
		ScoreType Score() const
		{
			return forward.minScorePerWordSlice.back() + forward.minScorePerWordSlice.back();
		}
		size_t sequenceSplitIndex;
		MatrixSlice forward;
		MatrixSlice backward;
	};
public:
	class AlignmentResult
	{
	public:
		AlignmentResult()
		{
		}
		AlignmentResult(vg::Alignment alignment, bool alignmentFailed, size_t cellsProcessed, size_t ms) :
		alignment(alignment),
		alignmentFailed(alignmentFailed),
		cellsProcessed(cellsProcessed),
		elapsedMilliseconds(ms)
		{
		}
		vg::Alignment alignment;
		bool alignmentFailed;
		size_t cellsProcessed;
		size_t elapsedMilliseconds;
		size_t firstPartStart;
		size_t firstPartEnd;
		size_t secondPartStart;
		size_t secondPartEnd;
		size_t secondPartFirstIndex;
	};

	GraphAligner(const AlignmentGraph& graph) :
	graph(graph)
	{
	}
	
	AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence, int dynamicWidth, LengthType dynamicRowStart) const
	{
		auto timeStart = std::chrono::system_clock::now();
		assert(graph.finalized);
		auto trace = getBacktraceFullStart(sequence, dynamicWidth);
		auto timeEnd = std::chrono::system_clock::now();
		size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
		//failed alignment, don't output
		if (std::get<0>(trace) == std::numeric_limits<ScoreType>::max()) return emptyAlignment(time, std::get<2>(trace));
		if (std::get<1>(trace).size() == 0) return emptyAlignment(time, std::get<2>(trace));
		auto result = traceToAlignment(seq_id, sequence, std::get<0>(trace), std::get<1>(trace), std::get<2>(trace));
		result.firstPartStart = std::get<1>(trace)[0].second;
		result.firstPartEnd = std::get<1>(trace).back().second;
		result.secondPartStart = 0;
		result.secondPartEnd = 0;
		timeEnd = std::chrono::system_clock::now();
		time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
		result.elapsedMilliseconds = time;
		return result;
	}

	AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence, int dynamicWidth, LengthType dynamicRowStart, const std::vector<std::pair<int, size_t>>& seedHits, int startBandwidth) const
	{
		auto timeStart = std::chrono::system_clock::now();
		assert(graph.finalized);
		assert(seedHits.size() > 0);
		TwoDirectionalSplitAlignment bestAlignment;
		std::pair<int, size_t> bestSeed;
		bool hasAlignment = false;
		for (size_t i = 0; i < seedHits.size(); i++)
		{
			std::cerr << "seed " << i << "/" << seedHits.size() << " " << seedHits[i].first << "," << seedHits[i].second << std::endl;
			auto result = getSplitAlignment(sequence, dynamicWidth, startBandwidth, seedHits[i].first, seedHits[i].second, hasAlignment ? bestAlignment.Score() : sequence.size() * 0.4);
			if (result.Score() > sequence.size() * 0.4) continue;
			if (!hasAlignment)
			{
				bestAlignment = std::move(result);
				bestSeed = seedHits[i];
				hasAlignment = true;
			}
			else
			{
				if (result.Score() < bestAlignment.Score())
				{
					bestAlignment = std::move(result);
					bestSeed = seedHits[i];
				}
			}
		}
		auto timeEnd = std::chrono::system_clock::now();
		size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
		//failed alignment, don't output
		if (!hasAlignment)
		{
			return emptyAlignment(time, 0);
		}
		auto bestTrace = getPiecewiseTracesFromSplit(bestAlignment, sequence);

		auto fwresult = traceToAlignment(seq_id, sequence, std::get<0>(bestTrace.first), std::get<1>(bestTrace.first), 0);
		auto bwresult = traceToAlignment(seq_id, sequence, std::get<0>(bestTrace.second), std::get<1>(bestTrace.second), 0);
		//failed alignment, don't output
		if (fwresult.alignmentFailed && bwresult.alignmentFailed)
		{
			return emptyAlignment(time, 0);
		}
		auto result = mergeAlignments(bwresult, fwresult);
		result.secondPartFirstIndex = bwresult.alignment.path().mapping_size();
		if (std::get<1>(bestTrace.second).size() == 0)
		{
			//this way around. 
			//bestTrace.second is backward and bestTrace.first is forward
			//but result.firstPart is backward and result.secondPart is forward
			result.firstPartStart = 0;
			result.firstPartEnd = 0;
		}
		else
		{
			result.firstPartStart = std::get<1>(bestTrace.second)[0].second;
			result.firstPartEnd = std::get<1>(bestTrace.second).back().second;
		}
		if (std::get<1>(bestTrace.first).size() == 0)
		{
			result.secondPartStart = 0;
			result.secondPartEnd = 0;
		}
		else
		{
			result.secondPartStart = std::get<1>(bestTrace.first)[0].second;
			result.secondPartEnd = std::get<1>(bestTrace.first).back().second;
		}
		timeEnd = std::chrono::system_clock::now();
		time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
		result.elapsedMilliseconds = time;
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

	std::pair<size_t, size_t> getLargestContiguousBlock(const std::vector<bool>& vec) const
	{
		size_t thisBlock = 0;
		size_t maxBlockSize = 0;
		size_t maxBlockEnd = 0;
		for (size_t i = 0; i < vec.size(); i++)
		{
			if (vec[i])
			{
				thisBlock++;
			}
			else
			{
				if (thisBlock > maxBlockSize)
				{
					assert(i > 0);
					assert(i >= thisBlock);
					maxBlockEnd = i-1;
					maxBlockSize = thisBlock - 1;
				}
				thisBlock = 0;
			}
		}
		if (thisBlock > maxBlockSize)
		{
			maxBlockEnd = vec.size()-1;
			maxBlockSize = thisBlock - 1;
		}
		assert(maxBlockEnd >= maxBlockSize);
		return std::make_pair(maxBlockEnd - maxBlockSize, maxBlockEnd);
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

	AlignmentResult emptyAlignment(size_t elapsedMilliseconds, size_t cellsProcessed) const
	{
		vg::Alignment result;
		result.set_score(std::numeric_limits<decltype(result.score())>::max());
		return AlignmentResult { result, true, cellsProcessed, elapsedMilliseconds };
	}

	bool posEqual(const vg::Position& pos1, const vg::Position& pos2) const
	{
		return pos1.node_id() == pos2.node_id() && pos1.is_reverse() == pos2.is_reverse();
	}

	AlignmentResult mergeAlignments(const AlignmentResult& first, const AlignmentResult& second) const
	{
		assert(!first.alignmentFailed || !second.alignmentFailed);
		if (first.alignmentFailed) return second;
		if (second.alignmentFailed) return first;
		assert(!first.alignmentFailed);
		assert(!second.alignmentFailed);
		AlignmentResult finalResult;
		finalResult.alignmentFailed = false;
		finalResult.cellsProcessed = first.cellsProcessed + second.cellsProcessed;
		finalResult.elapsedMilliseconds = first.elapsedMilliseconds + second.elapsedMilliseconds;
		finalResult.alignment = first.alignment;
		finalResult.alignment.set_score(first.alignment.score() + second.alignment.score());
		int start = 0;
		auto firstEndPos = first.alignment.path().mapping(first.alignment.path().mapping_size()-1).position();;
		auto secondStartPos = second.alignment.path().mapping(0).position();
		if (posEqual(firstEndPos, secondStartPos))
		{
			start = 1;
		}
		else if (graph.outNeighbors[graph.nodeLookup.at(firstEndPos.node_id())].count(graph.nodeLookup.at(secondStartPos.node_id())) == 1)
		{
			start = 0;
		}
		else
		{
			std::cerr << "Piecewise alignments can't be merged!";
			std::cerr << " first end: " << firstEndPos.node_id() << " " << (firstEndPos.is_reverse() ? "-" : "+");
			std::cerr << " second start: " << secondStartPos.node_id() << " " << (secondStartPos.is_reverse() ? "-" : "+") << std::endl;
		}
		for (int i = start; i < second.alignment.path().mapping_size(); i++)
		{
			auto mapping = finalResult.alignment.mutable_path()->add_mapping();
			*mapping = second.alignment.path().mapping(i);
		}
		return finalResult;
	}

	AlignmentResult traceToAlignment(const std::string& seq_id, const std::string& sequence, ScoreType score, const std::vector<MatrixPosition>& trace, size_t cellsProcessed) const
	{
		vg::Alignment result;
		result.set_name(seq_id);
		result.set_score(score);
		result.set_sequence(sequence);
		auto path = new vg::Path;
		result.set_allocated_path(path);
		if (trace.size() == 0) return AlignmentResult { result, true, cellsProcessed, std::numeric_limits<size_t>::max() };
		size_t pos = 0;
		size_t oldNode = graph.indexToNode[trace[0].first];
		while (oldNode == graph.dummyNodeStart)
		{
			pos++;
			if (pos == trace.size()) return emptyAlignment(std::numeric_limits<size_t>::max(), cellsProcessed);
			assert(pos < trace.size());
			assert(trace[pos].second >= trace[pos-1].second);
			oldNode = graph.indexToNode[trace[pos].first];
			assert(oldNode < graph.nodeIDs.size());
		}
		if (oldNode == graph.dummyNodeEnd) return emptyAlignment(std::numeric_limits<size_t>::max(), cellsProcessed);
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
			assert(trace[pos].second >= trace[pos-1].second);
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
		return AlignmentResult { result, false, cellsProcessed, std::numeric_limits<size_t>::max() };
	}

	class NodePosWithDistance
	{
	public:
		NodePosWithDistance(LengthType node, bool end, int distance) : node(node), end(end), distance(distance) {};
		LengthType node;
		bool end;
		int distance;
		bool operator<(const NodePosWithDistance& other) const
		{
			return distance < other.distance;
		}
		bool operator>(const NodePosWithDistance& other) const
		{
			return distance > other.distance;
		}
	};

	template <typename Container>
	void expandBandFromPositions(const Container& startpositions, LengthType dynamicWidth, std::unordered_map<size_t, size_t>& distanceAtNodeStart, std::unordered_map<size_t, size_t>& distanceAtNodeEnd, std::set<size_t>& bandOrder, std::set<size_t>& bandOrderOutOfOrder) const
	{
		std::priority_queue<NodePosWithDistance, std::vector<NodePosWithDistance>, std::greater<NodePosWithDistance>> queue;
		for (auto startpos : startpositions)
		{
			auto nodeIndex = graph.indexToNode[startpos];
			if (nodeIndex < graph.firstInOrder)
			{
				bandOrderOutOfOrder.insert(nodeIndex);
			}
			else if (nodeIndex >= graph.firstInOrder)
			{
				bandOrder.insert(nodeIndex);
			}
			auto start = graph.nodeStart[nodeIndex];
			auto end = graph.nodeEnd[nodeIndex];
			assert(end > startpos);
			assert(startpos >= start);
			assert(startpos - start >= 0);
			assert(end - startpos - 1 >= 0);
			queue.emplace(nodeIndex, false, startpos - start);
			queue.emplace(nodeIndex, true, end - startpos - 1);
		}
		int oldDistance = 0;
		while (queue.size() > 0)
		{
			NodePosWithDistance top = queue.top();
			assert(top.distance >= oldDistance);
			oldDistance = top.distance;
			assert(top.node < graph.nodeStart.size());
			queue.pop();
			assert(top.node < graph.nodeStart.size());
			if (top.distance > dynamicWidth) continue;
			if (top.end)
			{
				auto found = distanceAtNodeEnd.find(top.node);
				if (found != distanceAtNodeEnd.end() && found->second <= top.distance) continue;
				distanceAtNodeEnd[top.node] = top.distance;
			}
			else
			{
				auto found = distanceAtNodeStart.find(top.node);
				if (found != distanceAtNodeStart.end() && found->second <= top.distance) continue;
				distanceAtNodeStart[top.node] = top.distance;
			}
			auto nodeIndex = top.node;
			if (nodeIndex < graph.firstInOrder)
			{
				bandOrderOutOfOrder.insert(nodeIndex);
			}
			else if (nodeIndex >= graph.firstInOrder)
			{
				bandOrder.insert(nodeIndex);
			}
			assert(nodeIndex < graph.nodeEnd.size());
			assert(nodeIndex < graph.nodeStart.size());
			auto size = graph.nodeEnd[nodeIndex] - graph.nodeStart[nodeIndex];
			if (top.end)
			{
				assert(top.distance + size - 1 >= top.distance);
				queue.emplace(nodeIndex, false, top.distance + size - 1);
				assert(nodeIndex < graph.outNeighbors.size());
				for (auto neighbor : graph.outNeighbors[nodeIndex])
				{
					assert(top.distance + 1 >= top.distance);
					queue.emplace(neighbor, false, top.distance + 1);
				}
			}
			else
			{
				assert(top.distance + size - 1 >= top.distance);
				queue.emplace(nodeIndex, true, top.distance + size - 1);
				assert(nodeIndex < graph.inNeighbors.size());
				for (auto neighbor : graph.inNeighbors[nodeIndex])
				{
					assert(top.distance + 1 >= top.distance);
					queue.emplace(neighbor, true, top.distance + 1);
				}
			}
		}
	}

	std::pair<std::set<LengthType>, std::set<LengthType>> projectForwardAndExpandBand(LengthType previousMinimumIndex, LengthType dynamicWidth) const
	{
		std::set<LengthType> bandOrder;
		std::set<LengthType> bandOrderOutOfOrder;
		assert(previousMinimumIndex < graph.nodeSequences.size());
		auto nodeIndex = graph.indexToNode[previousMinimumIndex];
		auto end = graph.nodeEnd[nodeIndex];
		std::set<size_t> positions;
		positions.insert(previousMinimumIndex);
		positions = graph.ProjectForward(positions, WordConfiguration<Word>::WordSize);
		positions.insert(previousMinimumIndex);
		assert(positions.size() >= 1);
		if (nodeIndex < graph.firstInOrder)
		{
			bandOrderOutOfOrder.insert(nodeIndex);
		}
		else if (nodeIndex >= graph.firstInOrder)
		{
			bandOrder.insert(nodeIndex);
		}
		std::unordered_map<size_t, size_t> distanceAtNodeEnd;
		std::unordered_map<size_t, size_t> distanceAtNodeStart;
		expandBandFromPositions(positions, dynamicWidth, distanceAtNodeStart, distanceAtNodeEnd, bandOrder, bandOrderOutOfOrder);
		return std::make_pair(bandOrder, bandOrderOutOfOrder);
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

#ifdef EXTRABITVECTORASSERTIONS

	WordSlice getWordSliceCellByCell(size_t j, size_t w, const std::string& sequence, const NodeSlice<WordSlice>& currentSlice, const NodeSlice<WordSlice>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand) const
	{
		const auto lastBitMask = ((Word)1) << (WordConfiguration<Word>::WordSize-1);
		WordSlice result;
		auto nodeIndex = graph.indexToNode[w];
		assert(currentBand[nodeIndex]);
		const std::vector<WordSlice>& oldNode = previousBand[nodeIndex] ? previousSlice.node(nodeIndex) : currentSlice.node(nodeIndex);
		assert(currentBand[nodeIndex]);
		ScoreType current[66];
		current[0] = j+1;
		current[1] = j;
		if (j > 0 && previousBand[nodeIndex]) current[1] = std::min(current[1], oldNode[w-graph.nodeStart[nodeIndex]].scoreEnd);
		if (j > 0 && previousBand[nodeIndex]) current[0] = std::min(current[0], oldNode[w-graph.nodeStart[nodeIndex]].scoreEnd - ((oldNode[w-graph.nodeStart[nodeIndex]].VP & lastBitMask) ? 1 : 0) + ((oldNode[w-graph.nodeStart[nodeIndex]].VN & lastBitMask) ? 1 : 0));
		for (int i = 1; i < 65; i++)
		{
			current[i+1] = current[i]+1;
		}
		if (w == graph.nodeStart[nodeIndex])
		{
			for (auto neighbor : graph.inNeighbors[nodeIndex])
			{
				if (!previousBand[neighbor] && !currentBand[neighbor]) continue;
				const std::vector<WordSlice>& neighborSlice = currentBand[neighbor] ? currentSlice.node(neighbor) : previousSlice.node(neighbor);
				const std::vector<WordSlice>& oldNeighborSlice = previousBand[neighbor] ? previousSlice.node(neighbor) : currentSlice.node(neighbor);
				auto u = graph.nodeEnd[neighbor]-1;
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
			const std::vector<WordSlice>& slice = currentSlice.node(nodeIndex);
			const std::vector<WordSlice>& oldSlice = previousBand[nodeIndex] ? previousSlice.node(nodeIndex) : slice;
			auto u = w-1;
			ScoreType previous[66];
			previous[0] = slice[u-graph.nodeStart[nodeIndex]].scoreBeforeStart+1;
			previous[1] = slice[u-graph.nodeStart[nodeIndex]].scoreBeforeStart;
			if (previousBand[nodeIndex]) previous[0] = std::min(previous[0], oldSlice[u-graph.nodeStart[nodeIndex]].scoreEnd - ((oldSlice[u-graph.nodeStart[nodeIndex]].VP & lastBitMask) ? 1 : 0) + ((oldSlice[u-graph.nodeStart[nodeIndex]].VN & lastBitMask) ? 1 : 0));
			if (previousBand[nodeIndex]) previous[1] = std::min(previous[1], oldSlice[u-graph.nodeStart[nodeIndex]].scoreEnd);
			for (int i = 1; i < 65; i++)
			{
				previous[i+1] = previous[i];
				previous[i+1] += (slice[u-graph.nodeStart[nodeIndex]].VP & (((Word)1) << (i-1)) ? 1 : 0);
				previous[i+1] -= (slice[u-graph.nodeStart[nodeIndex]].VN & (((Word)1) << (i-1)) ? 1 : 0);
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

#ifdef EXTRABITVECTORASSERTIONS
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
#ifdef EXTRABITVECTORASSERTIONS
		auto correctValue = differenceMasksCellByCell(leftVP, leftVN, rightVP, rightVN, scoreDifference);
#endif
		assert(scoreDifference >= 0);
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
		if (scoreDifference == 128 && rightVN == allones && leftVP == allones)
		{
			return std::make_pair(allones ^ ((Word)1 << (WordConfiguration<Word>::WordSize-1)), allzeros);
		}
		else if (scoreDifference == 0 && rightVN == allones && leftVP == allones)
		{
			return std::make_pair(0, allones);
		}
		assert(scoreDifference >= 0);
		assert(scoreDifference < 128);
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
#ifdef EXTRABITVECTORASSERTIONS
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
#ifdef EXTRABITVECTORASSERTIONS
		auto correctValue = mergeTwoSlicesCellByCell(left, right);
#endif
		if (left.scoreBeforeStart > right.scoreBeforeStart) std::swap(left, right);
		WordSlice result;
		assert((left.VP & left.VN) == WordConfiguration<Word>::AllZeros);
		assert((right.VP & right.VN) == WordConfiguration<Word>::AllZeros);
		auto masks = differenceMasks(left.VP, left.VN, right.VP, right.VN, right.scoreBeforeStart - left.scoreBeforeStart);
		auto leftSmaller = masks.first;
		auto rightSmaller = masks.second;
		assert((leftSmaller & rightSmaller) == 0);
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
#ifdef EXTRABITVECTORASSERTIONS
		assert(result.VP == correctValue.VP);
		assert(result.VN == correctValue.VN);
		assert(result.scoreBeforeStart == correctValue.scoreBeforeStart);
		assert(result.scoreEnd == correctValue.scoreEnd);
#endif
		return result;
	}

#ifdef EXTRABITVECTORASSERTIONS
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

	WordSlice getNodeStartSlice(Word Eq, size_t nodeIndex, const NodeSlice<WordSlice>& previousSlice, const NodeSlice<WordSlice>& currentSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand, bool previousEq) const
	{
		//todo optimization: 10% time inclusive 4% exclusive. can this be improved?
		WordSlice previous;
		WordSlice previousUp;
		bool foundOne = false;
		bool foundOneUp = false;
		bool hasRealNeighbor = false;
		for (auto neighbor : graph.inNeighbors[nodeIndex])
		{
			if (currentBand[neighbor] && previousBand[neighbor]) assertSliceCorrectness(currentSlice.node(neighbor).back(), previousSlice.node(neighbor).back(), previousBand[neighbor]);
			if (previousBand[neighbor])
			{
				if (!foundOneUp)
				{
					previousUp = previousSlice.node(neighbor).back();
					foundOneUp = true;
				}
				else
				{
					auto competitor = previousSlice.node(neighbor).back();
					previousUp = mergeTwoSlices(previousUp, competitor);
				}
			}
			if (previousBand[neighbor] && !currentBand[neighbor])
			{
				if (!foundOne)
				{
					previous = getSourceSliceFromScore(previousSlice.node(neighbor).back().scoreEnd);
					foundOne = true;
				}
				else
				{
					auto competitor = getSourceSliceFromScore(previousSlice.node(neighbor).back().scoreEnd);
					previous = mergeTwoSlices(previous, competitor);
				}
			}
			if (!currentBand[neighbor]) continue;
			if (!foundOne)
			{
				previous = currentSlice.node(neighbor).back();
				foundOne = true;
				hasRealNeighbor = true;
			}
			else
			{
				auto competitor = currentSlice.node(neighbor).back();
				previous = mergeTwoSlices(previous, competitor);
				hasRealNeighbor = true;
			}
		}
		assert(foundOne);
		assertSliceCorrectness(previous, previousUp, foundOneUp);
		if (!hasRealNeighbor) Eq &= 1;
		auto result = getNextSlice(Eq, previous, previousBand[nodeIndex] && foundOneUp, foundOneUp, previousEq, previousUp);
		return result;
	}

	WordSlice getSourceSliceWithoutBefore(size_t row) const
	{
		return { WordConfiguration<Word>::AllOnes & ~((Word)1), WordConfiguration<Word>::AllZeros, row+WordConfiguration<Word>::WordSize, row+1 };
	}

	WordSlice getSourceSliceFromScore(ScoreType previousScore) const
	{
		return { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, previousScore+WordConfiguration<Word>::WordSize, previousScore };
	}

	WordSlice getSourceSlice(size_t nodeIndex, const NodeSlice<WordSlice>& previousSlice) const
	{
		return getSourceSliceFromScore(previousSlice.node(nodeIndex)[0].scoreEnd);
	}

	bool isSource(size_t nodeIndex, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand) const
	{
		for (auto neighbor : graph.inNeighbors[nodeIndex])
		{
			if (currentBand[neighbor]) return false;
			if (previousBand[neighbor]) return false;
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

	WordSlice getNextSlice(Word Eq, WordSlice slice, bool previousInsideBand, bool diagonalInsideBand, bool previousEq, WordSlice previous) const
	{
		//optimization: 13% of time. probably can't be improved easily.
		//http://www.gersteinlab.org/courses/452/09-spring/pdf/Myers.pdf
		//pages 405 and 408

		auto oldValue = slice.scoreBeforeStart;
		if (!diagonalInsideBand)
		{
			Eq &= ~((Word)1);
		}
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
		assert(slice.scoreBeforeStart >= debugLastRowMinScore);
		assert(slice.scoreEnd >= debugLastRowMinScore);
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

	bool firstZeroForced(const std::vector<bool>& previousBand, const std::vector<bool>& currentBand, LengthType neighborNodeIndex, WordSlice neighborSlice, Word currentEq) const
	{
		if (previousBand[neighborNodeIndex] && currentBand[neighborNodeIndex])
		{
			if (neighborSlice.VN & 1)
			{
				return true;
			}
			if (!(neighborSlice.VP & 1) && !(neighborSlice.VN & 1) && !(currentEq & 1))
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

	void assertSliceCorrectness(const WordSlice& current, const WordSlice& up, bool previousBand) const
	{
#ifndef NDEBUG
			auto wcvp = WordConfiguration<Word>::popcount(current.VP);
			auto wcvn = WordConfiguration<Word>::popcount(current.VN);
			assert(current.scoreEnd == current.scoreBeforeStart + wcvp - wcvn);

			assert(current.scoreBeforeStart >= 0);
			assert(current.scoreEnd >= 0);
			assert(current.scoreBeforeStart <= current.scoreEnd + WordConfiguration<Word>::WordSize);
			assert(current.scoreEnd <= current.scoreBeforeStart + WordConfiguration<Word>::WordSize);
			assert((current.VP & current.VN) == WordConfiguration<Word>::AllZeros);

			assert(!previousBand || current.scoreBeforeStart <= up.scoreEnd);
			assert(current.scoreBeforeStart >= 0);
			assert(current.scoreEnd >= debugLastRowMinScore);
			assert(current.scoreBeforeStart >= debugLastRowMinScore);
#endif
	}

	NodeCalculationResult calculateNode(size_t i, size_t j, const std::string& sequence, Word BA, Word BT, Word BC, Word BG, NodeSlice<WordSlice>& currentSlice, const NodeSlice<WordSlice>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand) const
	{
		//todo optimization: 42% inclusive 15% exclusive. can this be improved?
		NodeCalculationResult result;
		result.minScore = std::numeric_limits<ScoreType>::max();
		result.minScoreIndex = 0;
		result.cellsProcessed = 0;
		std::vector<WordSlice>& slice = currentSlice.node(i);
		const std::vector<WordSlice>& oldSlice = previousBand[i] ? previousSlice.node(i) : slice;
		assert(slice.size() == graph.nodeEnd[i] - graph.nodeStart[i]);
		auto nodeStart = graph.nodeStart[i];

#ifdef EXTRABITVECTORASSERTIONS
		WordSlice correctstart;
		correctstart = getWordSliceCellByCell(j, graph.nodeStart[i], sequence, currentSlice, previousSlice, currentBand, previousBand);
#endif

		if (isSource(i, currentBand, previousBand))
		{
			if (previousBand[i])
			{
				slice[0] = getSourceSlice(i, previousSlice);
			}
			else
			{
				slice[0] = getSourceSliceWithoutBefore(sequence.size());
			}
			if (slice[0].scoreEnd < result.minScore)
			{
				result.minScore = slice[0].scoreEnd;
				result.minScoreIndex = nodeStart;
			}
			assertSliceCorrectness(slice[0], oldSlice[0], previousBand[i]);
		}
		else
		{
			Word Eq = getEq(BA, BT, BC, BG, nodeStart);
			slice[0] = getNodeStartSlice(Eq, i, previousSlice, currentSlice, currentBand, previousBand, j > 0 && graph.nodeSequences[graph.nodeStart[i]] == sequence[j-1]);
			if (previousBand[i] && slice[0].scoreBeforeStart > oldSlice[0].scoreEnd)
			{
				slice[0] = mergeTwoSlices(getSourceSliceFromScore(oldSlice[0].scoreEnd), slice[0]);
			}
			// if (slice[0].scoreBeforeStart > j)
			// {
			// 	slice[0] = mergeTwoSlices(getSourceSliceWithoutBefore(j), slice[0]);
			// }
			if (slice[0].scoreEnd < result.minScore)
			{
				result.minScore = slice[0].scoreEnd;
				result.minScoreIndex = nodeStart;
			}
			assertSliceCorrectness(slice[0], oldSlice[0], previousBand[i]);
			//note: currentSlice[start].score - optimalInNeighborEndScore IS NOT within {-1, 0, 1} always because of the band
		}

#ifdef EXTRABITVECTORASSERTIONS
		assert(slice[0].scoreBeforeStart == correctstart.scoreBeforeStart);
		assert(slice[0].scoreEnd == correctstart.scoreEnd);
		assert(slice[0].VP == correctstart.VP);
		assert(slice[0].VN == correctstart.VN);
#endif

		for (LengthType w = 1; w < graph.nodeEnd[i] - graph.nodeStart[i]; w++)
		{
			Word Eq = getEq(BA, BT, BC, BG, nodeStart+w);
			bool forceFirstHorizontalPositive = firstZeroForced(previousBand, currentBand, i, slice[w-1], Eq);

			slice[w] = getNextSlice(Eq, slice[w-1], previousBand[i], previousBand[i], j > 0 && graph.nodeSequences[nodeStart+w] == sequence[j-1], oldSlice[w-1]);

			if (previousBand[i] && slice[w].scoreBeforeStart > oldSlice[w].scoreEnd)
			{
				slice[w] = mergeTwoSlices(getSourceSliceFromScore(oldSlice[w].scoreEnd), slice[w]);
			}
			// if (slice[w].scoreBeforeStart > j)
			// {
			// 	slice[w] = mergeTwoSlices(getSourceSliceWithoutBefore(j), slice[w]);
			// }

			assert(previousBand[i] || slice[w].scoreBeforeStart == j || slice[w].scoreBeforeStart == slice[w-1].scoreBeforeStart + 1);
			assertSliceCorrectness(slice[w], oldSlice[w], previousBand[i]);

			if (slice[w].scoreEnd <= result.minScore)
			{
				result.minScore = slice[w].scoreEnd;
				result.minScoreIndex = nodeStart + w;
			}

#ifdef EXTRABITVECTORASSERTIONS
			auto correctslice = getWordSliceCellByCell(j, nodeStart+w, sequence, currentSlice, previousSlice, currentBand, previousBand);
			assert(slice[w].scoreBeforeStart == correctslice.scoreBeforeStart);
			assert(slice[w].scoreEnd == correctslice.scoreEnd);
			assert(slice[w].VP == correctslice.VP);
			assert(slice[w].VN == correctslice.VN);
#endif
		}
		result.cellsProcessed = (graph.nodeEnd[i] - graph.nodeStart[i]) * WordConfiguration<Word>::WordSize;
		return result;
	}

	void getCycleCutReachability(size_t j, size_t cycleCut, size_t index, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand, std::vector<bool>& reachable, std::vector<bool>& source) const
	{
		assert(index < reachable.size());
		if (reachable[index]) return;
		reachable[index] = true;
		// assert(graph.notInOrder[cycleCut]);
		assert(currentBand[graph.cuts[cycleCut].nodes[index]]);
		if (graph.cuts[cycleCut].previousCut[index]) return;
		source[index] = true;
		for (auto otherIndex : graph.cuts[cycleCut].predecessors[index])
		{
			assert(otherIndex > index);
			if (previousBand[graph.cuts[cycleCut].nodes[otherIndex]])
			{
				source[index] = false;
			}
			if (currentBand[graph.cuts[cycleCut].nodes[otherIndex]])
			{
				getCycleCutReachability(j, cycleCut, otherIndex, currentBand, previousBand, reachable, source);
				source[index] = false;
			}
		}
	}

	std::unordered_map<LengthType, WordSlice> getCycleCutEndValues(size_t j, const std::string& sequence, Word BA, Word BT, Word BC, Word BG, const NodeSlice<WordSlice>& previousSlice, const std::vector<bool>& currentBand, const std::vector<bool>& previousBand, const std::set<size_t>& bandOrder, const std::set<size_t>& bandOrderOutOfOrder) const
	{
		if (graph.firstInOrder == 0) return {};
		NodeSlice<WordSlice> currentSlice;
		for (auto node : bandOrder)
		{
			currentSlice.addNode(node, graph.nodeEnd[node] - graph.nodeStart[node]);
			auto& slice = currentSlice.node(node);
			if (previousBand[node])
			{
				slice.back() = getSourceSliceFromScore(previousSlice.node(node).back().scoreEnd);
			}
			else
			{
				slice.back() = getSourceSliceWithoutBefore(sequence.size());
			}
		}
		for (auto node : bandOrderOutOfOrder)
		{
			currentSlice.addNode(node, graph.nodeEnd[node] - graph.nodeStart[node]);
			auto& slice = currentSlice.node(node);
			if (previousBand[node])
			{
				slice.back() = getSourceSliceFromScore(previousSlice.node(node).back().scoreEnd);
			}
			else
			{
				slice.back() = getSourceSliceWithoutBefore(sequence.size());
			}
		}
		//if there are cycles within 2*w of eachothers, calculating a latter slice may overwrite the earlier slice's value
		//store the correct values here and then merge them at the end
		std::unordered_map<size_t, WordSlice> correctEndValues;
		for (auto order : bandOrderOutOfOrder)
		{
			correctEndValues[order] = {WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::AllZeros, sequence.size(), sequence.size()};
		}
		for (auto i : bandOrderOutOfOrder)
		{
			if (i == 0) continue;
#ifndef NDEBUG
			std::vector<size_t> debugLengthAtNodeEnd;
			debugLengthAtNodeEnd.resize(graph.cuts[i].nodes.size(), 0);
#endif
			assert(currentBand[i]);
			assert(i > 0);
			assert(i < graph.firstInOrder);
			// assert(graph.notInOrder[i]);
			assert(graph.cuts[i].nodes.size() > 0);
			assert(graph.cuts[i].nodes[0] == i);
			std::vector<bool> reachable;
			std::vector<bool> source;
			reachable.resize(graph.cuts[i].nodes.size(), false);
			source.resize(graph.cuts[i].nodes.size(), false);
			getCycleCutReachability(j, i, 0, currentBand, previousBand, reachable, source);
			for (size_t index = graph.cuts[i].nodes.size()-1; index < graph.cuts[i].nodes.size(); index--)
			{
				if (!reachable[index]) continue;
				if (graph.cuts[i].previousCut[index])
				{
					assert(correctEndValues.find(graph.cuts[i].nodes[index]) != correctEndValues.end());
					assert(correctEndValues[graph.cuts[i].nodes[index]].scoreBeforeStart != sequence.size());
					currentSlice.node(graph.cuts[i].nodes[index]).back() = correctEndValues[graph.cuts[i].nodes[index]];
					if (previousBand[graph.cuts[i].nodes[index]]) assertSliceCorrectness(currentSlice.node(graph.cuts[i].nodes[index]).back(), previousSlice.node(graph.cuts[i].nodes[index]).back(), previousBand[graph.cuts[i].nodes[index]]);
#ifndef NDEBUG
					debugLengthAtNodeEnd[index] = WordConfiguration<Word>::WordSize*2;
#endif
				}
				else
				{
					calculateNode(graph.cuts[i].nodes[index], j, sequence, BA, BT, BC, BG, currentSlice, previousSlice, currentBand, previousBand);
#ifndef NDEBUG
					size_t debugSmallestPredecessor = WordConfiguration<Word>::WordSize * 2;
					bool debugFoundPredecessor = false;
					for (auto neighbor : graph.cuts[i].predecessors[index])
					{
						if (!reachable[neighbor]) continue;
						debugSmallestPredecessor = std::min(debugSmallestPredecessor, debugLengthAtNodeEnd[neighbor]);
						debugFoundPredecessor = true;
					}
					if (!debugFoundPredecessor) debugSmallestPredecessor = 0;
					debugLengthAtNodeEnd[index] = debugSmallestPredecessor + graph.nodeEnd[graph.cuts[i].nodes[index]] - graph.nodeStart[graph.cuts[i].nodes[index]];
#endif
					if (previousBand[graph.cuts[i].nodes[index]]) assertSliceCorrectness(currentSlice.node(graph.cuts[i].nodes[index]).back(), previousSlice.node(graph.cuts[i].nodes[index]).back(), previousBand[graph.cuts[i].nodes[index]]);
				}
			}
			assert(debugLengthAtNodeEnd[0] >= WordConfiguration<Word>::WordSize * 2);
			correctEndValues[i] = currentSlice.node(i).back();
			assert(graph.cuts[i].nodes[0] == i);
			for (size_t index = 1; index < graph.cuts[i].nodes.size(); index++)
			{
				auto node = graph.cuts[i].nodes[index];
				if (!currentBand[node]) continue;
				if (previousBand[node])
				{
					currentSlice.node(node).back() = getSourceSliceFromScore(previousSlice.node(node).back().scoreEnd);
				}
				else
				{
					currentSlice.node(node).back() = getSourceSliceWithoutBefore(sequence.size());
				}
			}
			currentSlice.node(i).back() = correctEndValues[i];
			if (previousBand[i]) assertSliceCorrectness(currentSlice.node(i).back(), previousSlice.node(i).back(), previousBand[i]);
		}
		return correctEndValues;
	}

	void getBandOrder(const std::vector<bool>& currentBand, std::set<size_t>& bandOrder, std::set<size_t>& bandOrderOutOfOrder) const
	{
		assert(currentBand.size() == graph.notInOrder.size());
		for (size_t i = 0; i < graph.firstInOrder; i++)
		{
			if (currentBand[i]) bandOrderOutOfOrder.insert(i);
		}
		for (size_t i = graph.firstInOrder; i < currentBand.size(); i++)
		{
			if (currentBand[i]) bandOrder.insert(i);
		}
	}

	std::tuple<std::set<LengthType>, std::set<LengthType>> getBandOrderFromSlice(const NodeSlice<WordSlice>& slice) const
	{
		std::set<LengthType> bandOrder;
		std::set<LengthType> bandOrderOutOfOrder;
		for (const auto& node : slice)
		{
			if (node.first < graph.firstInOrder)
			{
				bandOrderOutOfOrder.insert(node.first);
			}
			else
			{
				bandOrder.insert(node.first);
			}
		}
		return std::make_tuple(bandOrder, bandOrderOutOfOrder);
	}

	MatrixSlice getBacktraceSlice(const std::string& forwardSequence, const MatrixSlice& forward) const
	{
		assert(forward.scoreSlices.size() > 0);
		auto sequence = CommonUtils::ReverseComplement(forwardSequence);
		NodeSlice<WordSlice> startBand;
		auto minScore = forward.minScorePerWordSlice.back();
		for (const auto& pair : forward.scoreSlices.back())
		{
			auto node = pair.first;
			auto reverseNode = graph.GetReverseNode(node);
			assert(graph.nodeEnd[reverseNode] - graph.nodeStart[reverseNode] == graph.nodeEnd[node] - graph.nodeStart[node]);
			startBand.addNode(reverseNode, graph.nodeEnd[reverseNode] - graph.nodeStart[reverseNode]);
			auto& slice = startBand.node(reverseNode);
			assert(slice.size() == pair.second.size());
			for (size_t i = 0; i < pair.second.size(); i++)
			{
				assert(pair.second[i].scoreEnd >= minScore);
				slice[slice.size()-1-i] = {0, 0, pair.second[i].scoreEnd - minScore, pair.second[i].scoreEnd - minScore};
			}
		}
		assert(forwardSequence.size() % WordConfiguration<Word>::WordSize == 0);
		auto rowBandFunction = [this, &forward](size_t j, size_t previousMinimumIndex, const std::vector<bool>& previousBand)
		{
			std::set<LengthType> bandOrder;
			std::set<LengthType> bandOrderOutOfOrder;
			auto index = forward.scoreSlices.size() - (j / WordConfiguration<Word>::WordSize) - 1;
			for (const auto& pair : forward.scoreSlices[index])
			{
				auto reverseNode = graph.GetReverseNode(pair.first);
				if (reverseNode < graph.firstInOrder)
				{
					bandOrderOutOfOrder.insert(reverseNode);
				}
				else
				{
					bandOrder.insert(reverseNode);
				}
			}
			return std::make_tuple(bandOrder, bandOrderOutOfOrder);
		};
		auto result = getBitvectorSliceScoresAndFinalPosition(sequence, startBand, forward.minScorePerWordSlice.back()+1, rowBandFunction);
		assert(result.minScorePerWordSlice.size() == forward.minScorePerWordSlice.size());
#ifndef NDEBUG
		for (size_t i = 0; i < result.minScorePerWordSlice.size() - 1; i++)
		{
			//can be lower because the minimum score is not necessarily in the actual backtrace
			assert(result.minScorePerWordSlice[i] + forward.minScorePerWordSlice[result.minScorePerWordSlice.size() - 2 - i] <= forward.minScorePerWordSlice.back());
		}
#endif
		assert(result.minScorePerWordSlice.back() == forward.minScorePerWordSlice.back());
		return result;
	}

	std::tuple<std::set<LengthType>, std::set<LengthType>> defaultForwardRowBandFunction(LengthType j, LengthType previousMinimumIndex, const std::vector<bool>& previousBand, LengthType dynamicWidth, const NodeSlice<WordSlice>& beforeSlice) const
	{
		if (j == 0)
		{
			return getBandOrderFromSlice(beforeSlice);
		}
		assert(previousMinimumIndex != std::numeric_limits<LengthType>::max());
		return projectForwardAndExpandBand(previousMinimumIndex, dynamicWidth);
	}

#ifdef EXTRACORRECTNESSASSERTIONS

	void verifySlice(const std::string& sequence, const MatrixSlice& band, const NodeSlice<WordSlice>& initialSlice, ScoreType maxScore) const
	{
		auto slice = getMatrixSliceCellByCell(sequence, band, initialSlice, maxScore);
		assert(band.minScorePerWordSlice.size() == slice.minScorePerWordSlice.size());
		assert(slice.scoreSlices.size() == band.scoreSlices.size());
		assert(band.scoreSlices.size() == band.minScorePerWordSlice.size());
		for (size_t i = 0; i < band.minScorePerWordSlice.size(); i++)
		{
			auto pos = band.minScoreIndexPerWordSlice[i];
			auto node = graph.indexToNode[pos];
			assert(band.scoreSlices[i].node(node)[pos-graph.nodeStart[node]].scoreEnd == band.minScorePerWordSlice[i]);
			assert(band.minScorePerWordSlice[i] == slice.minScorePerWordSlice[i]);
			// assert(band.minScoreIndexPerWordSlice[i] == slice.minScoreIndexPerWordSlice[i]);
		}
		for (size_t i = 0; i < slice.scoreSlices.size(); i++)
		{
			std::vector<std::pair<size_t, std::vector<WordSlice>>> sortedList;
			sortedList.insert(sortedList.end(), slice.scoreSlices[i].begin(), slice.scoreSlices[i].end());
			std::sort(sortedList.begin(), sortedList.end(), [](auto& left, auto& right) { return left.first < right.first; });
			for (auto pair : sortedList)
			{
				volatile auto node = pair.first; //make sure this won't get optimized away
				assert(band.scoreSlices[i].hasNode(pair.first));
				auto compare = band.scoreSlices[i].node(pair.first);
				assert(compare.size() == pair.second.size());
				for (size_t ii = 0; ii < compare.size(); ii++)
				{
					assert(compare[ii].scoreEnd == pair.second[ii].scoreEnd);
					assert((compare[ii].VP | 3) == (pair.second[ii].VP | 3));
					assert((compare[ii].VN | 3) == (pair.second[ii].VN | 3));
				}
			}
		}
	}

	ScoreType getValue(const MatrixSlice& band, LengthType j, LengthType w) const
	{
		auto node = graph.indexToNode[w];
		auto slice = j / 64;
		auto word = band.scoreSlices[slice].node(node)[w - graph.nodeStart[node]];
		auto off = j % 64;
		auto mask = WordConfiguration<Word>::AllOnes;
		if (off < 63) mask = ~(WordConfiguration<Word>::AllOnes << (off + 1));
		auto value = word.scoreBeforeStart + WordConfiguration<Word>::popcount(word.VP & mask) - WordConfiguration<Word>::popcount(word.VN & mask);
		return value;
	}

	MatrixSlice getMatrixSliceCellByCell(const std::string& sequence, const MatrixSlice& band, const NodeSlice<WordSlice>& initialSlice, ScoreType maxScore) const
	{
		MatrixSlice result;
		result.cellsProcessed = 0;
		std::vector<ScoreType> currentRowScores;
		std::vector<ScoreType> previousRowScores;
		currentRowScores.resize(graph.nodeSequences.size(), sequence.size());
		previousRowScores.resize(graph.nodeSequences.size(), sequence.size());
		for (auto pair : initialSlice)
		{
			auto node = pair.first;
			auto nodestart = graph.nodeStart[node];
			for (size_t i = 0; i < pair.second.size(); i++)
			{
				previousRowScores[nodestart+i] = pair.second[i].scoreEnd;
			}
		}
		std::vector<bool> previousBand;
		std::vector<bool> currentBand;
		previousBand.resize(graph.nodeStart.size(), false);
		currentBand.resize(graph.nodeStart.size(), false);
		std::set<size_t> bandOrder;
		std::set<size_t> bandOrderOutOfOrder;
		std::tie(bandOrder, bandOrderOutOfOrder) = getBandOrderFromSlice(initialSlice);
		for (auto i : bandOrder)
		{
			currentBand[i] = true;
		}
		for (auto i : bandOrderOutOfOrder)
		{
			currentBand[i] = true;
		}
		std::vector<NodeSlice<ScoreType>> sliceScores;
		assert(band.scoreSlices.size() * WordConfiguration<Word>::WordSize == sequence.size());
		for (size_t j = 0; j < sequence.size(); j++)
		{
			if (j % 64 == 0)
			{
				sliceScores.clear();
				sliceScores.resize(64);
				previousBand = std::move(currentBand);
				currentBand.assign(graph.nodeStart.size(), false);
				std::tie(bandOrder, bandOrderOutOfOrder) = getBandOrderFromSlice(band.scoreSlices[j/64]);
				for (auto i : bandOrder)
				{
					currentBand[i] = true;
					for (size_t k = 0; k < 64; k++)
					{
						sliceScores[k].addNode(i, graph.nodeEnd[i] - graph.nodeStart[i]);
					}
				}
				for (auto i : bandOrderOutOfOrder)
				{
					currentBand[i] = true;
					for (size_t k = 0; k < 64; k++)
					{
						sliceScores[k].addNode(i, graph.nodeEnd[i] - graph.nodeStart[i]);
					}
				}
			}
			bool repeat = false;
			bool oneRepeatDone = false;
			ScoreType minscore = std::numeric_limits<ScoreType>::max();
			LengthType minindex = 0;
			do
			{
				repeat = false;
				for (auto node : bandOrder)
				{
					for (auto neighbor : graph.inNeighbors[node])
					{
						auto i = graph.nodeStart[node];
						auto oldscore = currentRowScores[i];
						auto previous = graph.nodeEnd[neighbor]-1;
						currentRowScores[i] = std::min(previousRowScores[i] + 1, currentRowScores[i]);
						currentRowScores[i] = std::min(currentRowScores[previous] + 1, currentRowScores[i]);
						if (graph.nodeSequences[i] == sequence[j] || sequence[j] == 'N')
						{
							currentRowScores[i] = std::min(previousRowScores[previous], currentRowScores[i]);
						}
						else
						{
							currentRowScores[i] = std::min(previousRowScores[previous] + 1, currentRowScores[i]);
						}
						// if (j % 64 == 63) currentRowScores[i] = std::min(currentRowScores[i], (ScoreType)j+1);
						if (currentRowScores[i] < oldscore) repeat = true;
						if (currentRowScores[i] <= minscore)
						{
							minscore = currentRowScores[i];
							minindex = i;
						}
						sliceScores[j%64].node(node)[i-graph.nodeStart[node]] = currentRowScores[i];
					}
					if (oneRepeatDone && currentRowScores[graph.nodeStart[node]] == sequence.size())
					{
						currentRowScores[graph.nodeStart[node]] = j+1;
						repeat = true;
					}
					for (size_t i = graph.nodeStart[node]+1; i < graph.nodeEnd[node]; i++)
					{
						auto oldscore = currentRowScores[i];
						currentRowScores[i] = std::min(previousRowScores[i] + 1, currentRowScores[i]);
						currentRowScores[i] = std::min(currentRowScores[i-1] + 1, currentRowScores[i]);
						if (graph.nodeSequences[i] == sequence[j] || sequence[j] == 'N')
						{
							currentRowScores[i] = std::min(previousRowScores[i-1], currentRowScores[i]);
						}
						else
						{
							currentRowScores[i] = std::min(previousRowScores[i-1] + 1, currentRowScores[i]);
						}
						// if (j % 64 == 63) currentRowScores[i] = std::min(currentRowScores[i], (ScoreType)j+1);
						if (currentRowScores[i] < oldscore) repeat = true;
						if (currentRowScores[i] <= minscore)
						{
							minscore = currentRowScores[i];
							minindex = i;
						}
						sliceScores[j%64].node(node)[i-graph.nodeStart[node]] = currentRowScores[i];
					}
					// if (j % 64 == 63)
					// {
					// 	for (auto neighbor : graph.outNeighbors[node])
					// 	{
					// 		auto i = graph.nodeStart[neighbor];
					// 		auto oldscore = currentRowScores[i];
					// 		auto previous = graph.nodeEnd[node]-1;
					// 		currentRowScores[i] = std::min(previousRowScores[i] + 1, currentRowScores[i]);
					// 		currentRowScores[i] = std::min(currentRowScores[previous] + 1, currentRowScores[i]);
					// 		if (graph.nodeSequences[i] == sequence[j] || sequence[j] == 'N')
					// 		{
					// 			currentRowScores[i] = std::min(previousRowScores[previous], currentRowScores[i]);
					// 		}
					// 		else
					// 		{
					// 			currentRowScores[i] = std::min(previousRowScores[previous] + 1, currentRowScores[i]);
					// 		}
					// 		// if (j % 64 == 63) currentRowScores[i] = std::min(currentRowScores[i], (ScoreType)j+1);
					// 		if (currentRowScores[i] < oldscore) repeat = true;
					// 		if (currentRowScores[i] <= minscore)
					// 		{
					// 			minscore = currentRowScores[i];
					// 			minindex = i;
					// 		}
					// 	}
					// }
				}
				for (auto node : bandOrderOutOfOrder)
				{
					for (auto neighbor : graph.inNeighbors[node])
					{
						auto i = graph.nodeStart[node];
						auto oldscore = currentRowScores[i];
						auto previous = graph.nodeEnd[neighbor]-1;
						currentRowScores[i] = std::min(previousRowScores[i] + 1, currentRowScores[i]);
						currentRowScores[i] = std::min(currentRowScores[previous] + 1, currentRowScores[i]);
						if (graph.nodeSequences[i] == sequence[j] || sequence[j] == 'N')
						{
							currentRowScores[i] = std::min(previousRowScores[previous], currentRowScores[i]);
						}
						else
						{
							currentRowScores[i] = std::min(previousRowScores[previous] + 1, currentRowScores[i]);
						}
						// if (j % 64 == 63) currentRowScores[i] = std::min(currentRowScores[i], (ScoreType)j+1);
						if (currentRowScores[i] < oldscore) repeat = true;
						if (currentRowScores[i] <= minscore)
						{
							minscore = currentRowScores[i];
							minindex = i;
						}
						sliceScores[j%64].node(node)[i-graph.nodeStart[node]] = currentRowScores[i];
					}
					if (oneRepeatDone && currentRowScores[graph.nodeStart[node]] == sequence.size())
					{
						currentRowScores[graph.nodeStart[node]] = j+1;
						repeat = true;
					}
					for (size_t i = graph.nodeStart[node]+1; i < graph.nodeEnd[node]; i++)
					{
						auto oldscore = currentRowScores[i];
						currentRowScores[i] = std::min(previousRowScores[i] + 1, currentRowScores[i]);
						currentRowScores[i] = std::min(currentRowScores[i-1] + 1, currentRowScores[i]);
						if (graph.nodeSequences[i] == sequence[j] || sequence[j] == 'N')
						{
							currentRowScores[i] = std::min(previousRowScores[i-1], currentRowScores[i]);
						}
						else
						{
							currentRowScores[i] = std::min(previousRowScores[i-1] + 1, currentRowScores[i]);
						}
						// if (j % 64 == 63) currentRowScores[i] = std::min(currentRowScores[i], (ScoreType)j+1);
						if (currentRowScores[i] < oldscore) repeat = true;
						if (currentRowScores[i] <= minscore)
						{
							minscore = currentRowScores[i];
							minindex = i;
						}
						sliceScores[j%64].node(node)[i-graph.nodeStart[node]] = currentRowScores[i];
					}
					// if (j % 64 == 63)
					// {
					// 	for (auto neighbor : graph.outNeighbors[node])
					// 	{
					// 		auto i = graph.nodeStart[neighbor];
					// 		auto oldscore = currentRowScores[i];
					// 		auto previous = graph.nodeEnd[node]-1;
					// 		currentRowScores[i] = std::min(previousRowScores[i] + 1, currentRowScores[i]);
					// 		currentRowScores[i] = std::min(currentRowScores[previous] + 1, currentRowScores[i]);
					// 		if (graph.nodeSequences[i] == sequence[j] || sequence[j] == 'N')
					// 		{
					// 			currentRowScores[i] = std::min(previousRowScores[previous], currentRowScores[i]);
					// 		}
					// 		else
					// 		{
					// 			currentRowScores[i] = std::min(previousRowScores[previous] + 1, currentRowScores[i]);
					// 		}
					// 		// if (j % 64 == 63) currentRowScores[i] = std::min(currentRowScores[i], (ScoreType)j+1);
					// 		if (currentRowScores[i] < oldscore) repeat = true;
					// 		if (currentRowScores[i] <= minscore)
					// 		{
					// 			minscore = currentRowScores[i];
					// 			minindex = i;
					// 		}
					// 	}
					// }
				}
				oneRepeatDone = true;
			} while (repeat);
			if (j % 64 == 63)
			{
				result.minScorePerWordSlice.push_back(minscore);
				result.minScoreIndexPerWordSlice.push_back(minindex);
				if (minscore > maxScore)
				{
					for (size_t i = j + 64; i < sequence.size(); i += 64)
					{
						result.minScorePerWordSlice.push_back(maxScore);
						result.minScoreIndexPerWordSlice.push_back(0);
					}
					break;
				}
				result.scoreSlices.emplace_back();
				for (auto node : bandOrder)
				{
					result.scoreSlices.back().addNode(node, graph.nodeEnd[node] - graph.nodeStart[node]);
					auto& slice = result.scoreSlices.back().node(node);
					for (size_t i = graph.nodeStart[node]; i < graph.nodeEnd[node]; i++)
					{
						slice[i - graph.nodeStart[node]] = {0, 0, currentRowScores[i], 0};
						auto oldscore = sliceScores[0].node(node)[i - graph.nodeStart[node]];
						for (size_t k = 1; k < 64; k++)
						{
							auto currentscore = sliceScores[k].node(node)[i-graph.nodeStart[node]];
							if (currentscore > oldscore) slice[i - graph.nodeStart[node]].VP |= ((Word)1) << k;
							if (currentscore < oldscore) slice[i - graph.nodeStart[node]].VN |= ((Word)1) << k;
							oldscore = currentscore;
						}
					}
				}
				for (auto node : bandOrderOutOfOrder)
				{
					result.scoreSlices.back().addNode(node, graph.nodeEnd[node] - graph.nodeStart[node]);
					auto& slice = result.scoreSlices.back().node(node);
					for (size_t i = graph.nodeStart[node]; i < graph.nodeEnd[node]; i++)
					{
						slice[i - graph.nodeStart[node]] = {0, 0, currentRowScores[i], 0};
						auto oldscore = sliceScores[0].node(node)[i - graph.nodeStart[node]];
						for (size_t k = 1; k < 64; k++)
						{
							auto currentscore = sliceScores[k].node(node)[i-graph.nodeStart[node]];
							if (currentscore > oldscore) slice[i - graph.nodeStart[node]].VP |= ((Word)1) << k;
							if (currentscore < oldscore) slice[i - graph.nodeStart[node]].VN |= ((Word)1) << k;
							oldscore = currentscore;
						}
					}
				}
			}
			for (auto node : bandOrder)
			{
				for (size_t i = graph.nodeStart[node]; i < graph.nodeEnd[node]; i++)
				{
					assert(currentRowScores[i] == getValue(band, j, i));
				}
			}
			for (auto node : bandOrderOutOfOrder)
			{
				for (size_t i = graph.nodeStart[node]; i < graph.nodeEnd[node]; i++)
				{
					assert(currentRowScores[i] == getValue(band, j, i));
				}
			}
			if (j % 64 == 0) previousRowScores = currentRowScores;
			for (auto node : bandOrder)
			{
				for (size_t i = graph.nodeStart[node]; i < graph.nodeEnd[node]; i++)
				{
					previousRowScores[i] = currentRowScores[i];
				}
				// if (j % 64 == 63)
				// {
				// 	for (auto neighbor : graph.outNeighbors[node])
				// 	{
				// 		previousRowScores[graph.nodeStart[neighbor]] = currentRowScores[graph.nodeStart[neighbor]];
				// 	}
				// }
			}
			for (auto node : bandOrderOutOfOrder)
			{
				for (size_t i = graph.nodeStart[node]; i < graph.nodeEnd[node]; i++)
				{
					previousRowScores[i] = currentRowScores[i];
				}
				// if (j % 64 == 63)
				// {
				// 	for (auto neighbor : graph.outNeighbors[node])
				// 	{
				// 		previousRowScores[graph.nodeStart[neighbor]] = currentRowScores[graph.nodeStart[neighbor]];
				// 	}
				// }
			}
			for (auto node : bandOrder)
			{
				for (size_t i = graph.nodeStart[node]; i < graph.nodeEnd[node]; i++)
				{
					currentRowScores[i] = sequence.size();
				}
				// if (j % 64 == 63)
				// {
				// 	for (auto neighbor : graph.outNeighbors[node])
				// 	{
				// 		currentRowScores[graph.nodeStart[neighbor]] = sequence.size();
				// 	}
				// }
			}
			for (auto node : bandOrderOutOfOrder)
			{
				for (size_t i = graph.nodeStart[node]; i < graph.nodeEnd[node]; i++)
				{
					currentRowScores[i] = sequence.size();
				}
				// if (j % 64 == 63)
				// {
				// 	for (auto neighbor : graph.outNeighbors[node])
				// 	{
				// 		currentRowScores[graph.nodeStart[neighbor]] = sequence.size();
				// 	}
				// }
			}
		}
		return result;
	}

#endif

	template <typename RowBandFunction>
	MatrixSlice getBitvectorSliceScoresAndFinalPosition(const std::string& sequence, const NodeSlice<WordSlice>& initialSlice, ScoreType maxScore, RowBandFunction rowBandFunction) const
	{
		//todo optimization: 82% inclusive 17% exclusive. can this be improved?
		MatrixSlice result;
		result.cellsProcessed = 0;

		NodeSlice<WordSlice> previousSlice;

		LengthType previousMinimumIndex = std::numeric_limits<LengthType>::max();
		std::vector<bool> currentBand;
		std::vector<bool> previousBand;
		assert(initialSlice.size() > 0);
		currentBand.resize(graph.nodeStart.size(), false);
		previousBand.resize(graph.nodeStart.size(), false);

		std::set<size_t> previousBandOrder;
		std::set<size_t> previousBandOrderOutOfOrder;

#ifndef NDEBUG
		debugLastRowMinScore = 0;
#endif

		for (size_t j = 0; j < sequence.size(); j += WordConfiguration<Word>::WordSize)
		{
			ScoreType currentMinimumScore = std::numeric_limits<ScoreType>::max();
			LengthType currentMinimumIndex = std::numeric_limits<LengthType>::max();
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
					case 'R':
					case 'r':
					BA |= mask;
					BG |= mask;
					break;
					case 'Y':
					case 'y':
					BC |= mask;
					BT |= mask;
					break;
					case 'K':
					case 'k':
					BG |= mask;
					BT |= mask;
					break;
					case 'M':
					case 'm':
					BA |= mask;
					BC |= mask;
					break;
					case 'S':
					case 's':
					BC |= mask;
					BG |= mask;
					break;
					case 'W':
					case 'w':
					BA |= mask;
					BT |= mask;
					break;
					case 'B':
					case 'b':
					BC |= mask;
					BG |= mask;
					BT |= mask;
					break;
					case 'D':
					case 'd':
					BA |= mask;
					BG |= mask;
					BT |= mask;
					break;
					case 'H':
					case 'h':
					BA |= mask;
					BC |= mask;
					BT |= mask;
					break;
					case 'V':
					case 'v':
					BA |= mask;
					BC |= mask;
					BG |= mask;
					break;
					default:
					assert(false);
					break;
				}
			}
			size_t slice = j / WordConfiguration<Word>::WordSize;
			NodeSlice<WordSlice> currentSlice;
			std::set<size_t> bandOrder;
			std::set<size_t> bandOrderOutOfOrder;
			std::swap(currentBand, previousBand);
			if (slice == 0)
			{
				previousSlice = initialSlice;
				std::tie(previousBandOrder, previousBandOrderOutOfOrder) = getBandOrderFromSlice(previousSlice);
				for (auto i : previousBandOrder)
				{
					previousBand[i] = true;
				}
				for (auto i : previousBandOrderOutOfOrder)
				{
					previousBand[i] = true;
				}
			}
			std::tie(bandOrder, bandOrderOutOfOrder) = rowBandFunction(j, previousMinimumIndex, previousBand);
			assert(bandOrder.size() > 0 || bandOrderOutOfOrder.size() > 0);
			for (auto i : bandOrder)
			{
				currentSlice.addNode(i, graph.nodeEnd[i] - graph.nodeStart[i]);
				currentBand[i] = true;
			}
			for (auto i : bandOrderOutOfOrder)
			{
				currentSlice.addNode(i, graph.nodeEnd[i] - graph.nodeStart[i]);
				currentBand[i] = true;
			}
			auto correctCycleCutEndValues = getCycleCutEndValues(j, sequence, BA, BT, BC, BG, previousSlice, currentBand, previousBand, bandOrder, bandOrderOutOfOrder);
			for (auto pair : correctCycleCutEndValues)
			{
				currentSlice.node(pair.first).back() = pair.second;
			}
			for (auto i : bandOrder)
			{
				assert(currentBand[i]);
				auto nodeCalc = calculateNode(i, j, sequence, BA, BT, BC, BG, currentSlice, previousSlice, currentBand, previousBand);
				assert(result.minScorePerWordSlice.size() == 0 || nodeCalc.minScore >= result.minScorePerWordSlice.back());
#ifndef NDEBUG
				auto debugslice = currentSlice.node(i);
				assert(nodeCalc.minScoreIndex >= graph.nodeStart[i]);
				assert(nodeCalc.minScoreIndex < graph.nodeEnd[i]);
				assert(debugslice[nodeCalc.minScoreIndex - graph.nodeStart[i]].scoreEnd == nodeCalc.minScore);
#endif
				if (nodeCalc.minScore < currentMinimumScore)
				{
					currentMinimumScore = nodeCalc.minScore;
					currentMinimumIndex = nodeCalc.minScoreIndex;
				}
				// if (nodeCalc.minScore <= currentMinimumScore)
				// {
				// 	if (nodeCalc.minScoreIndex == graph.nodeEnd[i]-1)
				// 	{
				// 		if (currentSlice.node(i).back().VP & ((Word)1 << (WordConfiguration<Word>::WordSize - 1)))
				// 		{
				// 			for (auto neighbor : graph.outNeighbors[i])
				// 			{
				// 				if (sequence[j + WordConfiguration<Word>::WordSize - 1] == graph.nodeSequences[graph.nodeStart[neighbor]])
				// 				{
				// 					assert(nodeCalc.minScore > 0);
				// 					currentMinimumScore = nodeCalc.minScore - 1;
				// 					currentMinimumIndex = graph.nodeStart[neighbor];
				// 				}
				// 			}
				// 		}
				// 	}
				// }
				result.cellsProcessed += nodeCalc.cellsProcessed;
			}
			for (auto i : bandOrderOutOfOrder)
			{
				assert(currentBand[i]);
				auto nodeCalc = calculateNode(i, j, sequence, BA, BT, BC, BG, currentSlice, previousSlice, currentBand, previousBand);
				assert(result.minScorePerWordSlice.size() == 0 || nodeCalc.minScore >= result.minScorePerWordSlice.back());
#ifndef NDEBUG
				auto debugslice = currentSlice.node(i);
				assert(nodeCalc.minScoreIndex >= graph.nodeStart[i]);
				assert(nodeCalc.minScoreIndex < graph.nodeEnd[i]);
				assert(debugslice[nodeCalc.minScoreIndex - graph.nodeStart[i]].scoreEnd == nodeCalc.minScore);
				auto debugcyclecut = correctCycleCutEndValues[i];
				assert(debugslice.back().VP == debugcyclecut.VP);
				assert(debugslice.back().VN == debugcyclecut.VN);
				assert(debugslice.back().scoreBeforeStart == debugcyclecut.scoreBeforeStart);
				assert(debugslice.back().scoreEnd == debugcyclecut.scoreEnd);
#endif
				if (nodeCalc.minScore < currentMinimumScore)
				{
					currentMinimumScore = nodeCalc.minScore;
					currentMinimumIndex = nodeCalc.minScoreIndex;
				}
				// if (nodeCalc.minScore <= currentMinimumScore)
				// {
				// 	if (nodeCalc.minScoreIndex == graph.nodeEnd[i]-1)
				// 	{
				// 		if (currentSlice.node(i).back().VP & ((Word)1 << (WordConfiguration<Word>::WordSize - 1)))
				// 		{
				// 			for (auto neighbor : graph.outNeighbors[i])
				// 			{
				// 				if (sequence[j + WordConfiguration<Word>::WordSize - 1] == graph.nodeSequences[graph.nodeStart[neighbor]])
				// 				{
				// 					assert(nodeCalc.minScore > 0);
				// 					currentMinimumScore = nodeCalc.minScore - 1;
				// 					currentMinimumIndex = graph.nodeStart[neighbor];
				// 				}
				// 			}
				// 		}
				// 	}
				// }
				result.cellsProcessed += nodeCalc.cellsProcessed;
			}
			for (auto node : previousBandOrder)
			{
				assert(previousBand[node]);
				previousBand[node] = false;
			}
			for (auto node : previousBandOrderOutOfOrder)
			{
				assert(previousBand[node]);
				previousBand[node] = false;
			}
			assert(currentMinimumIndex != std::numeric_limits<LengthType>::max());
			assert(result.minScorePerWordSlice.size() == 0 || currentMinimumScore >= result.minScorePerWordSlice.back());
#ifndef NDEBUG
			auto debugMinimumNode = graph.indexToNode[currentMinimumIndex];
			assert(currentSlice.hasNode(debugMinimumNode));
			auto debugslice = currentSlice.node(debugMinimumNode);
			assert(currentMinimumIndex >= graph.nodeStart[debugMinimumNode]);
			assert(currentMinimumIndex < graph.nodeEnd[debugMinimumNode]);
			assert(debugslice[currentMinimumIndex - graph.nodeStart[debugMinimumNode]].scoreEnd == currentMinimumScore);
#endif
			result.scoreSlices.push_back(currentSlice);
			previousSlice = std::move(currentSlice);
			previousMinimumIndex = currentMinimumIndex;
			result.minScorePerWordSlice.emplace_back(currentMinimumScore);
			result.minScoreIndexPerWordSlice.emplace_back(currentMinimumIndex);
			previousBandOrder = std::move(bandOrder);
			previousBandOrderOutOfOrder = std::move(bandOrderOutOfOrder);
#ifndef NDEBUG
			debugLastRowMinScore = currentMinimumScore;
#endif
			if (currentMinimumScore > maxScore)
			{
				for (int i = j + WordConfiguration<Word>::WordSize; i < sequence.size(); i += WordConfiguration<Word>::WordSize)
				{
					result.minScorePerWordSlice.push_back(sequence.size());
					result.minScoreIndexPerWordSlice.push_back(0);
				}
				break;
			}
		}
#ifndef NDEBUG
		for (size_t i = 1; i < result.minScorePerWordSlice.size(); i++)
		{
			assert(result.minScorePerWordSlice[i] >= result.minScorePerWordSlice[i-1]);
		}
#endif
#ifdef EXTRACORRECTNESSASSERTIONS
		verifySlice(sequence, result, initialSlice, maxScore);
#endif
		return result;
	}

	NodeSlice<WordSlice> getOnlyOneNodeBand(LengthType nodeIndex) const
	{
		NodeSlice<WordSlice> result;
		result.addNode(nodeIndex, graph.nodeEnd[nodeIndex] - graph.nodeStart[nodeIndex]);
		auto& slice = result.node(nodeIndex);
		for (size_t i = 0; i < slice.size(); i++)
		{
			slice[i] = {0, 0, 0, 0};
		}
		return result;
	}

	NodeSlice<WordSlice> getSeededNodeBandForward(LengthType nodeIndex, LengthType startExtensionWidth, size_t sequenceLen) const
	{
		NodeSlice<WordSlice> result;
		std::set<size_t> visited;
		std::priority_queue<NodePosWithDistance, std::vector<NodePosWithDistance>, std::greater<NodePosWithDistance>> queue;
		queue.emplace(nodeIndex, true, 0);
		while (queue.size() > 0)
		{
			auto top = queue.top();
			queue.pop();
			if (top.distance > startExtensionWidth) continue;
			if (visited.count(top.node) == 1) continue;
			visited.insert(top.node);
			assert(!result.hasNode(top.node));
			result.addNode(top.node, graph.nodeEnd[top.node] - graph.nodeStart[top.node]);
			auto& slice = result.node(top.node);
			for (size_t i = 0; i < slice.size(); i++)
			{
				if (top.node == nodeIndex)
				{
					slice[i] = {0, 0, 0, 0};
				}
				else
				{
					slice[i] = {0, 0, sequenceLen, sequenceLen};
				}
			}
			auto newDistance = top.distance + graph.nodeEnd[top.node] - graph.nodeStart[top.node];
			for (auto neighbor : graph.outNeighbors[top.node])
			{
				queue.emplace(neighbor, true, newDistance);
			}
		}
		return result;
	}

	TwoDirectionalSplitAlignment getSplitAlignment(const std::string& sequence, int dynamicWidth, int startExtensionWidth, LengthType matchBigraphNodeId, LengthType matchSequencePosition, ScoreType maxScore) const
	{
		assert(matchSequencePosition > 0);
		assert(matchSequencePosition < sequence.size() - 1);
		auto backwardPart = CommonUtils::ReverseComplement(sequence.substr(0, matchSequencePosition));
		auto forwardPart = sequence.substr(matchSequencePosition);
		int backwardpadding = (WordConfiguration<Word>::WordSize - (backwardPart.size() % WordConfiguration<Word>::WordSize)) % WordConfiguration<Word>::WordSize;
		assert(backwardpadding < WordConfiguration<Word>::WordSize);
		for (int i = 0; i < backwardpadding; i++)
		{
			backwardPart += 'N';
		}
		int forwardpadding = (WordConfiguration<Word>::WordSize - (forwardPart.size() % WordConfiguration<Word>::WordSize)) % WordConfiguration<Word>::WordSize;
		assert(forwardpadding < WordConfiguration<Word>::WordSize);
		for (int i = 0; i < forwardpadding; i++)
		{
			forwardPart += 'N';
		}
		assert(backwardPart.size() + forwardPart.size() <= sequence.size() + 2 * WordConfiguration<Word>::WordSize);
		auto forwardNode = graph.nodeLookup.at(matchBigraphNodeId * 2);
		auto backwardNode = graph.nodeLookup.at(matchBigraphNodeId * 2 + 1);
		assert(graph.nodeSequences.substr(graph.nodeStart[forwardNode], graph.nodeEnd[forwardNode] - graph.nodeStart[forwardNode]) == CommonUtils::ReverseComplement(graph.nodeSequences.substr(graph.nodeStart[backwardNode], graph.nodeEnd[backwardNode] - graph.nodeStart[backwardNode])));
		assert(graph.nodeEnd[forwardNode] - graph.nodeStart[forwardNode] == graph.nodeEnd[backwardNode] - graph.nodeStart[backwardNode]);
		auto forwardInitialBand = getOnlyOneNodeBand(forwardNode);
		auto forwardBand = getSeededNodeBandForward(forwardNode, startExtensionWidth, sequence.size());
		auto forwardRowBandFunction = [this, &forwardBand, dynamicWidth](LengthType j, LengthType previousMinimumIndex, const std::vector<bool>& previousBand) { return defaultForwardRowBandFunction(j, previousMinimumIndex, previousBand, dynamicWidth, forwardBand); };
		auto forwardSlice = getBitvectorSliceScoresAndFinalPosition(forwardPart, forwardInitialBand, maxScore, forwardRowBandFunction);
		auto backwardInitialBand = getOnlyOneNodeBand(backwardNode);
		auto backwardBand = getSeededNodeBandForward(backwardNode, startExtensionWidth, sequence.size());
		auto backwardRowBandFunction = [this, &backwardBand, dynamicWidth](LengthType j, LengthType previousMinimumIndex, const std::vector<bool>& previousBand) { return defaultForwardRowBandFunction(j, previousMinimumIndex, previousBand, dynamicWidth, backwardBand); };
		auto backwardSlice = getBitvectorSliceScoresAndFinalPosition(backwardPart, backwardInitialBand, maxScore, backwardRowBandFunction);
		auto reverseForwardSlice = getBitvectorSliceScoresAndFinalPosition(forwardPart, backwardInitialBand, maxScore, backwardRowBandFunction);
		auto reverseBackwardSlice = getBitvectorSliceScoresAndFinalPosition(backwardPart, forwardInitialBand, maxScore, forwardRowBandFunction);
		auto firstscore = forwardSlice.minScorePerWordSlice.back() + backwardSlice.minScorePerWordSlice.back();
		auto secondscore = reverseForwardSlice.minScorePerWordSlice.back() + reverseBackwardSlice.minScorePerWordSlice.back();
		assert(firstscore <= backwardPart.size() + forwardPart.size());
		assert(secondscore <= backwardPart.size() + forwardPart.size());
		std::cerr << "first direction score: " << firstscore << std::endl;
		std::cerr << "other direction score: " << secondscore << std::endl;
		if (firstscore < secondscore)
		{
			TwoDirectionalSplitAlignment result;
			result.sequenceSplitIndex = matchSequencePosition;
			result.forward = std::move(forwardSlice);
			result.backward = std::move(backwardSlice);
			return result;
		}
		else
		{
			TwoDirectionalSplitAlignment result;
			result.sequenceSplitIndex = matchSequencePosition;
			result.forward = std::move(reverseForwardSlice);
			result.backward = std::move(reverseBackwardSlice);
			return result;
		}
	}

	std::vector<ScoreType> getMergedSplitScores(const std::vector<ScoreType>& backward, const std::vector<ScoreType>& forward, size_t nodeSize, size_t startExtensionWidth) const
	{
		std::vector<ScoreType> partialScores;
		ScoreType endScore = 0;
		for (size_t i = backward.size()-1; i < backward.size(); i--)
		{
			endScore += backward[i];
			partialScores.push_back(endScore);
		}
		for (size_t i = 0; i < forward.size(); i++)
		{
			endScore += forward[i];
			partialScores.push_back(endScore);
		}
		partialScores.back() += nodeSize + startExtensionWidth * 2;
		return partialScores;
	}

	std::vector<MatrixPosition> reverseTrace(std::vector<MatrixPosition> trace, LengthType end) const
	{
		if (trace.size() == 0) return trace;
		std::reverse(trace.begin(), trace.end());
		for (size_t i = 0; i < trace.size(); i++)
		{
			trace[i].first = graph.GetReversePosition(trace[i].first);
			assert(trace[i].second <= end);
			trace[i].second = end - trace[i].second;
		}
		return trace;
	}

	template <typename T>
	T factorial(T n) const
	{
		T result {1};
		for (int i = 2; i <= n; i += 1)
		{
			result *= i;
		}
		return result;
	}

	template <typename T>
	T choose(T n, T k) const
	{
		return factorial<T>(n) / factorial<T>(k) / factorial<T>(n - k);
	}

	template <typename T>
	T powr(T base, int exponent) const
	{
		if (exponent == 0) return T{1};
		if (exponent == 1) return base;
		if (exponent % 2 == 0)
		{
			auto part = powr(base, exponent / 2);
			assert(part > 0);
			return part * part;
		}
		if (exponent % 2 == 1)
		{
			auto part = powr(base, exponent / 2);
			assert(part > 0);
			return part * part * base;
		}
		assert(false);
	}

	std::vector<bool> estimateCorrectAlignmentViterbi(const std::vector<ScoreType>& scores) const
	{
		const boost::rational<boost::multiprecision::cpp_int> correctMismatchProbability {15, 100}; //15% from pacbio error rate
		const boost::rational<boost::multiprecision::cpp_int> falseMismatchProbability {50, 100}; //50% empirically
		const boost::rational<boost::multiprecision::cpp_int> falseToCorrectTransitionProbability {1, 100000}; //10^-5. one slice at <=12 mismatches or two slices at <=15 mismatchs
		const boost::rational<boost::multiprecision::cpp_int> correctToFalseTransitionProbability {1LL, 1000000000000000LL}; //10^-15. three slices at >=25 mismatches or two slices at >=30 mismatches
		boost::rational<boost::multiprecision::cpp_int> correctProbability {80, 100}; //80% arbitrarily
		boost::rational<boost::multiprecision::cpp_int> falseProbability {20, 100}; //20% arbitrarily
		std::vector<bool> falseFromCorrectBacktrace;
		std::vector<bool> correctFromCorrectBacktrace;
		for (size_t i = 1; i < scores.size(); i++)
		{
			assert(scores[i] >= scores[i-1]);
			auto scorediff = scores[i] - scores[i-1];
			correctFromCorrectBacktrace.push_back(correctProbability * (1 - correctToFalseTransitionProbability) >= falseProbability * falseToCorrectTransitionProbability);
			falseFromCorrectBacktrace.push_back(correctProbability * correctToFalseTransitionProbability >= falseProbability * (1 - falseToCorrectTransitionProbability));
			boost::rational<boost::multiprecision::cpp_int> newCorrectProbability = std::max(correctProbability * (1 - correctToFalseTransitionProbability), falseProbability * falseToCorrectTransitionProbability);
			boost::rational<boost::multiprecision::cpp_int> newFalseProbability = std::max(correctProbability * correctToFalseTransitionProbability, falseProbability * (1 - falseToCorrectTransitionProbability));
			auto chooseresult = choose<boost::multiprecision::cpp_int>(WordConfiguration<Word>::WordSize, scorediff);
			auto correctMultiplier = chooseresult * powr(correctMismatchProbability, scorediff) * powr(1 - correctMismatchProbability, WordConfiguration<Word>::WordSize - scorediff);
			auto falseMultiplier = chooseresult * powr(falseMismatchProbability, scorediff) * powr(1 - falseMismatchProbability, WordConfiguration<Word>::WordSize - scorediff);
			newCorrectProbability *= correctMultiplier;
			newFalseProbability *= falseMultiplier;
			correctProbability = newCorrectProbability;
			falseProbability = newFalseProbability;
			auto normalizer = correctProbability + falseProbability;
			correctProbability /= normalizer;
			falseProbability /= normalizer;
		}
		assert(falseFromCorrectBacktrace.size() == scores.size() - 1);
		assert(correctFromCorrectBacktrace.size() == scores.size() - 1);
		bool currentCorrect = correctProbability > falseProbability;
		std::vector<bool> result;
		result.resize(scores.size() - 1);
		for (size_t i = scores.size()-2; i < scores.size(); i--)
		{
			result[i] = currentCorrect;
			if (currentCorrect)
			{
				currentCorrect = correctFromCorrectBacktrace[i];
			}
			else
			{
				currentCorrect = falseFromCorrectBacktrace[i];
			}
		}
		return result;
	}

	std::pair<std::tuple<ScoreType, std::vector<MatrixPosition>>, std::tuple<ScoreType, std::vector<MatrixPosition>>> getPiecewiseTracesFromSplit(const TwoDirectionalSplitAlignment& split, const std::string& sequence) const
	{
		std::string backtraceSequence;
		std::string backwardBacktraceSequence;
		auto startpartsize = split.sequenceSplitIndex;
		auto endpartsize = sequence.size() - split.sequenceSplitIndex;
		int startpadding = (WordConfiguration<Word>::WordSize - (startpartsize % WordConfiguration<Word>::WordSize)) % WordConfiguration<Word>::WordSize;
		int endpadding = (WordConfiguration<Word>::WordSize - (endpartsize % WordConfiguration<Word>::WordSize)) % WordConfiguration<Word>::WordSize;
		backtraceSequence = sequence.substr(split.sequenceSplitIndex);
		backwardBacktraceSequence = CommonUtils::ReverseComplement(sequence.substr(0, split.sequenceSplitIndex));
		backtraceSequence.reserve(sequence.size() + endpadding);
		backwardBacktraceSequence.reserve(sequence.size() + startpadding);
		for (int i = 0; i < startpadding; i++)
		{
			backwardBacktraceSequence += 'N';
		}
		for (int i = 0; i < endpadding; i++)
		{
			backtraceSequence += 'N';
		}
		assert(backtraceSequence.size() % WordConfiguration<Word>::WordSize == 0);
		assert(backwardBacktraceSequence.size() % WordConfiguration<Word>::WordSize == 0);
		auto backtraceSlice = getBacktraceSlice(backtraceSequence, split.forward);

		// auto backtraceresult = getTraceFromSlices(split.forward, backtraceSlice);
		std::pair<ScoreType, std::vector<MatrixPosition>> backtraceresult;

		std::cerr << "fw score: " << std::get<0>(backtraceresult) << std::endl;
		auto reverseBacktraceSlice = getBacktraceSlice(backwardBacktraceSequence, split.backward);

		// auto reverseBacktraceResult = getTraceFromSlices(split.backward, reverseBacktraceSlice);
		std::pair<ScoreType, std::vector<MatrixPosition>> reverseBacktraceResult;

		std::cerr << "bw score: " << std::get<0>(reverseBacktraceResult) << std::endl;

		while (backtraceresult.second.size() > 0 && backtraceresult.second.back().second >= backtraceSequence.size() - endpadding)
		{
			assert(backtraceresult.second.back().second >= backtraceSequence.size() - endpadding);
			backtraceresult.second.pop_back();
		}
		while (reverseBacktraceResult.second.size() > 0 && reverseBacktraceResult.second.back().second >= backwardBacktraceSequence.size() - startpadding)
		{
			assert(reverseBacktraceResult.second.back().second >= backwardBacktraceSequence.size() - startpadding);
			reverseBacktraceResult.second.pop_back();
		}

		reverseBacktraceResult.second = reverseTrace(reverseBacktraceResult.second, split.sequenceSplitIndex - 1);

		for (size_t i = 0; i < backtraceresult.second.size(); i++)
		{
			backtraceresult.second[i].second += split.sequenceSplitIndex;
		}

		return std::make_pair(backtraceresult, reverseBacktraceResult);
	}

	std::tuple<ScoreType, std::vector<MatrixPosition>, size_t> getBacktraceFullStart(std::string sequence, int dynamicWidth) const
	{
		int padding = (WordConfiguration<Word>::WordSize - (sequence.size() % WordConfiguration<Word>::WordSize)) % WordConfiguration<Word>::WordSize;
		for (int i = 0; i < padding; i++)
		{
			sequence += 'N';
		}
		NodeSlice<WordSlice> startBand;
		for (size_t i = 0; i < graph.nodeStart.size(); i++)
		{
			startBand.addNode(i, graph.nodeEnd[i] - graph.nodeStart[i]);
			auto& slice = startBand.node(i);
			for (size_t ii = 0; ii < slice.size(); ii++)
			{
				slice[ii] = {0, 0, 0, 0};
			}
		}
		auto rowBandFunction = [this, &startBand, dynamicWidth](LengthType j, LengthType previousMinimumIndex, const std::vector<bool>& previousBand) { return defaultForwardRowBandFunction(j, previousMinimumIndex, previousBand, dynamicWidth, startBand); };
		auto slice = getBitvectorSliceScoresAndFinalPosition(sequence, startBand, sequence.size() * 0.4, rowBandFunction);
		std::cerr << "score: " << slice.minScorePerWordSlice.back() << std::endl;
		if (slice.minScorePerWordSlice.back() > sequence.size() * 0.4)
		{
			return std::make_tuple(std::numeric_limits<ScoreType>::max(), std::vector<MatrixPosition>{}, slice.cellsProcessed);
		}
		auto backtraceSlice = getBacktraceSlice(sequence, slice);

		std::pair<ScoreType, std::vector<MatrixPosition>> backtraceresult;
		return std::make_tuple(backtraceresult.first, backtraceresult.second, slice.cellsProcessed);
		// auto backtraceresult = getTraceFromSlices(slice, backtraceSlice);
		// assert(backtraceresult.first <= slice.minScorePerWordSlice.back());
		// while (backtraceresult.second.back().second >= sequence.size() - padding)
		// {
		// 	backtraceresult.second.pop_back();
		// }
		// assert(backtraceresult.second[0].second == 0);
		// assert(backtraceresult.second.back().second == sequence.size() - padding - 1);
		// return std::make_tuple(backtraceresult.first, backtraceresult.second, slice.cellsProcessed);
	}

	const AlignmentGraph& graph;
};

#endif