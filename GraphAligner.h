#ifndef GraphAligner_H
#define GraphAligner_H

//http://biorxiv.org/content/early/2017/04/06/124941
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
#include "vg.pb.h"
#include "SliceRow.h"
#include "SparseBoolMatrix.h"
#include "SparseMatrix.h"
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
	    x = (x & 0x5555555555555555) + ((x >>  1) & 0x5555555555555555); //put count of each  2 bits into those  2 bits 
	    x = (x & 0x3333333333333333) + ((x >>  2) & 0x3333333333333333); //put count of each  4 bits into those  4 bits 
	    x = (x & 0x0F0F0F0F0F0F0F0F) + ((x >>  4) & 0x0F0F0F0F0F0F0F0F); //put count of each  8 bits into those  8 bits 
	    x = (x & 0x00FF00FF00FF00FF) + ((x >>  8) & 0x00FF00FF00FF00FF); //put count of each 16 bits into those 16 bits 
	    x = (x & 0x0000FFFF0000FFFF) + ((x >> 16) & 0x0000FFFF0000FFFF); //put count of each 32 bits into those 32 bits 
	    x = (x & 0x00000000FFFFFFFF) + ((x >> 32) & 0x00000000FFFFFFFF); //put count of each 64 bits into those 64 bits 
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
	class SeedHit
	{
	public:
		SeedHit(size_t seqPos, int nodeId, size_t nodePos) : sequencePosition(seqPos), nodeId(nodeId), nodePos(nodePos) {};
		size_t sequencePosition;
		int nodeId;
		size_t nodePos;
	};
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

	GraphAligner() :
	nodeStart(),
	indexToNode(),
	nodeLookup(),
	nodeIDs(),
	inNeighbors(),
	nodeSequences(),
	gapStartPenalty(1),
	gapContinuePenalty(1),
	finalized(false)
	{
		//add the start dummy node as the first node
		dummyNodeStart = nodeSequences.size();
		nodeIDs.push_back(0);
		nodeStart.push_back(nodeSequences.size());
		inNeighbors.emplace_back();
		outNeighbors.emplace_back();
		reverse.push_back(false);
		nodeSequences.push_back('N');
		indexToNode.resize(nodeSequences.size(), nodeStart.size()-1);
		nodeEnd.emplace_back(nodeSequences.size());
		notInOrder.push_back(false);
	}
	
	void AddNode(int nodeId, std::string sequence, bool reverseNode)
	{
		assert(!finalized);
		//subgraph extraction might produce different subgraphs with common nodes
		//don't add duplicate nodes
		if (nodeLookup.count(nodeId) != 0) return;

		assert(std::numeric_limits<LengthType>::max() - sequence.size() > nodeSequences.size());
		nodeLookup[nodeId] = nodeStart.size();
		nodeIDs.push_back(nodeId);
		nodeStart.push_back(nodeSequences.size());
		inNeighbors.emplace_back();
		outNeighbors.emplace_back();
		reverse.push_back(reverseNode);
		nodeSequences.insert(nodeSequences.end(), sequence.begin(), sequence.end());
		indexToNode.resize(nodeSequences.size(), nodeStart.size()-1);
		nodeEnd.emplace_back(nodeSequences.size());
		notInOrder.push_back(false);
		assert(nodeIDs.size() == nodeStart.size());
		assert(nodeStart.size() == inNeighbors.size());
		assert(inNeighbors.size() == nodeEnd.size());
		assert(nodeEnd.size() == notInOrder.size());
		assert(nodeSequences.size() == indexToNode.size());
		assert(inNeighbors.size() == outNeighbors.size());
	}
	
	void AddEdgeNodeId(int node_id_from, int node_id_to)
	{
		assert(!finalized);
		assert(nodeLookup.count(node_id_from) > 0);
		assert(nodeLookup.count(node_id_to) > 0);
		auto from = nodeLookup[node_id_from];
		auto to = nodeLookup[node_id_to];
		assert(to >= 0);
		assert(from >= 0);
		assert(to < inNeighbors.size());
		assert(from < nodeStart.size());

		//subgraph extraction might produce different subgraphs with common edges
		//don't add duplicate edges
		if (std::find(inNeighbors[to].begin(), inNeighbors[to].end(), from) != inNeighbors[to].end()) return;

		inNeighbors[to].push_back(from);
		outNeighbors[from].push_back(to);
		if (from >= to)
		{
			//todo fix for cyclic graphs
			assert(false);
			notInOrder[to] = true;
		}
	}

	void Finalize()
	{
		//add the end dummy node as the last node
		dummyNodeEnd = nodeSequences.size();
		nodeIDs.push_back(0);
		nodeStart.push_back(nodeSequences.size());
		reverse.push_back(false);
		inNeighbors.emplace_back();
		outNeighbors.emplace_back();
		nodeSequences.push_back('N');
		indexToNode.resize(nodeSequences.size(), nodeStart.size()-1);
		nodeEnd.emplace_back(nodeSequences.size());
		notInOrder.push_back(false);
		assert(nodeSequences.size() >= nodeStart.size());
		assert(nodeEnd.size() == nodeStart.size());
		assert(notInOrder.size() == nodeStart.size());
		assert(inNeighbors.size() == nodeStart.size());
		assert(outNeighbors.size() == nodeStart.size());
		assert(notInOrder.size() == nodeStart.size());
		assert(reverse.size() == nodeStart.size());
		assert(nodeIDs.size() == nodeStart.size());
		assert(indexToNode.size() == nodeSequences.size());
		finalized = true;
	}

	AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence, int dynamicWidth, LengthType dynamicRowStart) const
	{
		assert(finalized);
		auto band = getFullBand(sequence.size(), dynamicRowStart);
		auto trace = getBacktrace(sequence, dynamicWidth, dynamicRowStart, band);
		//failed alignment, don't output
		if (std::get<0>(trace) == std::numeric_limits<ScoreType>::min()) return emptyAlignment();
		auto result = traceToAlignment(seq_id, sequence, std::get<0>(trace), std::get<2>(trace), std::get<1>(trace), std::get<3>(trace));
		return result;
	}

	AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence, int dynamicWidth, LengthType dynamicRowStart, const std::vector<SeedHit>& seedHits, int startBandwidth) const
	{
		assert(finalized);
		auto band = getSeededStartBand(seedHits, dynamicRowStart, startBandwidth, sequence);
		auto trace = getBacktrace(sequence, dynamicWidth, dynamicRowStart, band);
		//failed alignment, don't output
		if (std::get<0>(trace) == std::numeric_limits<ScoreType>::min()) return emptyAlignment();
		auto result = traceToAlignment(seq_id, sequence, std::get<0>(trace), std::get<2>(trace), std::get<1>(trace), std::get<3>(trace));
		return result;
	}

	size_t SizeInBp()
	{
		return nodeSequences.size();
	}

private:

	std::vector<std::vector<bool>> getFullBand(size_t sequenceSize, LengthType dynamicRowStart) const
	{
		std::vector<std::vector<bool>> result;
		result.resize(dynamicRowStart/WordConfiguration<Word>::WordSize);
		for (size_t i = 0; i < dynamicRowStart/WordConfiguration<Word>::WordSize; i++)
		{
			result[i].resize(nodeStart.size(), true);
		}
		return result;
	}

	std::vector<std::vector<bool>> getSeededStartBand(const std::vector<SeedHit>& originalSeedHits, int dynamicRowStart, int startBandwidth, const std::string& sequence) const
	{
		auto seedHitsInMatrix = getSeedHitPositionsInMatrix(sequence, originalSeedHits);
		auto bandLocations = getBandLocations(sequence.size(), seedHitsInMatrix, dynamicRowStart);
		std::vector<std::vector<bool>> result;
		result.resize(dynamicRowStart/WordConfiguration<Word>::WordSize);
		for (LengthType j = 0; j < dynamicRowStart && j < sequence.size()+1; j += WordConfiguration<Word>::WordSize)
		{
			auto index = j/WordConfiguration<Word>::WordSize;
			result[index].resize(nodeStart.size(), false);
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

	std::vector<MatrixPosition> getSeedHitPositionsInMatrix(const std::string& sequence, const std::vector<SeedHit>& seedHits) const
	{
		std::vector<MatrixPosition> result;
		for (size_t i = 0; i < seedHits.size(); i++)
		{
			assert(nodeLookup.count(seedHits[i].nodeId) > 0);
			result.emplace_back(nodeStart[nodeLookup.at(seedHits[i].nodeId)] + seedHits[i].nodePos, seedHits[i].sequencePosition);
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
		size_t oldNode = indexToNode[trace[0].first];
		while (oldNode == dummyNodeStart)
		{
			pos++;
			if (pos == trace.size()) return emptyAlignment();
			assert(pos < trace.size());
			oldNode = indexToNode[trace[pos].first];
			assert(oldNode < nodeIDs.size());
		}
		if (oldNode == dummyNodeEnd) return emptyAlignment();
		int rank = 0;
		auto vgmapping = path->add_mapping();
		auto position = new vg::Position;
		vgmapping->set_allocated_position(position);
		vgmapping->set_rank(rank);
		position->set_node_id(nodeIDs[oldNode]);
		position->set_is_reverse(reverse[oldNode]);
		for (; pos < trace.size(); pos++)
		{
			if (indexToNode[trace[pos].first] == dummyNodeEnd) break;
			if (indexToNode[trace[pos].first] == oldNode) continue;
			oldNode = indexToNode[trace[pos].first];
			rank++;
			vgmapping = path->add_mapping();
			position = new vg::Position;
			vgmapping->set_allocated_position(position);
			vgmapping->set_rank(rank);
			position->set_node_id(nodeIDs[oldNode]);
			position->set_is_reverse(reverse[oldNode]);
		}
		result.set_score(score);
		result.set_sequence(sequence);
		return AlignmentResult { result, maxDistanceFromBand, false, cellsProcessed };
	}

	template <typename SparseMatrixType, typename MatrixType>
	std::tuple<ScoreType, int, std::vector<MatrixPosition>> backtrace(const std::vector<ScoreType>& Mslice, const SparseMatrixType& backtraceMatrix, const MatrixType& band, int sequenceLength, const std::vector<LengthType>& maxScorePositionPerRow, LengthType dynamicRowStart) const
	{
		assert(backtraceMatrix.sizeRows() == sequenceLength+1);
		assert(backtraceMatrix.sizeColumns() == nodeSequences.size());
		std::vector<MatrixPosition> trace;
		bool foundStart = false;
		MatrixPosition currentPosition = std::make_pair(0, sequenceLength);
		//start at the highest value at end of read
		for (size_t i = 0; i < Mslice.size(); i++)
		{
			if (!band(i, sequenceLength)) continue;
			MatrixPosition candidatePosition = std::make_pair(i, sequenceLength);
			if (!foundStart)
			{
				currentPosition = candidatePosition;
				foundStart = true;
			}
			if (Mslice[candidatePosition.first] > Mslice[currentPosition.first])
			{
				currentPosition = candidatePosition;
			}
		}
		assert(band(currentPosition.first, currentPosition.second));
		auto score = Mslice[currentPosition.first];
		trace.push_back(currentPosition);
		LengthType maxMinDistance = 0;
		while (currentPosition.second > 0)
		{
			assert(band(currentPosition.first, currentPosition.second));
			//the rows 0-dynamicRowStart don't use the dynamic band, don't include them here
			if (currentPosition.second > dynamicRowStart)
			{
				maxMinDistance = std::max(maxMinDistance, bandDistanceFromSeqToSeq(currentPosition.first, maxScorePositionPerRow[currentPosition.second]));
			}
			assert(currentPosition.second >= 0);
			assert(currentPosition.second < sequenceLength+1);
			assert(currentPosition.first >= 0);
			assert(currentPosition.first < nodeSequences.size());
			//If we're at the dummy node, we have to stay there
			if (currentPosition.first == 0) break;
			assert(backtraceMatrix.exists(currentPosition.first, currentPosition.second));
			auto newPos = backtraceMatrix(currentPosition.first, currentPosition.second);
			assert(newPos.second < currentPosition.second || (newPos.second == currentPosition.second && newPos.first < currentPosition.first));
			currentPosition = newPos;
			trace.push_back(currentPosition);
		}
		std::reverse(trace.begin(), trace.end());
		return std::make_tuple(score, maxMinDistance, trace);
	}

	void expandBandForwards(std::vector<std::vector<LengthType>>& result, LengthType w, LengthType j, size_t sequenceLength, LengthType maxRow) const
	{
		if (std::find(result[j].begin(), result[j].end(), w) != result[j].end()) return;
		auto nodeIndex = indexToNode[w];
		auto end = nodeEnd[nodeIndex];
		while (w != end && j < sequenceLength+1 && j < maxRow)
		{
			result[j].emplace_back(w);
			w++;
			j++;
		}
		if (w == end && j < sequenceLength+1 && j < maxRow)
		{
			for (size_t i = 0; i < outNeighbors[nodeIndex].size(); i++)
			{
				expandBandForwards(result, nodeStart[outNeighbors[nodeIndex][i]], j, sequenceLength, maxRow);
			}
		}
	}

	void expandBandBackwards(std::vector<std::vector<LengthType>>& result, LengthType w, LengthType j, size_t sequenceLength) const
	{
		if (std::find(result[j].begin(), result[j].end(), w) != result[j].end()) return;
		auto nodeIndex = indexToNode[w];
		auto start = nodeStart[nodeIndex];
		while (w != start && j > 0)
		{
			result[j].emplace_back(w);
			w--;
			j--;
		}
		result[j].emplace_back(w);
		if (w == start && j > 0)
		{
			for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
			{
				expandBandBackwards(result, nodeEnd[inNeighbors[nodeIndex][i]] - 1, j-1, sequenceLength);
			}
		}
	}

	void expandBandDynamically(std::vector<bool>& band, LengthType previousMinimumIndex, LengthType dynamicWidth) const
	{
		assert(previousMinimumIndex < nodeSequences.size());
		auto nodeIndex = indexToNode[previousMinimumIndex];
		band[nodeIndex] = true;
		LengthType start = nodeStart[nodeIndex];
		LengthType end = nodeEnd[nodeIndex];
		if (dynamicWidth > previousMinimumIndex - start)
		{
			for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
			{
				expandDynamicBandBackward(band, inNeighbors[nodeIndex][i], dynamicWidth - (previousMinimumIndex - start));
			}
		}
		if (dynamicWidth > end - previousMinimumIndex)
		{
			for (size_t i = 0; i < outNeighbors[nodeIndex].size(); i++)
			{
				expandDynamicBandForward(band, outNeighbors[nodeIndex][i], dynamicWidth - (end - previousMinimumIndex));
			}
		}
	}

	void expandDynamicBandBackward(std::vector<bool>& band, LengthType nodeIndex, LengthType dynamicWidth) const
	{
		if (dynamicWidth == 0) return;
		assert(nodeIndex < nodeStart.size());
		//todo fix: currently only the first path that reaches the node is considered
		//this means that it might not reach some nodes with distance < dynamicwidth
		if (band[nodeIndex]) return;
		band[nodeIndex] = true;
		for (size_t i = 0; i < outNeighbors[nodeIndex].size(); i++)
		{
			expandDynamicBandForward(band, outNeighbors[nodeIndex][i], dynamicWidth - 1);
		}
		auto nodeSize = nodeEnd[nodeIndex] - nodeStart[nodeIndex];
		if (dynamicWidth > nodeSize)
		{
			for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
			{
				expandDynamicBandBackward(band, inNeighbors[nodeIndex][i], dynamicWidth - nodeSize);
			}
		}
	}

	template <typename MatrixType>
	void expandDynamicBandForward(MatrixType& band, LengthType nodeIndex, LengthType dynamicWidth) const
	{
		if (dynamicWidth == 0) return;
		assert(nodeIndex < nodeStart.size());
		//todo fix: currently only the first path that reaches the node is considered
		//this means that it might not reach some nodes with distance < dynamicwidth if there's multiple paths and it arbitrarily picks the longest one first
		if (band[nodeIndex]) return;
		band[nodeIndex] = true;
		for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
		{
			expandDynamicBandBackward(band, inNeighbors[nodeIndex][i], dynamicWidth - 1);
		}
		auto nodeSize = nodeEnd[nodeIndex] - nodeStart[nodeIndex];
		if (dynamicWidth > nodeSize)
		{
			for (size_t i = 0; i < outNeighbors[nodeIndex].size(); i++)
			{
				expandDynamicBandForward(band, outNeighbors[nodeIndex][i], dynamicWidth - nodeSize);
			}
		}
	}

	WordSlice mergeTwoSlices(WordSlice left, WordSlice right) const
	{
		//currently O(w), O(log^2 w) is possible
		//todo figure out the details and implement
		//todo is O(log w) possible? an expected O(log w), worst case O(w) seems possible
		assert((left.VP & left.VN) == WordConfiguration<Word>::AllZeros);
		assert((right.VP & right.VN) == WordConfiguration<Word>::AllZeros);
		ScoreType leftScore = left.scoreStart;
		WordSlice merged;
		merged.scoreBeforeStart = std::min(left.scoreBeforeStart, right.scoreBeforeStart);
		merged.VP = WordConfiguration<Word>::AllZeros;
		merged.VN = WordConfiguration<Word>::AllZeros;
		ScoreType rightScore = right.scoreStart;
		merged.scoreStart = std::min(rightScore, leftScore);
		assert(merged.scoreStart >= merged.scoreBeforeStart - 1);
		assert(merged.scoreStart <= merged.scoreBeforeStart + 1);
		if (merged.scoreStart == merged.scoreBeforeStart - 1)
		{
			merged.VN |= 1;
		}
		else if (merged.scoreStart == merged.scoreBeforeStart + 1)
		{
			merged.VP |= 1;
		}
		ScoreType previousScore = merged.scoreStart;
		for (size_t j = 1; j < WordConfiguration<Word>::WordSize; j++)
		{
			Word mask = ((Word)1) << j;
			if (left.VP & mask) leftScore++;
			if (left.VN & mask) leftScore--;
			if (right.VP & mask) rightScore++;
			if (right.VN & mask) rightScore--;
			ScoreType betterScore = std::min(leftScore, rightScore);
			if (betterScore == previousScore+1) merged.VP |= mask;
			if (betterScore == previousScore-1) merged.VN |= mask;
			assert((merged.VP & merged.VN) == WordConfiguration<Word>::AllZeros);
			assert(betterScore >= previousScore-1);
			assert(betterScore <= previousScore+1);
			previousScore = betterScore;
		}
		merged.scoreEnd = previousScore;
		assert((merged.VP & merged.VN) == WordConfiguration<Word>::AllZeros);
		assert(merged.scoreEnd <= left.scoreEnd);
		assert(merged.scoreEnd <= right.scoreEnd);
		assert(merged.scoreStart <= left.scoreStart);
		assert(merged.scoreStart <= right.scoreStart);
		return merged;
	}

	WordSlice getPreviousSlice(size_t nodeIndex, const std::vector<WordSlice>& previousSlice, const std::vector<WordSlice>& currentSlice, const std::vector<bool>& band) const
	{
		auto start = nodeStart[nodeIndex];
		WordSlice result;
		bool foundOne = false;
		for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
		{
			if (!band[inNeighbors[nodeIndex][i]]) continue;
			if (!foundOne)
			{
				result = currentSlice[nodeEnd[inNeighbors[nodeIndex][i]]-1];
				foundOne = true;
			}
			else
			{
				result = mergeTwoSlices(result, currentSlice[nodeEnd[inNeighbors[nodeIndex][i]]-1]);
			}
		}
		assert(foundOne);
		return result;
	}

	WordSlice getSourceSlice(size_t nodeIndex, const std::vector<WordSlice>& previousSlice) const
	{
		auto start = nodeStart[nodeIndex];
		return { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, previousSlice[start].scoreEnd+1, previousSlice[start].scoreEnd+WordConfiguration<Word>::WordSize, previousSlice[start].scoreEnd };
	}

	bool isSource(size_t nodeIndex, const std::vector<bool>& band) const
	{
		for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
		{
			if (band[inNeighbors[nodeIndex][i]]) return false;
		}
		return true;
	}

	MatrixSlice getBitvectorSliceScoresAndFinalPosition(std::string sequence, int dynamicWidth, std::vector<std::vector<bool>>& startBand, LengthType dynamicRowStart) const
	{
		MatrixSlice result;
		result.cellsProcessed = 0;
		result.finalMinScore = 0;
		result.finalMinScoreColumn = 0;

		std::vector<WordSlice> slice1;
		std::vector<WordSlice> slice2;
		slice1.resize(nodeSequences.size());
		slice2.resize(nodeSequences.size());

		int padding = (WordConfiguration<Word>::WordSize - (sequence.size() % WordConfiguration<Word>::WordSize)) % WordConfiguration<Word>::WordSize;
		for (int i = 0; i < padding; i++)
		{
			sequence += 'N';
		}

		std::vector<WordSlice>& currentSlice = slice1;
		std::vector<WordSlice>& previousSlice = slice2;

		LengthType previousMinimumIndex;
		std::vector<bool> band;

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
			std::vector<bool>& currentBand = band;
			if (startBand.size() > slice)
			{
				currentBand = startBand[slice];
			}
			else
			{
				band.assign(band.size(), false);
				expandBandDynamically(band, previousMinimumIndex, dynamicWidth);
			}
			for (size_t i = 0; i < nodeStart.size(); i++)
			{
				if (!currentBand[i]) continue;
				WordSlice previous;
				LengthType start = nodeStart[i];
				if (isSource(i, currentBand))
				{
					previous = getSourceSlice(i, previousSlice);
					currentSlice[start] = previous;
					if (currentSlice[start].scoreEnd < currentMinimumScore)
					{
						currentMinimumScore = currentSlice[start].scoreEnd;
						currentMinimumIndex = start;
					}
					start++;
				}
				else
				{
					previous = getPreviousSlice(i, previousSlice, currentSlice, band);
				}
				ScoreType scoreStart = previous.scoreStart;
				ScoreType scoreEnd = previous.scoreEnd;
				Word VP = previous.VP;
				Word VN = previous.VN;

				for (LengthType w = start; w < nodeEnd[i]; w++)
				{
					Word Eq;
					switch(nodeSequences[w])
					{
						case 'A':
							Eq = BA;
							break;
						case 'T':
							Eq = BT;
							break;
						case 'C':
							Eq = BC;
							break;
						case 'G':
							Eq = BG;
							break;
						default:
							assert(false);
							break;
					}
					//http://www.gersteinlab.org/courses/452/09-spring/pdf/Myers.pdf
					//pages 405 and 408
					int hin = previousSlice[w].hout;
					if (j == 0) hin = 0;
					Word Xv = Eq | VN;
					//between 7 and 8
					if (hin < 0) Eq |= 1;
					Word Xh = (((Eq & VP) + VP) ^ VP) | Eq;
					Word Ph = VN | ~(Xh | VP);
					Word Mh = VP & Xh;
					if (Ph & ((Word)1))
					{
						scoreStart += 1;
					}
					else if (Mh & ((Word)1))
					{
						scoreStart -= 1;
					}
					Word lastBitMask = (((Word)1) << (WordConfiguration<Word>::WordSize - 1));
					if (Ph & lastBitMask)
					{
						scoreEnd += 1;
						currentSlice[w].hout = 1;
					}
					else if (Mh & lastBitMask)
					{
						scoreEnd -= 1;
						currentSlice[w].hout = -1;
					}
					else
					{
						currentSlice[w].hout = 0;
					}
					Ph <<= 1;
					Mh <<= 1;
					//between 16 and 17
					if (hin < 0) Mh |= 1; else if (hin > 0) Ph |= 1;
					VP = Mh | ~(Xv | Ph);
					VN = Ph & Xv;

					assert(scoreStart >= previousSlice[w].scoreEnd - 1);
					assert(scoreStart <= previousSlice[w].scoreEnd + 1);
					auto wcvp = WordConfiguration<Word>::popcount(VP);
					auto wcvn = WordConfiguration<Word>::popcount(VN);
					assert(scoreEnd == previousSlice[w].scoreEnd + wcvp - wcvn);
					auto wcvpExceptFirst = WordConfiguration<Word>::popcount(VP & ~((Word)1));
					auto wcvnExceptFirst = WordConfiguration<Word>::popcount(VN & ~((Word)1));
					assert(scoreEnd == scoreStart + wcvpExceptFirst - wcvnExceptFirst);
					assert(scoreStart >= 0);
					assert(scoreEnd >= 0);
					assert(scoreEnd <= j + WordConfiguration<Word>::WordSize);
					assert(scoreStart <= j + 1);
					assert(scoreStart <= scoreEnd + WordConfiguration<Word>::WordSize);
					assert(scoreEnd <= scoreStart + WordConfiguration<Word>::WordSize);
					assert(j == 0 || (VP & 1) == 0 || scoreStart == previousSlice[w].scoreEnd+1);
					assert(j == 0 || (VN & 1) == 0 || scoreStart == previousSlice[w].scoreEnd-1);
					assert((VP & VN) == WordConfiguration<Word>::AllZeros);
					currentSlice[w].scoreStart = scoreStart;
					currentSlice[w].scoreEnd = scoreEnd;
					currentSlice[w].VP = VP;
					currentSlice[w].VN = VN;
					currentSlice[w].scoreBeforeStart = scoreStart;
					if (currentSlice[w].VP & 1) currentSlice[w].scoreBeforeStart -= 1;
					if (currentSlice[w].VN & 1) currentSlice[w].scoreBeforeStart += 1;
					assert(currentSlice[w].scoreBeforeStart == previousSlice[w].scoreEnd);
					assert(currentSlice[w].scoreBeforeStart >= 0);
					assert(currentSlice[w].scoreBeforeStart <= j);
					if (scoreEnd < currentMinimumScore)
					{
						currentMinimumScore = scoreEnd;
						currentMinimumIndex = w;
					}
				}
				result.cellsProcessed += (nodeEnd[i] - nodeStart[i]) * WordConfiguration<Word>::WordSize;
			}
			std::swap(currentSlice, previousSlice);
			previousMinimumIndex = currentMinimumIndex;
			result.minScorePerWordSlice.emplace_back(currentMinimumScore);
		}
		result.finalMinScoreColumn = previousMinimumIndex - padding;
		result.finalMinScore = result.minScorePerWordSlice.back() - padding;
		return result;
	}

	std::tuple<ScoreType, int, std::vector<MatrixPosition>, size_t> getBacktrace(const std::string& sequence, int dynamicWidth, LengthType dynamicRowStart, std::vector<std::vector<bool>>& startBand) const
	{
		auto slice = getBitvectorSliceScoresAndFinalPosition(sequence, dynamicWidth, startBand, dynamicRowStart);
		return std::make_tuple(slice.finalMinScore, 0, std::vector<MatrixPosition>{ std::make_pair(slice.finalMinScoreColumn, 0) }, slice.cellsProcessed);
		//todo fix
		// auto backtraceresult = backtrace(slice.M, backtraceMatrix, startBand, sequence.size(), slice.maxScorePositionPerRow, dynamicRowStart);
		// return std::make_tuple(std::get<0>(backtraceresult), std::get<1>(backtraceresult), std::get<2>(backtraceresult), slice.cellsProcessed);
	}

	ScoreType gapPenalty(LengthType length) const
	{
		if (length == 0) return 0;
		return gapStartPenalty + gapContinuePenalty * (length - 1);
	}

	ScoreType matchScore(char graph, char sequence) const
	{
		return graph == sequence ? 1 : -1;
	}

	std::vector<bool> notInOrder;
	std::vector<LengthType> nodeStart;
	std::vector<LengthType> nodeEnd;
	std::vector<LengthType> indexToNode;
	std::map<int, LengthType> nodeLookup;
	std::vector<int> nodeIDs;
	std::vector<std::vector<LengthType>> inNeighbors;
	std::vector<std::vector<LengthType>> outNeighbors;
	std::vector<bool> reverse;
	std::string nodeSequences;
	ScoreType gapStartPenalty;
	ScoreType gapContinuePenalty;
	LengthType dummyNodeStart = 0;
	LengthType dummyNodeEnd = 1;
	bool finalized;
};

#endif