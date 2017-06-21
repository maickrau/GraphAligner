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
#include "2dArray.h"
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

template <typename LengthType, typename ScoreType>
class GraphAligner
{
public:
	class AlignmentResult
	{
	public:
		AlignmentResult(vg::Alignment alignment, int maxDistanceFromBand, bool alignmentFailed) :
		alignment(alignment),
		maxDistanceFromBand(maxDistanceFromBand),
		alignmentFailed(alignmentFailed)
		{
		}
		vg::Alignment alignment;
		int maxDistanceFromBand;
		bool alignmentFailed;
	};
	typedef std::pair<LengthType, LengthType> MatrixPosition;
	class MatrixSlice
	{
	public:
		std::vector<ScoreType> M;
		std::vector<ScoreType> Q;
		std::vector<ScoreType> R;
		std::vector<MatrixPosition> Rbacktrace;
		std::vector<MatrixPosition> Qbacktrace;
		std::vector<LengthType> maxScorePositionPerRow;
	};
	class SeedHit
	{
	public:
		SeedHit(size_t seqPos, int nodeId, size_t nodePos) : sequencePosition(seqPos), nodeId(nodeId), nodePos(nodePos) {};
		size_t sequencePosition;
		int nodeId;
		size_t nodePos;
	};

	GraphAligner() :
	nodeStart(),
	indexToNode(),
	nodeLookup(),
	nodeIDs(),
	inNeighbors(),
	nodeSequences(),
	gapStartPenalty(1),
	gapContinuePenalty(1)
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
		finalized = true;
	}

	AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence, int startBandWidth, int dynamicWidth, const std::vector<SeedHit>& seedHits, bool initialFullBand) const
	{
		assert(finalized);
		auto seedHitsInMatrix = getSeedHitPositionsInMatrix(sequence, seedHits);
		auto trace = getBacktrace(sequence, startBandWidth, dynamicWidth, seedHitsInMatrix, initialFullBand);
		//failed alignment, don't output
		if (std::get<0>(trace) == std::numeric_limits<ScoreType>::min()) return emptyAlignment();
		auto result = traceToAlignment(seq_id, sequence, std::get<0>(trace), std::get<2>(trace), std::get<1>(trace));
		return result;
	}

	size_t SizeInBp()
	{
		return nodeSequences.size();
	}

private:

	AlignmentResult emptyAlignment() const
	{
		vg::Alignment result;
		result.set_score(std::numeric_limits<decltype(result.score())>::min());
		return AlignmentResult { result, 0, true };
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

	AlignmentResult traceToAlignment(const std::string& seq_id, const std::string& sequence, ScoreType score, const std::vector<MatrixPosition>& trace, int maxDistanceFromBand) const
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
		return AlignmentResult { result, maxDistanceFromBand, false };
	}

	template <bool distanceMatrixOrder, typename SparseMatrixType, typename MatrixType>
	std::tuple<ScoreType, int, std::vector<MatrixPosition>> backtrace(const std::vector<ScoreType>& Mslice, const SparseMatrixType& backtraceMatrix, const MatrixType& band, int sequenceLength, const Array2D<LengthType, distanceMatrixOrder>& distanceMatrix, const std::vector<MatrixPosition>& seedHits, const std::vector<LengthType>& maxScorePositionPerRow) const
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
			//the rows 0-100 don't use the dynamic band, don't include them here
			if (currentPosition.second > 100)
			{
				maxMinDistance = std::max(maxMinDistance, bandDistanceFromSeqToSeq(currentPosition.first, maxScorePositionPerRow[currentPosition.second], distanceMatrix));
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
		result.resize(sequenceLength+1);
		for (size_t j = 0; j < maxRow; j++)
		{
			std::set<LengthType> rowResult;
			rowResult.insert(forwardResult[j].begin(), forwardResult[j].end());
			rowResult.insert(backwardResult[j].begin(), backwardResult[j].end());
			result[j].insert(result[j].end(), rowResult.begin(), rowResult.end());
		}
		return result;
	}

	template <typename MatrixType>
	std::pair<bool, std::vector<LengthType>> getProcessableColumns(const MatrixType& matrix, LengthType j, std::vector<bool>& rowBand) const
	{
		std::vector<LengthType> result;
		std::vector<LengthType> inOrder;
		bool hasWrongOrders = false;
		result.reserve(matrix.rowSize(j));
		inOrder.reserve(matrix.rowSize(j));
		for (auto iter = matrix.rowStart(j); iter != matrix.rowEnd(j); ++iter)
		{
			auto w = *iter;
			if (w == dummyNodeStart || w == dummyNodeEnd) continue;
			rowBand[w] = true;
			auto nodeIndex = indexToNode[w];
			if (nodeStart[nodeIndex] == w && notInOrder[nodeIndex])
			{
				result.push_back(w);
				hasWrongOrders = true;
			}
			else
			{
				inOrder.push_back(w);
			}
		}
		std::sort(inOrder.begin(), inOrder.end());
		result.insert(result.end(), inOrder.begin(), inOrder.end());
		return std::make_pair(hasWrongOrders, result);
	}

	template <typename MatrixType>
	void expandBandDynamically(MatrixType& band, LengthType previousMaximumIndex, LengthType j, LengthType dynamicWidth) const
	{
		assert(j < band.sizeRows());
		assert(previousMaximumIndex < nodeSequences.size());
		auto nodeIndex = indexToNode[previousMaximumIndex];
		LengthType end = nodeEnd[nodeIndex];
		LengthType start = nodeStart[nodeIndex];
		assert(end > previousMaximumIndex);
		assert(start <= previousMaximumIndex);
		LengthType blockStart = start;
		LengthType blockEnd = end-1; //-1 because block indices are inclusive but node ends exclusive
		//previousMaximumIndex - dynamicWidth > blockStart, change order to prevent underflow
		if (previousMaximumIndex > dynamicWidth + blockStart) blockStart = previousMaximumIndex - dynamicWidth;
		if (previousMaximumIndex + dynamicWidth < blockEnd) blockEnd = previousMaximumIndex + dynamicWidth;
		band.setBlock(blockStart, blockEnd, j);
		if (dynamicWidth > end - previousMaximumIndex)
		{
			for (size_t i = 0; i < outNeighbors[nodeIndex].size(); i++)
			{
				assert(nodeStart[outNeighbors[nodeIndex][i]] < nodeSequences.size());
				if (band(nodeStart[outNeighbors[nodeIndex][i]], j)) continue;
				expandDynamicBandForward(band, outNeighbors[nodeIndex][i], j, dynamicWidth - (end - previousMaximumIndex));
			}
		}
		if (dynamicWidth > previousMaximumIndex - start)
		{
			for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
			{
				assert(nodeEnd[inNeighbors[nodeIndex][i]]-1 < nodeSequences.size());
				if (band(nodeEnd[inNeighbors[nodeIndex][i]]-1, j)) continue;
				expandDynamicBandBackward(band, inNeighbors[nodeIndex][i], j, dynamicWidth - (previousMaximumIndex - start));
			}
		}
	}

	template <typename MatrixType>
	void expandDynamicBandBackward(MatrixType& band, LengthType nodeIndex, LengthType j, LengthType dynamicWidth) const
	{
		if (dynamicWidth == 0) return;
		assert(nodeIndex < nodeStart.size());
		auto start = nodeStart[nodeIndex];
		auto end = nodeEnd[nodeIndex];
		auto blockEnd = end - 1;
		auto blockStart = start;
		if (blockEnd - blockStart > dynamicWidth) blockStart = blockEnd - dynamicWidth;
		band.setBlock(blockStart, blockEnd, j);
		for (size_t i = 0; i < outNeighbors[nodeIndex].size(); i++)
		{
			assert(nodeStart[outNeighbors[nodeIndex][i]] < nodeSequences.size());
			if (band(nodeStart[outNeighbors[nodeIndex][i]], j)) continue;
			expandDynamicBandForward(band, outNeighbors[nodeIndex][i], j, dynamicWidth - 1);
		}
		if (dynamicWidth > end - start)
		{
			for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
			{
				assert(nodeEnd[inNeighbors[nodeIndex][i]]-1 < nodeSequences.size());
				if (band(nodeEnd[inNeighbors[nodeIndex][i]]-1, j)) continue;
				expandDynamicBandBackward(band, inNeighbors[nodeIndex][i], j, dynamicWidth - (end - start));
			}
		}
	}

	template <typename MatrixType>
	void expandDynamicBandForward(MatrixType& band, LengthType nodeIndex, LengthType j, LengthType dynamicWidth) const
	{
		if (dynamicWidth == 0) return;
		assert(nodeIndex < nodeStart.size());
		auto end = nodeEnd[nodeIndex];
		auto start = nodeStart[nodeIndex];
		LengthType blockEnd = end - 1;
		if (end - start > dynamicWidth) end = start + dynamicWidth;
		band.setBlock(start, blockEnd, j);
		for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
		{
			assert(nodeEnd[inNeighbors[nodeIndex][i]]-1 < nodeSequences.size());
			if (band(nodeEnd[inNeighbors[nodeIndex][i]]-1, j)) continue;
			expandDynamicBandBackward(band, inNeighbors[nodeIndex][i], j, dynamicWidth - 1);
		}
		if (dynamicWidth > end - start)
		{
			for (size_t i = 0; i < outNeighbors[nodeIndex].size(); i++)
			{
				assert(nodeStart[outNeighbors[nodeIndex][i]] < nodeSequences.size());
				if (band(nodeStart[outNeighbors[nodeIndex][i]], j)) continue;
				expandDynamicBandForward(band, outNeighbors[nodeIndex][i], j, dynamicWidth - (nodeEnd[nodeIndex] - nodeStart[nodeIndex]));
			}
		}
	}

	template<bool distanceMatrixOrder, typename SparseMatrixType, typename MatrixType>
	MatrixSlice getScoreAndBacktraceMatrix(const std::string& sequence, bool hasWrongOrders, const Array2D<LengthType, distanceMatrixOrder>& distanceMatrix, MatrixSlice& previous, int dynamicWidth, MatrixType& band, SparseMatrixType& backtrace) const
	{
		std::vector<bool> band1;
		std::vector<bool> band2;
		band1.resize(nodeSequences.size(), false);
		band2.resize(nodeSequences.size(), false);
		std::vector<ScoreType> M1;
		std::vector<ScoreType> M2;
		std::vector<ScoreType> Q1;
		std::vector<ScoreType> Q2;
		std::vector<ScoreType> R1;
		std::vector<ScoreType> R2;
		std::vector<MatrixPosition> Rbacktrace1;
		std::vector<MatrixPosition> Rbacktrace2;
		std::vector<LengthType> maxScorePositionPerRow;
		maxScorePositionPerRow.resize(sequence.size()+1, 0);
		assert(previous.M.size() == nodeSequences.size());
		assert(previous.R.size() == nodeSequences.size());
		assert(previous.Q.size() == nodeSequences.size());
		assert(previous.Rbacktrace.size() == nodeSequences.size());
		assert(previous.Qbacktrace.size() == nodeSequences.size());
		MatrixSlice result;
		std::vector<MatrixPosition> Qbacktrace;
		M1.resize(nodeSequences.size());
		Q1.resize(nodeSequences.size());
		R1.resize(nodeSequences.size());
		Rbacktrace1.resize(nodeSequences.size());
		std::vector<ScoreType>& currentM = M1;
		std::vector<ScoreType>& previousM = M2;
		std::vector<ScoreType>& currentQ = Q1;
		std::vector<ScoreType>& previousQ = Q2;
		std::vector<ScoreType>& currentR = R1;
		std::vector<ScoreType>& previousR = R2;
		std::vector<MatrixPosition>& currentRbacktrace = Rbacktrace1;
		std::vector<MatrixPosition>& previousRbacktrace = Rbacktrace2;
		std::vector<bool>& currentBand = band1;
		std::vector<bool>& previousBand = band2;
		previousM = std::move(previous.M);
		previousQ = std::move(previous.Q);
		previousR = std::move(previous.R);
		Qbacktrace = std::move(previous.Qbacktrace);
		previousRbacktrace = std::move(previous.Rbacktrace);
		currentR[dummyNodeStart] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
		previousR[dummyNodeStart] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
		currentM[dummyNodeStart] = -gapPenalty(1);
		previousM[dummyNodeStart] = 0;
		currentR[dummyNodeEnd] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
		previousR[dummyNodeEnd] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
		currentM[dummyNodeEnd] = -gapPenalty(sequence.size() - 1);
		previousM[dummyNodeEnd] = -gapPenalty(sequence.size());
		auto previousProcessableColumnsAndOrder = getProcessableColumns(band, 0, previousBand);
		LengthType previousRowMaximumIndex = 0;

		for (LengthType j = 1; j < sequence.size()+1; j++)
		{
			if (j >= 100)
			{
				assert(band(previousRowMaximumIndex, j-1));
				expandBandDynamically(band, previousRowMaximumIndex, j, dynamicWidth);
			}
			auto currentProcessableColumnsAndOrder = getProcessableColumns(band, j, currentBand);
			auto& previousProcessableColumns = previousProcessableColumnsAndOrder.second;
			auto& currentProcessableColumns = currentProcessableColumnsAndOrder.second;
			bool hasWrongOrders = currentProcessableColumnsAndOrder.first;
			currentM[dummyNodeStart] = -gapPenalty(j);
			currentR[dummyNodeStart] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
			backtrace.set(dummyNodeStart, j, std::make_pair(dummyNodeStart, j-1));
			LengthType maxScorePosition = dummyNodeStart;
			ScoreType maxScore = currentM[dummyNodeStart];
			std::vector<std::pair<LengthType, ScoreType>> Rhelper;
			if (hasWrongOrders) Rhelper = getRHelper(j, previousM, sequence, previousBand, previousProcessableColumns);

			for (LengthType w : currentProcessableColumns)
			{
				assert(band(w, j));
				assert(currentBand[w]);
				bool neighborInsideBand = hasInNeighborInsideBand(w, j, currentBand);
				auto nodeIndex = indexToNode[w];
				currentQ[w] = previousQ[w] - gapContinuePenalty;
				bool rCalculated = false;
				if (previousM[w] - gapPenalty(1) > currentQ[w])
				{
					currentQ[w] = previousM[w] - gapPenalty(1);
					Qbacktrace[w] = std::make_pair(w, j-1);
				}
				if (w == nodeStart[nodeIndex] && notInOrder[nodeIndex])
				{
					if (std::any_of(Rhelper.begin(), Rhelper.end(), [w](auto& x) { return x.first != w; }))
					{
						rCalculated = true;
						assert(hasWrongOrders);
						auto rr = fullR(w, j, Rhelper, distanceMatrix);
						currentR[w] = rr.first;
						currentRbacktrace[w] = rr.second;
						assert(currentRbacktrace[w].second < j || (currentRbacktrace[w].second == j && currentRbacktrace[w].first < w));
					}
				}
				else
				{
					if (neighborInsideBand && previousProcessableColumns.size() > 2)
					{
						rCalculated = true;
						auto rr = recurrenceR(w, j, currentM, currentR, currentRbacktrace, currentBand);
						currentR[w] = rr.first;
						currentRbacktrace[w] = rr.second;
						assert(currentRbacktrace[w].second < j || (currentRbacktrace[w].second == j && currentRbacktrace[w].first < w));
					}
				}
				//implicitly handle edges from dummy node by initializing M as coming from the dummy node
				currentM[w] = previousM[dummyNodeStart] + matchScore(nodeSequences[w], sequence[j-1]);
				MatrixPosition foundBacktrace = std::make_pair(dummyNodeStart, j-1);
				if (previousBand[w] && currentQ[w] > currentM[w])
				{
					foundBacktrace = Qbacktrace[w];
					assert(foundBacktrace.second < j || (foundBacktrace.second == j && foundBacktrace.first < w));
					currentM[w] = currentQ[w];
				}
				if (rCalculated)
				{
					if (currentR[w] > currentM[w])
					{
						currentM[w] = currentR[w];
						foundBacktrace = currentRbacktrace[w];
						assert(foundBacktrace.second < j || (foundBacktrace.second == j && foundBacktrace.first < w));
					}
				}
				if (w == nodeStart[nodeIndex])
				{
					for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
					{
						auto u = nodeEnd[inNeighbors[nodeIndex][i]]-1;
						if (!previousBand[u]) continue;
						//-1 because the rows in the DP matrix are one-based, eg. M[w][1] is the _first_ nucleotide of the read (sequence[0])
						if (previousM[u]+matchScore(nodeSequences[w], sequence[j - 1]) > currentM[w])
						{
							currentM[w] = previousM[u]+matchScore(nodeSequences[w], sequence[j - 1]);
							foundBacktrace = std::make_pair(u, j-1);
							assert(foundBacktrace.second < j || (foundBacktrace.second == j && foundBacktrace.first < w));
						}
					}
				}
				else
				{
					LengthType u = w-1;
					if (previousBand[u])
					{
						//-1 because the rows in the DP matrix are one-based, eg. M[w][1] is the _first_ nucleotide of the read (sequence[0])
						if (previousM[u]+matchScore(nodeSequences[w], sequence[j - 1]) > currentM[w])
						{
							currentM[w] = previousM[u]+matchScore(nodeSequences[w], sequence[j - 1]);
							foundBacktrace = std::make_pair(u, j-1);
							assert(foundBacktrace.second < j || (foundBacktrace.second == j && foundBacktrace.first < w));
						}
					}
				}
				//if the previous row was not inside the band, initialize Q as the current M
				if (!previousBand[w])
				{
					currentQ[w] = currentM[w];
					Qbacktrace[w] = std::make_pair(w, j);
				}
				//if R was unavaliable, initialize it as current M
				if (!rCalculated)
				{
					currentR[w] = currentM[w];
					currentRbacktrace[w] = std::make_pair(w, j);
				}
				assert(currentM[w] >= -std::numeric_limits<ScoreType>::min() + 100);
				assert(currentM[w] <= std::numeric_limits<ScoreType>::max() - 100);
				backtrace.set(w, j, foundBacktrace);
				assert(foundBacktrace.second < j || (foundBacktrace.second == j && foundBacktrace.first < w));
				if (currentM[w] > maxScore)
				{
					maxScore = currentM[w];
					maxScorePosition = w;
				}
			}

			currentM[dummyNodeEnd] = maxScore - gapPenalty(sequence.size() - j);
			backtrace.set(dummyNodeEnd, j, std::make_pair(maxScorePosition, j));

			previousRowMaximumIndex = maxScorePosition;
			maxScorePositionPerRow[j] = maxScorePosition;

			for (auto w : previousProcessableColumns)
			{
				previousBand[w] = false;
			}

			std::swap(currentBand, previousBand);
			std::swap(currentM, previousM);
			std::swap(currentQ, previousQ);
			std::swap(currentR, previousR);
			std::swap(currentRbacktrace, previousRbacktrace);
			previousProcessableColumnsAndOrder = std::move(currentProcessableColumnsAndOrder);
		}
		result.Qbacktrace = Qbacktrace;
		//use previous instead of current because the last line swapped them
		result.M = std::move(previousM);
		result.Q = std::move(previousQ);
		result.R = std::move(previousR);
		result.Rbacktrace = std::move(previousRbacktrace);
		result.maxScorePositionPerRow = std::move(maxScorePositionPerRow);
		return result;
	}

	SparseBoolMatrix<SliceRow<LengthType>> getBandedRows(const std::vector<MatrixPosition>& seedHits, int bandWidth, size_t sequenceLength) const
	{
		auto bandLocations = getBandLocations(sequenceLength, seedHits, 100);
		SparseBoolMatrix<SliceRow<LengthType>> result {nodeSequences.size(), sequenceLength+1};
		for (LengthType j = 0; j < 100 && j < sequenceLength+1; j++)
		{
			for (size_t i = 0; i < bandLocations[j].size(); i++)
			{
				expandBandDynamically(result, bandLocations[j][i], j, bandWidth);
			}
		}
		for (LengthType j = 0; j < sequenceLength+1; j++)
		{
			result.set(dummyNodeStart, j);
			result.set(dummyNodeEnd, j);
		}
		return result;
	}

	SparseBoolMatrix<SliceRow<LengthType>> getFullBand(size_t sequenceLength) const
	{
		SparseBoolMatrix<SliceRow<LengthType>> result {nodeSequences.size(), sequenceLength+1};
		for (LengthType j = 0; j < 100 && j < sequenceLength+1; j++)
		{
			result.setBlock(0, nodeSequences.size(), j);
		}
		return result;
	}

	std::tuple<ScoreType, int, std::vector<MatrixPosition>> getBacktrace(const std::string& sequence, int startBandWidth, int dynamicWidth, const std::vector<MatrixPosition>& seedHits, bool initialFullBand) const
	{
		SparseBoolMatrix<SliceRow<LengthType>> band {1, 1};
		if (initialFullBand)
		{
			band = getFullBand(sequence.size());
		}
		else
		{
			band = getBandedRows(seedHits, startBandWidth, sequence.size());
		}
		auto distanceMatrix = getDistanceMatrixBoostJohnson();
		bool hasWrongOrders = false;
		SparseMatrix<MatrixPosition, decltype(band)> backtraceMatrix {nodeSequences.size(), sequence.size() + 1, band};
		MatrixSlice lastRow = getFirstSlice(backtraceMatrix);
		auto slice = getScoreAndBacktraceMatrix(sequence, hasWrongOrders, distanceMatrix, lastRow, dynamicWidth, band, backtraceMatrix);
		auto result = backtrace(slice.M, backtraceMatrix, band, sequence.size(), distanceMatrix, seedHits, slice.maxScorePositionPerRow);
		return result;
	}

	template <typename SparseMatrixType>
	MatrixSlice getFirstSlice(SparseMatrixType& backtrace) const
	{
		MatrixSlice result;
		result.M.resize(nodeSequences.size(), 0);
		result.R.resize(nodeSequences.size(), 0);
		result.Q.resize(nodeSequences.size(), 0);
		result.Rbacktrace.reserve(nodeSequences.size());
		result.Qbacktrace.reserve(nodeSequences.size());
		for (LengthType i = 0; i < nodeSequences.size(); i++)
		{
			result.Qbacktrace.emplace_back(i, 0);
			result.Rbacktrace.emplace_back(i, 0);
		}
		result.R[0] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
		assert(result.M.size() == nodeSequences.size());
		assert(result.R.size() == nodeSequences.size());
		assert(result.Q.size() == nodeSequences.size());
		assert(result.Rbacktrace.size() == nodeSequences.size());
		assert(result.Qbacktrace.size() == nodeSequences.size());
		return result;
	}

	std::vector<std::pair<LengthType, ScoreType>> getRHelperZero() const
	{
		std::vector<std::pair<LengthType, ScoreType>> result;
		for (LengthType v = 0; v < nodeSequences.size(); v++)
		{
			result.emplace_back(v, 0);
		}
		return result;
	}

	std::vector<std::pair<LengthType, ScoreType>> getRHelperOne() const
	{
		std::vector<std::pair<LengthType, ScoreType>> result;
		for (LengthType v = 0; v < nodeSequences.size(); v++)
		{
			result.emplace_back(v, 0);
		}
		return result;
	}

	std::vector<std::pair<LengthType, ScoreType>> getRHelper(LengthType j, const std::vector<ScoreType>& previousM, const std::string& sequence, const std::vector<bool>& previousBand, const std::vector<LengthType>& previousProcessableColumns) const
	{
		if (j == 0) return getRHelperZero();
		if (j == 1) return getRHelperOne();
		std::vector<std::tuple<LengthType, ScoreType, ScoreType>> bestPerNode;
		bestPerNode.resize(nodeStart.size(), std::make_tuple(0, std::numeric_limits<ScoreType>::min() + 99, 0));
		for (auto v : previousProcessableColumns)
		{
			auto nodeIndex = indexToNode[v];
			if (nodeStart[nodeIndex] == v)
			{
				for (size_t neighbori = 0; neighbori < inNeighbors[nodeIndex].size(); neighbori++)
				{
					LengthType u = nodeEnd[inNeighbors[nodeIndex][neighbori]]-1;
					if (!previousBand[u]) continue;
					auto scoreHere = previousM[u] + matchScore(nodeSequences[v], sequence[j-1]);
					if (scoreHere - (ScoreType)(nodeEnd[nodeIndex] - v) * gapContinuePenalty > std::get<1>(bestPerNode[nodeIndex]) - std::get<2>(bestPerNode[nodeIndex]))
					{
						bestPerNode[nodeIndex] = std::make_tuple(v, scoreHere, (ScoreType)(nodeEnd[nodeIndex] - v) * gapContinuePenalty);
					}
				}
			}
			else
			{
				LengthType u = v-1;
				if (!previousBand[u]) continue;
				auto scoreHere = previousM[u] + matchScore(nodeSequences[v], sequence[j-1]);
				if (scoreHere - (ScoreType)(nodeEnd[nodeIndex] - v) * gapContinuePenalty > std::get<1>(bestPerNode[nodeIndex]) - std::get<2>(bestPerNode[nodeIndex]))
				{
					bestPerNode[nodeIndex] = std::make_tuple(v, scoreHere, (ScoreType)(nodeEnd[nodeIndex] - v) * gapContinuePenalty);
				}
			}
		}
		std::vector<std::pair<LengthType, ScoreType>> result;
		for (size_t i = 0; i < bestPerNode.size(); i++)
		{
			if (std::get<1>(bestPerNode[i]) > std::numeric_limits<ScoreType>::min() + 100)
			{
				result.emplace_back(std::get<0>(bestPerNode[i]), std::get<1>(bestPerNode[i]));
			}
		}
		assert(result.size() >= 1);
		return result;
	}

	bool hasInNeighborInsideBand(LengthType w, LengthType j, const std::vector<bool>& currentBand) const
	{
		auto nodeIndex = indexToNode[w];
		if (nodeStart[nodeIndex] == w)
		{
			for (size_t neighborI = 0; neighborI < inNeighbors[nodeIndex].size(); neighborI++)
			{
				if (currentBand[nodeEnd[inNeighbors[nodeIndex][neighborI]] - 1]) return true;
			}
		}
		else
		{
			return currentBand[w-1];
		}
		return false;
	}

	//compute R using the recurrence on page 3
	std::pair<ScoreType, MatrixPosition> recurrenceR(LengthType w, LengthType j, const std::vector<ScoreType>& currentM, const std::vector<ScoreType>& currentR, const std::vector<MatrixPosition>& currentRbacktrace, const std::vector<bool>& currentBand) const
	{
		assert(currentBand[w]);
		auto nodeIndex = indexToNode[w];
		assert(nodeStart[nodeIndex] != w || !notInOrder[nodeIndex]);
		MatrixPosition pos;
		ScoreType maxValue = std::numeric_limits<ScoreType>::min() + 99;
		if (nodeStart[nodeIndex] == w)
		{
			for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
			{
				auto u = nodeEnd[inNeighbors[nodeIndex][i]]-1;
				if (!currentBand[u]) continue;
				assert(u < w);
				if (currentM[u] - gapPenalty(1) > maxValue)
				{
					maxValue = currentM[u] - gapPenalty(1);
					pos = std::make_pair(u, j);
				}
				if (currentR[u] - gapContinuePenalty > maxValue)
				{
					maxValue = currentR[u] - gapContinuePenalty;
					pos = currentRbacktrace[u];
				}
			}
		}
		else
		{
			auto u = w-1;
			if (currentBand[u])
			{
				pos = currentRbacktrace[u];
				maxValue = currentR[u] - gapContinuePenalty;
				if (currentM[u] - gapPenalty(1) > maxValue)
				{
					pos = std::make_pair(u, j);
					maxValue = currentM[u] - gapPenalty(1);
				}
			}
		}
		assert(maxValue >= -std::numeric_limits<ScoreType>::min() + 100);
		assert(maxValue <= std::numeric_limits<ScoreType>::max() - 100);
		return std::make_pair(maxValue, pos);
	}

	//compute R using the slow, full definition on page 3
	template <bool distanceMatrixOrder>
	std::pair<ScoreType, MatrixPosition> fullR(LengthType w, LengthType j, const std::vector<std::pair<LengthType, ScoreType>>& RHelper, const Array2D<LengthType, distanceMatrixOrder>& distanceMatrix) const
	{
		assert(j > 0);
		assert(w > 0);
		auto nodeIndex = indexToNode[w];
		assert(nodeStart[nodeIndex] == w && notInOrder[nodeIndex]);
		MatrixPosition pos;
		ScoreType maxValue = std::numeric_limits<ScoreType>::min() + 99;
		for (auto pair : RHelper)
		{
			auto v = pair.first;
			if (v == w) continue;
			auto scoreHere = pair.second - gapPenalty(distanceFromSeqToSeq(v, w, distanceMatrix));
			if (scoreHere > maxValue)
			{
				maxValue = scoreHere;
				pos = std::make_pair(v, j-1);
			}
		}
		assert(maxValue >= -std::numeric_limits<ScoreType>::min() + 100);
		assert(maxValue <= std::numeric_limits<ScoreType>::max() - 100);
		return std::make_pair(maxValue, pos);
	}

	template <bool distanceMatrixOrder>
	LengthType bandDistanceFromSeqToSeq(LengthType start, LengthType end, const Array2D<LengthType, distanceMatrixOrder>& distanceMatrix) const
	{
		if (start == end) return 0;
		if (start == dummyNodeStart || start == dummyNodeEnd || end == dummyNodeStart || end == dummyNodeEnd) return 1;
		auto startNode = indexToNode[start];
		auto endNode = indexToNode[end];
		if (startNode == endNode) return std::min(end - start, start - end);
		if (distanceMatrix(startNode, endNode) == nodeEnd[startNode]-nodeStart[startNode])
		{
			return nodeEnd[startNode] - start + end - nodeStart[endNode];
		}
		if (distanceMatrix(endNode, startNode) == nodeEnd[endNode]-nodeStart[endNode])
		{
			return nodeEnd[endNode] - end + start - nodeStart[startNode];
		}
		LengthType minDistance = nodeSequences.size();
		for (size_t i = 0; i < distanceMatrix.sizeRows(); i++)
		{
			LengthType distanceFromStartToMid = distanceMatrix(startNode, i) + nodeStart[startNode] - start;
			if (distanceMatrix(i, startNode) + start - nodeStart[startNode] < distanceFromStartToMid)
			{
				distanceFromStartToMid = distanceMatrix(i, startNode) + start - nodeStart[startNode];
			}
			LengthType distanceFromMidToEnd = distanceMatrix(i, endNode) + end - nodeStart[endNode];
			if (distanceMatrix(endNode, i) + nodeStart[endNode] - end < distanceFromMidToEnd)
			{
				distanceFromMidToEnd = distanceMatrix(endNode, i) + nodeStart[endNode] - end;
			}
			minDistance = std::min(minDistance, distanceFromStartToMid + distanceFromMidToEnd);

			LengthType distanceFromStartToMidEnd = distanceMatrix(startNode, i) + nodeStart[startNode] - start + nodeEnd[i] - nodeStart[i];
			if (distanceMatrix(i, startNode) - (nodeEnd[i] - nodeStart[i]) + start - nodeStart[startNode] < distanceFromStartToMidEnd)
			{
				distanceFromStartToMidEnd = distanceMatrix(i, startNode) - (nodeEnd[i] - nodeStart[i]) + start - nodeStart[startNode];
			}
			LengthType distanceFromMidEndToEnd = distanceMatrix(i, endNode) - (nodeEnd[i] - nodeStart[i]) + end - nodeStart[endNode];
			if (distanceMatrix(endNode, i) + (nodeEnd[i] - nodeStart[i]) + nodeStart[endNode] - end < distanceFromMidToEnd)
			{
				distanceFromMidToEnd = distanceMatrix(endNode, i) + (nodeEnd[i] - nodeStart[i]) + nodeStart[endNode] - end;
			}
			minDistance = std::min(minDistance, distanceFromStartToMidEnd + distanceFromMidEndToEnd);
		}
		return minDistance;
	}

	template <bool distanceMatrixOrder>
	LengthType distanceFromSeqToSeq(LengthType start, LengthType end, const Array2D<LengthType, distanceMatrixOrder>& distanceMatrix) const
	{
		if (start == end) return 0;
		if (start == dummyNodeStart || start == dummyNodeEnd || end == dummyNodeStart || end == dummyNodeEnd) return 1;
		auto startNode = indexToNode[start];
		auto endNode = indexToNode[end];
		if (startNode == endNode && end >= start) return end - start;
		return distanceMatrix(startNode, endNode) + nodeStart[startNode] + end - nodeStart[endNode] - start;
	}

	void fillDistanceMatrixChain(std::vector<size_t>& chainStart, std::vector<LengthType>& distanceAlongChain, std::vector<std::tuple<LengthType, LengthType, LengthType>>& edges, size_t index) const
	{
		//starts inside a bubble
		if (inNeighbors[index].size() == 1 && outNeighbors[index].size() == 1)
		{
			auto before = inNeighbors[index][0];
			auto after = outNeighbors[index][0];
			if (outNeighbors[before].size() == 2 && inNeighbors[after].size() == 2)
			{
				auto first = outNeighbors[before][0];
				auto second = outNeighbors[before][1];
				if (inNeighbors[first].size() == 1 && outNeighbors[first].size() == 1)
				{
					if (inNeighbors[second].size() == 1 && outNeighbors[second].size() == 1)
					{
						if (inNeighbors[first][0] == inNeighbors[second][0] && outNeighbors[first][0] == outNeighbors[second][0])
						{
							index = inNeighbors[index][0];
						}
					}
				}
			}
		}
		while (true)
		{
			//just a chain
			if (inNeighbors[index].size() == 1 && outNeighbors[inNeighbors[index][0]].size() == 1)
			{
				index = inNeighbors[index][0];
				continue;
			}
			//a simple bubble
			// if (inNeighbors[index].size() == 2 && inNeighbors[inNeighbors[index][0]].size() == 1 && outNeighbors[inNeighbors[index][0]].size() == 1 && inNeighbors[inNeighbors[index][1]].size() == 1 && outNeighbors[inNeighbors[index][1]].size() == 1 && inNeighbors[inNeighbors[index][0]][0] == inNeighbors[inNeighbors[index][1]][0] && outNeighbors[inNeighbors[inNeighbors[index][0]][0]].size() == 2)
			if (inNeighbors[index].size() == 2)
			{
				auto first = inNeighbors[index][0];
				auto second = inNeighbors[index][1];
				if (inNeighbors[first].size() == 1 && outNeighbors[first].size() == 1)
				{
					if (inNeighbors[second].size() == 1 && outNeighbors[second].size() == 1)
					{
						assert(outNeighbors[first][0] == index);
						assert(outNeighbors[second][0] == index);
						if (inNeighbors[first][0] == inNeighbors[second][0])
						{
							auto before = inNeighbors[first][0];
							if (outNeighbors[before].size() == 2)
							{
								index = before;
								continue;
							}
						}
					}
				}
				// index = inNeighbors[inNeighbors[index][0]][0];
				// continue;
			}
			break;
		}
		auto start = index;
		assert(chainStart[index] == std::numeric_limits<size_t>::max());
		distanceAlongChain[index] = 0;
		chainStart[index] = start;
		LengthType pathLength = 0;
		pathLength += nodeEnd[index] - nodeStart[index];

		while (true)
		{
			//just a chain
			if (outNeighbors[index].size() == 1 && inNeighbors[outNeighbors[index][0]].size() == 1)
			{
				index = outNeighbors[index][0];
				assert(chainStart[index] == std::numeric_limits<size_t>::max());
				distanceAlongChain[index] = pathLength;
				chainStart[index] = start;
				pathLength += nodeEnd[index] - nodeStart[index];
				continue;
			}
			//a simple bubble
			// if (outNeighbors[index].size() == 2 && outNeighbors[outNeighbors[index][0]].size() == 1 && inNeighbors[outNeighbors[index][0]].size() == 1 && outNeighbors[outNeighbors[index][1]].size() == 1 && inNeighbors[outNeighbors[index][1]].size() == 1 && outNeighbors[outNeighbors[index][0]][0] == outNeighbors[outNeighbors[index][1]][0] && inNeighbors[outNeighbors[outNeighbors[index][0]][0]].size() == 2)
			if (outNeighbors[index].size() == 2)
			{
				auto first = outNeighbors[index][0];
				auto second = outNeighbors[index][1];
				if (inNeighbors[first].size() == 1 && outNeighbors[first].size() == 1)
				{
					if (inNeighbors[second].size() == 1 && outNeighbors[second].size() == 1)
					{
						assert(inNeighbors[first][0] == index);
						assert(inNeighbors[second][0] == index);
						if (outNeighbors[first][0] == outNeighbors[second][0])
						{
							auto after = outNeighbors[first][0];
							if (inNeighbors[after].size() == 2)
							{
								assert(chainStart[first] == std::numeric_limits<size_t>::max());
								assert(chainStart[second] == std::numeric_limits<size_t>::max());
								distanceAlongChain[first] = pathLength;
								distanceAlongChain[second] = pathLength;
								chainStart[first] = start;
								chainStart[second] = start;
								chainStart[after] = start;
								assert(nodeEnd[first] > nodeStart[first]);
								assert(nodeEnd[second] > nodeStart[second]);
								auto increase = std::min(nodeEnd[first] - nodeStart[first], nodeEnd[second] - nodeStart[second]);
								assert(increase > 0);
								assert(increase < nodeSequences.size());
								pathLength += increase;
								distanceAlongChain[after] = pathLength;
								pathLength += nodeEnd[after] - nodeStart[after];
								index = after;
								continue;
							}
						}
					}
				}
				// assert(chainStart[first] == std::numeric_limits<size_t>::max());
				// assert(chainStart[second] == std::numeric_limits<size_t>::max());
				// distanceAlongChain[first] = pathLength;
				// distanceAlongChain[second] = pathLength;
				// chainStart[first] = start;
				// chainStart[second] = start;
				// pathLength += std::min(nodeEnd[first] - nodeStart[first], nodeEnd[second] - nodeStart[second]);
				// index = outNeighbors[outNeighbors[index][0]][0];
				// continue;
			}
			break;
		}
		while (outNeighbors[index].size() == 1 && inNeighbors[outNeighbors[index][0]].size() == 1)
		{
			index = outNeighbors[index][0];
			assert(chainStart[index] == std::numeric_limits<size_t>::max());
			distanceAlongChain[index] = pathLength;
			chainStart[index] = start;
			pathLength += nodeEnd[index] - nodeStart[index];
		}

		for (size_t i = 0; i < outNeighbors[index].size(); i++)
		{
			edges.emplace_back(start, outNeighbors[index][i], pathLength);
		}
	}

	//http://www.boost.org/doc/libs/1_40_0/libs/graph/example/johnson-eg.cpp
	Array2D<LengthType, false> getDistanceMatrixBoostJohnson() const
	{
		//g++ can't optimize repeated calls into one call for some reason and it's expensive, so store them like this
		auto inneighborsSize = inNeighbors.size();
		auto nodesequencesSizePlusOne = nodeSequences.size() + 1;
		//Merge chains of nodes into one node, so that boost-johnson has a lower node and edge count
		std::vector<size_t> chainStart;
		std::vector<LengthType> distanceOnChain;
		chainStart.resize(inneighborsSize, std::numeric_limits<size_t>::max());
		distanceOnChain.resize(inneighborsSize, std::numeric_limits<LengthType>::max());
		std::vector<std::tuple<LengthType, LengthType, LengthType>> graphedges;
		Array2D<LengthType, false> distances {inneighborsSize, inneighborsSize, nodesequencesSizePlusOne};
		for (size_t i = 0; i < inneighborsSize; i++)
		{
			if (chainStart[i] == std::numeric_limits<size_t>::max())
			{
				fillDistanceMatrixChain(chainStart, distanceOnChain, graphedges, i);
			}
		}
		std::vector<size_t> actualCalculables;
		std::vector<size_t> helperLookup;
		helperLookup.resize(inneighborsSize, std::numeric_limits<size_t>::max());
		for (size_t i = 0; i < inneighborsSize; i++)
		{
			assert(chainStart[i] != std::numeric_limits<size_t>::max());
			assert(distanceOnChain[i] != std::numeric_limits<LengthType>::max());
			if (chainStart[i] == i) 
			{
				helperLookup[i] = actualCalculables.size();
				actualCalculables.push_back(i);
			}
		}
		auto V = actualCalculables.size();
		adjacency_list<vecS, vecS, directedS, no_property, property<edge_weight_t, int, property<edge_weight2_t, int>>> graph { V };

		std::sort(graphedges.begin(), graphedges.end(), [](auto& left, auto& right) { return std::get<0>(left) < std::get<0>(right) || (std::get<0>(left) == std::get<0>(right) && std::get<1>(left) < std::get<1>(right)); });

		for (size_t i = 0; i < graphedges.size(); i++)
		{
			boost::add_edge(helperLookup[std::get<0>(graphedges[i])], helperLookup[std::get<1>(graphedges[i])], graph);
		}

		property_map<adjacency_list<vecS, vecS, directedS, no_property, property<edge_weight_t, int, property<edge_weight2_t, int>>>, edge_weight_t>::type w = get(edge_weight, graph);
		graph_traits<adjacency_list<vecS, vecS, directedS, no_property, property<edge_weight_t, int, property<edge_weight2_t, int>>>>::edge_iterator e, e_end;
		int edgeindex = 0;
		for (boost::tie(e, e_end) = edges(graph); e != e_end; ++e)
		{
			assert((*e).m_source == helperLookup[std::get<0>(graphedges[edgeindex])]);
			assert((*e).m_target == helperLookup[std::get<1>(graphedges[edgeindex])]);
			w[*e] = std::get<2>(graphedges[edgeindex]);
			edgeindex++;
		}
		std::vector<int> d(V, nodesequencesSizePlusOne);
		int** D;
		D = new int*[V];
		for (size_t i = 0; i < V; i++)
		{
			D[i] = new int[V];
			for (size_t j = 0; j < V; j++)
			{
				D[i][j] = nodesequencesSizePlusOne;
			}
		}
		johnson_all_pairs_shortest_paths(graph, D, distance_map(&d[0]));
		for (size_t i = 0; i < V; i++)
		{
			for (size_t j = 0; j < V; j++)
			{
				//distances have to be positive
				assert(D[i][j] > 0 || (i == j && D[i][j] == 0));
				//distances are either reasonable or infinity
				assert(D[i][j] <= nodesequencesSizePlusOne || D[i][j] == std::numeric_limits<int>::max());
				if (D[i][j] == std::numeric_limits<int>::max())
				{
					D[i][j] = nodesequencesSizePlusOne;
				}
			}
		}
		//make sure that the distance to itself is not 0
		//we need to do this so distance calculation from a later point in the node to an earlier point in the node works correctly
		for (size_t i = 0; i < V; i++)
		{
			D[i][i] = nodesequencesSizePlusOne;
			for (size_t j = 0; j < V; j++)
			{
				if (j == i) continue;
				D[i][i] = std::min(D[i][i], D[i][j] + D[j][i]);
			}
		}

		for (size_t ii = 0; ii < inneighborsSize; ii++)
		{
			assert(distanceOnChain[ii] != std::numeric_limits<size_t>::max());
			assert(helperLookup[chainStart[ii]] != std::numeric_limits<size_t>::max());
			size_t i = helperLookup[chainStart[ii]];
			for (size_t jj = 0; jj < inneighborsSize; jj++)
			{
				assert(distanceOnChain[jj] != std::numeric_limits<size_t>::max());
				assert(helperLookup[chainStart[jj]] != std::numeric_limits<size_t>::max());
				size_t j = helperLookup[chainStart[jj]];
				if (i == j && distanceOnChain[jj] > distanceOnChain[ii])
				{
					distances(ii, jj) = distanceOnChain[jj] - distanceOnChain[ii];
				}
				else if (i == j && distanceOnChain[jj] < distanceOnChain[ii])
				{
					if (D[i][i] >= nodesequencesSizePlusOne)
					{
						distances(ii, jj) = nodesequencesSizePlusOne;
					}
					else
					{
						distances(ii, jj) = D[i][i] + distanceOnChain[jj] - distanceOnChain[ii];
					}
				}
				else if (i == j)
				{
					distances(ii, jj) = D[i][i];
				}
				else if (D[i][j] >= nodesequencesSizePlusOne)
				{
					distances(ii, jj) = nodesequencesSizePlusOne;
				}
				else
				{
					assert(distanceOnChain[ii] < D[i][j]);
					distances(ii, jj) = D[i][j] + distanceOnChain[jj] - distanceOnChain[ii];
				}
				assert(distances(ii, jj) <= nodesequencesSizePlusOne * 2);
				assert(distances(ii, jj) > 0);
				if (distances(ii, jj) > nodesequencesSizePlusOne) distances(ii, jj) = nodesequencesSizePlusOne;
			}
		}
		for (size_t i = 0; i < V; i++)
		{
			delete [] D[i];
		}
		delete [] D;

		return distances;

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