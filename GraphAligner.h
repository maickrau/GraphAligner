#ifndef GraphAligner_H
#define GraphAligner_H

//http://biorxiv.org/content/early/2017/04/06/124941
#include <chrono>
#include <algorithm>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include <boost/config.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include "vg.pb.h"
#include "2dArray.h"

using namespace boost;

template <typename LengthType, typename ScoreType>
class GraphAligner
{
public:
	typedef std::pair<LengthType, LengthType> MatrixPosition;
	class MatrixSlice
	{
	public:
		std::vector<ScoreType> M;
		std::vector<ScoreType> Q;
		std::vector<ScoreType> R;
		std::vector<MatrixPosition> Rbacktrace;
		std::vector<MatrixPosition> Qbacktrace;
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
		//add the dummy node as the first node
		nodeIDs.push_back(0);
		nodeStart.push_back(nodeSequences.size());
		inNeighbors.emplace_back();
		outNeighbors.emplace_back();
		nodeSequences.push_back('N');
		indexToNode.resize(nodeSequences.size(), nodeStart.size()-1);
		nodeEnd.emplace_back(nodeSequences.size());
		notInOrder.push_back(false);
		finalized = false;
	}
	
	void AddNode(int nodeId, std::string sequence)
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
			// std::cerr << "edge from " << from << " to " << to << std::endl;
			// std::cerr << "node ids " << node_id_from << " -> " << node_id_to << std::endl;
		}
	}

	void Finalize()
	{
		//add an edge from the dummy node to all nodes without an in-edge
		for (LengthType i = 0; i < nodeSequences.size(); i++)
		{
			if (nodeStart[indexToNode[i]] == i)
			{
				if (inNeighbors[indexToNode[i]].size() == 0)
				{
					inNeighbors[indexToNode[i]].push_back(0);
				}
			}
		}
		finalized = true;
	}

	vg::Alignment AlignOneWay(const std::string& seq_id, const std::string& sequence, bool reverse, int bandWidth, const std::vector<SeedHit>& seedHits) const
	{
		assert(finalized);
		auto seedHitsInMatrix = getSeedHitPositionsInMatrix(sequence, seedHits);
		auto trace = backtrackWithSquareRootSlices(sequence, bandWidth, seedHitsInMatrix);
		auto result = traceToAlignment(seq_id, trace, reverse);
		return result;
	}

	size_t SizeInBp()
	{
		return nodeSequences.size();
	}

private:

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

	vg::Alignment traceToAlignment(const std::string& seq_id, const std::pair<ScoreType, std::vector<MatrixPosition>>& traceWithScore, bool reverse) const
	{
		auto& trace = traceWithScore.second;
		vg::Alignment result;
		result.set_name(seq_id);
		auto path = new vg::Path;
		result.set_allocated_path(path);
		size_t pos = 0;
		size_t oldNode = indexToNode[trace[0].first];
		while (nodeIDs[oldNode] == 0)
		{
			pos++;
			assert(pos < trace.size());
			oldNode = indexToNode[trace[pos].first];
			assert(oldNode < nodeIDs.size());
		}
		int rank = 0;
		auto vgmapping = path->add_mapping();
		auto position = new vg::Position;
		vgmapping->set_allocated_position(position);
		vgmapping->set_rank(rank);
		position->set_node_id(nodeIDs[oldNode]);
		if (reverse) position->set_is_reverse(true);
		for (; pos < trace.size(); pos++)
		{
			if (indexToNode[trace[pos].first] == oldNode) continue;
			oldNode = indexToNode[trace[pos].first];
			if (nodeIDs[oldNode] == 0) break;
			rank++;
			vgmapping = path->add_mapping();
			position = new vg::Position;
			vgmapping->set_allocated_position(position);
			vgmapping->set_rank(rank);
			position->set_node_id(nodeIDs[oldNode]);
			if (reverse) position->set_is_reverse(true);
		}
		result.set_score(traceWithScore.first);
		return result;
	}

	template <bool backtraceOrder, bool bandOrder, bool distanceMatrixOrder>
	std::pair<ScoreType, std::vector<MatrixPosition>> backtrace(const std::vector<ScoreType>& Mslice, const Array2D<MatrixPosition, backtraceOrder>& backtraceMatrix, const Array2D<bool, bandOrder>& band, int sequenceLength, const Array2D<LengthType, distanceMatrixOrder>& distanceMatrix, const std::vector<MatrixPosition>& seedHits) const
	{
		auto bandLocations = getBandLocations(sequenceLength, seedHits);
		assert(backtraceMatrix.sizeColumns() == sequenceLength+1);
		assert(backtraceMatrix.sizeRows() == nodeSequences.size());
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
			// std::cerr << Mslice[candidatePosition.first] << "\t";
			if (Mslice[candidatePosition.first] > Mslice[currentPosition.first])
			{
				currentPosition = candidatePosition;
			}
		}
		// std::cerr << std::endl;
		assert(band(currentPosition.first, currentPosition.second));
		auto score = Mslice[currentPosition.first];
		trace.push_back(currentPosition);
		LengthType maxMinDistance = 0;
		while (currentPosition.second > 0)
		{
			assert(band(currentPosition.first, currentPosition.second));
			LengthType minDistance = nodeSequences.size();
			// std::cerr << "band locations: ";
			for (size_t i = 0; i < bandLocations[currentPosition.second].size(); i++)
			{
				// std::cerr << bandLocations[currentPosition.second][i] << " ";
				minDistance = std::min(minDistance, distanceFromSeqToSeq(currentPosition.first, bandLocations[currentPosition.second][i], distanceMatrix));
				minDistance = std::min(minDistance, distanceFromSeqToSeq(bandLocations[currentPosition.second][i], currentPosition.first, distanceMatrix));
			}
			// std::cerr << std::endl;
			maxMinDistance = std::max(maxMinDistance, minDistance);
			std::cerr << currentPosition.first << ", " << currentPosition.second << std::endl;
			// std::cerr << "row " << currentPosition.second << " distance from band: " << minDistance << std::endl;
			assert(currentPosition.second >= 0);
			assert(currentPosition.second < sequenceLength+1);
			assert(currentPosition.first >= 0);
			assert(currentPosition.first < nodeSequences.size());
			auto newPos = backtraceMatrix(currentPosition.first, currentPosition.second);
			//If we're at the dummy node, we have to stay there
			assert(currentPosition.first == 0 ? (newPos.first == 0) : true);
			assert(newPos.second < currentPosition.second || (newPos.second == currentPosition.second && newPos.first < currentPosition.first));
			currentPosition = newPos;
			trace.push_back(currentPosition);
		}
		std::cerr << "max distance from band: " << maxMinDistance << std::endl;
		std::reverse(trace.begin(), trace.end());
		return std::make_pair(score, trace);
	}

	void expandBandForwards(std::vector<std::vector<LengthType>>& result, LengthType w, LengthType j, size_t sequenceLength) const
	{
		if (std::find(result[j].begin(), result[j].end(), w) != result[j].end()) return;
		auto nodeIndex = indexToNode[w];
		auto end = nodeEnd[nodeIndex];
		while (w != end && j < sequenceLength+1)
		{
			result[j].emplace_back(w);
			w++;
			j++;
		}
		if (w == end && j < sequenceLength+1)
		{
			for (size_t i = 0; i < outNeighbors[nodeIndex].size(); i++)
			{
				expandBandForwards(result, nodeStart[outNeighbors[nodeIndex][i]], j, sequenceLength);
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
		if (j > 0) 
		{
			result[j].emplace_back(w);
		}
		if (w == start && j > 0)
		{
			for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
			{
				expandBandBackwards(result, nodeEnd[inNeighbors[nodeIndex][i]] - 1, j-1, sequenceLength);
			}
		}
	}

	std::vector<std::vector<LengthType>> getBandLocations(int sequenceLength, const std::vector<MatrixPosition>& seedHits) const
	{
		std::vector<std::vector<LengthType>> forwardResult;
		std::vector<std::vector<LengthType>> backwardResult;
		backwardResult.resize(sequenceLength+1);
		forwardResult.resize(sequenceLength+1);
		backwardResult[0].emplace_back(0);
		forwardResult[0].emplace_back(0);
		for (auto hit : seedHits)
		{
			expandBandForwards(forwardResult, hit.first, hit.second, sequenceLength);
			expandBandBackwards(backwardResult, hit.first, hit.second, sequenceLength);
		}
		std::vector<std::vector<LengthType>> result;
		result.resize(sequenceLength+1);
		for (size_t j = 0; j < sequenceLength+1; j++)
		{
			std::set<LengthType> rowResult;
			rowResult.insert(forwardResult[j].begin(), forwardResult[j].end());
			rowResult.insert(backwardResult[j].begin(), backwardResult[j].end());
			result[j].insert(result[j].end(), rowResult.begin(), rowResult.end());
		}
		return result;
	}

	template<bool distanceMatrixOrder, bool bandOrder, bool backtraceOrder>
	MatrixSlice getScoreAndBacktraceMatrixSlice(const std::string& sequence, bool hasWrongOrders, const std::vector<LengthType>& nodeOrdering, const Array2D<LengthType, distanceMatrixOrder>& distanceMatrix, MatrixSlice& previous, LengthType start, LengthType end, int bandWidth, const Array2D<bool, bandOrder>& band, Array2D<MatrixPosition, backtraceOrder>& backtrace) const
	{
		std::vector<ScoreType> M1;
		std::vector<ScoreType> M2;
		std::vector<ScoreType> Q1;
		std::vector<ScoreType> Q2;
		std::vector<ScoreType> R1;
		std::vector<ScoreType> R2;
		std::vector<LengthType> processableColumns1;
		std::vector<LengthType> processableColumns2;
		std::vector<MatrixPosition> Rbacktrace1;
		std::vector<MatrixPosition> Rbacktrace2;
		assert(previous.M.size() == nodeSequences.size());
		assert(previous.R.size() == nodeSequences.size());
		assert(previous.Q.size() == nodeSequences.size());
		assert(previous.Rbacktrace.size() == nodeSequences.size());
		assert(previous.Qbacktrace.size() == nodeSequences.size());
		MatrixSlice result;
		std::vector<MatrixPosition> Qbacktrace;
		processableColumns1.reserve(nodeSequences.size());
		processableColumns2.reserve(nodeSequences.size());
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
		previousM = std::move(previous.M);
		previousQ = std::move(previous.Q);
		previousR = std::move(previous.R);
		Qbacktrace = std::move(previous.Qbacktrace);
		previousRbacktrace = std::move(previous.Rbacktrace);
		currentR[0] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
		previousR[0] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
		currentM[0] = -gapPenalty(start + 1);
		previousM[0] = -gapPenalty(start);
		std::vector<LengthType>& currentProcessableColumns = processableColumns1;
		std::vector<LengthType>& previousProcessableColumns = processableColumns2;
		for (auto w : nodeOrdering)
		{
			if (band(w, start)) previousProcessableColumns.push_back(w);
		}
		for (LengthType j = 1; j < end - start; j++)
		{
			currentProcessableColumns.clear();
			int insideBand = 0;
			for (LengthType w : nodeOrdering)
			{
				if (band(w, start+j)) 
				{
					currentProcessableColumns.push_back(w);
					insideBand++;
				}
			}
			std::cerr << "inside band: " << insideBand << std::endl;
			std::cerr << "current processable columns: " << currentProcessableColumns.size() << std::endl;
			std::cerr << "previous processable columns: " << previousProcessableColumns.size() << std::endl;
			currentM[0] = -gapPenalty(start + j);
			currentR[0] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
			std::vector<std::pair<LengthType, ScoreType>> Rhelper;
			if (hasWrongOrders) Rhelper = getRHelper(j, start, previousM, sequence, band, previousProcessableColumns);

			for (LengthType w : currentProcessableColumns)
			{
				assert(band(w, j+start));
				bool neighborInsideBand = hasInNeighborInsideBand(w, j, start, band);
				auto nodeIndex = indexToNode[w];
				currentQ[w] = previousQ[w] - gapContinuePenalty;
				if (previousM[w] - gapPenalty(1) > currentQ[w])
				{
					currentQ[w] = previousM[w] - gapPenalty(1);
					Qbacktrace[w] = std::make_pair(w, j-1 + start);
				}
				if (w == nodeStart[nodeIndex] && notInOrder[nodeIndex])
				{
					assert(hasWrongOrders);
					auto rr = fullR(w, j, Rhelper, distanceMatrix, start);
					currentR[w] = rr.first;
					currentRbacktrace[w] = rr.second;
					assert(currentRbacktrace[w].second < (j + start) || (currentRbacktrace[w].second == (j + start) && currentRbacktrace[w].first < w));
				}
				else
				{
					if (neighborInsideBand)
					{
						auto rr = recurrenceR(w, j, start, currentM, currentR, currentRbacktrace, band);
						currentR[w] = rr.first;
						currentRbacktrace[w] = rr.second;
						assert(currentRbacktrace[w].second < (j + start) || (currentRbacktrace[w].second == (j + start) && currentRbacktrace[w].first < w));
					}
				}
				currentM[w] = std::numeric_limits<ScoreType>::min() + 99;
				if (band(w, start+j-1))
				{
					backtrace(w, j) = Qbacktrace[w];
					assert(backtrace(w, j).second < (j + start) || (backtrace(w, j).second == (j + start) && backtrace(w, j).first < w));
					currentM[w] = currentQ[w];
				}
				//allow this only if R has been computed, so only if fullR condition is true or fullR condition is false and has inneighbors inside the band
				if ((w == nodeStart[nodeIndex] && notInOrder[nodeIndex]) || (!(w == nodeStart[nodeIndex] && notInOrder[nodeIndex]) && neighborInsideBand))
				{
					if (currentR[w] > currentM[w])
					{
						currentM[w] = currentR[w];
						backtrace(w, j) = currentRbacktrace[w];
						assert(backtrace(w, j).second < (j + start) || (backtrace(w, j).second == (j + start) && backtrace(w, j).first < w));
					}
				}
				if (w == nodeStart[nodeIndex])
				{
					for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
					{
						auto u = nodeEnd[inNeighbors[nodeIndex][i]]-1;
						if (!band(u, start+j-1)) continue;
						//-1 because the rows in the DP matrix are one-based, eg. M[w][1] is the _first_ nucleotide of the read (sequence[0])
						if (previousM[u]+matchScore(nodeSequences[w], sequence[j + start - 1]) > currentM[w])
						{
							currentM[w] = previousM[u]+matchScore(nodeSequences[w], sequence[j + start - 1]);
							backtrace(w, j) = std::make_pair(u, j-1 + start);
							assert(backtrace(w, j).second < (j + start) || (backtrace(w, j).second == (j + start) && backtrace(w, j).first < w));
						}
					}
				}
				else
				{
					LengthType u = w-1;
					if (band(u, start+j-1))
					{
						//-1 because the rows in the DP matrix are one-based, eg. M[w][1] is the _first_ nucleotide of the read (sequence[0])
						if (previousM[u]+matchScore(nodeSequences[w], sequence[j + start - 1]) > currentM[w])
						{
							currentM[w] = previousM[u]+matchScore(nodeSequences[w], sequence[j + start - 1]);
							backtrace(w, j) = std::make_pair(u, j-1 + start);
							assert(backtrace(w, j).second < (j + start) || (backtrace(w, j).second == (j + start) && backtrace(w, j).first < w));
						}
					}
				}
				//if the previous row was not inside the band, initialize Q as the current M
				if (!band(w, start+j-1))
				{
					currentQ[w] = currentM[w];
					Qbacktrace[w] = std::make_pair(w, j + start);
				}
				//if we calculated R using the recurrence but previous columns are not inside the band, initialize it as current M
				if ((!(w == nodeStart[nodeIndex] && notInOrder[nodeIndex]) && !neighborInsideBand))
				{
					currentR[w] = currentM[w];
					currentRbacktrace[w] = std::make_pair(w, j + start);
				}
				assert(currentM[w] >= -std::numeric_limits<ScoreType>::min() + 100);
				assert(currentM[w] <= std::numeric_limits<ScoreType>::max() - 100);
				assert(backtrace(w, j).second < (j + start) || (backtrace(w, j).second == (j + start) && backtrace(w, j).first < w));
			}

			std::swap(currentM, previousM);
			std::swap(currentQ, previousQ);
			std::swap(currentR, previousR);
			std::swap(currentRbacktrace, previousRbacktrace);
			std::swap(currentProcessableColumns, previousProcessableColumns);
		}
		result.Qbacktrace = Qbacktrace;
		//use previous instead of current because the last line swapped them
		result.M = std::move(previousM);
		result.Q = std::move(previousQ);
		result.R = std::move(previousR);
		result.Rbacktrace = std::move(previousRbacktrace);
		return result;
	}

	template <bool matrixOrder>
	void expandBandRightwards(Array2D<int, matrixOrder>& matrix, LengthType w, LengthType j, int bandWidth) const
	{
		auto nodeIndex = indexToNode[w];
		auto end = nodeEnd[nodeIndex];
		while (w != end && bandWidth > 0)
		{
			matrix(w, j) = bandWidth;
			w++;
			bandWidth--;
			if (w != end && matrix(w, j) >= bandWidth) return;
		}
		if (w == end && bandWidth > 0)
		{
			for (size_t i = 0; i < outNeighbors[nodeIndex].size(); i++)
			{
				expandBandRightwards(matrix, nodeStart[outNeighbors[nodeIndex][i]], j, bandWidth);
			}
		}
	}

	template <bool matrixOrder>
	void expandBandLeftwards(Array2D<int, matrixOrder>& matrix, LengthType w, LengthType j, int bandWidth) const
	{
		auto nodeIndex = indexToNode[w];
		auto start = nodeStart[nodeIndex];
		while (w != start && bandWidth > 0)
		{
			matrix(w, j) = bandWidth;
			w--;
			bandWidth--;
			if (w != start && matrix(w, j) >= bandWidth) return;
		}
		if (w == start && bandWidth > 0)
		{
			matrix(w, j) = bandWidth;
			for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
			{
				expandBandLeftwards(matrix, nodeEnd[inNeighbors[nodeIndex][i]] - 1, j, bandWidth-1);
			}
		}
	}

	Array2D<bool, false> getBandedRows(const std::vector<MatrixPosition>& seedHits, int bandWidth, size_t sequenceLength) const
	{
		Array2D<int, false> forward {nodeSequences.size(), sequenceLength+1, 0};
		Array2D<int, false> backward {nodeSequences.size(), sequenceLength+1, 0};
		Array2D<bool, false> result {nodeSequences.size(), sequenceLength+1, false};
		for (auto pos : seedHits)
		{
			std::cerr << "seed hit: " << pos.first << ", " << pos.second << std::endl;
			forward(pos.first, pos.second) = bandWidth;
			backward(pos.first, pos.second) = bandWidth;
			expandBandRightwards(forward, pos.first, pos.second, bandWidth);
			expandBandRightwards(backward, pos.first, pos.second, bandWidth);
			expandBandLeftwards(forward, pos.first, pos.second, bandWidth);
			expandBandLeftwards(backward, pos.first, pos.second, bandWidth);
		}
		for (size_t j = 0; j < sequenceLength+1; j++)
		{
			for (size_t w = 0; w < nodeSequences.size(); w++)
			{
				if (forward(w, j) > 0) continue;
				auto nodeIndex = indexToNode[w];
				if (nodeStart[nodeIndex] == w)
				{
					for (size_t neighbori = 0; neighbori < inNeighbors[nodeIndex].size(); neighbori++)
					{
						assert(inNeighbors[nodeIndex][neighbori] < nodeEnd.size());
						auto u = nodeEnd[inNeighbors[nodeIndex][neighbori]]-1;
						assert(u < nodeSequences.size());
						if (forward(u, j-1) > 0) forward(w, j) = 1;
					}
				}
				else
				{
					if (forward(w-1, j-1) > 0) forward(w, j) = 1;
				}
			}
		}
		for (size_t j = sequenceLength; j < sequenceLength+1; j--)
		{
			for (size_t w = nodeSequences.size()-1; w < nodeSequences.size(); w--)
			{
				if (backward(w, j) > 0) continue;
				auto nodeIndex = indexToNode[w];
				if (nodeEnd[nodeIndex]-1 == w)
				{
					for (size_t neighbori = 0; neighbori < outNeighbors[nodeIndex].size(); neighbori++)
					{
						assert(outNeighbors[nodeIndex][neighbori] < nodeEnd.size());
						auto u = nodeStart[outNeighbors[nodeIndex][neighbori]];
						assert(u < nodeSequences.size());
						if (backward(u, j+1) > 0) backward(w, j) = 1;
					}
				}
				else
				{
					if (backward(w+1, j+1) > 0) backward(w, j) = 1;
				}
			}
			result(0, j) = true;
		}
		for (size_t j = 0; j < sequenceLength+1; j++)
		{
			for (size_t w = 1; w < nodeSequences.size(); w++)
			{
				result(w, j) = forward(w, j) > 0 || backward(w, j) > 0;
			}
		}
		return result;
	}

	std::pair<ScoreType, std::vector<MatrixPosition>> backtrackWithSquareRootSlices(const std::string& sequence, int bandWidth, const std::vector<MatrixPosition>& seedHits) const
	{
		auto band = getBandedRows(seedHits, bandWidth, sequence.size());
		auto distanceMatrix = getDistanceMatrixBoostJohnson();
		bool hasWrongOrders = false;
		std::vector<LengthType> nodeOrdering;
		std::vector<LengthType> nodeNotOrdering;
		nodeOrdering.reserve(nodeSequences.size());
		for (LengthType i = 1; i < nodeSequences.size(); i++)
		{
			auto nodeIndex = indexToNode[i];
			if (i == nodeStart[nodeIndex] && notInOrder[nodeIndex])
			{
				nodeOrdering.emplace_back(i);
				hasWrongOrders = true;
			}
			else{
				nodeNotOrdering.emplace_back(i);
			}
		}
		nodeOrdering.insert(nodeOrdering.end(), nodeNotOrdering.begin(), nodeNotOrdering.end());
		nodeNotOrdering.clear();
		assert(nodeOrdering.size() == nodeSequences.size() - 1);
		Array2D<MatrixPosition, true> backtraceMatrix {nodeSequences.size(), sequence.size() + 1, std::make_pair(0, 0)};
		MatrixSlice lastRow = getFirstSlice(bandWidth, backtraceMatrix);
		int sliceSize = sequence.size();
		std::vector<ScoreType> lastRowScore;
		LengthType start = 1;
		//size+1 because the rows in the DP matrix are one-based, eg. M[w][1] is the _first_ nucleotide of the read (sequence[0])
		while (start < sequence.size()+1)
		{
			LengthType end = start + sliceSize;
			if (end > sequence.size()+1) end = sequence.size();
			auto slice = getScoreAndBacktraceMatrixSlice(sequence, hasWrongOrders, nodeOrdering, distanceMatrix, lastRow, start-1, end, bandWidth, band, backtraceMatrix);
			lastRowScore = slice.M;
			lastRow = std::move(slice);
			start = end;
		}
		auto result = backtrace(lastRow.M, backtraceMatrix, band, sequence.size(), distanceMatrix, seedHits);
		return result;
	}

	template <bool backtraceOrder>
	MatrixSlice getFirstSlice(int bandWidth, Array2D<MatrixPosition, backtraceOrder>& backtrace) const
	{
		MatrixSlice result;
		result.M.resize(nodeSequences.size(), 0);
		result.R.resize(nodeSequences.size(), 0);
		result.Q.resize(nodeSequences.size(), 0);
		result.Rbacktrace.reserve(nodeSequences.size());
		result.Qbacktrace.reserve(nodeSequences.size());
		for (LengthType i = 0; i < nodeSequences.size(); i++)
		{
			backtrace(i, 0) = std::make_pair(i, 0);
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

	template <bool bandOrder>
	std::vector<std::pair<LengthType, ScoreType>> getRHelper(LengthType j, LengthType start, const std::vector<ScoreType>& previousM, const std::string& sequence, const Array2D<bool, bandOrder>& band, const std::vector<LengthType>& previousProcessableColumns) const
	{
		if (j == 0) return getRHelperZero();
		if (j == 1 && start == 0) return getRHelperOne();
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
					if (!band(u, start+j-1)) continue;
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
				if (!band(u, start+j-1)) continue;
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

	template <bool bandOrder>
	bool hasInNeighborInsideBand(LengthType w, LengthType j, LengthType start, const Array2D<bool, bandOrder>& band) const
	{
		auto nodeIndex = indexToNode[w];
		if (nodeStart[nodeIndex] == w)
		{
			for (size_t neighborI = 0; neighborI < inNeighbors[nodeIndex].size(); neighborI++)
			{
				if (band(nodeEnd[inNeighbors[nodeIndex][neighborI]] - 1, j+start)) return true;
			}
		}
		else
		{
			return band(w-1, start+j);
		}
		return false;
	}

	//compute R using the recurrence on page 3
	template <bool bandOrder>
	std::pair<ScoreType, MatrixPosition> recurrenceR(LengthType w, LengthType j, LengthType start, const std::vector<ScoreType>& currentM, const std::vector<ScoreType>& currentR, const std::vector<MatrixPosition>& currentRbacktrace, const Array2D<bool, bandOrder>& band) const
	{
		assert(band(w, start+j));
		auto nodeIndex = indexToNode[w];
		assert(nodeStart[nodeIndex] != w || !notInOrder[nodeIndex]);
		MatrixPosition pos;
		ScoreType maxValue = std::numeric_limits<ScoreType>::min() + 99;
		if (nodeStart[nodeIndex] == w)
		{
			for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
			{
				auto u = nodeEnd[inNeighbors[nodeIndex][i]]-1;
				if (!band(u, start+j)) continue;
				assert(u < w);
				if (currentM[u] - gapPenalty(1) > maxValue)
				{
					maxValue = currentM[u] - gapPenalty(1);
					pos = std::make_pair(u, j + start);
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
			if (band(u, start+j))
			{
				pos = currentRbacktrace[u];
				maxValue = currentR[u] - gapContinuePenalty;
				if (currentM[u] - gapPenalty(1) > maxValue)
				{
					pos = std::make_pair(u, j + start);
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
	std::pair<ScoreType, MatrixPosition> fullR(LengthType w, LengthType j, const std::vector<std::pair<LengthType, ScoreType>>& RHelper, const Array2D<LengthType, distanceMatrixOrder>& distanceMatrix, LengthType start) const
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
				pos = std::make_pair(v, j-1 + start);
			}
		}
		assert(maxValue >= -std::numeric_limits<ScoreType>::min() + 100);
		assert(maxValue <= std::numeric_limits<ScoreType>::max() - 100);
		return std::make_pair(maxValue, pos);
	}

	template <bool distanceMatrixOrder>
	LengthType distanceFromSeqToSeq(LengthType start, LengthType end, const Array2D<LengthType, distanceMatrixOrder>& distanceMatrix) const
	{
		auto startNode = indexToNode[start];
		auto endNode = indexToNode[end];
		if (startNode == endNode && end >= start) return end - start;
		return distanceMatrix(startNode, endNode) + nodeStart[startNode] + end - nodeStart[endNode] - start;
	}

	Array2D<LengthType, false> getDistanceMatrixBoostJohnson() const
	{
		//http://www.boost.org/doc/libs/1_40_0/libs/graph/example/johnson-eg.cpp
		auto V = inNeighbors.size();
		adjacency_list<vecS, vecS, directedS, no_property, property<edge_weight_t, int, property<edge_weight2_t, int>>> graph { inNeighbors.size() };
		for (size_t i = 0; i < inNeighbors.size(); i++)
		{
			for (size_t j = 0; j < inNeighbors[i].size(); j++)
			{
				boost::add_edge(inNeighbors[i][j], i, graph);
			}
		}
		property_map<adjacency_list<vecS, vecS, directedS, no_property, property<edge_weight_t, int, property<edge_weight2_t, int>>>, edge_weight_t>::type w = get(edge_weight, graph);
		graph_traits<adjacency_list<vecS, vecS, directedS, no_property, property<edge_weight_t, int, property<edge_weight2_t, int>>>>::edge_iterator e, e_end;
		for (boost::tie(e, e_end) = edges(graph); e != e_end; ++e)
		{
			auto startIndex = (*e).m_source;
			w[*e] = nodeEnd[startIndex] - nodeStart[startIndex];
		}
		std::vector<int> d(V, nodeSequences.size()+1);
		int** D;
		D = new int*[inNeighbors.size()];
		for (size_t i = 0; i < inNeighbors.size(); i++)
		{
			D[i] = new int[inNeighbors.size()];
			for (size_t j = 0; j < inNeighbors.size(); j++)
			{
				D[i][j] = nodeSequences.size()+1;
			}
		}
		johnson_all_pairs_shortest_paths(graph, D, distance_map(&d[0]));
		Array2D<LengthType, false> result {inNeighbors.size(), inNeighbors.size(), std::numeric_limits<LengthType>::max()};
		for (size_t i = 0; i < inNeighbors.size(); i++)
		{
			for (size_t j = 0; j < inNeighbors.size(); j++)
			{
				if (D[i][j] == std::numeric_limits<int>::max())
				{
					result(i, j) = nodeSequences.size()+1;
				}
				else
				{
					result(i, j) = D[i][j];
				}
			}
		}
		for (size_t i = 0; i < inNeighbors.size(); i++)
		{
			delete [] D[i];
		}
		delete [] D;
		//make sure that the distance to itself is not 0
		//we need to do this so distance calculation from a later point in the node to an earlier point in the node works correctly
		for (size_t i = 0; i < inNeighbors.size(); i++)
		{
			result(i, i) = nodeSequences.size()+1;
			for (size_t j = 0; j < inNeighbors.size(); j++)
			{
				if (j == i) continue;
				result(i, i) = std::min(result(i, i), result(i, j) + result(j, i));
			}
		}
		return result;
	}

	ScoreType gapPenalty(LengthType length) const
	{
		if (length == 0) return 0;
		return gapStartPenalty + gapContinuePenalty * (length - 1);
	}

	ScoreType matchScore(char graph, char sequence) const
	{
		return graph == sequence ? 1 : -4;
	}

	std::vector<bool> notInOrder;
	std::vector<LengthType> nodeStart;
	std::vector<LengthType> nodeEnd;
	std::vector<LengthType> indexToNode;
	std::map<int, LengthType> nodeLookup;
	std::vector<int> nodeIDs;
	std::vector<std::vector<LengthType>> inNeighbors;
	std::vector<std::vector<LengthType>> outNeighbors;
	std::string nodeSequences;
	ScoreType gapStartPenalty;
	ScoreType gapContinuePenalty;
	bool finalized;
};

#endif