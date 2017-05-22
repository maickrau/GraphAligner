#ifndef GraphAligner_H
#define GraphAligner_H

//http://biorxiv.org/content/early/2017/04/06/124941
#include <algorithm>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include "vg.pb.h"

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
		std::vector<std::vector<MatrixPosition>> backtrace;
		std::vector<bool> insideBand;
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
		nodeSequences.insert(nodeSequences.end(), sequence.begin(), sequence.end());
		indexToNode.resize(nodeSequences.size(), nodeStart.size()-1);
		nodeEnd.emplace_back(nodeSequences.size());
		notInOrder.push_back(false);
		assert(nodeIDs.size() == nodeStart.size());
		assert(nodeStart.size() == inNeighbors.size());
		assert(inNeighbors.size() == nodeEnd.size());
		assert(nodeEnd.size() == notInOrder.size());
		assert(nodeSequences.size() == indexToNode.size());
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

	vg::Alignment AlignOneWay(const std::string& seq_id, const std::string& sequence, bool reverse, int bandWidth) const
	{
		assert(finalized);
		auto trace = backtrackWithSquareRootSlices(sequence, bandWidth);
		auto result = traceToAlignment(seq_id, trace, reverse);
		return result;
	}

	size_t SizeInBp()
	{
		return nodeSequences.size();
	}

private:

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

	std::pair<ScoreType, std::vector<MatrixPosition>> backtrace(const std::vector<ScoreType>& Mslice, const std::vector<std::vector<MatrixPosition>>& backtrace, const std::vector<bool>& insideBand, int sequenceLength, const std::vector<std::vector<LengthType>>& distanceMatrix) const
	{
		auto bandLocations = getBandLocations(sequenceLength);
		assert(backtrace.size() == nodeSequences.size());
		assert(backtrace[0].size() > 0);
		std::vector<MatrixPosition> trace;
		bool foundStart = false;
		MatrixPosition currentPosition = std::make_pair(0, backtrace[0].size()-1);
		//start at the highest value at end of read
		for (size_t i = 0; i < Mslice.size(); i++)
		{
			if (!insideBand[i]) continue;
			MatrixPosition candidatePosition = std::make_pair(i, backtrace[0].size()-1);
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
		assert(insideBand[currentPosition.first]);
		auto score = Mslice[currentPosition.first];
		trace.push_back(currentPosition);
		LengthType maxMinDistance = 0;
		while (currentPosition.second > 0)
		{
			LengthType minDistance = nodeSequences.size();
			for (size_t i = 0; i < bandLocations[currentPosition.second].size(); i++)
			{
				minDistance = std::min(minDistance, distanceFromSeqToSeq(currentPosition.first, bandLocations[currentPosition.second][i], distanceMatrix));
				minDistance = std::min(minDistance, distanceFromSeqToSeq(bandLocations[currentPosition.second][i], currentPosition.first, distanceMatrix));
			}
			maxMinDistance = std::max(maxMinDistance, minDistance);
			std::cerr << currentPosition.first << ", " << currentPosition.second << std::endl;
			std::cerr << "row " << currentPosition.second << " distance from band: " << minDistance << std::endl;
			assert(currentPosition.second >= 0);
			assert(currentPosition.second < backtrace[0].size());
			assert(currentPosition.first >= 0);
			assert(currentPosition.first < nodeSequences.size());
			auto newPos = backtrace[currentPosition.first][currentPosition.second];
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

	std::vector<std::vector<LengthType>> getBandLocations(int sequenceLength) const
	{
		std::vector<std::vector<LengthType>> result;
		result.resize(sequenceLength+1);
		result[0].emplace_back(0);
		for (LengthType j = 1; j < sequenceLength+1; j++)
		{
			for (LengthType w = 1; w < nodeSequences.size(); w++)
			{
				auto nodeIndex = indexToNode[w];
				if (nodeStart[nodeIndex] == w)
				{
					for (size_t neighborI = 0; neighborI < inNeighbors[nodeIndex].size(); neighborI++)
					{
						if (std::find(result[j-1].begin(), result[j-1].end(), nodeEnd[inNeighbors[nodeIndex][neighborI]] - 1) != result[j-1].end())
						{
							result[j].emplace_back(w);
							break;
						}
					}
				}
				else
				{
					if (std::find(result[j-1].begin(), result[j-1].end(), w-1) != result[j-1].end())
					{
						result[j].emplace_back(w);
					}
				}
			}
		}
		return result;
	}

	MatrixSlice getScoreAndBacktraceMatrixSlice(const std::string& sequence, bool hasWrongOrders, const std::vector<LengthType>& nodeOrdering, const std::vector<std::vector<LengthType>>& distanceMatrix, MatrixSlice& previous, LengthType start, LengthType end, int bandWidth) const
	{
		std::vector<ScoreType> M1;
		std::vector<ScoreType> M2;
		std::vector<ScoreType> Q1;
		std::vector<ScoreType> Q2;
		std::vector<ScoreType> R1;
		std::vector<ScoreType> R2;
		std::vector<MatrixPosition> Rbacktrace1;
		std::vector<MatrixPosition> Rbacktrace2;
		std::vector<bool> insideBand1;
		std::vector<bool> insideBand2;
		assert(previous.M.size() == nodeSequences.size());
		assert(previous.R.size() == nodeSequences.size());
		assert(previous.Q.size() == nodeSequences.size());
		assert(previous.Rbacktrace.size() == nodeSequences.size());
		assert(previous.Qbacktrace.size() == nodeSequences.size());
		assert(previous.backtrace.size() == nodeSequences.size());
		MatrixSlice result;
		std::vector<std::vector<MatrixPosition>> backtrace;
		std::vector<MatrixPosition> Qbacktrace;
		insideBand1.resize(nodeSequences.size(), false);
		insideBand2.resize(nodeSequences.size(), false);
		M1.resize(nodeSequences.size());
		Q1.resize(nodeSequences.size());
		R1.resize(nodeSequences.size());
		Rbacktrace1.resize(nodeSequences.size());
		backtrace.resize(nodeSequences.size());
		std::vector<ScoreType>& currentM = M1;
		std::vector<ScoreType>& previousM = M2;
		std::vector<ScoreType>& currentQ = Q1;
		std::vector<ScoreType>& previousQ = Q2;
		std::vector<ScoreType>& currentR = R1;
		std::vector<ScoreType>& previousR = R2;
		std::vector<MatrixPosition>& currentRbacktrace = Rbacktrace1;
		std::vector<MatrixPosition>& previousRbacktrace = Rbacktrace2;
		std::vector<bool>& currentInsideBand = insideBand1;
		std::vector<bool>& previousInsideBand = insideBand2;
		previousInsideBand = std::move(previous.insideBand);
		previousM = std::move(previous.M);
		previousQ = std::move(previous.Q);
		previousR = std::move(previous.R);
		Qbacktrace = std::move(previous.Qbacktrace);
		previousRbacktrace = std::move(previous.Rbacktrace);
		for (LengthType w = 0; w < nodeSequences.size(); w++)
		{
			backtrace[w].resize(end-start);
			backtrace[w][0] = previous.backtrace[w].back();
		}
		currentR[0] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
		previousR[0] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
		currentM[0] = -gapPenalty(start + 1);
		previousM[0] = -gapPenalty(start);
		std::vector<LengthType> processableColumns;
		processableColumns.reserve(nodeSequences.size());
		for (LengthType j = 1; j < end - start; j++)
		{
			if (j + start < bandWidth) 
			{
				currentInsideBand[0] = true;
			}
			else
			{
				currentInsideBand[0] = false;
			}
			int insideBand = 0;
			for (LengthType i = 1; i < nodeSequences.size(); i++)
			{
				auto nodeIndex = indexToNode[i];
				if (nodeStart[nodeIndex] == i)
				{
					currentInsideBand[i] = false;
					for (size_t neighborI = 0; neighborI < inNeighbors[nodeIndex].size(); neighborI++)
					{
						if (previousInsideBand[nodeEnd[inNeighbors[nodeIndex][neighborI]] - 1])
						{
							currentInsideBand[i] = true;
							break;
						}
					}
				}
				else
				{
					currentInsideBand[i] = previousInsideBand[i-1];
				}
				if (currentInsideBand[i]) insideBand++;
			}
			processableColumns.clear();
			for (LengthType w : nodeOrdering)
			{
				if (currentInsideBand[w]) processableColumns.push_back(w);
			}
			std::cerr << "inside band: " << insideBand << std::endl;
			currentM[0] = -gapPenalty(start + j);
			currentR[0] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
			std::vector<std::pair<LengthType, ScoreType>> Rhelper;
			if (hasWrongOrders) Rhelper = getRHelper(j, start, previousM, sequence, previousInsideBand);

			for (LengthType w : processableColumns)
			{
				if (!currentInsideBand[w]) continue;
				bool neighborInsideBand = hasInNeighborInsideBand(w, currentInsideBand);
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
					assert(currentRbacktrace[w].second < (j + start) || (currentRbacktrace[w].second == (j + start) && backtrace[w][j].first < w));
				}
				else
				{
					if (neighborInsideBand)
					{
						auto rr = recurrenceR(w, j, currentM, currentR, currentRbacktrace, start, currentInsideBand);
						currentR[w] = rr.first;
						currentRbacktrace[w] = rr.second;
						assert(currentRbacktrace[w].second < (j + start) || (currentRbacktrace[w].second == (j + start) && backtrace[w][j].first < w));
					}
				}
				currentM[w] = std::numeric_limits<ScoreType>::min();
				if (previousInsideBand[w])
				{
					backtrace[w][j] = Qbacktrace[w];
					assert(backtrace[w][j].second < (j + start) || (backtrace[w][j].second == (j + start) && backtrace[w][j].first < w));
					currentM[w] = currentQ[w];
				}
				//allow this only if R has been computed, so only if fullR condition is true or fullR condition is false and has inneighbors inside the band
				if ((w == nodeStart[nodeIndex] && notInOrder[nodeIndex]) || (!(w == nodeStart[nodeIndex] && notInOrder[nodeIndex]) && neighborInsideBand))
				{
					if (currentR[w] > currentM[w])
					{
						currentM[w] = currentR[w];
						backtrace[w][j] = currentRbacktrace[w];
						assert(backtrace[w][j].second < (j + start) || (backtrace[w][j].second == (j + start) && backtrace[w][j].first < w));
					}
				}
				if (w == nodeStart[nodeIndex])
				{
					for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
					{
						auto u = nodeEnd[inNeighbors[nodeIndex][i]]-1;
						if (!previousInsideBand[u]) continue;
						//-1 because the rows in the DP matrix are one-based, eg. M[w][1] is the _first_ nucleotide of the read (sequence[0])
						if (previousM[u]+matchScore(nodeSequences[w], sequence[j + start - 1]) > currentM[w])
						{
							currentM[w] = previousM[u]+matchScore(nodeSequences[w], sequence[j + start - 1]);
							backtrace[w][j] = std::make_pair(u, j-1 + start);
							assert(backtrace[w][j].second < (j + start) || (backtrace[w][j].second == (j + start) && backtrace[w][j].first < w));
						}
					}
				}
				else
				{
					LengthType u = w-1;
					if (previousInsideBand[u])
					{
						//-1 because the rows in the DP matrix are one-based, eg. M[w][1] is the _first_ nucleotide of the read (sequence[0])
						if (previousM[u]+matchScore(nodeSequences[w], sequence[j + start - 1]) > currentM[w])
						{
							currentM[w] = previousM[u]+matchScore(nodeSequences[w], sequence[j + start - 1]);
							backtrace[w][j] = std::make_pair(u, j-1 + start);
							assert(backtrace[w][j].second < (j + start) || (backtrace[w][j].second == (j + start) && backtrace[w][j].first < w));
						}
					}
				}
				//if the previous row was not inside the band, initialize Q as the current M
				if (!previousInsideBand[w])
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
				assert(backtrace[w][j].second < (j + start) || (backtrace[w][j].second == (j + start) && backtrace[w][j].first < w));
			}

			std::swap(currentM, previousM);
			std::swap(currentQ, previousQ);
			std::swap(currentR, previousR);
			std::swap(currentRbacktrace, previousRbacktrace);
			std::swap(currentInsideBand, previousInsideBand);
		}
		result.backtrace = backtrace;
		result.Qbacktrace = Qbacktrace;
		//use previous instead of current because the last line swapped them
		result.M = std::move(previousM);
		result.Q = std::move(previousQ);
		result.R = std::move(previousR);
		result.Rbacktrace = std::move(previousRbacktrace);
		result.insideBand = std::move(previousInsideBand);
		return result;
	}

	void addBacktraceMatrix(std::vector<std::vector<MatrixPosition>>& backtrace, const std::vector<std::vector<MatrixPosition>>& addThese) const
	{
		assert(backtrace.size() == nodeSequences.size());
		assert(backtrace.size() == addThese.size());
		assert(backtrace.size() > 0);
		assert(addThese[0].size() > 1);
		for (LengthType w = 0; w < backtrace.size(); w++)
		{
			backtrace[w].reserve(backtrace[w].size() + addThese[w].size());
			assert(w == 0 || addThese[w].size() == addThese[w-1].size());
			for (LengthType j = 1; j < addThese[w].size(); j++)
			{
				backtrace[w].emplace_back(addThese[w][j]);
			}
		}
	}

	std::pair<ScoreType, std::vector<MatrixPosition>> backtrackWithSquareRootSlices(const std::string& sequence, int bandWidth) const
	{
		std::vector<std::vector<LengthType>> distanceMatrix = getDistanceMatrix();
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
		MatrixSlice lastRow = getFirstSlice(distanceMatrix, bandWidth);
		int sliceSize = sequence.size();
		// int sliceSize = sqrt(sequence.size());
		std::vector<std::vector<MatrixPosition>> backtraceMatrix;
		backtraceMatrix.resize(nodeSequences.size());
		assert(lastRow.backtrace.size() == nodeSequences.size());
		for (LengthType w = 0; w < nodeSequences.size(); w++)
		{
			assert(lastRow.backtrace[w].size() == 1);
			backtraceMatrix[w].emplace_back(lastRow.backtrace[w][0]);
		}
		std::vector<ScoreType> lastRowScore;
		LengthType start = 1;
		//size+1 because the rows in the DP matrix are one-based, eg. M[w][1] is the _first_ nucleotide of the read (sequence[0])
		while (start < sequence.size()+1)
		{
			LengthType end = start + sliceSize;
			if (end > sequence.size()+1) end = sequence.size();
			auto slice = getScoreAndBacktraceMatrixSlice(sequence, hasWrongOrders, nodeOrdering, distanceMatrix, lastRow, start-1, end, bandWidth);
			addBacktraceMatrix(backtraceMatrix, slice.backtrace);
			lastRowScore = slice.M;
			lastRow = std::move(slice);
			start = end;
		}
		auto result = backtrace(lastRow.M, backtraceMatrix, lastRow.insideBand, sequence.size(), distanceMatrix);
		return result;
	}

	MatrixSlice getFirstSlice(const std::vector<std::vector<LengthType>>& distanceMatrix, int bandWidth) const
	{
		MatrixSlice result;
		result.M.resize(nodeSequences.size(), 0);
		result.R.resize(nodeSequences.size(), 0);
		result.Q.resize(nodeSequences.size(), 0);
		result.Rbacktrace.reserve(nodeSequences.size());
		result.Qbacktrace.reserve(nodeSequences.size());
		result.backtrace.resize(nodeSequences.size());
		result.insideBand.resize(nodeSequences.size(), false);
		result.insideBand[0] = true;
		for (LengthType i = 1; i < nodeSequences.size(); i++)
		{
			if (distanceFromSeqToSeq(0, i, distanceMatrix) < bandWidth) result.insideBand[i] = true;
		}
		for (LengthType i = 0; i < nodeSequences.size(); i++)
		{
			result.backtrace[i].emplace_back(i, 0);
			result.Qbacktrace.emplace_back(i, 0);
			result.Rbacktrace.emplace_back(i, 0);
		}
		result.R[0] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
		assert(result.M.size() == nodeSequences.size());
		assert(result.R.size() == nodeSequences.size());
		assert(result.Q.size() == nodeSequences.size());
		assert(result.Rbacktrace.size() == nodeSequences.size());
		assert(result.Qbacktrace.size() == nodeSequences.size());
		assert(result.backtrace.size() == nodeSequences.size());
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

	std::vector<std::pair<LengthType, ScoreType>> getRHelper(LengthType j, LengthType start, const std::vector<ScoreType>& previousM, const std::string& sequence, const std::vector<bool>& previousInsideBand) const
	{
		if (j == 0) return getRHelperZero();
		if (j == 1 && start == 0) return getRHelperOne();
		std::vector<std::pair<LengthType, ScoreType>> result;
		for (size_t i = 1; i < nodeStart.size(); i++)
		{
			assert(nodeEnd[i] > nodeStart[i]);
			assert(inNeighbors[i].size() > 0);
			ScoreType maxValue = std::numeric_limits<ScoreType>::min();
			ScoreType distancePenaltyAtMaxValue = 0;
			LengthType resultv = -1;
			LengthType v = nodeStart[i];
			bool insideBand = false;
			for (size_t neighbori = 0; neighbori < inNeighbors[i].size(); neighbori++)
			{
				LengthType u = nodeEnd[inNeighbors[i][neighbori]]-1;
				if (!previousInsideBand[u]) continue;
				if (!previousInsideBand[v]) continue;
				insideBand = true;
				assert(u < nodeSequences.size());
				//j-1 because the rows in the DP matrix are one-based, eg. M[w][1] is the _first_ nucleotide of the read (sequence[0])
				auto scoreHere = previousM[u] + matchScore(nodeSequences[v], sequence[j-1]);
				assert(previousM[u] <= std::numeric_limits<ScoreType>::max() - 100);
				assert(previousM[u] >= -std::numeric_limits<ScoreType>::min() + 100);
				if (scoreHere - (ScoreType)(nodeEnd[i] - v) * gapContinuePenalty > maxValue - distancePenaltyAtMaxValue)
				{
					resultv = v;
					maxValue = scoreHere;
					distancePenaltyAtMaxValue = (ScoreType)(nodeEnd[i] - v) * gapContinuePenalty;
				}
			}
			v++;
			while (v < nodeEnd[i])
			{
				LengthType u = v-1;
				if (!previousInsideBand[u])
				{
					v++;
					continue;
				}
				if (!previousInsideBand[v])
				{
					v++;
					continue;
				}
				insideBand = true;
				auto scoreHere = previousM[u] + matchScore(nodeSequences[v], sequence[j-1]);
				if (scoreHere - (ScoreType)(nodeEnd[i] - v) * gapContinuePenalty > maxValue - distancePenaltyAtMaxValue)
				{
					resultv = v;
					maxValue = scoreHere;
					distancePenaltyAtMaxValue = (ScoreType)(nodeEnd[i] - v) * gapContinuePenalty;
				}
				v++;
			}
			if (!insideBand)
			{
				continue;
			}
			assert(maxValue <= std::numeric_limits<ScoreType>::max() - 100);
			assert(maxValue >= -std::numeric_limits<ScoreType>::min() + 100);
			assert(resultv != -1);
			result.emplace_back(resultv, maxValue);
		}
		assert(result.size() >= 1);
		return result;
	}

	bool hasInNeighborInsideBand(LengthType w, const std::vector<bool>& currentInsideBand) const
	{
		auto nodeIndex = indexToNode[w];
		if (nodeStart[nodeIndex] == w)
		{
			for (size_t neighborI = 0; neighborI < inNeighbors[nodeIndex].size(); neighborI++)
			{
				if (currentInsideBand[nodeEnd[inNeighbors[nodeIndex][neighborI]] - 1]) return true;
			}
		}
		else
		{
			return currentInsideBand[w-1];
		}
		return false;
	}

	//compute R using the recurrence on page 3
	std::pair<ScoreType, MatrixPosition> recurrenceR(LengthType w, LengthType j, const std::vector<ScoreType>& currentM, const std::vector<ScoreType>& currentR, const std::vector<MatrixPosition>& currentRbacktrace, LengthType start, const std::vector<bool>& currentInsideBand) const
	{
		assert(currentInsideBand[w]);
		auto nodeIndex = indexToNode[w];
		assert(nodeStart[nodeIndex] != w || !notInOrder[nodeIndex]);
		MatrixPosition pos;
		ScoreType maxValue = std::numeric_limits<ScoreType>::min();
		if (nodeStart[nodeIndex] == w)
		{
			for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
			{
				auto u = nodeEnd[inNeighbors[nodeIndex][i]]-1;
				if (!currentInsideBand[u]) continue;
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
			if (currentInsideBand[u])
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
	std::pair<ScoreType, MatrixPosition> fullR(LengthType w, LengthType j, const std::vector<std::pair<LengthType, ScoreType>>& RHelper, const std::vector<std::vector<LengthType>>& distanceMatrix, LengthType start) const
	{
		assert(j > 0);
		assert(w > 0);
		auto nodeIndex = indexToNode[w];
		assert(nodeStart[nodeIndex] == w && notInOrder[nodeIndex]);
		MatrixPosition pos;
		ScoreType maxValue = std::numeric_limits<ScoreType>::min();
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

	LengthType distanceFromSeqToSeq(LengthType start, LengthType end, const std::vector<std::vector<LengthType>>& distanceMatrix) const
	{
		auto startNode = indexToNode[start];
		auto endNode = indexToNode[end];
		if (startNode == endNode && end >= start) return end - start;
		return distanceMatrix[startNode][endNode] + nodeStart[startNode] + end - nodeStart[endNode] - start;
	}

	std::vector<std::vector<LengthType>> getDistanceMatrix() const
	{
		//https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm
		std::vector<std::vector<LengthType>> result;
		result.resize(inNeighbors.size());
		for (size_t i = 0; i < inNeighbors.size(); i++)
		{
			result[i].resize(inNeighbors.size(), nodeSequences.size()+1);
		}
		for (size_t i = 0; i < inNeighbors.size(); i++)
		{
			result[i][i] = 0;
			for (size_t j = 0; j < inNeighbors[i].size(); j++)
			{
				LengthType neighbor = inNeighbors[i][j];
				LengthType neighborLength = nodeEnd[neighbor] - nodeStart[neighbor];
				result[neighbor][i] = neighborLength;
			}
		}
		//put the distance from self to self as the maximum distance.
		//this is for calculating the distance from a later point inside the node to an earlier point inside the node
		for (size_t i = 0; i < inNeighbors.size(); i++)
		{
			result[i][i] = nodeSequences.size()+1;
		}
		for (size_t k = 0; k < inNeighbors.size(); k++)
		{
			for (size_t i = 0; i < inNeighbors.size(); i++)
			{
				for (size_t j = 0; j < inNeighbors.size(); j++)
				{
					if (result[i][j] > result[i][k] + result[k][j]) 
					{
						result[i][j] = result[i][k] + result[k][j];
					}
				}
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
	std::string nodeSequences;
	ScoreType gapStartPenalty;
	ScoreType gapContinuePenalty;
	bool finalized;
};

#endif