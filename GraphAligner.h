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
		ScoreType localMaximumScore;
		MatrixPosition localMaximum;
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

	vg::Alignment AlignOneWay(const std::string& seq_id, const std::string& sequence, bool reverse) const
	{
		assert(finalized);
		auto trace = backtrackWithSquareRootSlices(sequence, false);
		auto result = traceToAlignment(seq_id, trace, reverse);
		return result;
	}

	std::tuple<ScoreType, LengthType, LengthType> GetLocalAlignmentSequencePosition(const std::string& sequence)
	{
		assert(finalized);
		auto trace = backtrackWithSquareRootSlices(sequence, true);
		LengthType minpos = trace.second[0].second;
		LengthType maxpos = trace.second[0].second;
		for (size_t i = 1; i < trace.second.size(); i++)
		{
			minpos = std::min(minpos, trace.second[i].second);
			maxpos = std::max(maxpos, trace.second[i].second);
		}
		return std::make_tuple(trace.first, minpos, maxpos);
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

	std::pair<ScoreType, std::vector<MatrixPosition>> backtraceLocal(ScoreType score, MatrixPosition start, const std::vector<std::vector<MatrixPosition>>& backtrace) const
	{
		std::vector<MatrixPosition> positions;
		while (start.first > 0 && start.second > 0)
		{
			positions.push_back(start);
			start = backtrace[start.first][start.second];
		}
		return std::make_pair(score, positions);
	}

	std::pair<ScoreType, std::vector<MatrixPosition>> backtrace(const std::vector<ScoreType>& Mslice, const std::vector<std::vector<MatrixPosition>>& backtrace) const
	{
		assert(backtrace.size() == nodeSequences.size());
		assert(backtrace[0].size() > 0);
		std::vector<MatrixPosition> trace;
		MatrixPosition currentPosition = std::make_pair(0, backtrace[0].size()-1);
		//start at the highest value at end of read
		for (size_t i = 0; i < Mslice.size(); i++)
		{
			MatrixPosition candidatePosition = std::make_pair(i, backtrace[0].size()-1);
			if (Mslice[candidatePosition.first] > Mslice[currentPosition.first])
			{
				currentPosition = candidatePosition;
			}
		}
		auto score = Mslice[currentPosition.first];
		trace.push_back(currentPosition);
		while (currentPosition.second > 0)
		{
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
		std::reverse(trace.begin(), trace.end());
		return std::make_pair(score, trace);
	}

	MatrixSlice getScoreAndBacktraceMatrixSlice(const std::string& sequence, bool hasWrongOrders, const std::vector<LengthType>& nodeOrdering, const std::vector<std::vector<LengthType>>& distanceMatrix, const MatrixSlice& previous, LengthType start, LengthType end, bool local) const
	{
		ScoreType localMaximumScore = std::numeric_limits<ScoreType>::min();
		MatrixPosition localMaximum = std::make_pair(0, 0);
		std::vector<ScoreType> M1;
		std::vector<ScoreType> M2;
		std::vector<ScoreType> Q1;
		std::vector<ScoreType> Q2;
		std::vector<ScoreType> R1;
		std::vector<ScoreType> R2;
		std::vector<MatrixPosition> Rbacktrace1;
		std::vector<MatrixPosition> Rbacktrace2;
		assert(previous.M.size() == nodeSequences.size());
		assert(previous.R.size() == nodeSequences.size());
		assert(previous.Q.size() == nodeSequences.size());
		assert(previous.Rbacktrace.size() == nodeSequences.size());
		assert(previous.Qbacktrace.size() == nodeSequences.size());
		assert(previous.backtrace.size() == nodeSequences.size());
		MatrixSlice result;
		std::vector<std::vector<MatrixPosition>> backtrace;
		std::vector<MatrixPosition> Qbacktrace;
		M1.reserve(nodeSequences.size());
		Q1.reserve(nodeSequences.size());
		R1.reserve(nodeSequences.size());
		Rbacktrace1.reserve(nodeSequences.size());
		M2.reserve(nodeSequences.size());
		Q2.reserve(nodeSequences.size());
		R2.reserve(nodeSequences.size());
		Rbacktrace2.reserve(nodeSequences.size());
		Qbacktrace.reserve(nodeSequences.size());
		backtrace.resize(nodeSequences.size());
		std::vector<ScoreType>& currentM = M1;
		std::vector<ScoreType>& previousM = M2;
		std::vector<ScoreType>& currentQ = Q1;
		std::vector<ScoreType>& previousQ = Q2;
		std::vector<ScoreType>& currentR = R1;
		std::vector<ScoreType>& previousR = R2;
		std::vector<MatrixPosition>& currentRbacktrace = Rbacktrace1;
		std::vector<MatrixPosition>& previousRbacktrace = Rbacktrace2;
		for (LengthType w = 0; w < nodeSequences.size(); w++)
		{
			previousM[w] = previous.M[w];
			previousQ[w] = previous.Q[w];
			previousR[w] = previous.R[w];
			Qbacktrace[w] = previous.Qbacktrace[w];
			previousRbacktrace[w] = previous.Rbacktrace[w];
			backtrace[w].resize(end-start);
			backtrace[w][0] = previous.backtrace[w].back();
		}
		currentR[0] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
		previousR[0] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
		currentM[0] = -gapPenalty(start + 1);
		previousM[0] = -gapPenalty(start);
		for (LengthType j = 1; j < end - start; j++)
		{
			currentM[0] = -gapPenalty(start + j);
			currentR[0] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
			if (local) currentM[0] = 0;
			std::vector<std::pair<LengthType, ScoreType>> Rhelper;
			if (hasWrongOrders) Rhelper = getRHelper(j, previousM, sequence);

			for (LengthType w : nodeOrdering)
			{
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
					auto rr = recurrenceR(w, j, currentM, currentR, currentRbacktrace, start);
					currentR[w] = rr.first;
					currentRbacktrace[w] = rr.second;
					assert(currentRbacktrace[w].second < (j + start) || (currentRbacktrace[w].second == (j + start) && backtrace[w][j].first < w));
				}
				backtrace[w][j] = Qbacktrace[w];
				assert(backtrace[w][j].second < (j + start) || (backtrace[w][j].second == (j + start) && backtrace[w][j].first < w));
				currentM[w] = currentQ[w];
				if (currentR[w] > currentM[w])
				{
					currentM[w] = currentR[w];
					backtrace[w][j] = currentRbacktrace[w];
					assert(backtrace[w][j].second < (j + start) || (backtrace[w][j].second == (j + start) && backtrace[w][j].first < w));
				}
				if (w == nodeStart[nodeIndex])
				{
					for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
					{
						auto u = nodeEnd[inNeighbors[nodeIndex][i]]-1;
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
					//-1 because the rows in the DP matrix are one-based, eg. M[w][1] is the _first_ nucleotide of the read (sequence[0])
					if (previousM[u]+matchScore(nodeSequences[w], sequence[j + start - 1]) > currentM[w])
					{
						currentM[w] = previousM[u]+matchScore(nodeSequences[w], sequence[j + start - 1]);
						backtrace[w][j] = std::make_pair(u, j-1 + start);
						assert(backtrace[w][j].second < (j + start) || (backtrace[w][j].second == (j + start) && backtrace[w][j].first < w));
					}
				}
				if (local)
				{
					if (currentM[w] < 0)
					{
						currentM[w] = 0;
						backtrace[w][j] = std::make_pair(0, 0);
					}
				}
				if (currentM[w] > localMaximumScore)
				{
					localMaximumScore = currentM[w];
					localMaximum = std::make_pair(w, j + start);
				}
				assert(currentM[w] >= -std::numeric_limits<ScoreType>::min() + 100);
				assert(currentM[w] <= std::numeric_limits<ScoreType>::max() - 100);
				assert(backtrace[w][j].second < (j + start) || (backtrace[w][j].second == (j + start) && backtrace[w][j].first < w));
			}

			std::swap(currentM, previousM);
			std::swap(currentQ, previousQ);
			std::swap(currentR, previousR);
			std::swap(currentRbacktrace, previousRbacktrace);
		}
		result.M.reserve(nodeSequences.size());
		result.Q.reserve(nodeSequences.size());
		result.R.reserve(nodeSequences.size());
		result.Rbacktrace.reserve(nodeSequences.size());
		result.backtrace = backtrace;
		result.Qbacktrace = Qbacktrace;
		result.localMaximumScore = localMaximumScore;
		result.localMaximum = localMaximum;
		for (LengthType w = 0; w < nodeSequences.size(); w++)
		{
			//use previous instead of current because the last line swapped them
			result.M.emplace_back(previousM[w]);
			result.Q.emplace_back(previousQ[w]);
			result.R.emplace_back(previousR[w]);
			result.Rbacktrace.emplace_back(previousRbacktrace[w]);
		}
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

	std::pair<ScoreType, std::vector<MatrixPosition>> backtrackWithSquareRootSlices(const std::string& sequence, bool local) const
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
		MatrixSlice lastRow = getFirstSlice(sequence, hasWrongOrders, nodeOrdering);
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
		MatrixPosition localMaximum = std::make_pair(0, 0);
		ScoreType localMaximumScore = std::numeric_limits<ScoreType>::min();
		//size+1 because the rows in the DP matrix are one-based, eg. M[w][1] is the _first_ nucleotide of the read (sequence[0])
		while (start < sequence.size()+1)
		{
			LengthType end = start + sliceSize;
			if (end > sequence.size()+1) end = sequence.size();
			auto slice = getScoreAndBacktraceMatrixSlice(sequence, hasWrongOrders, nodeOrdering, distanceMatrix, lastRow, start-1, end, local);
			addBacktraceMatrix(backtraceMatrix, slice.backtrace);
			lastRowScore = slice.M;
			lastRow = std::move(slice);
			start = end;
			if (lastRow.localMaximumScore > localMaximumScore)
			{
				localMaximumScore = lastRow.localMaximumScore;
				localMaximum = lastRow.localMaximum;
			}
		}
		if (local)
		{
			assert(localMaximumScore > 0);
			assert(localMaximum.first != 0 || localMaximum.second != 0);
			auto localresult = backtraceLocal(localMaximumScore, localMaximum, backtraceMatrix);
			return localresult;
		}
		auto result = backtrace(lastRow.M, backtraceMatrix);
		return result;
	}

	MatrixSlice getFirstSlice(const std::string& sequence, bool hasWrongOrders, const std::vector<LengthType>& nodeOrdering) const
	{
		MatrixSlice result;
		result.M.resize(nodeSequences.size(), 0);
		result.R.resize(nodeSequences.size(), 0);
		result.Q.resize(nodeSequences.size(), 0);
		result.Rbacktrace.reserve(nodeSequences.size());
		result.Qbacktrace.reserve(nodeSequences.size());
		result.backtrace.resize(nodeSequences.size());
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

	std::vector<std::pair<LengthType, ScoreType>> getRHelper(LengthType j, const std::vector<ScoreType>& previousM, const std::string& sequence) const
	{
		if (j == 0) return getRHelperZero();
		std::vector<std::pair<LengthType, ScoreType>> result;
		for (LengthType v = 1; v < nodeSequences.size(); v++)
		{
			ScoreType maxValue = std::numeric_limits<ScoreType>::min();
			LengthType resultv = -1;
			auto otherNode = indexToNode[v];
			if (v == nodeStart[otherNode])
			{
				for (size_t neighbori = 0; neighbori < inNeighbors[otherNode].size(); neighbori++)
				{
					LengthType u = nodeEnd[inNeighbors[otherNode][neighbori]]-1;
					//j-1 because the rows in the DP matrix are one-based, eg. M[w][1] is the _first_ nucleotide of the read (sequence[0])
					auto scoreHere = previousM[u] + matchScore(nodeSequences[v], sequence[j-1]);
					if (scoreHere > maxValue)
					{
						resultv = v;
						maxValue = scoreHere;
					}
				}
			}
			else
			{
				LengthType u = v-1;
				auto scoreHere = previousM[u]+matchScore(nodeSequences[v], sequence[j-1]);
				if (scoreHere > maxValue)
				{
					resultv = v;
					maxValue = scoreHere;
				}
			}
			assert(maxValue <= std::numeric_limits<ScoreType>::max() - 100);
			assert(maxValue >= -std::numeric_limits<ScoreType>::min() + 100);
			assert(resultv != -1);
			result.emplace_back(resultv, maxValue);
		}
		return result;
	}

	//compute R using the recurrence on page 3
	std::pair<ScoreType, MatrixPosition> recurrenceR(LengthType w, LengthType j, const std::vector<ScoreType>& currentM, const std::vector<ScoreType>& currentR, const std::vector<MatrixPosition>& currentRbacktrace, LengthType start) const
	{
		auto nodeIndex = indexToNode[w];
		assert(nodeStart[nodeIndex] != w || !notInOrder[nodeIndex]);
		MatrixPosition pos;
		ScoreType maxValue = std::numeric_limits<ScoreType>::min();
		if (nodeStart[nodeIndex] == w)
		{
			for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
			{
				auto u = nodeEnd[inNeighbors[nodeIndex][i]]-1;
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
			pos = currentRbacktrace[u];
			maxValue = currentR[u] - gapContinuePenalty;
			if (currentM[u] - gapPenalty(1) > maxValue)
			{
				pos = std::make_pair(w-1, j + start);
				maxValue = currentM[u] - gapPenalty(1);
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
		return distanceMatrix[startNode][endNode] + end + nodeEnd[startNode] - nodeStart[endNode] - start;
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