//http://biorxiv.org/content/early/2017/04/06/124941
#include <algorithm>
#include <string>
#include <vector>
#include <cassert>
#include "vg.pb.h"

template <typename LengthType, typename ScoreType>
class GraphAligner
{
public:
	typedef std::pair<LengthType, LengthType> MatrixPosition;

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
	}
	
	void AddNode(int nodeId, std::string sequence)
	{
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
	
	void AddEdge(int node_id_from, int node_id_to)
	{
		assert(nodeLookup.count(node_id_from) > 0);
		assert(nodeLookup.count(node_id_to) > 0);
		auto from = nodeLookup[node_id_from];
		auto to = nodeLookup[node_id_to];
		assert(to >= 0);
		assert(from >= 0);
		assert(to < inNeighbors.size());
		assert(from < nodeStart.size());
		inNeighbors[to].push_back(from);
		if (from > to)
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
	}

	vg::Alignment AlignOneWay(std::string seq_id, std::string sequence, bool reverse) const
	{
		auto matrix = getScoreAndBacktraceMatrix(seq_id, sequence);
		auto trace = backtrace(matrix.first, matrix.second);
		auto result = traceToAlignment(seq_id, trace, reverse);
		return result;
	}

private:

	vg::Alignment traceToAlignment(std::string seq_id, std::pair<ScoreType, std::vector<MatrixPosition>> traceWithScore, bool reverse) const
	{
		auto& trace = traceWithScore.second;
		vg::Alignment result;
		result.set_name(seq_id);
		auto path = new vg::Path;
		result.set_allocated_path(path);
		size_t oldNode = indexToNode[trace[0].first];
		int rank = 0;
		auto vgmapping = path->add_mapping();
		auto position = new vg::Position;
		vgmapping->set_allocated_position(position);
		vgmapping->set_rank(rank);
		position->set_node_id(nodeIDs[oldNode]);
		if (reverse) position->set_is_reverse(true);
		std::cerr << nodeIDs[oldNode] << " ";
		for (size_t i = 1; i < trace.size(); i++)
		{
			if (indexToNode[trace[i].first] == oldNode) continue;
			oldNode = indexToNode[trace[i].first];
			rank++;
			vgmapping = path->add_mapping();
			position = new vg::Position;
			vgmapping->set_allocated_position(position);
			vgmapping->set_rank(rank);
			position->set_node_id(nodeIDs[oldNode]);
			if (reverse) position->set_is_reverse(true);
			std::cerr << nodeIDs[oldNode] << " ";
		}
		std::cerr << "\n";
		result.set_score(traceWithScore.first);
		return result;
	}

	std::pair<ScoreType, std::vector<MatrixPosition>> backtrace(const std::vector<std::vector<ScoreType>>& M, const std::vector<std::vector<MatrixPosition>>& backtrace) const
	{
		std::vector<MatrixPosition> trace;
		MatrixPosition currentPosition = std::make_pair(0, M[0].size()-1);
		//start at the highest value at end of read
		for (size_t i = 0; i < M.size(); i++)
		{
			MatrixPosition candidatePosition = std::make_pair(i, M[0].size()-1);
			if (M[candidatePosition.first][candidatePosition.second] > M[currentPosition.first][currentPosition.second])
			{
				currentPosition = candidatePosition;
			}
		}
		auto score = M[currentPosition.first][currentPosition.second];
		trace.push_back(currentPosition);
		std::cerr << currentPosition.first << ", " << currentPosition.second << " ";
		std::cerr << M[currentPosition.first][currentPosition.second] << std::endl;
		while (currentPosition.second > 0)
		{
			auto newPos = backtrace[currentPosition.first][currentPosition.second];
			assert(newPos.second < currentPosition.second || (newPos.second == currentPosition.second && newPos.first < currentPosition.first));
			currentPosition = newPos;
			trace.push_back(currentPosition);
			std::cerr << currentPosition.first << ", " << currentPosition.second << " ";
			std::cerr << M[currentPosition.first][currentPosition.second] << std::endl;
		}
		std::reverse(trace.begin(), trace.end());
		return std::make_pair(score, trace);
	}

	std::pair<std::vector<std::vector<ScoreType>>, std::vector<std::vector<MatrixPosition>>> getScoreAndBacktraceMatrix(std::string seq_id, std::string sequence) const
	{
		std::vector<std::vector<ScoreType>> M;
		std::vector<std::vector<ScoreType>> Q;
		std::vector<std::vector<ScoreType>> R;
		std::vector<std::vector<MatrixPosition>> backtrace;
		std::vector<std::vector<MatrixPosition>> Rbacktrace;
		std::vector<MatrixPosition> Qbacktrace;
		M.resize(nodeSequences.size());
		Q.resize(nodeSequences.size());
		R.resize(nodeSequences.size());
		Rbacktrace.resize(nodeSequences.size());
		backtrace.resize(nodeSequences.size());
		for (LengthType i = 0; i < nodeSequences.size(); i++)
		{
			M[i].resize(sequence.size());
			Q[i].resize(sequence.size());
			R[i].resize(sequence.size());
			Rbacktrace[i].resize(sequence.size());
			backtrace[i].resize(sequence.size());
			Qbacktrace.emplace_back(i, 0);
			Rbacktrace[i][0] = std::make_pair(i, 0);
		}
		for (LengthType i = 1; i < sequence.size(); i++)
		{
			R[0][i] = std::numeric_limits<ScoreType>::min() + gapContinuePenalty + 100;
		}
		std::vector<std::vector<LengthType>> distanceMatrix = getDistanceMatrix();
		for (LengthType w = 0; w < nodeSequences.size(); w++)
		{
			M[w][0] = 0;
			Q[w][0] = 0;
			R[w][0] = 0;
		}
		std::vector<LengthType> nodeOrdering;
		for (LengthType i = 1; i < nodeSequences.size(); i++)
		{
			auto nodeIndex = indexToNode[i];
			if (i == nodeStart[nodeIndex] && notInOrder[nodeIndex])
			{
				nodeOrdering.push_back(i);
			}
		}
		for (LengthType i = 1; i < nodeSequences.size(); i++)
		{
			auto nodeIndex = indexToNode[i];
			if (!(i == nodeStart[nodeIndex] && notInOrder[nodeIndex]))
			{
				nodeOrdering.push_back(i);
			}
		}
		assert(nodeOrdering.size() == nodeSequences.size() - 1);
		for (LengthType j = 1; j < sequence.size(); j++)
		{
			for (LengthType w : nodeOrdering)
			{
				auto nodeIndex = indexToNode[w];
				Q[w][j] = Q[w][j-1] - gapContinuePenalty;
				if (M[w][j-1] - gapPenalty(1) > Q[w][j])
				{
					Q[w][j] = M[w][j-1] - gapPenalty(1);
					Qbacktrace[w] = std::make_pair(w, j-1);
				}
				if (w == nodeStart[nodeIndex] && notInOrder[nodeIndex])
				{
					auto rr = fullR(w, j, M, sequence, distanceMatrix);
					R[w][j] = rr.first;
					Rbacktrace[w][j] = rr.second;
				}
				else
				{
					auto rr = recurrenceR(w, j, M, R, Rbacktrace);
					R[w][j] = rr.first;
					Rbacktrace[w][j] = rr.second;
				}
				if (w == nodeStart[nodeIndex])
				{
					backtrace[w][j] = Qbacktrace[w];
					M[w][j] = Q[w][j];
					if (R[w][j] > M[w][j])
					{
						M[w][j] = R[w][j];
						backtrace[w][j] = Rbacktrace[w][j];
					}
					for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
					{
						auto u = nodeEnd[inNeighbors[nodeIndex][i]]-1;
						if (M[u][j-1]+matchScore(nodeSequences[w], sequence[j]) > M[w][j])
						{
							M[w][j] = M[u][j-1]+matchScore(nodeSequences[w], sequence[j]);
							backtrace[w][j] = std::make_pair(u, j-1);
						}
					}
				}
				else
				{
					LengthType u = w-1;
					backtrace[w][j] = Qbacktrace[w];
					M[w][j] = Q[w][j];
					if (R[w][j] > M[w][j])
					{
						M[w][j] = R[w][j];
						backtrace[w][j] = Rbacktrace[w][j];
					}
					if (M[u][j-1]+matchScore(nodeSequences[w], sequence[j]) > M[w][j])
					{
						M[w][j] = M[u][j-1]+matchScore(nodeSequences[w], sequence[j]);
						backtrace[w][j] = std::make_pair(u, j-1);
					}
				}
				assert(M[w][j] >= -10000);
				assert(M[w][j] <= 10000);
			}
		}
		return std::make_pair(M, backtrace);
	}

	//compute R using the recurrence on page 3
	std::pair<ScoreType, MatrixPosition> recurrenceR(LengthType w, LengthType j, const std::vector<std::vector<ScoreType>>& M, const std::vector<std::vector<ScoreType>>& R, const std::vector<std::vector<MatrixPosition>>& Rbacktrace) const
	{
		auto nodeIndex = indexToNode[w];
		MatrixPosition pos;
		ScoreType maxValue = std::numeric_limits<ScoreType>::min();
		if (nodeStart[nodeIndex] == w)
		{
			for (size_t i = 0; i < inNeighbors[nodeIndex].size(); i++)
			{
				auto neighborEnd = nodeEnd[inNeighbors[nodeIndex][i]]-1;
				if (M[neighborEnd][j] - gapPenalty(1) > maxValue)
				{
					maxValue = M[neighborEnd][j] - gapPenalty(1);
					pos = std::make_pair(neighborEnd, j);
				}
				if (R[neighborEnd][j] - gapContinuePenalty > maxValue)
				{
					maxValue = R[neighborEnd][j] - gapContinuePenalty;
					pos = Rbacktrace[neighborEnd][j];
				}
			}
		}
		else
		{
			pos = Rbacktrace[w-1][j];
			maxValue = R[w-1][j] - gapContinuePenalty;
			if (M[w-1][j] - gapPenalty(1) > maxValue)
			{
				pos = std::make_pair(w-1, j);
				maxValue = M[w-1][j] - gapPenalty(1);
			}
		}
		assert(maxValue >= -10000);
		assert(maxValue <= 10000);
		return std::make_pair(maxValue, pos);
	}

	//compute R using the slow, full definition on page 3
	std::pair<ScoreType, MatrixPosition> fullR(LengthType w, LengthType j, const std::vector<std::vector<ScoreType>>& M, const std::string& sequence, const std::vector<std::vector<LengthType>>& distanceMatrix) const
	{
		MatrixPosition pos;
		ScoreType maxValue = std::numeric_limits<ScoreType>::min();
		for (LengthType v = 0; v < nodeSequences.size(); v++)
		{
			if (v == w) continue;
			auto otherNode = indexToNode[v];
			if (v == nodeStart[otherNode])
			{
				for (size_t neighbori = 0; neighbori < inNeighbors[otherNode].size(); neighbori++)
				{
					LengthType u = inNeighbors[otherNode][neighbori];
					auto scoreHere = M[u][j-1] + matchScore(nodeSequences[v], sequence[j]) - gapPenalty(distanceFromSeqToSeq(v, w, distanceMatrix));
					if (scoreHere > maxValue)
					{
						pos = std::make_pair(v, j-1);
						maxValue = scoreHere;
					}
				}
			}
			else
			{
				LengthType u = v-1;
				auto scoreHere = M[u][j-1]+matchScore(nodeSequences[v], sequence[j]) - gapPenalty(distanceFromSeqToSeq(v, w, distanceMatrix));
				if (scoreHere > maxValue)
				{
					pos = std::make_pair(u, j-1);
					maxValue = scoreHere;
				}
			}
		}
		assert(maxValue >= -10000);
		assert(maxValue <= 10000);
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
};
