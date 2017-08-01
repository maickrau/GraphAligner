#include <iostream>
#include <limits>
#include <algorithm>
#include "AlignmentGraph.h"


AlignmentGraph::AlignmentGraph() :
nodeStart(),
indexToNode(),
nodeLookup(),
nodeIDs(),
inNeighbors(),
nodeSequences(),
finalized(false),
firstInOrder(0)
{
	//add the start dummy node as the first node
	dummyNodeStart = nodeSequences.size();
	nodeIDs.push_back(0);
	nodeStart.push_back(nodeSequences.size());
	inNeighbors.emplace_back();
	outNeighbors.emplace_back();
	reverse.push_back(false);
	nodeSequences.push_back('-');
	indexToNode.resize(nodeSequences.size(), nodeStart.size()-1);
	nodeEnd.emplace_back(nodeSequences.size());
	notInOrder.push_back(false);
}

void AlignmentGraph::AddNode(int nodeId, std::string sequence, bool reverseNode)
{
	assert(!finalized);
		//subgraph extraction might produce different subgraphs with common nodes
		//don't add duplicate nodes
	if (nodeLookup.count(nodeId) != 0) return;

	assert(std::numeric_limits<size_t>::max() - sequence.size() > nodeSequences.size());
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

void AlignmentGraph::AddEdgeNodeId(int node_id_from, int node_id_to)
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
	//skip duplicate edges
	if (std::find(inNeighbors[to].begin(), inNeighbors[to].end(), from) != inNeighbors[to].end())
	{
		return;
	}

	inNeighbors[to].push_back(from);
	outNeighbors[from].push_back(to);
	if (from >= to)
	{
		notInOrder[to] = true;
	}
}

void AlignmentGraph::Finalize(int wordSize)
{
	//add the end dummy node as the last node
	dummyNodeEnd = nodeSequences.size();
	nodeIDs.push_back(0);
	nodeStart.push_back(nodeSequences.size());
	reverse.push_back(false);
	inNeighbors.emplace_back();
	outNeighbors.emplace_back();
	nodeSequences.push_back('-');
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
	int specialNodes = 0;
	for (size_t i = 0; i < inNeighbors.size(); i++)
	{
		if (inNeighbors[i].size() >= 2) specialNodes++;
	}
	firstInOrder = 0;
	for (size_t i = 1; i < notInOrder.size(); i++)
	{
		if (notInOrder[i]) firstInOrder = i+1;
		//all not-in-order nodes have to be at the start
		assert(i == 1 || !notInOrder[i] || notInOrder[i-1]);
	}
	cycleCuttingNodes.resize(firstInOrder);
	cycleCuttingNodePredecessors.resize(firstInOrder);
	for (size_t i = 1; i < firstInOrder; i++)
	{
		calculateCycleCutters(i, wordSize);
	}
	std::cerr << nodeStart.size() << " nodes" << std::endl;
	std::cerr << nodeSequences.size() << "bp" << std::endl;
	std::cerr << specialNodes << " nodes with in-degree >= 2" << std::endl;
	if (firstInOrder != 0)
	{
		std::cerr << (firstInOrder - 1) << " nodes out of order" << std::endl;
		std::cerr << "cycle cuts:" << std::endl;
		for (size_t i = 1; i < cycleCuttingNodes.size(); i++)
		{
			int cuttersbp = 0;
			for (size_t j = 0; j < cycleCuttingNodes[i].size(); j++)
			{
				cuttersbp += nodeEnd[cycleCuttingNodes[i][j]] - nodeStart[cycleCuttingNodes[i][j]];
			}
			std::cerr << i << ": id " << nodeIDs[i] << ", cutting nodes " << cycleCuttingNodes[i].size() << ", " << cuttersbp << "bp" << std::endl;
			if (cuttersbp > nodeSequences.size() * 0.1)
			{
				std::cerr << "WARNING: large cut!" << std::endl;
			}
		}
	}
	else
	{
		std::cerr << "0 nodes out of order" << std::endl;
	}
}

size_t AlignmentGraph::SizeInBp() const
{
	return nodeSequences.size();
}

std::vector<AlignmentGraph::MatrixPosition> AlignmentGraph::GetSeedHitPositionsInMatrix(const std::string& sequence, const std::vector<SeedHit>& seedHits) const
{
	std::vector<AlignmentGraph::MatrixPosition> result;
	for (size_t i = 0; i < seedHits.size(); i++)
	{
		assert(nodeLookup.count(seedHits[i].nodeId) > 0);
		result.emplace_back(nodeStart[nodeLookup.at(seedHits[i].nodeId)] + seedHits[i].nodePos, seedHits[i].sequencePosition);
	}
	return result;
}

void AlignmentGraph::calculateCycleCuttersRec(size_t cycleStart, size_t node, int wordSize, int lengthLeft, std::vector<size_t>& cutters, std::vector<std::vector<size_t>>& predecessors)
{
	lengthLeft -= nodeEnd[node] - nodeStart[node];
	if (lengthLeft > 0)
	{
		for (size_t i = 0; i < inNeighbors[node].size(); i++)
		{
			calculateCycleCuttersRec(cycleStart, inNeighbors[node][i], wordSize, lengthLeft, cutters, predecessors);
		}
	}
	predecessors.emplace_back();
	if (lengthLeft > 0)
	{
		for (size_t i = 0; i < inNeighbors[node].size(); i++)
		{
			predecessors.back().push_back(inNeighbors[node][i]);
		}
	}
	cutters.push_back(node);
}

void AlignmentGraph::calculateCycleCutters(size_t cycleStart, int wordSize)
{
	calculateCycleCuttersRec(cycleStart, cycleStart, wordSize, 2*wordSize, cycleCuttingNodes[cycleStart], cycleCuttingNodePredecessors[cycleStart]);
	assert(cycleCuttingNodes[cycleStart].size() > 0);
	assert(cycleCuttingNodes[cycleStart].back() == cycleStart);
	assert(cycleCuttingNodes[cycleStart].size() == cycleCuttingNodePredecessors[cycleStart].size());
}
