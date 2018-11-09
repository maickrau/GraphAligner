#include <iostream>
#include <limits>
#include <algorithm>
#include <queue>
#include "AlignmentGraph.h"
#include "CommonUtils.h"
#include "ThreadReadAssertion.h"

AlignmentGraph::AlignmentGraph() :
DBGOverlap(0),
nodeLength(),
nodeLookup(),
unitigStartNode(),
nodeIDs(),
inNeighbors(),
nodeSequences(),
finalized(false)
{
}

void AlignmentGraph::ReserveNodes(size_t numNodes, size_t numSplitNodes)
{
	nodeSequences.reserve(numSplitNodes);
	nodeLookup.reserve(numNodes);
	unitigStartNode.reserve(numNodes);
	nodeIDs.reserve(numSplitNodes);
	nodeLength.reserve(numSplitNodes);
	inNeighbors.reserve(numSplitNodes);
	outNeighbors.reserve(numSplitNodes);
	reverse.reserve(numSplitNodes);
	nodeOffset.reserve(numSplitNodes);
	nodeUnitigReid.reserve(numSplitNodes);
	unitigReidInNeighbors.reserve(numNodes);
	unitigReidOutNeighbors.reserve(numNodes);
}

void AlignmentGraph::AddNode(int nodeId, const std::string& sequence, bool reverseNode)
{
	assert(!finalized);
	//subgraph extraction might produce different subgraphs with common nodes
	//don't add duplicate nodes
	if (nodeLookup.count(nodeId) != 0) return;
	originalNodeSize[nodeId] = sequence.size();
	size_t reid = unitigReidInNeighbors.size();
	unitigReidLookup[nodeId] = reid;
	unitigReidInNeighbors.emplace_back();
	unitigReidOutNeighbors.emplace_back();
	unitigReidLength.push_back(sequence.size());
	assert(DBGOverlap >= 0);
	for (size_t i = 0; i < (size_t)DBGOverlap; i += SPLIT_NODE_SIZE)
	{
		size_t size = SPLIT_NODE_SIZE;
		if (DBGOverlap - i < size) size = DBGOverlap - i;
		if (size == 0) continue;
		AddNode(nodeId, i, sequence.substr(i, size), reverseNode, reid);
		if (i > 0)
		{
			assert(outNeighbors.size() >= 2);
			assert(outNeighbors.size() == inNeighbors.size());
			assert(nodeIDs.size() == outNeighbors.size());
			assert(nodeOffset.size() == outNeighbors.size());
			assert(nodeIDs[outNeighbors.size()-2] == nodeIDs[outNeighbors.size()-1]);
			assert(nodeOffset[outNeighbors.size()-2] + SPLIT_NODE_SIZE == nodeOffset[outNeighbors.size()-1]);
			outNeighbors[outNeighbors.size()-2].push_back(outNeighbors.size()-1);
			inNeighbors[inNeighbors.size()-1].push_back(inNeighbors.size()-2);
		}
	}
	unitigStartNode[nodeId] = outNeighbors.size();
	for (size_t i = DBGOverlap; i < sequence.size(); i += SPLIT_NODE_SIZE)
	{
		AddNode(nodeId, i, sequence.substr(i, SPLIT_NODE_SIZE), reverseNode, reid);
		if (i > 0)
		{
			assert(outNeighbors.size() >= 2);
			assert(outNeighbors.size() == inNeighbors.size());
			assert(nodeIDs.size() == outNeighbors.size());
			assert(nodeOffset.size() == outNeighbors.size());
			assert(nodeIDs[outNeighbors.size()-2] == nodeIDs[outNeighbors.size()-1]);
			assert(nodeOffset[outNeighbors.size()-2] + nodeLength[outNeighbors.size()-2] == nodeOffset[outNeighbors.size()-1]);
			outNeighbors[outNeighbors.size()-2].push_back(outNeighbors.size()-1);
			inNeighbors[inNeighbors.size()-1].push_back(inNeighbors.size()-2);
		}
	}
}

void AlignmentGraph::AddNode(int nodeId, int offset, const std::string& sequence, bool reverseNode, size_t reid)
{
	assert(!finalized);
	assert(sequence.size() <= SPLIT_NODE_SIZE);

	nodeUnitigReid.push_back(reid);
	nodeLookup[nodeId].push_back(nodeLength.size());
	nodeLength.push_back(sequence.size());
	nodeIDs.push_back(nodeId);
	inNeighbors.emplace_back();
	outNeighbors.emplace_back();
	reverse.push_back(reverseNode);
	nodeOffset.push_back(offset);
	nodeSequences.emplace_back();
	for (size_t i = 0; i < sequence.size(); i++)
	{
		size_t chunk = i / BP_IN_CHUNK;
		size_t offset = (i % BP_IN_CHUNK) * 2;
		size_t c = 0;
		switch(sequence[i])
		{
			case 'a':
			case 'A':
				c = 0;
				break;
			case 'c':
			case 'C':
				c = 1;
				break;
			case 'g':
			case 'G':
				c = 2;
				break;
			case 't':
			case 'T':
				c = 3;
				break;
			default:
				assert(false);
		}
		nodeSequences.back()[chunk] |= c << offset;
	}
	assert(nodeIDs.size() == nodeLength.size());
	assert(nodeLength.size() == inNeighbors.size());
	assert(inNeighbors.size() == outNeighbors.size());
}

void AlignmentGraph::AddEdgeNodeId(int node_id_from, int node_id_to)
{
	assert(!finalized);
	assert(nodeLookup.count(node_id_from) > 0);
	assert(nodeLookup.count(node_id_to) > 0);
	auto from = nodeLookup[node_id_from].back();
	auto to = unitigStartNode[node_id_to];
	assert(to >= 0);
	assert(from >= 0);
	assert(to < inNeighbors.size());
	assert(from < nodeLength.size());

	assert(unitigReidLookup.count(node_id_from) == 1);
	assert(unitigReidLookup.count(node_id_to) == 1);
	auto reidFrom = unitigReidLookup[node_id_from];
	auto reidTo = unitigReidLookup[node_id_to];
	assert(reidFrom < unitigReidOutNeighbors.size());
	assert(reidTo < unitigReidInNeighbors.size());
	unitigReidOutNeighbors[reidFrom].push_back(reidTo);
	unitigReidInNeighbors[reidTo].push_back(reidFrom);

	//don't add double edges
	if (std::find(inNeighbors[to].begin(), inNeighbors[to].end(), from) == inNeighbors[to].end()) inNeighbors[to].push_back(from);
	if (std::find(outNeighbors[from].begin(), outNeighbors[from].end(), to) == outNeighbors[from].end()) outNeighbors[from].push_back(to);
}

void AlignmentGraph::Finalize(int wordSize)
{
	assert(nodeSequences.size() == nodeLength.size());
	assert(inNeighbors.size() == nodeLength.size());
	assert(outNeighbors.size() == nodeLength.size());
	assert(reverse.size() == nodeLength.size());
	assert(nodeIDs.size() == nodeLength.size());
	assert(nodeUnitigReid.size() == nodeLength.size());
	assert(unitigReidOutNeighbors.size() == unitigReidInNeighbors.size());
	assert(unitigReidLookup.size() == unitigReidInNeighbors.size());
	assert(unitigReidLength.size() == unitigReidInNeighbors.size());
	std::cout << nodeLookup.size() << " original nodes" << std::endl;
	std::cout << nodeLength.size() << " split nodes" << std::endl;
	finalized = true;
	int specialNodes = 0;
	size_t edges = 0;
	for (size_t i = 0; i < inNeighbors.size(); i++)
	{
		inNeighbors[i].shrink_to_fit();
		outNeighbors[i].shrink_to_fit();
		if (inNeighbors[i].size() >= 2) specialNodes++;
		edges += inNeighbors[i].size();
	}
	std::cout << edges << " edges" << std::endl;
	std::cout << specialNodes << " nodes with in-degree >= 2" << std::endl;
	assert(nodeSequences.size() == nodeLength.size());
	assert(inNeighbors.size() == nodeLength.size());
	assert(outNeighbors.size() == nodeLength.size());
	assert(reverse.size() == nodeLength.size());
	assert(nodeIDs.size() == nodeLength.size());
	assert(nodeOffset.size() == nodeLength.size());
	nodeLength.shrink_to_fit();
	nodeIDs.shrink_to_fit();
	inNeighbors.shrink_to_fit();
	outNeighbors.shrink_to_fit();
	reverse.shrink_to_fit();
	nodeSequences.shrink_to_fit();
}

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
size_t AlignmentGraph::NodeLength(size_t index) const
{
	return nodeLength[index];
}

char AlignmentGraph::NodeSequences(size_t node, size_t pos) const
{
	assert(node < nodeSequences.size());
	assert(pos < nodeLength[node]);
	size_t chunk = pos / BP_IN_CHUNK;
	size_t offset = (pos % BP_IN_CHUNK) * 2;
	return "ACGT"[(nodeSequences[node][chunk] >> offset) & 3];
}

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
AlignmentGraph::NodeChunkSequence AlignmentGraph::NodeChunks(size_t index) const
{
	assert(index < nodeSequences.size());
	return nodeSequences[index];
}

size_t AlignmentGraph::NodeSize() const
{
	return nodeLength.size();
}

class NodeWithDistance
{
public:
	NodeWithDistance(size_t node, bool start, size_t distance) : node(node), start(start), distance(distance) {};
	bool operator>(const NodeWithDistance& other) const
	{
		return distance > other.distance;
	}
	size_t node;
	bool start;
	size_t distance;
};

size_t AlignmentGraph::GetUnitigNode(int nodeId, size_t offset) const
{
	auto nodes = nodeLookup.at(nodeId);
	size_t index = 0;
	while (index < nodes.size() && (nodeOffset[nodes[index]] > offset || nodeOffset[nodes[index]] + NodeLength(nodes[index]) <= offset)) index++;
	assert(index != nodes.size());
	size_t result = nodes[index];
	assert(nodeOffset[result] <= offset);
	assert(nodeOffset[result] + NodeLength(result) > offset);
	return result;
}

std::pair<int, size_t> AlignmentGraph::GetReversePosition(int nodeId, size_t offset) const
{
	assert(nodeLookup.count(nodeId) == 1);
	auto nodes = nodeLookup.at(nodeId);
	size_t originalSize = originalNodeSize.at(nodeId);
	assert(offset < originalSize);
	size_t newOffset = originalSize - offset - 1;
	assert(newOffset < originalSize);
	int reverseNodeId;
	if (nodeId % 2 == 0)
	{
		reverseNodeId = (nodeId / 2) * 2 + 1;
	}
	else
	{
		reverseNodeId = (nodeId / 2) * 2;
	}
	return std::make_pair(reverseNodeId, newOffset);
}


// size_t AlignmentGraph::GetReverseNode(size_t node) const
// {
// 	assert(node < nodeLength.size());

// 	size_t originalNodeSize = (nodeLookup.at(nodeIDs[node]).size() - 1) * SPLIT_NODE_SIZE + NodeLength(nodeLookup.at(nodeIDs[node]).back());
// 	size_t currentOffset = nodeOffset[node];
// 	assert(currentOffset < originalNodeSize);
// 	size_t reverseOffset = originalNodeSize - currentOffset - 1;
// 	assert(reverseOffset < originalNodeSize);
// 	size_t reverseNodeOriginalId;
// 	if (nodeIDs[node] % 2 == 0)
// 	{
// 		reverseNodeOriginalId = (nodeIDs[node] / 2) * 2 + 1;
// 	}
// 	else
// 	{
// 		reverseNodeOriginalId = (nodeIDs[node] / 2) * 2;
// 	}
// 	assert(nodeLookup.count(reverseNodeOriginalId) == 1);
// 	assert(nodeLookup.at(reverseNodeOriginalId).size() > reverseOffset / SPLIT_NODE_SIZE);
// 	size_t reverseNode = nodeLookup.at(reverseNodeOriginalId)[reverseOffset / SPLIT_NODE_SIZE];

// 	return reverseNode;
// }

AlignmentGraph::MatrixPosition::MatrixPosition(size_t node, size_t nodeOffset, size_t seqPos) :
	node(node),
	nodeOffset(nodeOffset),
	seqPos(seqPos)
{
}

bool AlignmentGraph::MatrixPosition::operator==(const AlignmentGraph::MatrixPosition& other) const
{
	return node == other.node && nodeOffset == other.nodeOffset && seqPos == other.seqPos;
}

bool AlignmentGraph::MatrixPosition::operator!=(const AlignmentGraph::MatrixPosition& other) const
{
	return !(*this == other);
}

AlignmentGraph AlignmentGraph::GetSubgraph(const std::unordered_map<size_t, size_t>& nodeMapping) const
{
	AlignmentGraph result;
	result.nodeLength.resize(nodeMapping.size());
	result.nodeOffset.resize(nodeMapping.size());
	result.nodeIDs.resize(nodeMapping.size());
	result.inNeighbors.resize(nodeMapping.size());
	result.outNeighbors.resize(nodeMapping.size());
	result.reverse.resize(nodeMapping.size());
	result.nodeSequences.resize(nodeMapping.size());

	for (auto pair : nodeMapping)
	{
		for (auto inNeighbor : inNeighbors[pair.first])
		{
			if (nodeMapping.count(inNeighbor) == 1)
			{
				result.inNeighbors[pair.second].push_back(nodeMapping.at(inNeighbor));
			}
		}
		for (auto outNeighbor : outNeighbors[pair.first])
		{
			if (nodeMapping.count(outNeighbor) == 1)
			{
				result.outNeighbors[pair.second].push_back(nodeMapping.at(outNeighbor));
			}
		}
		result.nodeLength[pair.second] = nodeLength[pair.first];
		result.nodeOffset[pair.second] = nodeOffset[pair.first];
		result.nodeIDs[pair.second] = nodeIDs[pair.first];
		result.reverse[pair.second] = reverse[pair.first];
		result.nodeSequences[pair.second] = nodeSequences[pair.first];
		result.nodeLookup[result.nodeIDs[pair.second]].push_back(pair.second);
		result.originalNodeSize[result.nodeIDs[pair.second]] = originalNodeSize.at(nodeIDs[pair.first]);
	}

	result.finalized = true;
	return result;
}

size_t AlignmentGraph::UnitigReidSize() const
{
	return unitigReidInNeighbors.size();
}

size_t AlignmentGraph::ReidLength(size_t reid) const
{
	return unitigReidLength[reid];
}
