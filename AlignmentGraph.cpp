#include <iostream>
#include <limits>
#include <algorithm>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <queue>
#include "AlignmentGraph.h"
#include "CommonUtils.h"

AlignmentGraph::AlignmentGraph() :
DBGOverlap(0),
nodeStart(),
nodeLookup(),
nodeIDs(),
inNeighbors(),
nodeSequencesATorCG(),
nodeSequencesACorTG(),
finalized(false)
{
	//add the start dummy node as the first node
	dummyNodeStart = 0;
	nodeIDs.push_back(0);
	nodeStart.push_back(0);
	inNeighbors.emplace_back();
	outNeighbors.emplace_back();
	reverse.push_back(false);
	nodeSequencesATorCG.push_back(false);
	nodeSequencesACorTG.push_back(false);
}

void AlignmentGraph::ReserveNodes(size_t numNodes, size_t sequenceLength)
{
	numNodes += 2; //dummy start and end nodes
	sequenceLength += 2; //dummy start and end nodes
	nodeSequencesATorCG.reserve(sequenceLength);
	nodeSequencesACorTG.reserve(sequenceLength);
	nodeLookup.reserve(numNodes);
	nodeIDs.reserve(numNodes);
	nodeStart.reserve(numNodes);
	inNeighbors.reserve(numNodes);
	outNeighbors.reserve(numNodes);
	reverse.reserve(numNodes);
}

void AlignmentGraph::AddNode(int nodeId, const std::string& sequence, bool reverseNode)
{
	assert(!finalized);
	//subgraph extraction might produce different subgraphs with common nodes
	//don't add duplicate nodes
	if (nodeLookup.count(nodeId) != 0) return;

	assert(std::numeric_limits<size_t>::max() - sequence.size() > nodeSequencesATorCG.size());
	nodeLookup[nodeId] = nodeStart.size();
	nodeIDs.push_back(nodeId);
	nodeStart.push_back(nodeSequencesATorCG.size());
	inNeighbors.emplace_back();
	outNeighbors.emplace_back();
	reverse.push_back(reverseNode);
	for (auto c : sequence)
	{
		switch(c)
		{
			case 'A':
				nodeSequencesATorCG.push_back(false);
				nodeSequencesACorTG.push_back(false);
				break;
			case 'T':
				nodeSequencesATorCG.push_back(false);
				nodeSequencesACorTG.push_back(true);
				break;
			case 'C':
				nodeSequencesATorCG.push_back(true);
				nodeSequencesACorTG.push_back(false);
				break;
			case 'G':
				nodeSequencesATorCG.push_back(true);
				nodeSequencesACorTG.push_back(true);
				break;
			default:
				assert(false);
				std::abort();
		}
	}
	assert(nodeIDs.size() == nodeStart.size());
	assert(nodeStart.size() == inNeighbors.size());
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

	//don't add double edges
	if (std::find(inNeighbors[to].begin(), inNeighbors[to].end(), from) == inNeighbors[to].end()) inNeighbors[to].push_back(from);
	if (std::find(outNeighbors[from].begin(), outNeighbors[from].end(), to) == outNeighbors[from].end()) outNeighbors[from].push_back(to);
}

void AlignmentGraph::Finalize(int wordSize)
{
	//add the end dummy node as the last node
	dummyNodeEnd = nodeSequencesATorCG.size();
	nodeIDs.push_back(0);
	nodeStart.push_back(nodeSequencesATorCG.size());
	reverse.push_back(false);
	inNeighbors.emplace_back();
	outNeighbors.emplace_back();
	nodeSequencesATorCG.push_back(false);
	nodeSequencesACorTG.push_back(false);
	assert(nodeSequencesATorCG.size() == nodeSequencesACorTG.size());
	assert(nodeSequencesATorCG.size() >= nodeStart.size());
	assert(inNeighbors.size() == nodeStart.size());
	assert(outNeighbors.size() == nodeStart.size());
	assert(reverse.size() == nodeStart.size());
	assert(nodeIDs.size() == nodeStart.size());
	std::cerr << nodeStart.size() << " nodes" << std::endl;
	std::cerr << nodeSequencesATorCG.size() << "bp" << std::endl;
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
	std::cerr << edges << " edges" << std::endl;
	std::cerr << specialNodes << " nodes with in-degree >= 2" << std::endl;
#ifndef NDEBUG
	assert(nodeSequencesATorCG.size() == nodeSequencesACorTG.size());
	assert(nodeSequencesATorCG.size() >= nodeStart.size());
	assert(inNeighbors.size() == nodeStart.size());
	assert(outNeighbors.size() == nodeStart.size());
	assert(reverse.size() == nodeStart.size());
	assert(nodeIDs.size() == nodeStart.size());
#endif
	nodeStart.shrink_to_fit();
	nodeIDs.shrink_to_fit();
	inNeighbors.shrink_to_fit();
	outNeighbors.shrink_to_fit();
	reverse.shrink_to_fit();
	nodeSequencesATorCG.shrink_to_fit();
	nodeSequencesACorTG.shrink_to_fit();
}

size_t AlignmentGraph::SizeInBp() const
{
	return nodeSequencesATorCG.size();
}

std::set<size_t> AlignmentGraph::ProjectForward(const std::set<size_t>& startpositions, size_t amount) const
{
	std::vector<std::set<size_t>> positions;
	positions.resize(amount+1);
	positions[0].insert(startpositions.begin(), startpositions.end());
	for (size_t i = 0; i < amount; i++)
	{
		auto left = amount - i;
		for (auto pos : positions[i])
		{
			auto nodeIndex = IndexToNode(pos);
			auto end = NodeEnd(nodeIndex);
			if (pos + left < end)
			{
				assert(i + end - pos > amount);
				positions.back().insert(pos + left);
			}
			else if (pos + left == end)
			{
				assert(i + end - pos == amount);
				for (auto neighbor : outNeighbors[nodeIndex])
				{
					positions.back().insert(nodeStart[neighbor]);
				}
			}
			else
			{
				assert(i + end - pos < amount);
				for (auto neighbor : outNeighbors[nodeIndex])
				{
					positions[i + end - pos].insert(nodeStart[neighbor]);
				}
			}
		}
	}
	return positions.back();
}

size_t AlignmentGraph::GetReverseNode(size_t nodeIndex) const
{
	auto bigraphNodeId = nodeIDs[nodeIndex] / 2;
	size_t otherNode;
	if (nodeIDs[nodeIndex] % 2 == 1)
	{
		otherNode = nodeLookup.at(bigraphNodeId * 2);
	}
	else
	{
		otherNode = nodeLookup.at(bigraphNodeId * 2 + 1);
	}
	assert(otherNode != nodeIndex);
	assert(NodeEnd(otherNode) - NodeStart(otherNode) == NodeEnd(nodeIndex) - NodeStart(nodeIndex));
	return otherNode;
}

size_t AlignmentGraph::GetReversePosition(size_t pos) const
{
	assert(pos < nodeSequencesATorCG.size());
	assert(pos > 0);
	auto originalNode = IndexToNode(pos);
	auto otherNode = GetReverseNode(originalNode);
	size_t newPos = (NodeEnd(otherNode) - 1) - (pos - nodeStart[originalNode]);
	return newPos;
}

size_t AlignmentGraph::IndexToNode(size_t index) const
{
	assert(index < nodeSequencesATorCG.size());
	auto nextnode = std::upper_bound(nodeStart.begin(), nodeStart.end(), index);
	auto nextindex = nextnode - nodeStart.begin();
	assert(nextindex > 0);
	assert(nextindex <= nodeStart.size());
	return nextindex-1;
}

size_t AlignmentGraph::NodeStart(size_t index) const
{
	return nodeStart[index];
}

size_t AlignmentGraph::NodeEnd(size_t index) const
{
	if (index == nodeStart.size()-1) return nodeSequencesATorCG.size();
	return nodeStart[index+1];
}

size_t AlignmentGraph::NodeLength(size_t index) const
{
	return NodeEnd(index) - NodeStart(index);
}

char AlignmentGraph::NodeSequences(size_t index) const
{
	assert(index < nodeSequencesATorCG.size());
	//dummy nodes
	if (index == 0 || index == nodeSequencesACorTG.size()-1) return '-';
	int first = nodeSequencesATorCG[index];
	int second = nodeSequencesACorTG[index];
	return "ATCG"[first*2+second];
}

size_t AlignmentGraph::NodeSequencesSize() const
{
	return nodeSequencesACorTG.size();
}

size_t AlignmentGraph::NodeSize() const
{
	return nodeStart.size();
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

size_t AlignmentGraph::MinDistance(size_t pos, const std::vector<size_t>& targets) const
{
	assert(targets.size() > 0);
	std::set<size_t> validNodes;
	for (auto target : targets)
	{
		validNodes.insert(IndexToNode(target));
	}
	std::unordered_map<size_t, size_t> distanceAtNodeEnd;
	std::unordered_map<size_t, size_t> distanceAtNodeStart;
	std::priority_queue<NodeWithDistance, std::vector<NodeWithDistance>, std::greater<NodeWithDistance>> queue;
	size_t mindist = std::numeric_limits<size_t>::max();
	{
		auto node = IndexToNode(pos);
		queue.emplace(node, true, pos - NodeStart(node));
		queue.emplace(node, false, NodeEnd(node) - 1 - pos);
	}
	size_t lastdist = 0;
	while (queue.size() > 0)
	{
		NodeWithDistance top = queue.top();
		assert(top.distance >= lastdist);
		lastdist = top.distance;
		if (top.distance >= mindist) break;
		queue.pop();
		if (top.start)
		{
			if (distanceAtNodeStart.count(top.node) > 0 && distanceAtNodeStart[top.node] <= top.distance) continue;
			distanceAtNodeStart[top.node] = top.distance;
		}
		else
		{
			if (distanceAtNodeEnd.count(top.node) > 0 && distanceAtNodeEnd[top.node] <= top.distance) continue;
			distanceAtNodeEnd[top.node] = top.distance;
		}
		if (validNodes.count(top.node) > 0)
		{
			for (auto target : targets)
			{
				if (IndexToNode(target) == top.node)
				{
					if (top.start)
					{
						mindist = std::min(mindist, top.distance + target - NodeStart(top.node));
					}
					else
					{
						mindist = std::min(mindist, top.distance + NodeEnd(top.node) - 1 - target);
					}
				}
			}
		}
		if (top.start)
		{
			queue.emplace(top.node, false, top.distance + NodeLength(top.node) - 1);
			for (auto neighbor : inNeighbors[top.node])
			{
				queue.emplace(neighbor, false, top.distance + 1);
			}
		}
		else
		{
			queue.emplace(top.node, true, top.distance + NodeLength(top.node) - 1);
			for (auto neighbor : outNeighbors[top.node])
			{
				queue.emplace(neighbor, true, top.distance + 1);
			}
		}
	}
	return mindist;
}