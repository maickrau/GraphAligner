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
#include "ThreadReadAssertion.h"

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
}

void AlignmentGraph::ReserveNodes(size_t numNodes, size_t numSplitNodes, size_t sequenceLength)
{
	nodeSequencesATorCG.reserve(sequenceLength);
	nodeSequencesACorTG.reserve(sequenceLength);
	nodeLookup.reserve(numNodes);
	nodeIDs.reserve(numSplitNodes);
	nodeStart.reserve(numSplitNodes);
	inNeighbors.reserve(numSplitNodes);
	outNeighbors.reserve(numSplitNodes);
	reverse.reserve(numSplitNodes);
	nodeOffset.reserve(numSplitNodes);
}

void AlignmentGraph::AddNode(int nodeId, const std::string& sequence, bool reverseNode)
{
	assert(!finalized);
	//subgraph extraction might produce different subgraphs with common nodes
	//don't add duplicate nodes
	if (nodeLookup.count(nodeId) != 0) return;
	for (size_t i = 0; i < sequence.size(); i += SPLIT_NODE_SIZE)
	{
		AddNode(nodeId, i, sequence.substr(i, SPLIT_NODE_SIZE), reverseNode);
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
}

void AlignmentGraph::AddNode(int nodeId, int offset, std::string sequence, bool reverseNode)
{
	assert(!finalized);

	assert(std::numeric_limits<size_t>::max() - sequence.size() > nodeSequencesATorCG.size());
	nodeLookup[nodeId].push_back(nodeStart.size());
	nodeIDs.push_back(nodeId);
	nodeStart.push_back(nodeSequencesATorCG.size());
	inNeighbors.emplace_back();
	outNeighbors.emplace_back();
	reverse.push_back(reverseNode);
	nodeOffset.push_back(offset);
	for (auto c : sequence)
	{
		switch(c)
		{
			case 'A':
			case 'a':
			nodeSequencesATorCG.push_back(false);
			nodeSequencesACorTG.push_back(false);
			break;
			case 't':
			case 'T':
			nodeSequencesATorCG.push_back(false);
			nodeSequencesACorTG.push_back(true);
			break;
			case 'c':
			case 'C':
			nodeSequencesATorCG.push_back(true);
			nodeSequencesACorTG.push_back(false);
			break;
			case 'g':
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
	auto from = nodeLookup[node_id_from].back();
	auto to = nodeLookup[node_id_to][0];
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
	assert(nodeSequencesATorCG.size() == nodeSequencesACorTG.size());
	assert(nodeSequencesATorCG.size() >= nodeStart.size());
	assert(inNeighbors.size() == nodeStart.size());
	assert(outNeighbors.size() == nodeStart.size());
	assert(reverse.size() == nodeStart.size());
	assert(nodeIDs.size() == nodeStart.size());
	std::cout << nodeLookup.size() << " original nodes" << std::endl;
	std::cout << nodeStart.size() << " split nodes" << std::endl;
	std::cout << nodeSequencesATorCG.size() << "bp" << std::endl;
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
	assert(nodeSequencesATorCG.size() == nodeSequencesACorTG.size());
	assert(nodeSequencesATorCG.size() >= nodeStart.size());
	assert(inNeighbors.size() == nodeStart.size());
	assert(outNeighbors.size() == nodeStart.size());
	assert(reverse.size() == nodeStart.size());
	assert(nodeIDs.size() == nodeStart.size());
	assert(nodeOffset.size() == nodeStart.size());
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

size_t AlignmentGraph::GetReversePosition(size_t pos) const
{
	assert(pos < nodeSequencesATorCG.size());
	assert(pos > 0);
	size_t forwardNode = IndexToNode(pos);
	size_t originalNodeSize = (nodeLookup.at(nodeIDs[forwardNode]).size() - 1) * SPLIT_NODE_SIZE + NodeLength(nodeLookup.at(nodeIDs[forwardNode]).back());
	size_t currentOffset = pos - NodeStart(forwardNode) + nodeOffset[forwardNode];
	assert(currentOffset < originalNodeSize);
	size_t reverseOffset = originalNodeSize - currentOffset - 1;
	assert(reverseOffset < originalNodeSize);
	size_t reverseNodeOriginalId;
	if (nodeIDs[forwardNode] % 2 == 0)
	{
		reverseNodeOriginalId = (nodeIDs[forwardNode] / 2) * 2 + 1;
	}
	else
	{
		reverseNodeOriginalId = (nodeIDs[forwardNode] / 2) * 2;
	}
	assert(nodeLookup.count(reverseNodeOriginalId) == 1);
	assert(nodeLookup.at(reverseNodeOriginalId).size() > reverseOffset / SPLIT_NODE_SIZE);
	size_t reverseNode = nodeLookup.at(reverseNodeOriginalId)[reverseOffset / SPLIT_NODE_SIZE];
	size_t reverseNodeOffset = reverseOffset % SPLIT_NODE_SIZE;
	size_t newPos = NodeStart(reverseNode) + reverseNodeOffset;
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
		if (validNodes.count(node) == 1)
		{
			for (auto target : targets)
			{
				if (IndexToNode(target) != node) continue;
				if (pos <= target) mindist = std::min(mindist, target - pos);
				if (target <= pos) mindist = std::min(mindist, pos - target);
			}
		}
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

void AlignmentGraph::connect(size_t node, std::vector<std::vector<size_t>>& result, size_t& indexnum, std::vector<size_t>& index, std::vector<size_t>& lowlink, std::vector<bool>& onStack, std::vector<size_t>& S) const
{
	index[node] = indexnum;
	lowlink[node] = indexnum;
	indexnum++;
	S.push_back(node);
	onStack[node] = true;
	for (auto neighbor : outNeighbors[node])
	{
		if (index[neighbor] == -1)
		{
			connect(neighbor, result, indexnum, index, lowlink, onStack, S);
			lowlink[node] = std::min(lowlink[neighbor], lowlink[node]);
		}
		else if (onStack[neighbor])
		{
			lowlink[node] = std::min(lowlink[node], index[neighbor]);
		}
	}
	if (lowlink[node] == index[node])
	{
		result.emplace_back();
		while (true)
		{
			auto w = S.back();
			S.pop_back();
			onStack[w] = false;
			result.back().push_back(w);
			if (w == node) break;
		}
	}
}

//https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
std::pair<std::vector<std::vector<size_t>>, std::vector<size_t>> AlignmentGraph::TopologicalOrderOfComponents() const
{
	std::vector<size_t> index;
	std::vector<size_t> lowlink;
	std::vector<bool> onStack;
	std::vector<size_t> S;
	index.resize(NodeSize(), -1);
	lowlink.resize(NodeSize(), -1);
	onStack.resize(NodeSize(), false);
	size_t indexnum = 0;
	std::vector<std::vector<size_t>> result;
	for (size_t i = 0; i < NodeSize(); i++)
	{
		if (index[i] == -1) connect(i, result, indexnum, index, lowlink, onStack, S);
	}
	std::reverse(result.begin(), result.end());
	assert(result.size() > 0);
	assert(result[0].size() > 0);
	assert(result.back().size() > 0);
	std::vector<size_t> belongsToComponent;
	belongsToComponent.resize(NodeSize(), -1);
	for (size_t component = 0; component < result.size(); component++)
	{
		assert(result[component].size() > 0);
		for (auto node : result[component])
		{
			belongsToComponent[node] = component;
		}
	}
#ifndef NDEBUG
	//verify that it's correct
	for (size_t i = 0; i < NodeSize(); i++)
	{
		assert(belongsToComponent[i] != -1);
		for (auto neighbor : inNeighbors[i])
		{
			assert(belongsToComponent[neighbor] <= belongsToComponent[i]);
		}
	}
#endif
	return std::make_pair(result, belongsToComponent);
}
