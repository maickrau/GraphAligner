#include <iostream>
#include <cassert>
#include <unordered_map>
#include "CommonUtils.h"
#include "vg.pb.h"
#include "fastqloader.h"
#include "BigraphToDigraph.h"
#include "ThreadReadAssertion.h"

DirectedGraph::SeedHit::SeedHit(int nodeId, size_t nodePos, size_t seqPos) :
nodeId(nodeId),
nodePos(nodePos),
seqPos(seqPos)
{
}

DirectedGraph::Node::Node(int nodeId, int originalNodeId, bool rightEnd, std::string sequence) :
nodeId(nodeId),
originalNodeId(originalNodeId),
rightEnd(rightEnd),
sequence(sequence)
{
}

DirectedGraph::Edge::Edge(size_t from, size_t to) :
fromIndex(from),
toIndex(to)
{
}

DirectedGraph::DirectedGraph() :
totalSequenceLength(0)
{
}

DirectedGraph::DirectedGraph(const vg::Graph& bigraph) :
totalSequenceLength(0)
{
	std::unordered_map<int, std::pair<size_t, size_t>> nodeMapping;
	nodeMapping.reserve(bigraph.node_size());
	for (int i = 0; i < bigraph.node_size(); i++)
	{
		assert(nodeMapping.count(bigraph.node(i).id()) == 0);
		nodeMapping[bigraph.node(i).id()] = std::make_pair(nodes.size()+1, nodes.size());
		nodes.emplace_back(bigraph.node(i).id() * 2, bigraph.node(i).id(), true, bigraph.node(i).sequence());
		nodes.emplace_back(bigraph.node(i).id() * 2 + 1, bigraph.node(i).id(), false, CommonUtils::ReverseComplement(bigraph.node(i).sequence()));
		totalSequenceLength += bigraph.node(i).sequence().size() * 2;
	}
	for (int i = 0; i < bigraph.edge_size(); i++)
	{
		assert(nodeMapping.count(bigraph.edge(i).from()) == 1);
		assert(nodeMapping.count(bigraph.edge(i).to()) == 1);
		size_t fromLeft, fromRight, toLeft, toRight;
		if (bigraph.edge(i).from_start())
		{
			fromLeft = nodeMapping[bigraph.edge(i).from()].second;
			fromRight = nodeMapping[bigraph.edge(i).from()].first;
		}
		else
		{
			fromLeft = nodeMapping[bigraph.edge(i).from()].first;
			fromRight = nodeMapping[bigraph.edge(i).from()].second;
		}
		if (bigraph.edge(i).to_end())
		{
			toLeft = nodeMapping[bigraph.edge(i).to()].second;
			toRight = nodeMapping[bigraph.edge(i).to()].first;
		}
		else
		{
			toLeft = nodeMapping[bigraph.edge(i).to()].first;
			toRight = nodeMapping[bigraph.edge(i).to()].second;
		}
		edges.emplace_back(fromRight, toRight);
		edges.emplace_back(toLeft, fromLeft);
	}
	assert(edgesPointToValidNodes());
	assert(nodeIdsAreValid());
}

DirectedGraph::DirectedGraph(const GfaGraph& deBruijn) :
totalSequenceLength(0)
{
	std::unordered_map<int, std::pair<size_t, size_t>> nodeMapping;
	nodeMapping.reserve(deBruijn.nodes.size());
	for (int i = 0; i < deBruijn.nodes.size(); i++)
	{
		assert(nodeMapping.count(deBruijn.nodes[i].id) == 0);
		nodeMapping[deBruijn.nodes[i].id] = std::make_pair(nodes.size()+1, nodes.size());
		nodes.emplace_back(deBruijn.nodes[i].id * 2, deBruijn.nodes[i].id, true, deBruijn.nodes[i].sequence.substr(0, deBruijn.nodes[i].sequence.size() - deBruijn.edgeOverlap));
		nodes.emplace_back(deBruijn.nodes[i].id * 2 + 1, deBruijn.nodes[i].id, false, CommonUtils::ReverseComplement(deBruijn.nodes[i].sequence).substr(0, deBruijn.nodes[i].sequence.size() - deBruijn.edgeOverlap));
		totalSequenceLength += deBruijn.nodes[i].sequence.size() * 2;
	}
	for (int i = 0; i < deBruijn.edges.size(); i++)
	{
		assert(nodeMapping.count(deBruijn.edges[i].from) == 1);
		assert(nodeMapping.count(deBruijn.edges[i].to) == 1);
		size_t fromLeft, fromRight, toLeft, toRight;
		if (deBruijn.edges[i].fromStart)
		{
			fromLeft = nodeMapping[deBruijn.edges[i].from].second;
			fromRight = nodeMapping[deBruijn.edges[i].from].first;
		}
		else
		{
			fromLeft = nodeMapping[deBruijn.edges[i].from].first;
			fromRight = nodeMapping[deBruijn.edges[i].from].second;
		}
		if (deBruijn.edges[i].toEnd)
		{
			toLeft = nodeMapping[deBruijn.edges[i].to].second;
			toRight = nodeMapping[deBruijn.edges[i].to].first;
		}
		else
		{
			toLeft = nodeMapping[deBruijn.edges[i].to].first;
			toRight = nodeMapping[deBruijn.edges[i].to].second;
		}
		edges.emplace_back(fromRight, toRight);
		edges.emplace_back(toLeft, fromLeft);
	}
	assert(edgesPointToValidNodes());
	assert(nodeIdsAreValid());
}

void DirectedGraph::ReorderByNodeIds(const std::vector<int>& nodeIdOrder)
{
	assert(edgesPointToValidNodes());
	assert(nodeIdsAreValid());
	assert(nodeIdOrder.size() == nodes.size());
	std::vector<size_t> indexOrder;
	indexOrder.resize(nodes.size(), -1);
	std::unordered_map<int, size_t> mapping;
	mapping.reserve(nodes.size());
	for (size_t i = 0; i < nodes.size(); i++)
	{
		mapping[nodes[i].nodeId] = i;
	}
	for (size_t i = 0; i < nodeIdOrder.size(); i++)
	{
		indexOrder[i] = mapping[nodeIdOrder[i]];
	}
	std::vector<Node> newNodes;
	std::vector<size_t> reverseIndexOrder;
	reverseIndexOrder.resize(nodes.size(), -1);
	for (size_t i = 0; i < indexOrder.size(); i++)
	{
		assert(indexOrder[i] != -1);
		newNodes.emplace_back(std::move(nodes[indexOrder[i]]));
		reverseIndexOrder[indexOrder[i]] = i;
	}
	for (size_t i = 0; i < edges.size(); i++)
	{
		assert(reverseIndexOrder[edges[i].fromIndex] != -1);
		assert(reverseIndexOrder[edges[i].toIndex] != -1);
		assert(reverseIndexOrder[edges[i].fromIndex] < nodes.size());
		assert(reverseIndexOrder[edges[i].toIndex] < nodes.size());
		edges[i].fromIndex = reverseIndexOrder[edges[i].fromIndex];
		edges[i].toIndex = reverseIndexOrder[edges[i].toIndex];
	}
	assert(newNodes.size() == nodes.size());
	nodes = std::move(newNodes);
	assert(edgesPointToValidNodes());
	assert(nodeIdsAreValid());
}

void DirectedGraph::RemoveNodes(const std::set<int>& nodeIndices)
{
	assert(edgesPointToValidNodes());
	assert(nodeIdsAreValid());
	size_t nodesBefore = nodes.size();
	size_t offset = 0;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		while (nodeIndices.count(i+offset) > 0)
		{
			totalSequenceLength	-= nodes[i+offset].sequence.size();
			offset++;
		}
		assert(i+offset <= nodes.size());
		if (i+offset == nodes.size()) break;
		if (offset > 0) nodes[i] = std::move(nodes[i+offset]);
	}
	assert(offset == nodeIndices.size());
	nodes.erase(nodes.end()-offset, nodes.end());
	offset = 0;
	std::vector<size_t> newIndex;
	newIndex.resize(nodesBefore, 1);
	newIndex[0] = 0;
	for (auto x : nodeIndices)
	{
		assert(x < newIndex.size());
		assert(x >= 0);
		assert((x == 0 && newIndex[x] == 0) || (x != 0 && newIndex[x] == 1));
		newIndex[x] -= 1;
	}
	for (size_t i = 1; i < nodesBefore; i++)
	{
		newIndex[i] += newIndex[i-1];
	}
	assert(newIndex[nodesBefore-1] == nodesBefore-1-nodeIndices.size());
	for (size_t i = 0; i < edges.size(); i++)
	{
		while (i+offset < edges.size() && (nodeIndices.count(edges[i+offset].fromIndex) > 0 || nodeIndices.count(edges[i+offset].toIndex) > 0)) offset++;
		if (i+offset == edges.size()) break;
		assert(newIndex[edges[i+offset].fromIndex] >= 0);
		assert(newIndex[edges[i+offset].toIndex] >= 0);
		assert(edges[i+offset].fromIndex == 0 || newIndex[edges[i+offset].fromIndex] == newIndex[edges[i+offset].fromIndex-1]+1);
		assert(edges[i+offset].toIndex == 0 || newIndex[edges[i+offset].toIndex] == newIndex[edges[i+offset].toIndex-1]+1);
		edges[i].fromIndex = newIndex[edges[i+offset].fromIndex];
		edges[i].toIndex = newIndex[edges[i+offset].toIndex];
	}
	edges.erase(edges.end()-offset, edges.end());
	assert(nodes.size() + nodeIndices.size() == nodesBefore);
	assert(edgesPointToValidNodes());
	assert(nodeIdsAreValid());
}

bool DirectedGraph::edgesPointToValidNodes()
{
	for (size_t i = 0; i < edges.size(); i++)
	{
		if (edges[i].fromIndex >= nodes.size()) return false;
		if (edges[i].toIndex >= nodes.size()) return false;
	}
	return true;
}

bool DirectedGraph::nodeIdsAreValid()
{
	for (size_t i = 0; i < nodes.size(); i++)
	{
		if (nodes[i].nodeId / 2 != nodes[i].originalNodeId) return false;
	}
	return true;
}

void DirectedGraph::AddSubgraph(const DirectedGraph& subgraph)
{
	std::set<int> existingNodes;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		existingNodes.insert(nodes[i].nodeId);
	}
	for (size_t i = 0; i < subgraph.nodes.size(); i++)
	{
		if (existingNodes.count(subgraph.nodes[i].nodeId) == 0)
		{
			nodes.emplace_back(subgraph.nodes[i]);
			totalSequenceLength += subgraph.nodes[i].sequence.size();
		}
	}
	std::map<int, size_t> idToIndex;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		idToIndex[nodes[i].nodeId] = i;
	}
	std::set<std::pair<size_t, size_t>> existingEdges;
	for (size_t i = 0; i < edges.size(); i++)
	{
		existingEdges.insert(std::make_pair(edges[i].fromIndex, edges[i].toIndex));
	}
	for (size_t i = 0; i < subgraph.edges.size(); i++)
	{
		auto newFromIndex = idToIndex[subgraph.nodes[subgraph.edges[i].fromIndex].nodeId];
		auto newToIndex = idToIndex[subgraph.nodes[subgraph.edges[i].toIndex].nodeId];
		if (existingEdges.count(std::make_pair(newFromIndex, newToIndex)) == 0)
		{
			edges.emplace_back(newFromIndex, newToIndex);
		}
	}
}

void DirectedGraph::ConnectComponents(const std::vector<int>& previousSinksIds, const std::vector<int>& nextSourcesIds)
{
	std::map<int, size_t> idToIndex;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		idToIndex[nodes[i].nodeId] = i;
	}
	for (size_t i = 0; i < previousSinksIds.size(); i++)
	{
		for (size_t j = 0; j < nextSourcesIds.size(); j++)
		{
			edges.emplace_back(idToIndex[previousSinksIds[i]], idToIndex[nextSourcesIds[j]]);
		}
	}
}

void DirectedGraph::addReachable(std::vector<bool>& reachable, const std::vector<std::vector<size_t>>& outNeighbors, size_t current)
{
	assert(current < nodes.size());
	reachable[current] = true;
	for (size_t i = 0; i < outNeighbors[current].size(); i++)
	{
		if (reachable[outNeighbors[current][i]]) continue;
		addReachable(reachable, outNeighbors, outNeighbors[current][i]);
	}
}
