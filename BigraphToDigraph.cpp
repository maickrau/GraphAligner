#include "vg.pb.h"
#include "fastqloader.h"
#include "BigraphToDigraph.h"

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

DirectedGraph::DirectedGraph(const vg::Graph& bigraph)
{
	std::map<int, std::pair<size_t, size_t>> nodeMapping;
	for (int i = 0; i < bigraph.node_size(); i++)
	{
		assert(nodeMapping.count(bigraph.node(i).id()) == 0);
		nodeMapping[bigraph.node(i).id()] = std::make_pair(nodes.size()+1, nodes.size());
		nodes.emplace_back(bigraph.node(i).id() * 2, bigraph.node(i).id(), true, bigraph.node(i).sequence());
		nodes.emplace_back(bigraph.node(i).id() * 2 + 1, bigraph.node(i).id(), false, FastQ::reverseComplement(bigraph.node(i).sequence()));
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

void DirectedGraph::ReorderByNodeIds(const std::vector<int>& nodeIdOrder)
{
	assert(edgesPointToValidNodes());
	assert(nodeIdsAreValid());
	assert(nodeIdOrder.size() == nodes.size());
	std::vector<size_t> indexOrder;
	indexOrder.resize(nodes.size(), -1);
	std::map<int, size_t> mapping;
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
		newNodes.emplace_back(nodes[indexOrder[i]]);
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
	nodes = newNodes;
	assert(edgesPointToValidNodes());
	assert(nodeIdsAreValid());
}

size_t smallerThan(const std::set<int>& set, int comparison)
{
	size_t result = 0;
	for (auto x : set)
	{
		if (x < comparison) result++;
	}
	return result;
}

void DirectedGraph::RemoveNodes(const std::set<int>& nodeIndices)
{
	assert(edgesPointToValidNodes());
	assert(nodeIdsAreValid());
	size_t nodesBefore = nodes.size();
	for (size_t i = nodes.size()-1; i >= 0 && i < nodes.size(); i--)
	{
		if (nodeIndices.count(i) > 0)
		{
			nodes.erase(nodes.begin() + i);
		}
	}
	for (size_t i = edges.size()-1; i >= 0 && i < edges.size(); i--)
	{
		if (nodeIndices.count(edges[i].fromIndex) > 0 || nodeIndices.count(edges[i].toIndex) > 0)
		{
			edges.erase(edges.begin() + i);
		}
		edges[i].fromIndex -= smallerThan(nodeIndices, edges[i].fromIndex);
		edges[i].toIndex -= smallerThan(nodeIndices, edges[i].toIndex);
	}
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