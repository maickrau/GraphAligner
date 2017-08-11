#include <iostream>
#include <cassert>
#include "CommonUtils.h"
#include "ssw_cpp.h"
#include "vg.pb.h"
#include "fastqloader.h"
#include "BigraphToDigraph.h"

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

DirectedGraph::DirectedGraph()
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
		nodes.emplace_back(bigraph.node(i).id() * 2 + 1, bigraph.node(i).id(), false, CommonUtils::ReverseComplement(bigraph.node(i).sequence()));
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

void DirectedGraph::RemoveNodes(const std::set<int>& nodeIndices)
{
	assert(edgesPointToValidNodes());
	assert(nodeIdsAreValid());
	size_t nodesBefore = nodes.size();
	size_t offset = 0;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		while (nodeIndices.count(i+offset) > 0) offset++;
		assert(i+offset <= nodes.size());
		if (i+offset == nodes.size()) break;
		if (offset > 0) nodes[i] = nodes[i+offset];
	}
	assert(offset == nodeIndices.size());
	nodes.erase(nodes.end()-offset, nodes.end());
	offset = 0;
	std::vector<size_t> newIndex;
	newIndex.resize(nodes.size(), 1);
	newIndex[0] = 0;
	for (auto x : nodeIndices)
	{
		assert((x == 0 && newIndex[x] == 0) || (x != 0 && newIndex[x] == 1));
		newIndex[x] -= 1;
	}
	for (size_t i = 1; i < nodes.size(); i++)
	{
		newIndex[i] += newIndex[i-1];
	}
	assert(newIndex[nodes.size()-1] == nodes.size()-1-nodeIndices.size());
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

std::tuple<size_t, size_t, size_t> DirectedGraph::getLocalAlignmentSsw(std::string sequence, size_t nodeIndex) const
{
	StripedSmithWaterman::Aligner aligner {2, 5, 2, 1};
	aligner.SetReferenceSequence(sequence.c_str(), sequence.size());
	StripedSmithWaterman::Filter filter;
	filter.report_begin_position = true;
	filter.report_cigar = false;
	StripedSmithWaterman::Alignment alignment;
	aligner.Align(nodes[nodeIndex].sequence.c_str(), filter, &alignment);
	return std::make_tuple(alignment.query_begin, alignment.ref_begin, alignment.sw_score);
}

std::vector<DirectedGraph::SeedHit> DirectedGraph::GetSeedHits(const std::string& sequence, const std::vector<std::pair<int, size_t>>& hitsOriginalNodeIds) const
{
	std::vector<SeedHit> result;
	for (size_t i = 0; i < hitsOriginalNodeIds.size(); i++)
	{
		std::pair<SeedHit, size_t> bestHit {{0, 0, 0}, 0};
		for (size_t j = 0; j < nodes.size(); j++)
		{
			if (nodes[j].originalNodeId == hitsOriginalNodeIds[i].first)
			{
				auto hit = getLocalAlignmentSsw(sequence.substr(hitsOriginalNodeIds[i].second, 200), j);
				bestHit = std::make_pair(SeedHit(nodes[j].nodeId, std::get<0>(hit), std::get<1>(hit)+hitsOriginalNodeIds[i].second), std::get<2>(hit));
				if (bestHit.second > 0)
				{
					result.push_back(bestHit.first);
				}
			}
		}
	}
	return result;
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

void DirectedGraph::PruneByReachability(const std::vector<int>& startNodeIds)
{
	assert(edgesPointToValidNodes());
	assert(nodeIdsAreValid());
	std::vector<std::vector<size_t>> outNeighbors;
	std::vector<std::vector<size_t>> inNeighbors;
	outNeighbors.resize(nodes.size());
	inNeighbors.resize(nodes.size());
	for (size_t i = 0; i < edges.size(); i++)
	{
		outNeighbors[edges[i].fromIndex].push_back(edges[i].toIndex);
		inNeighbors[edges[i].toIndex].push_back(edges[i].fromIndex);
	}
	std::vector<bool> reachableOut;
	std::vector<bool> reachableIn;
	reachableOut.resize(nodes.size(), false);
	reachableIn.resize(nodes.size(), false);
	std::set<int> start;
	start.insert(startNodeIds.begin(), startNodeIds.end());
	for (size_t i = 0; i < nodes.size(); i++)
	{
		if (start.count(nodes[i].nodeId) > 0)
		{
			addReachable(reachableOut, outNeighbors, i);
			addReachable(reachableIn, inNeighbors, i);
		}
	}
	std::set<int> removeThese;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		if (!reachableOut[i] && !reachableIn[i]) removeThese.insert(i);
	}
	RemoveNodes(removeThese);
	assert(edgesPointToValidNodes());
	assert(nodeIdsAreValid());
}
