#include <iostream>
#include <limits>
#include <algorithm>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "AlignmentGraph.h"
#include "CommonUtils.h"

AlignmentGraph::AlignmentGraph() :
nodeStart(),
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
	nodeEnd.emplace_back(nodeSequences.size());
	notInOrder.push_back(false);
}

void AlignmentGraph::ReserveNodes(size_t numNodes, size_t sequenceLength)
{
	numNodes += 2; //dummy start and end nodes
	sequenceLength += 2; //dummy start and end nodes
	nodeSequences.reserve(sequenceLength);
	nodeLookup.reserve(numNodes);
	nodeIDs.reserve(numNodes);
	nodeStart.reserve(numNodes);
	inNeighbors.reserve(numNodes);
	outNeighbors.reserve(numNodes);
	reverse.reserve(numNodes);
	nodeEnd.reserve(numNodes);
	notInOrder.reserve(numNodes);
}

void AlignmentGraph::AddNode(int nodeId, const std::string& sequence, bool reverseNode)
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
	nodeEnd.push_back(nodeSequences.size());
	notInOrder.push_back(false);
	assert(nodeIDs.size() == nodeStart.size());
	assert(nodeStart.size() == inNeighbors.size());
	assert(inNeighbors.size() == nodeEnd.size());
	assert(nodeEnd.size() == notInOrder.size());
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

	inNeighbors[to].insert(from);
	outNeighbors[from].insert(to);
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
	std::cerr << nodeStart.size() << " nodes" << std::endl;
	std::cerr << nodeSequences.size() << "bp" << std::endl;
	finalized = true;
	int specialNodes = 0;
	size_t edges = 0;
	for (size_t i = 0; i < inNeighbors.size(); i++)
	{
		if (inNeighbors[i].size() >= 2) specialNodes++;
		edges += inNeighbors[i].size();
	}
	std::cerr << edges << " edges" << std::endl;
	std::cerr << specialNodes << " nodes with in-degree >= 2" << std::endl;
#ifndef NDEBUG
	assert(nodeSequences.size() >= nodeStart.size());
	assert(nodeEnd.size() == nodeStart.size());
	assert(notInOrder.size() == nodeStart.size());
	assert(inNeighbors.size() == nodeStart.size());
	assert(outNeighbors.size() == nodeStart.size());
	assert(notInOrder.size() == nodeStart.size());
	assert(reverse.size() == nodeStart.size());
	assert(nodeIDs.size() == nodeStart.size());
	for (size_t i = 0; i < nodeStart.size(); i++)
	{
		assert(nodeEnd[i] > nodeStart[i]);
		assert(nodeEnd[i] <= nodeSequences.size());
		if (i > 0) assert(nodeStart[i] == nodeEnd[i-1]);
	}
#endif
	notInOrder.shrink_to_fit();
	nodeStart.shrink_to_fit();
	nodeEnd.shrink_to_fit();
	nodeIDs.shrink_to_fit();
	inNeighbors.shrink_to_fit();
	outNeighbors.shrink_to_fit();
	reverse.shrink_to_fit();
	nodeSequences.shrink_to_fit();
}

size_t AlignmentGraph::SizeInBp() const
{
	return nodeSequences.size();
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
			auto end = nodeEnd[nodeIndex];
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
	assert(nodeEnd[otherNode] - nodeStart[otherNode] == nodeEnd[nodeIndex] - nodeStart[nodeIndex]);
	// assert(nodeSequences.substr(nodeStart[nodeIndex], nodeEnd[nodeIndex] - nodeStart[nodeIndex]) == CommonUtils::ReverseComplement(nodeSequences.substr(nodeStart[otherNode], nodeEnd[otherNode] - nodeStart[otherNode])));
	return otherNode;
}

size_t AlignmentGraph::GetReversePosition(size_t pos) const
{
	assert(pos < nodeSequences.size());
	assert(pos > 0);
	auto originalNode = IndexToNode(pos);
	auto otherNode = GetReverseNode(originalNode);
	size_t newPos = (nodeEnd[otherNode] - 1) - (pos - nodeStart[originalNode]);
	// assert(nodeSequences.substr(pos, 1) == CommonUtils::ReverseComplement(nodeSequences.substr(newPos, 1)));
	return newPos;
}

size_t AlignmentGraph::IndexToNode(size_t index) const
{
	assert(index < nodeSequences.size());
	auto nextnode = std::upper_bound(nodeStart.begin(), nodeStart.end(), index);
	auto nextindex = nextnode - nodeStart.begin();
	assert(nextindex > 0);
	assert(nextindex <= nodeStart.size());
	return nextindex-1;
}
