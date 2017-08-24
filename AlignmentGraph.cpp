#include <fstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "AlignmentGraph.h"
#include "TopologicalSort.h"
#include "CycleCutCalculation.h"

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

void AlignmentGraph::ReserveNodes(size_t numNodes, size_t sequenceLength)
{
	nodeSequences.reserve(sequenceLength);
	nodeLookup.reserve(numNodes);
	nodeIDs.reserve(numNodes);
	nodeStart.reserve(numNodes);
	inNeighbors.reserve(numNodes);
	outNeighbors.reserve(numNodes);
	reverse.reserve(numNodes);
	indexToNode.reserve(sequenceLength);
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
	for (size_t i = indexToNode.size(); i < nodeSequences.size(); i++)
	{
		indexToNode.push_back(nodeStart.size()-1);
	}
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

	inNeighbors[to].insert(from);
	outNeighbors[from].insert(to);
	if (from >= to)
	{
		notInOrder[to] = true;
	}
}

void AlignmentGraph::Finalize(int wordSize, std::string cutFilename)
{
	//add the end dummy node as the last node
	dummyNodeEnd = nodeSequences.size();
	nodeIDs.push_back(0);
	nodeStart.push_back(nodeSequences.size());
	reverse.push_back(false);
	inNeighbors.emplace_back();
	outNeighbors.emplace_back();
	nodeSequences.push_back('-');
	indexToNode.push_back('-');
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
	std::cerr << nodeStart.size() << " nodes" << std::endl;
	std::cerr << nodeSequences.size() << "bp" << std::endl;
	finalized = true;
	int specialNodes = 0;
	for (size_t i = 0; i < inNeighbors.size(); i++)
	{
		if (inNeighbors[i].size() >= 2) specialNodes++;
	}
	std::cerr << specialNodes << " nodes with in-degree >= 2" << std::endl;
	firstInOrder = 0;
	size_t inOrdersWrongfullyClassified = 0;
	for (size_t i = 1; i < notInOrder.size(); i++)
	{
		if (i > 1 && notInOrder[i] && !notInOrder[i-1])
		{
			inOrdersWrongfullyClassified += i - firstInOrder;
		}
		if (notInOrder[i]) firstInOrder = i+1;
		//relax this constraint for now until we figure out what is going on with the inOrdersWrongfullyClassified nodes
		// //all not-in-order nodes have to be at the start
		// assert(i == 1 || !notInOrder[i] || notInOrder[i-1]);
	}
	if (inOrdersWrongfullyClassified > 0)
	{
		std::cerr << inOrdersWrongfullyClassified << " nodes wrongly(?) classified as MFVS vertices!!" << std::endl;
	}
	if (firstInOrder != 0)
	{
		std::cerr << (firstInOrder - 1) << " nodes out of order" << std::endl;
	}
	else
	{
		std::cerr << "0 nodes out of order" << std::endl;
	}
	if (cutFilename != "")
	{
		if (!loadCycleCut(cutFilename))
		{
			calculateCycleCuts(wordSize);
			saveCycleCut(cutFilename);
		}
	}
	else
	{
		calculateCycleCuts(wordSize);
	}
	if (firstInOrder != 0)
	{
		std::cerr << "cycle cuts:" << std::endl;
		size_t totalCuttersbp = 0;
		for (size_t i = 1; i < cuts.size(); i++)
		{
			size_t cuttersbp = 0;
			for (size_t j = 0; j < cuts[i].nodes.size(); j++)
			{
				cuttersbp += nodeEnd[cuts[i].nodes[j]] - nodeStart[cuts[i].nodes[j]];
			}
			std::cerr << i << ": id " << nodeIDs[i] << ", cutting nodes " << cuts[i].nodes.size() << ", " << cuttersbp << "bp" << std::endl;
			totalCuttersbp += cuttersbp;
		}
		std::cerr << "total cut: " << totalCuttersbp << "bp (" << (double)totalCuttersbp / (double)nodeSequences.size() * 100 << "%)" << std::endl;
	}
}

void AlignmentGraph::calculateCycleCuts(int wordSize)
{
	CycleCutCalculation cutCalculator(*this);
	std::cerr << "calculating cycle cuts" << std::endl;
	cuts.resize(firstInOrder);
	for (size_t i = 1; i < firstInOrder; i++)
	{
		std::cerr << "cut " << i << "/" << firstInOrder << " node id " << nodeIDs[i] << std::endl;
		calculateCycleCutters(cutCalculator, i, wordSize);
	}
}

void AlignmentGraph::saveCycleCut(std::string filename)
{
	std::ofstream file {filename};
	boost::archive::text_oarchive oa(file);
	oa << cuts;
}

bool AlignmentGraph::loadCycleCut(std::string filename)
{
	std::ifstream file {filename};
	if (!file.good()) return false;
	boost::archive::text_iarchive ia(file);
	ia >> cuts;
	return true;
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

void AlignmentGraph::calculateCycleCutters(const CycleCutCalculation& cutCalculation, size_t cycleStart, int wordSize)
{
	assert(cuts.size() > cycleStart);
	assert(cuts[cycleStart].nodes.size() == 0);
	assert(cuts[cycleStart].predecessors.size() == 0);
	assert(cuts[cycleStart].previousCut.size() == 0);

	// cuts[cycleStart] = cutCalculation.GetCycleCutByDumbWay(cycleStart, wordSize);
	// cuts[cycleStart] = cutCalculation.GetCycleCutBySupersequence(cycleStart, wordSize);
	// cuts[cycleStart] = cutCalculation.GetCycleCutByIndex(cycleStart, wordSize);
	// cuts[cycleStart] = cutCalculation.GetCycleCutByOrder(cycleStart, wordSize);
	cuts[cycleStart] = cutCalculation.GetCycleCutBySupersequenceOverEdgeCoveringPaths(cycleStart, wordSize);

	assert(cuts[cycleStart].nodes.size() > 0);
	assert(cuts[cycleStart].predecessors.size() == cuts[cycleStart].nodes.size());
	assert(cuts[cycleStart].previousCut.size() == cuts[cycleStart].nodes.size());
	assert(cuts[cycleStart].nodes[0] == cycleStart);
	assert(cuts[cycleStart].predecessors.size() == 1 || cuts[cycleStart].predecessors[0].size() > 0);
}
