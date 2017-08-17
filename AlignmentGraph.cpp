#include <fstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "2dArray.h"
#include "AlignmentGraph.h"
#include "TopologicalSort.h"

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
			cycleCuttingNodes.resize(firstInOrder);
			cycleCuttingNodePredecessor.resize(firstInOrder);
			cycleCutPreviousCut.resize(firstInOrder);
			for (size_t i = 1; i < firstInOrder; i++)
			{
				calculateCycleCutters(i, wordSize);
			}
			saveCycleCut(cutFilename);
		}
	}
	else
	{
		cycleCuttingNodes.resize(firstInOrder);
		cycleCuttingNodePredecessor.resize(firstInOrder);
		cycleCutPreviousCut.resize(firstInOrder);
		for (size_t i = 1; i < firstInOrder; i++)
		{
			calculateCycleCutters(i, wordSize);
		}
	}
	if (firstInOrder != 0)
	{
		std::cerr << "cycle cuts:" << std::endl;
		size_t totalCuttersbp = 0;
		for (size_t i = 1; i < cycleCuttingNodes.size(); i++)
		{
			size_t cuttersbp = 0;
			for (size_t j = 0; j < cycleCuttingNodes[i].size(); j++)
			{
				cuttersbp += nodeEnd[cycleCuttingNodes[i][j]] - nodeStart[cycleCuttingNodes[i][j]];
			}
			std::cerr << i << ": id " << nodeIDs[i] << ", cutting nodes " << cycleCuttingNodes[i].size() << ", " << cuttersbp << "bp" << std::endl;
			totalCuttersbp += cuttersbp;
		}
		std::cerr << "total cut: " << totalCuttersbp << "bp (" << (double)totalCuttersbp / (double)nodeSequences.size() * 100 << "%)" << std::endl;
	}
}

void AlignmentGraph::saveCycleCut(std::string filename)
{
	std::ofstream file {filename};
	boost::archive::text_oarchive oa(file);
	oa << cycleCuttingNodes;
	oa << cycleCuttingNodePredecessor;
	oa << cycleCutPreviousCut;
}

bool AlignmentGraph::loadCycleCut(std::string filename)
{
	std::ifstream file {filename};
	if (!file.good()) return false;
	boost::archive::text_iarchive ia(file);
	ia >> cycleCuttingNodes;
	ia >> cycleCuttingNodePredecessor;
	ia >> cycleCutPreviousCut;
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

void AlignmentGraph::iterateOverCycleCuttingTreeRec(size_t cycleStart, size_t node, int sizeLeft, std::vector<size_t>& currentStack, std::function<void(const std::vector<size_t>&)> f)
{
	currentStack.push_back(node);
	auto nodeSize = nodeEnd[node]-nodeStart[node];
	if (node >= cycleStart && nodeSize < sizeLeft)
	{
		sizeLeft -= nodeSize;
		for (auto neighbor : inNeighbors[node])
		{
			iterateOverCycleCuttingTreeRec(cycleStart, neighbor, sizeLeft, currentStack, f);
		}
		currentStack.pop_back();
		return;
	}
	f(currentStack);
	currentStack.pop_back();
}

void AlignmentGraph::iterateOverCycleCuttingTree(size_t cycleStart, int sizeLeft, std::function<void(const std::vector<size_t>&)> f)
{
	std::vector<size_t> currentStack;
	iterateOverCycleCuttingTreeRec(cycleStart, cycleStart, sizeLeft, currentStack, f);
}

void AlignmentGraph::getCycleCuttersSupersequence(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut)
{
	std::unordered_map<size_t, std::vector<size_t>> supersequenceIndex;
	size_t supersequenceSize = 0;
	iterateOverCycleCuttingTree(cycleStart, sizeLeft, [&supersequenceIndex, &supersequenceSize](const std::vector<size_t>& currentStack) {
		size_t currentPos = 0;
		size_t stackProcessed = 0;
		for (size_t i = 0; i < currentStack.size(); i++)
		{
			auto list = supersequenceIndex[currentStack[i]];
			assert(i != 0 || supersequenceSize == 0 || (list.size() > 0 && list[0] == 0));
			if (list.size() == 0 || list.back() <= currentPos) break;
			if (i == 0)
			{
				assert(list[0] == 0);
				stackProcessed++;
				continue;
			}
			auto pos = std::upper_bound(list.begin(), list.end(), currentPos);
			assert(pos != list.end());
			currentPos = *pos;
			stackProcessed++;
		}

		for (size_t i = stackProcessed; i < currentStack.size(); i++)
		{
			supersequenceIndex[currentStack[i]].push_back(supersequenceSize);
			supersequenceSize++;
		}
	});

	std::vector<size_t> toobigSupersequence;
	toobigSupersequence.resize(supersequenceSize, std::numeric_limits<size_t>::max());
	for (auto pair : supersequenceIndex)
	{
		for (auto pos : pair.second)
		{
			assert(pos < toobigSupersequence.size());
			assert(toobigSupersequence[pos] == std::numeric_limits<size_t>::max());
			toobigSupersequence[pos] = pair.first;
		}
	}

	#ifndef NDEBUG
	for (size_t i = 0; i < toobigSupersequence.size(); i++)
	{
		assert(toobigSupersequence[i] != std::numeric_limits<size_t>::max());
	}
	#endif

	supersequence = toobigSupersequence;
	// std::vector<bool> used;
	// used.resize(toobigSupersequence.size(), false);
	// used[0] = true;
	// iterateOverCycleCuttingTree(cycleStart, sizeLeft, [&toobigSupersequence, &used](const std::vector<size_t>& currentStack) {
	// 	size_t offset = 0;
	// 	assert(toobigSupersequence[0] == currentStack[0]);
	// 	assert(toobigSupersequence.size() >= currentStack.size());
	// 	for (size_t i = 1; i < currentStack.size(); i++)
	// 	{
	// 		while (toobigSupersequence[i+offset] != currentStack[i])
	// 		{
	// 			offset++;
	// 			assert(i+offset < toobigSupersequence.size());
	// 		}
	// 		used[i+offset] = true;
	// 	}
	// });

	// supersequence.reserve(toobigSupersequence.size());
	// for (size_t i = 0; i < toobigSupersequence.size(); i++)
	// {
	// 	if (used[i]) supersequence.push_back(toobigSupersequence[i]);
	// }

	supersequencePredecessors.resize(supersequence.size());

	iterateOverCycleCuttingTree(cycleStart, sizeLeft, [&supersequence, &supersequencePredecessors](const std::vector<size_t>& currentStack) {
		size_t offset = 0;
		size_t lastIndex = 0;
		assert(supersequence[0] == currentStack[0]);
		assert(supersequence.size() >= currentStack.size());
		for (size_t i = 1; i < currentStack.size(); i++)
		{
			while (supersequence[i+offset] != currentStack[i])
			{
				offset++;
				assert(i+offset < supersequence.size());
			}
			supersequencePredecessors[lastIndex].insert(i+offset);
			lastIndex = i+offset;
		}
	});

	for (size_t i = 0; i < supersequence.size(); i++)
	{
		previousCut.push_back(supersequence[i] < cycleStart);
	}
}

void AlignmentGraph::calculateCycleCutters(size_t cycleStart, int wordSize)
{
	assert(cycleCuttingNodes.size() > cycleStart);
	assert(cycleCuttingNodePredecessor.size() > cycleStart);
	assert(cycleCuttingNodes[cycleStart].size() == 0);
	assert(cycleCuttingNodePredecessor[cycleStart].size() == 0);

	std::vector<size_t> nodes;
	std::vector<size_t> parents;
	getCycleCuttersSupersequence(cycleStart, 2*wordSize, cycleCuttingNodes[cycleStart], cycleCuttingNodePredecessor[cycleStart], cycleCutPreviousCut[cycleStart]);

	assert(cycleCuttingNodes[cycleStart].size() > 0);
	assert(cycleCuttingNodes[cycleStart][0] == cycleStart);
	assert(cycleCuttingNodes[cycleStart].size() == cycleCuttingNodePredecessor[cycleStart].size());
	assert(cycleCuttingNodePredecessor[cycleStart].size() == 1 || cycleCuttingNodePredecessor[cycleStart][0].size() > 0);
	assert(cycleCutPreviousCut[cycleStart].size() == cycleCuttingNodes[cycleStart].size());
}
