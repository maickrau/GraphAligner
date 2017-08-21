#include <limits>
#include <algorithm>
#include "AlignmentGraph.h"
#include "CycleCutCalculation.h"

CycleCutCalculation::CycleCutCalculation(const AlignmentGraph& graph) : graph(graph)
{
}

void CycleCutCalculation::iterateOverCycleCuttingTreeRec(size_t cycleStart, size_t node, int sizeLeft, std::vector<size_t>& currentStack, std::function<void(const std::vector<size_t>&)> f) const
{
	currentStack.push_back(node);
	auto nodeSize = graph.nodeEnd[node]-graph.nodeStart[node];
	if (node >= cycleStart && nodeSize < sizeLeft)
	{
		sizeLeft -= nodeSize;
		for (auto neighbor : graph.inNeighbors[node])
		{
			iterateOverCycleCuttingTreeRec(cycleStart, neighbor, sizeLeft, currentStack, f);
		}
		currentStack.pop_back();
		return;
	}
	f(currentStack);
	currentStack.pop_back();
}

void CycleCutCalculation::iterateOverCycleCuttingTree(size_t cycleStart, int sizeLeft, std::function<void(const std::vector<size_t>&)> f) const
{
	std::vector<size_t> currentStack;
	iterateOverCycleCuttingTreeRec(cycleStart, cycleStart, sizeLeft, currentStack, f);
}

std::vector<size_t> CycleCutCalculation::getSupersequenceIndexingAndPredecessors(size_t cycleStart, int sizeLeft, std::vector<std::set<size_t>>& predecessors) const
{
	std::unordered_map<size_t, std::vector<size_t>> supersequenceIndex;
	size_t supersequenceSize = 0;
	iterateOverCycleCuttingTree(cycleStart, sizeLeft, [&supersequenceIndex, &supersequenceSize, &predecessors](const std::vector<size_t>& currentStack) {
		assert(currentStack.size() > 0);
		size_t currentPos = 0;
		if (supersequenceSize == 0)
		{
			supersequenceIndex[currentStack[0]].push_back(supersequenceSize);
			predecessors.emplace_back();
			supersequenceSize++;
		}
		else
		{
			auto& list = supersequenceIndex[currentStack[0]];
			assert(list.size() > 0 && list[0] == 0);
		}
		size_t stackProcessed = 1;
		for (size_t i = 1; i < currentStack.size();)
		{
			auto& list = supersequenceIndex[currentStack[i]];
			if (list.size() == 0 || list.back() <= currentPos) break;
			auto pos = std::upper_bound(list.begin(), list.end(), currentPos);
			assert(pos != list.end());
			do
			{
				predecessors[currentPos].insert(*pos);
				currentPos = *pos;
				stackProcessed++;
				++pos;
				i++;
			} while (pos != list.end() && i < currentStack.size() && currentStack[i-1] == currentStack[i]);
		}

		for (size_t i = stackProcessed; i < currentStack.size(); i++)
		{
			predecessors[currentPos].insert(supersequenceSize);
			currentPos = supersequenceSize;
			supersequenceIndex[currentStack[i]].push_back(supersequenceSize);
			supersequenceSize++;
			predecessors.emplace_back();
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
	assert(predecessors.size() == toobigSupersequence.size());

	#ifndef NDEBUG
	for (size_t i = 0; i < toobigSupersequence.size(); i++)
	{
		assert(toobigSupersequence[i] != std::numeric_limits<size_t>::max());
	}
	#endif

	return toobigSupersequence;
}

void CycleCutCalculation::getPredecessorsFromSupersequence(size_t cycleStart, int sizeLeft, const std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const
{
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
}

void CycleCutCalculation::getCycleCuttersSupersequence(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const
{
	supersequence = getSupersequenceIndexingAndPredecessors(cycleStart, sizeLeft, supersequencePredecessors);
	// getPredecessorsFromSupersequence(cycleStart, sizeLeft, supersequence, supersequencePredecessors, previousCut);
	for (size_t i = 0; i < supersequence.size(); i++)
	{
		previousCut.push_back(supersequence[i] < cycleStart);
	}
}

AlignmentGraph::CycleCut CycleCutCalculation::GetCycleCut(size_t startNode, int wordSize) const
{
	AlignmentGraph::CycleCut result;
	getCycleCuttersSupersequence(startNode, wordSize*2, result.nodes, result.predecessors, result.previousCut);
	return result;
}

