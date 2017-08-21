#include <limits>
#include <algorithm>
#include "AlignmentGraph.h"
#include "CycleCutCalculation.h"
#include "2dArray.h"

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

std::vector<size_t> CycleCutCalculation::getCycleCuttersSupersequence(size_t cycleStart, int sizeLeft) const
{
	std::vector<size_t> supersequence;
	iterateOverCycleCuttingTree(cycleStart, sizeLeft, [&supersequence](const std::vector<size_t>& currentStack) {
		if (supersequence.size() == 0)
		{
			assert(currentStack.size() > 0);
			supersequence = currentStack;
			return;
		}
		assert(supersequence.size() > 0);
		assert(currentStack.size() > 0);

		Array2D<size_t, false> scores { supersequence.size(), currentStack.size(), std::numeric_limits<size_t>::max() };
		Array2D<char, false> backtrace { supersequence.size(), currentStack.size(), '-' };

		for (size_t i = 0; i < supersequence.size(); i++)
		{
			scores(i, 0) = 0;
			backtrace(i, 0) = 'L';
		}
		for (size_t j = 0; j < currentStack.size(); j++)
		{
			scores(0, j) = j;
			backtrace(0, j) = 'U';
		}
		for (size_t i = 1; i < supersequence.size(); i++)
		{
			for (size_t j = 1; j < currentStack.size(); j++)
			{
				size_t value = scores(i-1, j);
				char source = 'L';
				if (scores(i, j-1)+1 < value)
				{
					value = scores(i, j-1)+1;
					source = 'U';
				}
				if (supersequence[i] == currentStack[j] && scores(i-1, j-1) < value)
				{
					value = scores(i-1, j-1);
					source = 'D';
				}
				scores(i, j) = value;
				backtrace(i, j) = source;
			}
		}
		std::vector<size_t> newSupersequence;
		size_t i = supersequence.size()-1;
		size_t j = currentStack.size()-1;
		while (i != 0 || j != 0)
		{
			assert(i < supersequence.size());
			assert(j < currentStack.size());
			char dir = backtrace(i, j);
			switch(dir)
			{
				case 'L':
					newSupersequence.push_back(supersequence[i]);
					i--;
					break;
				case 'U':
					newSupersequence.push_back(currentStack[j]);
					j--;
					break;
				case 'D':
					newSupersequence.push_back(supersequence[i]);
					assert(supersequence[i] == currentStack[j]);
					i--;
					j--;
					break;
				default:
					assert(false);
			}
		}
		assert(supersequence[0] == currentStack[0]);
		newSupersequence.push_back(supersequence[0]);
		std::reverse(newSupersequence.begin(), newSupersequence.end());
		assert(newSupersequence.size() >= supersequence.size());
		supersequence = std::move(newSupersequence);
	});
	return supersequence;
}

void CycleCutCalculation::getCycleCuttersBySupersequence(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const
{
	supersequence = getCycleCuttersSupersequence(cycleStart, sizeLeft);
	getPredecessorsFromSupersequence(cycleStart, sizeLeft, supersequence, supersequencePredecessors, previousCut);
	for (size_t i = 0; i < supersequence.size(); i++)
	{
		previousCut.push_back(supersequence[i] < cycleStart);
	}
}

void CycleCutCalculation::getCycleCuttersByIndex(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const
{
	supersequence = getSupersequenceIndexingAndPredecessors(cycleStart, sizeLeft, supersequencePredecessors);
	// getPredecessorsFromSupersequence(cycleStart, sizeLeft, supersequence, supersequencePredecessors, previousCut);
	for (size_t i = 0; i < supersequence.size(); i++)
	{
		previousCut.push_back(supersequence[i] < cycleStart);
	}
}

AlignmentGraph::CycleCut CycleCutCalculation::GetCycleCutByIndex(size_t startNode, int wordSize) const
{
	AlignmentGraph::CycleCut result;
	getCycleCuttersByIndex(startNode, wordSize*2, result.nodes, result.predecessors, result.previousCut);
	return result;
}

AlignmentGraph::CycleCut CycleCutCalculation::GetCycleCutBySupersequence(size_t startNode, int wordSize) const
{
	AlignmentGraph::CycleCut result;
	getCycleCuttersBySupersequence(startNode, wordSize*2, result.nodes, result.predecessors, result.previousCut);
	return result;
}

