#include <iostream>
#include <limits>
#include <algorithm>
#include <boost/config.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include "TopologicalSort.h"
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

void CycleCutCalculation::iterateOverCycleCuttingCyclesRec(size_t cycleStart, size_t node, std::vector<size_t>& currentStack, std::unordered_set<size_t>& cyclic, int sizeLeft, std::function<void(const std::vector<size_t>&)> cycleFunction) const
{
	currentStack.push_back(node);
	auto nodeSize = graph.nodeEnd[node]-graph.nodeStart[node];
	bool recursed = false;
	if (node >= cycleStart && nodeSize < sizeLeft)
	{
		sizeLeft -= nodeSize;
		for (auto neighbor : graph.inNeighbors[node])
		{
			if (cyclic.count(neighbor) > 0)
			{
				iterateOverCycleCuttingCyclesRec(cycleStart, neighbor, currentStack, cyclic, sizeLeft, cycleFunction);
				recursed = true;
			}
		}
		if (recursed)
		{
			currentStack.pop_back();
			return;
		}
	}
	cycleFunction(currentStack);
	currentStack.pop_back();
}

void CycleCutCalculation::getReachableRec(size_t node, int sizeLeft, std::unordered_set<size_t>& reachable) const
{
	reachable.insert(node);
	auto nodeSize = graph.nodeEnd[node] - graph.nodeStart[node];
	if (nodeSize < sizeLeft)
	{
		sizeLeft -= nodeSize;
		for (auto neighbor : graph.inNeighbors[node])
		{
			if (reachable.count(neighbor) > 0) continue;
			getReachableRec(neighbor, sizeLeft, reachable);
		}
	}
}

void CycleCutCalculation::DFSSplitCyclicAndNoncyclic(size_t node, std::vector<size_t>& currentStack, std::unordered_set<size_t>& visited, std::unordered_set<size_t>& cyclic, const std::unordered_set<size_t>& reachable) const
{
	if (visited.count(node) > 0)
	{
		if (cyclic.count(node) > 0)
		{
			for (size_t j = 0; j < currentStack.size(); j++)
			{
				cyclic.insert(currentStack[j]);
			}
			return;
		}
		for (size_t i = 0; i < currentStack.size(); i++)
		{
			if (currentStack[i] == node)
			{
				for (size_t j = 0; j < currentStack.size(); j++)
				{
					cyclic.insert(currentStack[j]);
				}
				break;
			}
		}
		return;
	}
	currentStack.push_back(node);
	visited.insert(node);
	for (auto neighbor : graph.inNeighbors[node])
	{
		if (reachable.count(neighbor) == 0) continue;
		DFSSplitCyclicAndNoncyclic(neighbor, currentStack, visited, cyclic, reachable);
	}
	currentStack.pop_back();
}

void CycleCutCalculation::iterateOverCycleCuttingCycles(size_t cycleStart, int sizeLeft, std::function<void(const std::vector<size_t>&)> cycleFunction, std::function<void(const std::unordered_set<size_t>&)> uncycleFunction) const
{
	std::unordered_set<size_t> visited;
	std::unordered_set<size_t> cyclic;
	std::vector<size_t> currentStack;
	std::unordered_set<size_t> reachable;
	getReachableRec(cycleStart, sizeLeft, reachable);
	DFSSplitCyclicAndNoncyclic(cycleStart, currentStack, visited, cyclic, reachable);
	std::cerr << "visited " << visited.size() << ", cyclic " << cyclic.size() << std::endl;
	assert(currentStack.size() == 0);
	if (cyclic.count(cycleStart) > 0)
	{
		iterateOverCycleCuttingCyclesRec(cycleStart, cycleStart, currentStack, cyclic, sizeLeft, cycleFunction);
	}
	for (auto node : cyclic)
	{
		visited.erase(node);
	}
	assert(cyclic.size() == 0 || cyclic.count(cycleStart) == 1);
	assert(visited.count(cycleStart) != cyclic.count(cycleStart));
	uncycleFunction(visited);
}

std::vector<size_t> CycleCutCalculation::getSupersequenceIndexingAndPredecessors(size_t cycleStart, int sizeLeft, std::vector<std::set<size_t>>& predecessors) const
{
	std::unordered_map<size_t, std::vector<size_t>> supersequenceIndex;
	size_t supersequenceSize = 0;
	std::unordered_set<size_t> uncyclic;
	// iterateOverCycleCuttingTree(cycleStart, sizeLeft, [&supersequenceIndex, &supersequenceSize, &predecessors](const std::vector<size_t>& currentStack) {
	iterateOverCycleCuttingCycles(cycleStart, sizeLeft, [&supersequenceIndex, &supersequenceSize, &predecessors](const std::vector<size_t>& currentStack) {
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
	}
	,[&uncyclic](const std::unordered_set<size_t>& set) {
		uncyclic = set;
	}
	);

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

	appendNonCyclicParts(uncyclic, toobigSupersequence, predecessors);

	return toobigSupersequence;
}

void CycleCutCalculation::filterUnnecessaryCharacters(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors) const
{
	std::vector<bool> isPredecessor;
	isPredecessor.resize(supersequence.size(), false);
	for (size_t i = 0; i < supersequencePredecessors.size(); i++)
	{
		for (auto predecessor : supersequencePredecessors[i])
		{
			isPredecessor[predecessor] = true;
		}
	}
	std::set<size_t> removeIndices;
	for (size_t i = 1; i < isPredecessor.size(); i++)
	{
		if (!isPredecessor[i])
		{
			removeIndices.insert(i);
		}
	}
	
	if (removeIndices.size() == 0) return;

	std::vector<size_t> newIndex;
	newIndex.resize(supersequence.size(), 1);
	newIndex[0] = 0;
	for (auto x : removeIndices)
	{
		newIndex[x]--;
	}
	for (size_t i = 1; i < newIndex.size(); i++)
	{
		newIndex[i] = newIndex[i-1] + newIndex[i];
	}
#ifndef NDEBUG
	size_t lastErase = supersequence.size();
#endif
	for (auto iter = removeIndices.rbegin(); iter != removeIndices.rend(); ++iter)
	{
		auto x = *iter;
		assert(x < lastErase);
#ifndef NDEBUG
		lastErase = x;
#endif
		supersequence.erase(supersequence.begin()+x);
		supersequencePredecessors.erase(supersequencePredecessors.begin()+x);
	}
	for (size_t i = 0; i < supersequencePredecessors.size(); i++)
	{
		std::set<size_t> olds;
		std::swap(olds, supersequencePredecessors[i]);
		for (auto x : olds)
		{
			supersequencePredecessors[i].insert(newIndex[x]);
		}
	}
}

void CycleCutCalculation::getPredecessorsFromSupersequence(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const
{
	supersequencePredecessors.resize(supersequence.size());
	std::unordered_set<size_t> uncyclic;

	iterateOverCycleCuttingCycles(cycleStart, sizeLeft, [&supersequence, &supersequencePredecessors](const std::vector<size_t>& currentStack) {
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
	},
	[&uncyclic](const std::unordered_set<size_t>& set) {
		uncyclic = set;
	});

	filterUnnecessaryCharacters(cycleStart, sizeLeft, supersequence, supersequencePredecessors);

#ifndef NDEBUG
	std::vector<bool> isPredecessor;
	isPredecessor.resize(supersequence.size(), false);
	for (size_t i = 0; i < supersequencePredecessors.size(); i++)
	{
		for (auto predecessor : supersequencePredecessors[i])
		{
			isPredecessor[predecessor] = true;
		}
	}
	for (size_t i = 1; i < isPredecessor.size(); i++)
	{
		assert(isPredecessor[i]);
	}
#endif

	appendNonCyclicParts(uncyclic, supersequence, supersequencePredecessors);
}

void CycleCutCalculation::getPredecessorsFromSupersequenceOverEdgeCoveringPaths(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const
{
	supersequencePredecessors.resize(supersequence.size());
	std::unordered_set<size_t> uncyclic;

	iterateOverEdgeCoveringPaths(cycleStart, sizeLeft, [&supersequence, &supersequencePredecessors](const std::vector<size_t>& currentStack) {
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

	filterUnnecessaryCharacters(cycleStart, sizeLeft, supersequence, supersequencePredecessors);

#ifndef NDEBUG
	std::vector<bool> isPredecessor;
	isPredecessor.resize(supersequence.size(), false);
	for (size_t i = 0; i < supersequencePredecessors.size(); i++)
	{
		for (auto predecessor : supersequencePredecessors[i])
		{
			isPredecessor[predecessor] = true;
		}
	}
	for (size_t i = 1; i < isPredecessor.size(); i++)
	{
		assert(isPredecessor[i]);
	}
#endif
}

void CycleCutCalculation::appendNonCyclicParts(const std::unordered_set<size_t>& uncyclic, std::vector<size_t>& nodes, std::vector<std::set<size_t>>& predecessors) const
{
	std::unordered_set<size_t> cyclic;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		cyclic.insert(nodes[i]);
		assert(uncyclic.count(nodes[i]) == 0);
	}
	std::vector<std::vector<size_t>> topologicalSortGraph;
	std::unordered_map<size_t, size_t> topologicalIndex;
	std::vector<size_t> reverseTopologicalIndex;
	for (auto node : uncyclic)
	{
		topologicalIndex[node] = topologicalSortGraph.size();
		reverseTopologicalIndex.push_back(node);
		topologicalSortGraph.emplace_back();
	}
	for (size_t i = 0; i < topologicalSortGraph.size(); i++)
	{
		for (auto neighbor : graph.inNeighbors[reverseTopologicalIndex[i]])
		{
			if (uncyclic.count(neighbor) == 0) continue;
			topologicalSortGraph[i].push_back(topologicalIndex[neighbor]);
		}
	}
	auto order = topologicalSort(topologicalSortGraph);
	assert(order.size() == reverseTopologicalIndex.size());
#ifndef NDEBUG
	for (size_t i = 0; i < order.size(); i++)
	{
		auto node = reverseTopologicalIndex[order[i]];
		for (size_t j = 0; j < i; j++)
		{
			assert(graph.inNeighbors[node].count(reverseTopologicalIndex[order[j]]) == 0);
		}
	}
#endif
	std::unordered_map<size_t, std::vector<size_t>> nodeIndex;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		nodeIndex[nodes[i]].push_back(i);
	}
	for (size_t i = 0; i < order.size(); i++)
	{
		assert(order[i] < reverseTopologicalIndex.size());
		nodes.push_back(reverseTopologicalIndex[order[i]]);
		nodeIndex[nodes.back()].push_back(nodes.size()-1);
		predecessors.emplace_back();
#ifndef NDEBUG
		bool inserted = false;
#endif
		for (auto neighbor : graph.outNeighbors[nodes.back()])
		{
			assert((cyclic.count(neighbor) == 0 && uncyclic.count(neighbor) == 0) || nodeIndex[neighbor].size() > 0);
			for (auto index : nodeIndex[neighbor])
			{
				assert(index < nodes.size()-1);
				assert(graph.inNeighbors[nodes[index]].count(nodes.back()) == 1);
				predecessors[index].insert(nodes.size()-1);
#ifndef NDEBUG
				inserted = true;
#endif
			}
		}
		assert(inserted || nodes.size() == 1);
	}
#ifndef NDEBUG
	for (size_t i = 0; i < nodes.size(); i++)
	{
		for (auto predecessor : predecessors[i])
		{
			assert(predecessor > i);
			assert(graph.inNeighbors[nodes[i]].count(nodes[predecessor]) == 1);
		}
		if (uncyclic.count(nodes[i]) == 1)
		{
			for (size_t j = 0; j < i; j++)
			{
				assert(predecessors[j].count(i) == graph.inNeighbors[nodes[j]].count(nodes[i]));
			}
		}
	}
	for (size_t i = 1; i < nodes.size(); i++)
	{
		bool foundOne = false;
		for (size_t j = 0; j < i; j++)
		{
			if (predecessors[j].count(i) == 1)
			{
				foundOne = true;
				break;
			}
		}
		assert(foundOne);
	}
#endif
}

std::vector<size_t> CycleCutCalculation::getCycleCuttersTheDumbWay(size_t cycleStart, int sizeLeft) const
{
	std::vector<size_t> supersequence;
	iterateOverCycleCuttingCycles(cycleStart, sizeLeft, [&supersequence](const std::vector<size_t>& currentStack) {
		supersequence.insert(supersequence.end(), currentStack.begin(), currentStack.end());
	},
	[](const std::unordered_set<size_t>& set) {}
	);
	return supersequence;
}

std::vector<size_t> getPairwiseSupersequenceByAlignment(const std::vector<size_t>& supersequence, const std::vector<size_t>& currentStack)
{
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
	return newSupersequence;
}

std::vector<size_t> CycleCutCalculation::getCycleCuttersSupersequenceOverEdgeCoveringPaths(size_t cycleStart, int sizeLeft) const
{
	std::vector<size_t> supersequence;
	iterateOverEdgeCoveringPaths(cycleStart, sizeLeft, [&supersequence](const std::vector<size_t>& currentStack) {
		if (supersequence.size() == 0)
		{
			assert(currentStack.size() > 0);
			supersequence = currentStack;
			return;
		}
		supersequence = getPairwiseSupersequenceByAlignment(supersequence, currentStack);
	});
	return supersequence;
}

std::vector<size_t> CycleCutCalculation::getCycleCuttersSupersequence(size_t cycleStart, int sizeLeft) const
{
	std::vector<size_t> supersequence;
	iterateOverCycleCuttingCycles(cycleStart, sizeLeft, [&supersequence](const std::vector<size_t>& currentStack) {
		if (supersequence.size() == 0)
		{
			assert(currentStack.size() > 0);
			supersequence = currentStack;
			return;
		}
		supersequence = getPairwiseSupersequenceByAlignment(supersequence, currentStack);
	},
	[](const std::unordered_set<size_t>& set) {}
	);
	return supersequence;
}

std::vector<size_t> CycleCutCalculation::getCycleCuttersOrder(size_t cycleStart, int sizeLeft, std::vector<std::set<size_t>>& predecessors) const
{
	std::vector<size_t> supersequence;
	std::map<std::pair<size_t, size_t>, size_t> positionInSupersequence;
	std::vector<std::set<size_t>> nodes;
	nodes.resize(sizeLeft);
	std::vector<std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t>>> edges;
	std::unordered_set<size_t> reachable;
	getReachableRec(cycleStart, sizeLeft, reachable);
	nodes[0].insert(cycleStart);
	for (size_t i = 0; i < sizeLeft; i++)
	{
		for (auto node : nodes[i])
		{
			positionInSupersequence[std::make_pair(i, node)] = supersequence.size();
			supersequence.push_back(node);
			auto nodeSize = graph.nodeEnd[node] - graph.nodeStart[node];
			if (i + nodeSize >= sizeLeft) continue;
			for (auto neighbor : graph.inNeighbors[node])
			{
				nodes[i + nodeSize].insert(neighbor);
				edges.emplace_back(std::make_pair(i, node), std::make_pair(i + nodeSize, neighbor));
			}
		}
	}
	predecessors.resize(supersequence.size());
	for (auto edge : edges)
	{
		auto from = edge.first;
		auto to = edge.second;
		predecessors[positionInSupersequence[from]].insert(positionInSupersequence[to]);
	}
	return supersequence;
}

std::map<std::pair<size_t, size_t>, size_t> findFeasibleFlow(const std::vector<size_t>& supersequence, const std::vector<std::set<size_t>>& predecessors)
{
	std::map<std::pair<size_t, size_t>, size_t> result;
	std::map<size_t, std::vector<std::pair<size_t, size_t>>> pathBack;
	std::map<size_t, std::vector<std::pair<size_t, size_t>>> pathForward;
	std::vector<std::set<size_t>> successors;
	successors.resize(supersequence.size());
	for (size_t i = 0; i < supersequence.size(); i++)
	{
		for (auto predecessor : predecessors[i])
		{
			pathBack[predecessor] = pathBack[i];
			pathBack[predecessor].emplace_back(i, predecessor);
			successors[predecessor].insert(i);
		}
	}
	for (size_t i = supersequence.size()-1; i < supersequence.size(); i--)
	{
		for (auto successor : successors[i])
		{
			pathForward[successor] = pathForward[i];
			pathForward[successor].emplace_back(successor, i);
		}
	}
	for (size_t i = 0; i < supersequence.size(); i++)
	{
		for (auto predecessor : predecessors[i])
		{
			auto thispart = std::make_pair(i, predecessor);
			if (result[thispart] > 0) continue;
			result[thispart]++;;
			for (auto part : pathBack[i])
			{
				result[part]++;
			}
			for (auto part : pathForward[predecessor])
			{
				result[part]++;
			}
		}
	}
	return result;
}

void getOneFlowPath(const std::vector<size_t>& supersequence, const std::vector<std::set<size_t>>& predecessors, const std::map<std::pair<size_t, size_t>, size_t>& flows, size_t node, std::vector<size_t>& result)
{
	result.push_back(node);
	for (auto predecessor : predecessors[node])
	{
		assert(predecessor > node);
		if (flows.at(std::make_pair(node, predecessor)) > 0)
		{
			getOneFlowPath(supersequence, predecessors, flows, predecessor, result);
			return;
		}
	}
}

std::vector<std::vector<size_t>> findFlowPaths(const std::vector<size_t>& supersequence, const std::vector<std::set<size_t>>& predecessors, std::map<std::pair<size_t, size_t>, size_t> flows)
{
	assert(supersequence.size() > 0);
	std::vector<std::vector<size_t>> result;
	if (supersequence.size() == 1)
	{
		result.emplace_back();
		result.back().push_back(0);
		return result;
	}
	while (true)
	{
		std::vector<size_t> path;
		getOneFlowPath(supersequence, predecessors, flows, 0, path);
		assert(path.size() > 0);
		if (path.size() == 1) break;
		for (size_t i = 1; i < path.size(); i++)
		{
			assert(flows[std::make_pair(path[i-1], path[i])] >= 1);
			flows[std::make_pair(path[i-1], path[i])] -= 1;
		}
		std::vector<size_t> pathID;
		pathID.reserve(path.size());
		for (size_t i = 0; i < path.size(); i++)
		{
			pathID.push_back(supersequence[path[i]]);
		}
		result.push_back(pathID);
	}
	assert(result.size() > 0);
#ifndef NDEBUG
	for (auto x : flows)
	{
		assert(x.second == 0);
	}
#endif
	return result;
}

void CycleCutCalculation::iterateOverEdgeCoveringPaths(size_t cycleStart, int sizeLeft, std::function<void(const std::vector<size_t>&)> function) const
{
	std::vector<std::set<size_t>> predecessors;
	auto supersequence = getCycleCuttersOrder(cycleStart, sizeLeft, predecessors);

	//https://stackoverflow.com/questions/18598399/what-algorithm-should-i-use-to-find-the-minimum-flow-on-a-digraph-where-there-ar
	//http://www.boost.org/doc/libs/1_41_0/libs/graph/example/max_flow.cpp
	//http://www.boost.org/doc/libs/1_41_0/boost/graph/read_dimacs.hpp
	typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> Traits;
	typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, 
		boost::property<boost::vertex_name_t, std::string>,
		boost::property<boost::edge_capacity_t, long,
		boost::property<boost::edge_residual_capacity_t, long,
		boost::property<boost::edge_reverse_t, Traits::edge_descriptor>>>
	> Graph;
	Graph g;
	boost::property_map<Graph, boost::edge_capacity_t>::type 
		capacity = get(boost::edge_capacity, g);
	boost::property_map<Graph, boost::edge_reverse_t>::type 
		rev = get(boost::edge_reverse, g);
	boost::property_map<Graph, boost::edge_residual_capacity_t>::type 
		residual_capacity = get(boost::edge_residual_capacity, g);
	Traits::vertex_descriptor s, t;
	std::vector<Traits::vertex_descriptor> verts;
	for (size_t i = 0; i < supersequence.size(); i++)
	{
		verts.push_back(boost::add_vertex(g));
	}
	auto src = boost::add_vertex(g);
	auto sink = boost::add_vertex(g);
	auto startFlow = findFeasibleFlow(supersequence, predecessors);
	size_t flowBeforeReduction = 0;
	for (auto predecessor : predecessors[0])
	{
		flowBeforeReduction += startFlow[std::make_pair(0, predecessor)];
	}
	std::cerr << "flow before reduction: " << flowBeforeReduction << std::endl;
	auto flowPathsBefore = findFlowPaths(supersequence, predecessors, startFlow).size();
	std::cerr << "flow paths before: " << flowPathsBefore << std::endl;
	{
		Traits::edge_descriptor e1, e2;
		bool in1, in2;
		boost::tie(e1, in1) = add_edge(src, verts[0], g);
		boost::tie(e2, in2) = add_edge(verts[0], src, g);
		if (!in1 || !in2) {
			std::cerr << "unable to add edge (src, 0)" << std::endl;
			std::abort();
		}
		capacity[e1] = supersequence.size() * supersequence.size();
		capacity[e2] = 0;
		rev[e1] = e2;
		rev[e2] = e1;
	}
	for (size_t i = 0; i < supersequence.size(); i++)
	{
		for (auto predecessor : predecessors[i])
		{
			auto tail = i;
			auto head = predecessor;
			Traits::edge_descriptor e1, e2;
			bool in1, in2;
			boost::tie(e1, in1) = add_edge(verts[tail], verts[head], g);
			boost::tie(e2, in2) = add_edge(verts[head], verts[tail], g);
			if (!in1 || !in2) {
				std::cerr << "unable to add edge (" << head << "," << tail << ")" << std::endl;
				std::abort();
			}
			assert(startFlow[std::make_pair(i, predecessor)] >= 1);
			capacity[e1] = startFlow[std::make_pair(i, predecessor)] - 1;
			capacity[e2] = 0;
			rev[e1] = e2;
			rev[e2] = e1;
		}
		if (predecessors[i].size() == 0)
		{
			Traits::edge_descriptor e1, e2;
			bool in1, in2;
			boost::tie(e1, in1) = add_edge(verts[i], sink, g);
			boost::tie(e2, in2) = add_edge(sink, verts[i], g);
			if (!in1 || !in2) {
				std::cerr << "unable to add edge (" << i << ", sink)" << std::endl;
				std::abort();
			}
			capacity[e1] = supersequence.size() * supersequence.size();
			capacity[e2] = 0;
			rev[e1] = e2;
			rev[e2] = e1;
		}
	}

	// auto boostFlow = edmonds_karp_max_flow(g, s, t);
	// auto boostFlow = boykov_kolmogorov_max_flow(g, s, t);
	// auto boostFlow = push_relabel_max_flow(g, s, t);

	// boost::graph_traits<Graph>::vertex_iterator u_iter, u_end;
	// boost::graph_traits<Graph>::out_edge_iterator ei, e_end;
	// for (boost::tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter)
	// {
	// 	for (boost::tie(ei, e_end) = out_edges(*u_iter, g); ei != e_end; ++ei)
	// 	{
	// 		if (capacity[*ei] > 0)
	// 		{
	// 			auto from = *u_iter;
	// 			auto to = target(*ei, g);
	// 			if (from >= supersequence.size() || to >= supersequence.size()) continue;
	// 			startFlow[std::make_pair(from, to)] -= (capacity[*ei] - residual_capacity[*ei]);
	// 		}
	// 	}
	// }
	// size_t flow = 0;
	// for (auto predecessor : predecessors[0])
	// {
	// 	flow += startFlow[std::make_pair(0, predecessor)];
	// }
	// auto flowPathsAfter = findFlowPaths(supersequence, predecessors, startFlow).size();
	// std::cerr << "flow paths after: " << flowPathsAfter << std::endl;
	// std::cerr << "flow: " << flow << std::endl;
	// std::cerr << "boost flow: " << boostFlow << std::endl;
	// assert(flow == flowBeforeReduction - boostFlow);

	for (auto path : findFlowPaths(supersequence, predecessors, startFlow))
	{
		function(path);
	}
}

size_t getNumberOfPaths(const std::vector<size_t>& supersequence, const std::vector<std::set<size_t>>& predecessors)
{
	std::vector<size_t> pathsEndingHere;
	pathsEndingHere.resize(supersequence.size(), 0);
	pathsEndingHere[0] = 1;
	size_t total = 0;
	for (size_t i = 0; i < supersequence.size(); i++)
	{
		for (auto predecessor : predecessors[i])
		{
			pathsEndingHere[predecessor] += pathsEndingHere[i];
		}
		if (predecessors[i].size() == 0)
		{
			total += pathsEndingHere[i];
		}
	}
	return total;
}

void CycleCutCalculation::getCycleCuttersByOrder(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const
{
	supersequence = getCycleCuttersOrder(cycleStart, sizeLeft, supersequencePredecessors);
	// getPredecessorsFromSupersequence(cycleStart, sizeLeft, supersequence, supersequencePredecessors, previousCut);
	for (size_t i = 0; i < supersequence.size(); i++)
	{
		previousCut.push_back(supersequence[i] < cycleStart);
	}
}

void CycleCutCalculation::getCycleCuttersBySupersequenceOverEdgeCoveringPaths(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const
{
	supersequence = getCycleCuttersSupersequenceOverEdgeCoveringPaths(cycleStart, sizeLeft);
	getPredecessorsFromSupersequenceOverEdgeCoveringPaths(cycleStart, sizeLeft, supersequence, supersequencePredecessors, previousCut);
	for (size_t i = 0; i < supersequence.size(); i++)
	{
		previousCut.push_back(supersequence[i] < cycleStart);
	}
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

void CycleCutCalculation::getCycleCuttersByDumbWay(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const
{
	supersequence = getCycleCuttersTheDumbWay(cycleStart, sizeLeft);
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

void CycleCutCalculation::PrintTreeSize(size_t startNode, int wordSize) const
{
	size_t total = 0;
	iterateOverCycleCuttingCycles(startNode, wordSize*2, [&total](const std::vector<size_t>& currentStack) {
		total++;
		if (total % 1000000 == 0)
		{
			std::cerr << total << " ";
		}
	}
	,[&total](const std::unordered_set<size_t>& uncyclic) {
		std::cerr << " uncyclic " << uncyclic.size() << std::endl;
		total += uncyclic.size();
	}
	);
	std::cerr << "total treesize: " << total << std::endl;
}

AlignmentGraph::CycleCut CycleCutCalculation::GetCycleCutByIndex(size_t startNode, int wordSize) const
{
	AlignmentGraph::CycleCut result;
	getCycleCuttersByIndex(startNode, wordSize*2, result.nodes, result.predecessors, result.previousCut);
	std::cerr << "paths: " << getNumberOfPaths(result.nodes, result.predecessors) << std::endl;
	return result;
}

AlignmentGraph::CycleCut CycleCutCalculation::GetCycleCutByOrder(size_t startNode, int wordSize) const
{
	AlignmentGraph::CycleCut result;
	getCycleCuttersByOrder(startNode, wordSize*2, result.nodes, result.predecessors, result.previousCut);
	std::cerr << "paths: " << getNumberOfPaths(result.nodes, result.predecessors) << std::endl;
	return result;
}

AlignmentGraph::CycleCut CycleCutCalculation::GetCycleCutBySupersequenceOverEdgeCoveringPaths(size_t startNode, int wordSize) const
{
	AlignmentGraph::CycleCut result;
	getCycleCuttersBySupersequenceOverEdgeCoveringPaths(startNode, wordSize*2, result.nodes, result.predecessors, result.previousCut);
	std::cerr << "paths: " << getNumberOfPaths(result.nodes, result.predecessors) << std::endl;
	return result;
}

AlignmentGraph::CycleCut CycleCutCalculation::GetCycleCutBySupersequence(size_t startNode, int wordSize) const
{
	AlignmentGraph::CycleCut result;
	getCycleCuttersBySupersequence(startNode, wordSize*2, result.nodes, result.predecessors, result.previousCut);
	std::cerr << "paths: " << getNumberOfPaths(result.nodes, result.predecessors) << std::endl;
	return result;
}

AlignmentGraph::CycleCut CycleCutCalculation::GetCycleCutByDumbWay(size_t startNode, int wordSize) const
{
	AlignmentGraph::CycleCut result;
	getCycleCuttersByDumbWay(startNode, wordSize*2, result.nodes, result.predecessors, result.previousCut);
	std::cerr << "paths: " << getNumberOfPaths(result.nodes, result.predecessors) << std::endl;
	return result;
}


