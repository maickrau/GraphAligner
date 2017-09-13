#include <queue>
#include <iostream>
#include <limits>
#include <algorithm>
#include <boost/config.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include "TopologicalSort.h"
#include "AlignmentGraph.h"
#include "CycleCutCalculation.h"
#include "2dArray.h"

CycleCutCalculation::CycleCutCalculation(const AlignmentGraph& graph) : graph(graph)
{
}

class PriorityNode
{
public:
	PriorityNode(size_t n, size_t d) : nodeIndex(n), distance(d) {};
	size_t nodeIndex;
	size_t distance;
	bool operator<(const PriorityNode& other) const
	{
		return distance < other.distance;
	}
	bool operator>(const PriorityNode& other) const
	{
		return distance > other.distance;
	}
};

std::unordered_set<size_t> CycleCutCalculation::getReachable(size_t cycleStart, size_t sizeLeft) const
{
	std::unordered_set<size_t> result;
	std::priority_queue<PriorityNode, std::vector<PriorityNode>, std::greater<PriorityNode>> queue;
	queue.emplace(cycleStart, 0);
	while (queue.size() > 0)
	{
		auto top = queue.top();
		queue.pop();
		if (top.distance > sizeLeft) continue;
		if (result.count(top.nodeIndex) == 1) continue;
		result.insert(top.nodeIndex);
		auto nextDistance = top.distance + graph.nodeEnd[top.nodeIndex] - graph.nodeStart[top.nodeIndex];
		assert(nextDistance > top.distance);
		for (auto neighbor : graph.inNeighbors[top.nodeIndex])
		{
			queue.emplace(neighbor, nextDistance);
		}
	}
	return result;
}

void CycleCutCalculation::splitCyclicAndNoncyclicRec(std::vector<size_t>& stack, size_t currentNode, const std::unordered_set<size_t>& reachable, std::unordered_set<size_t>& visited, std::unordered_set<size_t>& cyclic) const
{
	for (size_t i = 0; i < stack.size(); i++)
	{
		if (stack[i] == currentNode)
		{
			for (auto node : stack)
			{
				cyclic.insert(node);
			}
			visited.insert(currentNode);
			return;
		}
	}
	if (cyclic.count(currentNode) > 0)
	{
		for (auto node : stack)
		{
			cyclic.insert(node);
		}
		visited.insert(currentNode);
		return;
	}
	if (visited.count(currentNode) > 0) return;
	visited.insert(currentNode);
	stack.push_back(currentNode);
	for (auto neighbor : graph.inNeighbors[currentNode])
	{
		if (reachable.count(neighbor) == 0) continue;
		splitCyclicAndNoncyclicRec(stack, neighbor, reachable, visited, cyclic);
	}
	stack.pop_back();
}

std::pair<std::unordered_set<size_t>, std::unordered_set<size_t>> CycleCutCalculation::splitCyclicAndNoncyclic(size_t cycleStart, int sizeLeft) const
{
	const auto reachable = getReachable(cycleStart, sizeLeft);
	// return std::make_pair(reachable, std::unordered_set<size_t>{});
	std::unordered_set<size_t> cyclic;
	std::unordered_set<size_t> visited;
	std::vector<size_t> stack;
	splitCyclicAndNoncyclicRec(stack, cycleStart, reachable, visited, cyclic);
	assert(visited.size() == reachable.size());
	for (auto node : cyclic)
	{
		visited.erase(node);
	}
	assert(visited.size() + cyclic.size() == reachable.size());
	assert(cyclic.count(cycleStart) != visited.count(cycleStart));
	assert((cyclic.size() == 0 && cyclic.count(cycleStart) == 0) || (cyclic.size() > 0 && cyclic.count(cycleStart) == 1));
	return std::make_pair(cyclic, visited);
}

std::vector<size_t> CycleCutCalculation::getCycleCuttersOrder(size_t cycleStart, int sizeLeft, std::vector<std::set<size_t>>& predecessors) const
{
	auto cyclicUncyclic = splitCyclicAndNoncyclic(cycleStart, sizeLeft);
	auto& cyclic = cyclicUncyclic.first;
	auto& uncyclicset = cyclicUncyclic.second;
	assert((cyclic.size() == 0 && cyclic.count(cycleStart) == 0) || (cyclic.size() > 0 && cyclic.count(cycleStart) == 1));
	assert((cyclic.size() == 0 && uncyclicset.count(cycleStart) == 1) || (cyclic.size() > 0 && uncyclicset.count(cycleStart) == 0));
	std::vector<size_t> supersequence;
	std::vector<std::unordered_map<size_t, size_t>> positionInSupersequence;
	std::unordered_map<size_t, std::vector<size_t>> existingIndexesForNode;
	positionInSupersequence.resize(sizeLeft);
	std::vector<std::unordered_set<size_t>> nodes;
	nodes.resize(sizeLeft);
	std::vector<std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t>>> edges;
	if (cyclic.size() > 0)
	{
		assert(cyclic.count(cycleStart) > 0);
		nodes[0].insert(cycleStart);
		for (size_t i = 0; i < sizeLeft; i++)
		{
			for (auto node : nodes[i])
			{
				assert(cyclic.count(node) > 0);
				positionInSupersequence[i][node] = supersequence.size();
				existingIndexesForNode[node].push_back(i);
				supersequence.push_back(node);
				auto nodeSize = graph.nodeEnd[node] - graph.nodeStart[node];
				if (i + nodeSize >= sizeLeft) continue;
				for (auto neighbor : graph.inNeighbors[node])
				{
					if (cyclic.count(neighbor) == 0) continue;
					bool found = false;
					for (size_t j = i+1; j <= i+nodeSize; j++)
					{
						if (nodes[j].count(neighbor) > 0)
						{
							edges.emplace_back(std::make_pair(i, node), std::make_pair(j, neighbor));
							found = true;
							break;
						}
					}
					if (!found)
					{
						nodes[i + nodeSize].insert(neighbor);
						edges.emplace_back(std::make_pair(i, node), std::make_pair(i + nodeSize, neighbor));
					}
				}
			}
		}
	}
	if (uncyclicset.size() > 0)
	{
		std::vector<size_t> uncyclic;
		{
			std::unordered_map<size_t, size_t> nodeIndex;
			std::vector<size_t> topologicalSortNodes;
			std::vector<std::vector<size_t>> topologicalSortNeighbors;
			for (auto node : uncyclicset)
			{
				nodeIndex[node] = topologicalSortNodes.size();
				topologicalSortNodes.push_back(node);
			}
			topologicalSortNeighbors.resize(topologicalSortNodes.size());
			for (size_t i = 0; i < topologicalSortNodes.size(); i++)
			{
				for (auto neighbor : graph.inNeighbors[topologicalSortNodes[i]])
				{
					assert(neighbor != topologicalSortNodes[i]);
					if (uncyclicset.count(neighbor) == 0) continue;
					topologicalSortNeighbors[i].push_back(nodeIndex[neighbor]);
				}
			}
			auto order = topologicalSort(topologicalSortNeighbors);
			for (size_t i = 0; i < order.size(); i++)
			{
				uncyclic.push_back(topologicalSortNodes[order[i]]);
			}
		}
		for (auto node : uncyclic)
		{
			assert(nodes[0].count(node) == 0);
			assert(existingIndexesForNode.count(node) == 0);
			assert(positionInSupersequence[0].count(node) == 0);
			nodes[0].insert(node);
			existingIndexesForNode[node].push_back(0);
			positionInSupersequence[0][node] = supersequence.size();
			supersequence.push_back(node);
		}
		for (auto node : uncyclic)
		{
			for (auto neighbor : graph.outNeighbors[node])
			{
				if (cyclic.count(neighbor) == 0 && uncyclicset.count(neighbor) == 0) continue;
				assert(existingIndexesForNode.count(neighbor) > 0);
				for (auto pos : existingIndexesForNode[neighbor])
				{
					edges.emplace_back(std::make_pair(pos, neighbor), std::make_pair(0, node));
					assert(positionInSupersequence[pos].count(neighbor) == 1);
					assert(positionInSupersequence[pos][neighbor] < positionInSupersequence[0][node]);
				}
			}
#ifndef NDEBUG
			for (auto neighbor : graph.inNeighbors[node])
			{
				assert(cyclic.count(neighbor) == 0);
			}
#endif
		}
	}
	predecessors.resize(supersequence.size());
	for (auto edge : edges)
	{
		auto from = edge.first;
		auto to = edge.second;
		assert(positionInSupersequence[from.first][from.second] < positionInSupersequence[to.first][to.second]);
		predecessors[positionInSupersequence[from.first][from.second]].insert(positionInSupersequence[to.first][to.second]);
	}
	return supersequence;
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

void CycleCutCalculation::getPredecessorsFromSupersequenceOverEdgeCoveringPaths(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut, const std::vector<std::vector<size_t>>& paths) const
{
	supersequencePredecessors.resize(supersequence.size());
	std::unordered_set<size_t> uncyclic;

	for (const auto& path : paths)
	{
		size_t offset = 0;
		size_t lastIndex = 0;
		assert(supersequence[0] == path[0]);
		assert(supersequence.size() >= path.size());
		for (size_t i = 1; i < path.size(); i++)
		{
			while (supersequence[i+offset] != path[i])
			{
				offset++;
				assert(i+offset < supersequence.size());
			}
			supersequencePredecessors[lastIndex].insert(i+offset);
			lastIndex = i+offset;
		}
	}

	filterUnnecessaryCharacters(cycleStart, sizeLeft, supersequence, supersequencePredecessors);
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

std::vector<size_t> CycleCutCalculation::getCycleCuttersSupersequence(size_t cycleStart, int sizeLeft, const std::vector<std::vector<size_t>>& paths) const
{
	std::vector<size_t> supersequence;
	for (const auto& path : paths)
	{
		if (supersequence.size() == 0)
		{
			assert(path.size() > 0);
			supersequence = path;
			continue;
		}
		supersequence = getPairwiseSupersequenceByAlignment(supersequence, path);
	}
	return supersequence;
}

std::vector<std::map<size_t, size_t>> findFeasibleFlow(const std::vector<size_t>& supersequence, const std::vector<std::set<size_t>>& predecessors)
{
	std::vector<std::map<size_t, size_t>> result;
	std::vector<std::vector<std::pair<size_t, size_t>>> pathBack;
	std::vector<std::vector<std::pair<size_t, size_t>>> pathForward;
	std::vector<std::set<size_t>> successors;
	pathBack.resize(supersequence.size());
	pathForward.resize(supersequence.size());
	successors.resize(supersequence.size());
	result.resize(supersequence.size());
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
			if (result[i][predecessor] > 0) continue;
			result[i][predecessor]++;;
			for (auto part : pathBack[i])
			{
				result[part.first][part.second]++;
			}
			for (auto part : pathForward[predecessor])
			{
				result[part.first][part.second]++;
			}
		}
	}
	return result;
}

void getOneFlowPath(const std::vector<size_t>& supersequence, const std::vector<std::set<size_t>>& predecessors, const std::vector<std::map<size_t, size_t>>& flows, size_t node, std::vector<size_t>& result)
{
	result.push_back(node);
	for (auto predecessor : predecessors[node])
	{
		assert(predecessor > node);
		if (flows[node].at(predecessor) > 0)
		{
			getOneFlowPath(supersequence, predecessors, flows, predecessor, result);
			return;
		}
	}
}

std::vector<std::vector<size_t>> findFlowPaths(const std::vector<size_t>& supersequence, const std::vector<std::set<size_t>>& predecessors, std::vector<std::map<size_t, size_t>> flows)
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
			assert(flows[path[i-1]][path[i]] >= 1);
			flows[path[i-1]][path[i]] -= 1;
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
	for (size_t i = 0; i < flows.size(); i++)
	{
		for (auto x : flows[i])
		{
			assert(x.second == 0);
		}
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
	Traits::vertex_descriptor src, sink;
	std::vector<Traits::vertex_descriptor> verts;
	for (size_t i = 0; i < supersequence.size(); i++)
	{
		verts.push_back(boost::add_vertex(g));
	}
	src = boost::add_vertex(g);
	sink = boost::add_vertex(g);
	auto startFlow = findFeasibleFlow(supersequence, predecessors);
	{
		Traits::edge_descriptor e1, e2;
		bool in1, in2;
		boost::tie(e1, in1) = add_edge(src, verts[0], g);
		boost::tie(e2, in2) = add_edge(verts[0], src, g);
		assert(in1 && in2);
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
			assert(in1 && in2);
			assert(startFlow[i][predecessor] >= 1);
			capacity[e1] = startFlow[i][predecessor] - 1;
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
			assert(in1 && in2);
			capacity[e1] = supersequence.size() * supersequence.size();
			capacity[e2] = 0;
			rev[e1] = e2;
			rev[e2] = e1;
		}
	}

	auto boostFlow = push_relabel_max_flow(g, src, sink);

	boost::graph_traits<Graph>::vertex_iterator u_iter, u_end;
	boost::graph_traits<Graph>::out_edge_iterator ei, e_end;
	for (boost::tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter)
	{
		for (boost::tie(ei, e_end) = out_edges(*u_iter, g); ei != e_end; ++ei)
		{
			if (capacity[*ei] > 0)
			{
				auto from = *u_iter;
				auto to = target(*ei, g);
				if (from >= supersequence.size() || to >= supersequence.size()) continue;
				startFlow[from][to] -= (capacity[*ei] - residual_capacity[*ei]);
			}
		}
	}

	for (auto path : findFlowPaths(supersequence, predecessors, startFlow))
	{
		function(path);
	}
}

void CycleCutCalculation::getCycleCutters(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const
{
	std::vector<std::vector<size_t>> paths;
	iterateOverEdgeCoveringPaths(cycleStart, sizeLeft, [&paths](const std::vector<size_t>& path)
	{
		paths.push_back(path);
	});
	supersequence = getCycleCuttersSupersequence(cycleStart, sizeLeft, paths);
	getPredecessorsFromSupersequenceOverEdgeCoveringPaths(cycleStart, sizeLeft, supersequence, supersequencePredecessors, previousCut, paths);
	for (size_t i = 0; i < supersequence.size(); i++)
	{
		previousCut.push_back(supersequence[i] < cycleStart);
	}
}

void CycleCutCalculation::getCycleCuttersTooBig(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const
{
	supersequence = getCycleCuttersOrder(cycleStart, sizeLeft, supersequencePredecessors);
	for (size_t i = 0; i < supersequence.size(); i++)
	{
		previousCut.push_back(supersequence[i] < cycleStart);
	}
}

AlignmentGraph::CycleCut CycleCutCalculation::GetCycleCutTooBig(size_t startNode, int wordSize) const
{
	AlignmentGraph::CycleCut result;
	getCycleCuttersTooBig(startNode, wordSize*2+1, result.nodes, result.predecessors, result.previousCut);
	return result;
}

AlignmentGraph::CycleCut CycleCutCalculation::GetCycleCut(size_t startNode, int wordSize) const
{
	AlignmentGraph::CycleCut result;
	getCycleCutters(startNode, wordSize*2, result.nodes, result.predecessors, result.previousCut);
	return result;
}
