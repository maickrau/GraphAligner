#include <map>
#include <iostream>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <algorithm>
#include "GfaGraph.h"

std::unordered_map<NodePos, size_t> getNodeMapping(const GfaGraph& graph)
{
	std::unordered_map<NodePos, size_t> result;
	for (auto node : graph.nodes)
	{
		size_t id = result.size();
		result[NodePos { node.first, true }] = id;
		id = result.size();
		result[NodePos { node.first, false }] = id;
	}
	return result;
}

std::vector<size_t> getLengths(const std::unordered_map<NodePos, size_t>& nodeMapping, const GfaGraph& graph)
{
	std::vector<size_t> result;
	result.resize(nodeMapping.size(), 0);
	for (auto node : graph.nodes)
	{
		result[nodeMapping.at(NodePos{ node.first, true })] = node.second.size();
		result[nodeMapping.at(NodePos{ node.first, false })] = node.second.size();
	}
	return result;
}

std::vector<std::vector<size_t>> getOutEdges(const std::unordered_map<NodePos, size_t>& nodeMapping, const GfaGraph& graph)
{
	std::vector<std::vector<size_t>> result;
	result.resize(nodeMapping.size());
	for (auto edge : graph.edges)
	{
		NodePos source = edge.first;
		NodePos revSource = source.Reverse();
		for (auto target : edge.second)
		{
			NodePos revTarget = target.Reverse();
			assert(nodeMapping.at(source) < result.size());
			assert(nodeMapping.at(target) < result.size());
			assert(nodeMapping.at(revSource) < result.size());
			assert(nodeMapping.at(revTarget) < result.size());
			result[nodeMapping.at(source)].push_back(nodeMapping.at(target));
			result[nodeMapping.at(revTarget)].push_back(nodeMapping.at(revSource));
		}
	}
	return result;
}

std::vector<size_t> getNodeDepths(const std::vector<std::vector<size_t>>& componentNodes, const std::vector<size_t>& nodeLengths, const std::vector<std::vector<size_t>>& edges, const std::map<std::pair<size_t, size_t>, size_t>& overlaps, const size_t dbgOverlap)
{
	std::vector<size_t> result;
	result.resize(nodeLengths.size(), 0);
	for (size_t i = componentNodes.size()-1; i < componentNodes.size(); i--)
	{
		if (componentNodes[i].size() > 1)
		{
			for (auto node : componentNodes[i])
			{
				result[node] = std::numeric_limits<size_t>::max();
			}
		}
		else
		{
			auto node = componentNodes[i][0];
			result[node] = 0;
			for (auto neighbor : edges[node])
			{
				if (result[neighbor] == std::numeric_limits<size_t>::max())
				{
					result[node] = std::numeric_limits<size_t>::max();
					break;
				}
				if (neighbor == node)
				{
					result[node] = std::numeric_limits<size_t>::max();
					break;
				}
				size_t overlap = dbgOverlap;
				if (overlaps.count(std::make_pair(node, neighbor)) == 1)
				{
					overlap = overlaps.at(std::make_pair(node, neighbor));
				}
				assert(nodeLengths[neighbor] > overlap);
				result[node] = std::max(result[node], result[neighbor] + nodeLengths[neighbor] - overlap);
			}
		}
	}
	return result;
}

void removeRec(std::vector<bool>& keepers, size_t pos, const std::vector<std::vector<size_t>>& edges)
{
	if (!keepers[pos]) return;
	keepers[pos] = false;
	for (auto neighbor : edges[pos])
	{
		removeRec(keepers, neighbor, edges);
	}
}

std::vector<bool> getKeepers(const std::vector<size_t>& depths, const std::vector<std::vector<size_t>>& edges, const size_t maxRemovableLen, const size_t minSafeLen, const double fraction, const std::map<std::pair<size_t, size_t>, size_t>& overlaps, const size_t dbgOverlap, const std::vector<size_t>& nodeLengths)
{
	std::vector<bool> result;
	result.resize(depths.size(), true);
	for (size_t i = 0; i < depths.size(); i++)
	{
		if (!result[i]) continue;
		size_t bigLength = 0;
		for (auto neighbor : edges[i])
		{
			if (depths[neighbor] == std::numeric_limits<size_t>::max())
			{
				bigLength = std::numeric_limits<size_t>::max();
				break;
			}
			size_t overlap = dbgOverlap;
			if (overlaps.count(std::make_pair(i, neighbor)) == 1)
			{
				overlap = overlaps.at(std::make_pair(i, neighbor));
			}
			assert(nodeLengths[neighbor] > overlap);
			bigLength = std::max(bigLength, depths[neighbor] + nodeLengths[neighbor] - overlap);
		}
		if (bigLength < minSafeLen) continue;
		size_t removableLen = bigLength * fraction;
		removableLen = std::min(removableLen, maxRemovableLen);
		for (auto neighbor : edges[i])
		{
			if (depths[neighbor] == std::numeric_limits<size_t>::max()) continue;
			size_t overlap = dbgOverlap;
			if (overlaps.count(std::make_pair(i, neighbor)) == 1)
			{
				overlap = overlaps.at(std::make_pair(i, neighbor));
			}
			assert(nodeLengths[neighbor] > overlap);
			if (depths[neighbor] + nodeLengths[neighbor] - overlap <= removableLen)
			{
				removeRec(result, neighbor, edges);
			}
		}
	}
	return result;
}

void strongConnectIterative(size_t node, size_t& i, std::vector<size_t>& index, std::vector<size_t>& lowlink, std::vector<bool>& onStack, std::vector<size_t>& S, std::vector<std::vector<size_t>>& result, const std::vector<std::vector<size_t>>& edges)
{
	std::vector<std::tuple<int, size_t, size_t>> stack;
	stack.emplace_back(0, node, 0);
	while (stack.size() > 0)
	{
		auto top = stack.back();
		size_t node = std::get<1>(top);
		size_t neighborI = std::get<2>(top);
		stack.pop_back();
		switch(std::get<0>(top))
		{
			case 0:
				assert(!onStack[node]);
				assert(index[node] == -1);
				assert(lowlink[node] == -1);
				index[node] = i;
				lowlink[node] = i;
				i++;
				S.push_back(node);
				onStack[node] = true;
			START_LOOP:
			case 1:
				if (neighborI < edges[node].size())
				{
					auto neighbor = edges[node][neighborI];
					if (index[neighbor] == -1)
					{
						stack.emplace_back(2, node, neighborI);
						stack.emplace_back(0, edges[node][neighborI], 0);
						continue;
					}
					else if (onStack[neighbor])
					{
						assert(index[neighbor] != -1);
						lowlink[node] = std::min(lowlink[node], index[neighbor]);
					}
					neighborI++;
				}
				if (neighborI < edges[node].size()) goto START_LOOP;
				goto END_LOOP;
			case 2:
				{
					auto neighbor = edges[node][neighborI];
					assert(lowlink[neighbor] != -1);
					lowlink[node] = std::min(lowlink[node], lowlink[neighbor]);
					neighborI++;
					goto START_LOOP;
				}
			END_LOOP:
			case 3:
				assert(lowlink[node] != -1);
				assert(index[node] != -1);
				if (lowlink[node] == index[node])
				{
					result.emplace_back();
					size_t stacknode;
					do
					{
						assert(S.size() > 0);
						stacknode = S.back();
						S.pop_back();
						assert(onStack[stacknode]);
						onStack[stacknode] = false;
						result.back().push_back(stacknode);
					} while (stacknode != node);
				}
		}
	}
}

std::vector<std::vector<size_t>> topologicalSort(const std::vector<std::vector<size_t>>& edges)
{
	std::vector<size_t> index;
	std::vector<size_t> lowlink;
	std::vector<bool> onStack;
	index.resize(edges.size(), -1);
	lowlink.resize(edges.size(), -1);
	onStack.resize(edges.size(), false);
	std::vector<size_t> S;
	std::vector<std::vector<size_t>> result;
	size_t i = 0;
	for (size_t node = 0; node < edges.size(); node++)
	{
		if (index[node] == -1) strongConnectIterative(node, i, index, lowlink, onStack, S, result, edges);
		assert(S.size() == 0);
	}
	assert(i == edges.size());
	std::reverse(result.begin(), result.end());
	std::vector<size_t> belongsToComponent;
	belongsToComponent.resize(edges.size(), -1);
	for (size_t i = 0; i < result.size(); i++)
	{
		for (auto node : result[i])
		{
			belongsToComponent[node] = i;
		}
	}
	for (size_t i = 0; i < edges.size(); i++)
	{
		assert(belongsToComponent[i] != -1);
		for (auto edge : edges[i])
		{
			assert(belongsToComponent[edge] != -1);
			assert(belongsToComponent[edge] >= belongsToComponent[i]);
		}
	}
	return result;
}

std::map<std::pair<size_t, size_t>, size_t> getOverlaps(const std::unordered_map<NodePos, size_t>& nodeMapping, const GfaGraph& graph)
{
	std::map<std::pair<size_t, size_t>, size_t> result;
	for (auto pair : graph.varyingOverlaps)
	{
		auto fromGfa = pair.first.first;
		auto toGfa = pair.first.second;
		auto overlap = pair.second;
		auto from = nodeMapping.at(fromGfa);
		auto to = nodeMapping.at(toGfa);
		result[std::make_pair(from, to)] = overlap;
	}
	return result;
}

std::unordered_set<int> filterNodes(const GfaGraph& graph, const int maxRemovableLen, const int minSafeLen, const double fraction)
{
	auto nodeMapping = getNodeMapping(graph);
	auto lengths = getLengths(nodeMapping, graph);
	auto overlaps = getOverlaps(nodeMapping, graph);
	auto edges = getOutEdges(nodeMapping, graph);
	auto order = topologicalSort(edges);
	auto depths = getNodeDepths(order, lengths, edges, overlaps, graph.edgeOverlap);
	auto keepers = getKeepers(depths, edges, maxRemovableLen, minSafeLen, fraction, overlaps, graph.edgeOverlap, lengths);
	std::unordered_set<int> result;
	for (auto node : graph.nodes)
	{
		if (keepers[nodeMapping[NodePos { node.first, true }]] && keepers[nodeMapping[NodePos { node.first, false }]])
		{
			result.emplace(node.first);
		}
	}
	return result;
}

int main(int argc, char** argv)
{
	int maxRemovableLen = std::stoi(argv[1]);
	int minSafeLen = std::stoi(argv[2]);
	double fraction = std::stod(argv[3]);
	auto graph = GfaGraph::LoadFromStream(std::cin, true);
	//write to cout

	auto keptNodes = filterNodes(graph, maxRemovableLen, minSafeLen, fraction);
	auto filteredGraph = graph.GetSubgraph(keptNodes);
	filteredGraph.SaveToStream(std::cout);
}