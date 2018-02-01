#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "CommonUtils.h"
#include "GfaGraph.h"

int main(int argc, char** argv)
{
	std::string ingraphfilename {argv[1]};
	std::string inalignmentfilename {argv[2]};
	std::string outfilename {argv[3]};

	std::cerr << "load graph" << std::endl;
	auto graph = GfaGraph::LoadFromFile(ingraphfilename);
	std::cerr << "process graph" << std::endl;
	std::unordered_map<int, std::unordered_map<std::string, size_t>> baseCounts;
	std::unordered_map<int, std::vector<int>> outNeighbors;
	std::unordered_map<int, std::vector<int>> leftInneighbors;
	std::unordered_map<int, std::vector<int>> rightInneighbors;
	std::unordered_map<int, size_t> counts;
	{
		std::unordered_map<NodePos, std::unordered_set<NodePos>> edges;
		for (auto edge : graph.edges)
		{
			for (auto target : edge.second)
			{
				edges[edge.first].insert(target);
				edges[NodePos{target.id, !target.end}].insert(NodePos{edge.first.id, !edge.first.end});
			}
		}
		for (auto node : graph.nodes)
		{
			NodePos end {node.first, true};
			if (edges.count(end) == 1 && edges[end].size() == 1)
			{
				outNeighbors[node.first].push_back(edges[end].begin()->id);
				if (edges[end].begin()->end)
				{
					rightInneighbors[edges[end].begin()->id].push_back(node.first);
				}
				else
				{
					leftInneighbors[edges[end].begin()->id].push_back(node.first);
				}
			}
			NodePos start {node.first, false};
			if (edges.count(start) == 1 && edges[start].size() == 1)
			{
				outNeighbors[node.first].push_back(edges[start].begin()->id);
				if (edges[start].begin()->end)
				{
					rightInneighbors[edges[start].begin()->id].push_back(node.first);
				}
				else
				{
					leftInneighbors[edges[start].begin()->id].push_back(node.first);
				}
			}
			if (edges.count(start) == 1) counts[node.first] = std::max(counts[node.first], edges[start].size());
			if (edges.count(end) == 1) counts[node.first] = std::max(counts[node.first], edges[end].size());
		}
	}

	std::cerr << "load alignment" << std::endl;
	auto alignments = CommonUtils::LoadVGAlignments(inalignmentfilename);
	for (auto aln : alignments)
	{
		for (size_t i = 0; i < aln.path().mapping_size(); i++)
		{
			baseCounts[aln.path().mapping(i).position().node_id()][aln.name()] += 1;
		}
	}
	std::cerr << "init counts" << std::endl;
	for (auto pair : baseCounts)
	{
		for (auto count : pair.second)
		{
			counts[pair.first] = std::max(counts[pair.first], count.second);
		}
	}

	std::cerr << "iterate" << std::endl;
	std::vector<int> updateQueue;
	updateQueue.reserve(graph.nodes.size());
	for (const auto& node : graph.nodes)
	{
		updateQueue.push_back(node.first);
	}
	std::cerr << "numnodes " << updateQueue.size() << std::endl;
	size_t iterated = 0;
	int maxcount = 0;
	while (updateQueue.size() > 0)
	{
		auto node = updateQueue.back();
		updateQueue.pop_back();
		iterated++;
		if (iterated % 1000000 == 0) std::cerr << "iterated " << iterated << std::endl;
		int leftCountShouldBe = 0;
		if (leftInneighbors.count(node) == 1)
		{
			for (auto neighbor : leftInneighbors.at(node))
			{
				leftCountShouldBe += counts[neighbor];
			}
		}
		int rightCountShouldBe = 0;
		if (rightInneighbors.count(node) == 1)
		{
			for (auto neighbor : rightInneighbors.at(node))
			{
				rightCountShouldBe += counts[neighbor];
			}
		}
		if (counts[node] >= leftCountShouldBe && counts[node] >= rightCountShouldBe) continue;
		counts[node] = std::max(leftCountShouldBe, rightCountShouldBe);
		if (counts[node] > maxcount)
		{
			maxcount = counts[node];
			std::cerr << "node " << node << " iter " << iterated << " maxcount " << maxcount << std::endl;
		}
		if (outNeighbors.count(node) == 1)
		{
			for (auto neighbor : outNeighbors.at(node))
			{
				updateQueue.push_back(neighbor);
			}
		}
	}
	std::cerr << "iteration done with " << iterated << std::endl;

	std::cerr << "write result" << std::endl;
	std::ofstream out {outfilename};
	out << "node,_minalntoporepeatcount";
	out << std::endl;
	std::vector<int> nodevec;
	for (const auto& node : graph.nodes)
	{
		nodevec.push_back(node.first);
	}
	std::sort(nodevec.begin(), nodevec.end());
	for (auto node : nodevec)
	{
		out << node;
		out << "," << counts[node];
		out << std::endl;
	}
}
