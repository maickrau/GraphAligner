#include <cassert>
#include <string>
#include <unordered_map>
#include "CommonUtils.h"
#include "GfaGraph.h"

GfaGraph unitigify(const GfaGraph& graph)
{
	GfaGraph result;
	result.edgeOverlap = graph.edgeOverlap;
	std::unordered_map<NodePos, int> belongsInUnitig;
	std::unordered_set<int> nodesHandled;
	std::unordered_map<int, NodePos> unitigLeft;
	std::unordered_map<int, NodePos> unitigRight;
	std::vector<std::vector<NodePos>> nodesInUnitig;

	for (auto node : graph.nodes)
	{
		if (nodesHandled.count(node.first) == 1) continue;
		NodePos left { node.first, false };		
		NodePos right { node.first, true };
		bool leftBreaks = true;
		bool rightBreaks = true;
		if (graph.edges.count(left) == 1 && graph.edges.at(left).size() == 1)
		{
			auto neighbor = graph.edges.at(left)[0];
			assert(graph.edges.count(neighbor.Reverse()) == 1);
			if (graph.edges.at(neighbor.Reverse()).size() == 1) leftBreaks = false;
		}
		if (graph.edges.count(right) == 1 && graph.edges.at(right).size() == 1)
		{
			auto neighbor = graph.edges.at(right)[0];
			assert(graph.edges.count(neighbor.Reverse()) == 1);
			if (graph.edges.at(neighbor.Reverse()).size() == 1) rightBreaks = false;
		}
		if (leftBreaks && rightBreaks)
		{
			assert(nodesHandled.count(node.first) == 0);
			nodesHandled.insert(node.first);
			unitigLeft[nodesInUnitig.size()] = NodePos { node.first, true };
			unitigRight[nodesInUnitig.size()] = NodePos { node.first, true };
			assert(belongsInUnitig.count(right) == 0);
			belongsInUnitig[right] = nodesInUnitig.size();
			nodesInUnitig.emplace_back();
			nodesInUnitig.back().emplace_back(node.first, true);
			continue;
		}
		if (!leftBreaks && !rightBreaks)
		{
			continue;
		}
		assert((leftBreaks && !rightBreaks) || (rightBreaks && !leftBreaks));
		NodePos start;
		int id = nodesInUnitig.size();
		nodesInUnitig.emplace_back();
		start.id = node.first;
		start.end = leftBreaks;
		assert(belongsInUnitig.count(start) == 0);
		assert(belongsInUnitig.count(start.Reverse()) == 0);
		assert(nodesHandled.count(start.id) == 0);
		unitigLeft[id] = start;
		unitigRight[id] = start;
		nodesHandled.insert(start.id);
		belongsInUnitig[start] = id;
		nodesInUnitig.back().push_back(start);
		assert(graph.edges.count(start) == 1 && graph.edges.at(start).size() == 1);
		while (graph.edges.count(start) == 1 && graph.edges.at(start).size() == 1)
		{
			start = graph.edges.at(start)[0];
			assert(graph.edges.count(start.Reverse()) == 1);
			if (graph.edges.at(start.Reverse()).size() != 1) break;
			assert(belongsInUnitig.count(start) == 0);
			assert(belongsInUnitig.count(start.Reverse()) == 0);
			assert(nodesHandled.count(start.id) == 0);
			unitigRight[id] = start;
			nodesHandled.insert(start.id);
			belongsInUnitig[start] = id;
			nodesInUnitig.back().push_back(start);
		}
	}
	//circular separate components
	for (auto node : graph.nodes)
	{
		if (nodesHandled.count(node.first) == 1) continue;
		NodePos left { node.first, false };		
		NodePos right { node.first, true };
		assert(graph.edges.count(left) == 1 && graph.edges.at(left).size() == 1);
		assert(graph.edges.count(right) == 1 && graph.edges.at(right).size() == 1);
		NodePos start = right;
		int id = nodesInUnitig.size();
		nodesInUnitig.emplace_back();
		unitigLeft[id] = start;
		unitigRight[id] = start;
		do
		{
			nodesHandled.insert(node.first);
			belongsInUnitig[start] = id;
			nodesInUnitig.back().push_back(start);
			assert(graph.edges.count(start) == 1 && graph.edges.at(start).size() == 1);
			assert(graph.edges.count(start.Reverse()) == 1 && graph.edges.at(start.Reverse()).size() == 1);
			start = graph.edges.at(start)[0];
		} while (start.id != node.first);
		result.edges[NodePos { id, true }].emplace_back(id, true);
	}
	assert(nodesHandled.size() == graph.nodes.size());
	assert(belongsInUnitig.size() == graph.nodes.size());
	for (size_t i = 0; i < nodesInUnitig.size(); i++)
	{
		std::string seq;
		assert(nodesInUnitig[i].size() > 0);
		seq = graph.nodes.at(nodesInUnitig[i][0].id);
		if (!nodesInUnitig[i][0].end) seq = CommonUtils::ReverseComplement(seq);
		seq = seq.substr(0, graph.edgeOverlap);
		for (auto node : nodesInUnitig[i])
		{
			std::string add;
			add = graph.nodes.at(node.id);
			if (!node.end) add = CommonUtils::ReverseComplement(add);
			add = add.substr(graph.edgeOverlap);
			seq += add;
		}
		result.nodes[i] = seq;
	}
	for (auto edge : graph.edges)
	{
		NodePos src = edge.first;
		NodePos from;
		assert(belongsInUnitig.count(src) == 1 || belongsInUnitig.count(src.Reverse()) == 1);
		if (belongsInUnitig.count(src) == 1)
		{
			assert(belongsInUnitig.count(src) == 1);
			assert(belongsInUnitig.count(src.Reverse()) == 0);
			if (unitigRight[belongsInUnitig[src]] != src) continue;
			from = NodePos { belongsInUnitig[src], true };
		}
		else
		{
			assert(belongsInUnitig.count(src) == 0);
			assert(belongsInUnitig.count(src.Reverse()) == 1);
			if (unitigLeft[belongsInUnitig[src.Reverse()]] != src.Reverse()) continue;
			from = NodePos { belongsInUnitig[src.Reverse()], false };
		}
		for (auto dst : edge.second)
		{
			NodePos to;
			assert(belongsInUnitig.count(dst) == 1 || belongsInUnitig.count(dst.Reverse()) == 1);
			if (belongsInUnitig.count(dst) == 1)
			{
				assert(belongsInUnitig.count(dst) == 1);
				assert(belongsInUnitig.count(dst.Reverse()) == 0);
				if (unitigLeft[belongsInUnitig[dst]] != dst) continue;
				to = NodePos { belongsInUnitig[dst], true };
			}
			else
			{
				assert(belongsInUnitig.count(dst) == 0);
				assert(belongsInUnitig.count(dst.Reverse()) == 1);
				if (unitigRight[belongsInUnitig[dst.Reverse()]] != dst.Reverse()) continue;
				to = NodePos { belongsInUnitig[dst.Reverse()], false };
			}
			result.edges[from].push_back(to);
		}
	}
	return result;
}

int main(int argc, char** argv)
{
	std::string inputGraph { argv[1] };
	std::string outputGraph { argv[2] };

	auto graph = GfaGraph::LoadFromFile(inputGraph);
	graph.confirmDoublesidedEdges();
	auto result = unitigify(graph);
	result.SaveToFile(outputGraph);
}