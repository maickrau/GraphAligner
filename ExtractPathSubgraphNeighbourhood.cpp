#include <fstream>
#include <iostream>
#include <queue>
#include "GfaGraph.h"
#include "CommonUtils.h"

class PriorityNode
{
public:
	PriorityNode(NodePos pos, int priority) :
	pos(pos),
	priority(priority)
	{}
	NodePos pos;
	int priority;
	bool operator>(const PriorityNode& other) const
	{
		return priority > other.priority;
	}
};

int main(int argc, char** argv)
{
	std::string infile {argv[1]};
	std::string outfile {argv[2]};
	std::string alignmentfile {argv[3]};
	int length = std::stoi(argv[4]);
	std::cerr << "length: " << length << std::endl;
	auto alignment = CommonUtils::LoadVGAlignment(alignmentfile);
	auto graph = GfaGraph::LoadFromFile(infile);
	std::priority_queue<PriorityNode, std::vector<PriorityNode>, std::greater<PriorityNode>> queue;
	for (const auto& pos : alignment.path().mapping())
	{
		queue.emplace(NodePos {pos.position().node_id(), pos.position().is_reverse()}, 0);
	}
	std::unordered_map<NodePos, size_t> distance;
	while (queue.size() != 0)
	{
		auto top = queue.top();
		queue.pop();
		if (top.priority > length) break;
		if (distance.count(top.pos) == 1 && distance[top.pos] <= top.priority) continue;
		distance[top.pos] = top.priority;
		if (graph.edges.count(top.pos) == 1)
		{
			for (auto edge : graph.edges.at(top.pos))
			{
				assert(graph.nodes.at(top.pos.id).size() > graph.edgeOverlap);
				queue.emplace(edge, top.priority + graph.nodes.at(top.pos.id).size() - graph.edgeOverlap);
			}
		}
	}
	std::unordered_set<int> picked;
	for (auto pair : distance)
	{
		picked.insert(pair.first.id);
	}
	std::cerr << picked.size() << std::endl;
	auto result = graph.GetSubgraph(picked);
	result.SaveToFile(outfile);
}