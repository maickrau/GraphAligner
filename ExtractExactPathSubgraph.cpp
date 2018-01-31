#include <fstream>
#include <iostream>
#include <queue>
#include "GfaGraph.h"
#include "CommonUtils.h"

int main(int argc, char** argv)
{
	std::string infile {argv[1]};
	std::string outfile {argv[2]};
	std::string alignmentfile {argv[3]};
	auto alignments = CommonUtils::LoadVGAlignments(alignmentfile);
	auto graph = GfaGraph::LoadFromFile(infile);
	std::unordered_set<int> pickedNodes;
	std::unordered_set<std::pair<NodePos, NodePos>> pickedEdges;
	for (const auto& alignment : alignments)
	{
		pickedNodes.insert(alignment.path().mapping(0).position().node_id());
		for (int i = 1; i < alignment.path().mapping_size(); i++)
		{
			pickedNodes.insert(alignment.path().mapping(i).position().node_id());
			NodePos from {alignment.path().mapping(i-1).position().node_id(), alignment.path().mapping(i-1).position().is_reverse()};
			NodePos to {alignment.path().mapping(i).position().node_id(), alignment.path().mapping(i).position().is_reverse()};
			pickedEdges.emplace(from, to);
		}
	}
	std::cerr << pickedNodes.size() << " nodes, ~" << pickedEdges.size() << " edges" << std::endl;
	auto result = graph.GetSubgraph(pickedNodes, pickedEdges);
	result.SaveToFile(outfile);
}