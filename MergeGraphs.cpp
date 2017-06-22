#include <fstream>
#include <algorithm>
#include <set>
#include "vg.pb.h"
#include "stream.hpp"

// Separate script to merge augmented vg graphs, code from aligner.cpp
vg::Graph mergeGraphs(const std::vector<vg::Graph>& parts)
{
	vg::Graph newGraph;
	std::vector<const vg::Node*> allNodes;
	std::vector<const vg::Edge*> allEdges;
	for (size_t i = 0; i < parts.size(); i++)
	{
		for (int j = 0; j < parts[i].node_size(); j++)
		{
			allNodes.push_back(&parts[i].node(j));
		}
		for (int j = 0; j < parts[i].edge_size(); j++)
		{
			allEdges.push_back(&parts[i].edge(j));
		}
	}
	for (size_t i = 0; i < allNodes.size(); i++)
	{
		auto node = newGraph.add_node();
		node->set_id(allNodes[i]->id());
		node->set_sequence(allNodes[i]->sequence());
		node->set_name(allNodes[i]->name());
	}
	for (size_t i = 0; i < allEdges.size(); i++)
	{
		auto edge = newGraph.add_edge();
		edge->set_from(allEdges[i]->from());
		edge->set_to(allEdges[i]->to());
		edge->set_from_start(allEdges[i]->from_start());
		edge->set_to_end(allEdges[i]->to_end());
		edge->set_overlap(allEdges[i]->overlap());
	}
	return newGraph;
}

void outputGraph(std::string filename, const vg::Graph& graph)
{
	std::ofstream alignmentOut { filename, std::ios::out | std::ios::binary };
	std::vector<vg::Graph> writeVector {graph};
	stream::write_buffered(alignmentOut, writeVector, 0);
}

int main(int argc, char** argv)
{
	int n =4; // TODO: Hard-coded 4 vg graphs here, can be generalized.
	vg::Graph graph;
	std::vector<vg::Graph> parts;
	for (int i=0; i < n; i++){

		std::ifstream graphfile { argv[i], std::ios::in | std::ios::binary };
		std::function<void(vg::Graph&)> lambda = [&parts](vg::Graph& g) {
			parts.push_back(g);
		};
		stream::for_each(graphfile, lambda);
	}		
		graph = mergeGraphs(parts);
		outputGraph("finalaugmentedgraph.vg", graph);
}
