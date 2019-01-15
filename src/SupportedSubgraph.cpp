#include <algorithm>
#include <fstream>
#include <unordered_set>
#include "CommonUtils.h"
#include "vg.pb.h"
#include "stream.hpp"


int main(int argc, char** argv)
{
	std::string graphFile { argv[1] };
	std::string alnFile { argv[2] };
	std::string outputGraph { argv[3] };

	vg::Graph graph = CommonUtils::LoadVGGraph(graphFile);

	std::vector<vg::Alignment> alignments;
	{
		std::ifstream alignmentfile {alnFile, std::ios::in | std::ios::binary};
		std::function<void(vg::Alignment&)> lambda = [&alignments](vg::Alignment& g) {
			alignments.push_back(g);
		};
		stream::for_each(alignmentfile, lambda);
	}

	std::map<int, std::set<int>> existingEdges;
	for (size_t i = 0; i < graph.edge_size(); i++)
	{
		existingEdges[graph.edge(i).from()].insert(graph.edge(i).to());
	}

	std::map<int, std::set<int>> supportedEdges;
	std::unordered_set<int> supportedNodes;

	for (size_t i = 0; i < alignments.size(); i++)
	{
		std::cout << "alignment " << alignments[i].name() << std::endl;
		for (size_t j = 0; j < alignments[i].path().mapping_size()-1; j++)
		{
			auto from = alignments[i].path().mapping(j).position().node_id();
			auto to = alignments[i].path().mapping(j+1).position().node_id();
			supportedNodes.insert(from);
			supportedNodes.insert(to);
			if (existingEdges[from].count(to) == 0 && existingEdges[to].count(from) == 0)
			{
				std::cout << "nonexistant alignment from " << from << " to " << to << std::endl;
			}
			supportedEdges[from].insert(to);
		}
	}

	vg::Graph resultGraph;
	for (int i = 0 ; i < graph.node_size(); i++)
	{
		if (supportedNodes.count(graph.node(i).id()) == 0) continue;
		auto* node = resultGraph.add_node();
		node->set_sequence(graph.node(i).sequence());
		node->set_id(graph.node(i).id());
		node->set_name(graph.node(i).name());
	}
	for (int i = 0; i < graph.edge_size(); i++)
	{
		auto from = graph.edge(i).from();
		auto to = graph.edge(i).to();
		bool foundForward = supportedEdges[from].count(to) == 1;
		auto foundBackward = supportedEdges[to].count(from) == 1;
		if (!foundForward && !foundBackward)
		{
			continue;
		}
		auto* edge = resultGraph.add_edge();
		edge->set_from(graph.edge(i).from());
		edge->set_to(graph.edge(i).to());
		edge->set_from_start(graph.edge(i).from_start());
		edge->set_to_end(graph.edge(i).to_end());
		edge->set_overlap(graph.edge(i).overlap());
	}

	std::ofstream graphOut { outputGraph, std::ios::out | std::ios::binary };
	std::vector<vg::Graph> writeVector {resultGraph};
	stream::write_buffered(graphOut, writeVector, 0);
}