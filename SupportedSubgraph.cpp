#include <algorithm>
#include <fstream>
#include "vg.pb.h"
#include "stream.hpp"


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

int main(int argc, char** argv)
{
	vg::Graph graph;
	{
		std::ifstream graphfile { argv[1], std::ios::in | std::ios::binary };
		std::vector<vg::Graph> parts;
		std::function<void(vg::Graph&)> lambda = [&parts](vg::Graph& g) {
			parts.push_back(g);
		};
		stream::for_each(graphfile, lambda);
		graph = mergeGraphs(parts);
	}

	std::vector<vg::Alignment> alignments;
	{
		std::ifstream alignmentfile {argv[2], std::ios::in | std::ios::binary};
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

	for (size_t i = 0; i < alignments.size(); i++)
	{
		std::cout << "alignment " << alignments[i].name() << std::endl;
		for (size_t j = 0; j < alignments[i].path().mapping_size()-1; j++)
		{
			auto from = alignments[i].path().mapping(j).position().node_id();
			auto to = alignments[i].path().mapping(j+1).position().node_id();
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

	std::ofstream graphOut { argv[3], std::ios::out | std::ios::binary };
	std::vector<vg::Graph> writeVector {resultGraph};
	stream::write_buffered(graphOut, writeVector, 0);
}