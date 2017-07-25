#include "CommonUtils.h"
#include "stream.hpp"

namespace CommonUtils
{

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

	vg::Graph LoadVGGraph(std::string filename)
	{
		vg::Graph result;
		std::ifstream graphfile { filename, std::ios::in | std::ios::binary };
		std::vector<vg::Graph> parts;
		std::function<void(vg::Graph&)> lambda = [&parts](vg::Graph& g) {
			parts.push_back(g);
		};
		stream::for_each(graphfile, lambda);
		result = mergeGraphs(parts);
		return result;
	}
}