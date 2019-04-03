#include <algorithm>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <iostream>
#include <functional>
#include "GfaGraph.h"
#include "vg.pb.h"
#include "stream.hpp"
#include "CommonUtils.h"

void printPath(const std::unordered_map<int, int>& ids, std::function<std::string(int)> seqGetter, const vg::Alignment& v)
{
	std::cout << ">" << v.name() << "_" << v.query_position() << "_" << (v.query_position() + v.sequence().size()) << std::endl;
	for (int i = 0; i < v.path().mapping_size(); i++)
	{
		auto nodeid = v.path().mapping(i).position().node_id();
		auto sequence = seqGetter(ids.at(nodeid));
		int len = 0;
		for (int j = 0; j < v.path().mapping(i).edit_size(); j++)
		{
			len += v.path().mapping(i).edit(j).from_length();
		}
		if (v.path().mapping(i).position().is_reverse())
		{
			sequence = CommonUtils::ReverseComplement(sequence);
		}
		if (v.path().mapping(i).position().offset() > 0)
		{
			sequence = sequence.substr(v.path().mapping(i).position().offset());
		}
		sequence = sequence.substr(0, len);
		std::cout << sequence;
	}
	std::cout << std::endl;
}

void printPath(const vg::Graph& g, const std::unordered_map<int, int>& ids, const vg::Alignment& v)
{
	printPath(ids, [&g](int id) {return g.node(id).sequence();}, v);
}

void printPath(const GfaGraph& g, const std::unordered_map<int, int>& ids, const vg::Alignment& v)
{
	printPath(ids, [&g](int id) {return g.nodes.at(id);}, v);
}

int main(int argc, char** argv)
{
	std::string graphfilename {argv[1]};
	std::string alnfilename { argv[2] };
	std::unordered_map<int, int> ids;
	if (graphfilename.substr(graphfilename.size()-3) == ".vg")
	{
		vg::Graph graph = CommonUtils::LoadVGGraph(argv[1]);
		for (int i = 0; i < graph.node_size(); i++)
		{
			ids[graph.node(i).id()] = i;
		}
		{
			std::ifstream graphfile { argv[2], std::ios::in | std::ios::binary };
			std::function<void(vg::Alignment&)> lambda = [&graph, &ids](vg::Alignment& g) {
				std::cerr << g.name() << std::endl;
				printPath(graph, ids, g);
			};
			stream::for_each(graphfile, lambda);
		}
	}
	else if (graphfilename.substr(graphfilename.size() - 4) == ".gfa")
	{
		GfaGraph graph = GfaGraph::LoadFromFile(argv[1]);
		for (auto node : graph.nodes)
		{
			ids[node.first] = node.first;
		}
		{
			std::ifstream graphfile { argv[2], std::ios::in | std::ios::binary };
			std::function<void(vg::Alignment&)> lambda = [&graph, &ids](vg::Alignment& g) {
				std::cerr << g.name() << std::endl;
				printPath(graph, ids, g);
			};
			stream::for_each(graphfile, lambda);
		}
	}

}