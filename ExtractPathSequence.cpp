#include <algorithm>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <functional>
#include "GfaGraph.h"
#include "vg.pb.h"
#include "stream.hpp"
#include "CommonUtils.h"

void printPath(const std::map<int, int>& ids, std::function<std::string(int)> seqGetter, const vg::Alignment& v)
{
	std::cout << ">" << v.name() << std::endl;
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

void printPath(const vg::Graph& g, const vg::Alignment& v)
{
	std::map<int, int> ids;
	for (int i = 0; i < g.node_size(); i++)
	{
		ids[g.node(i).id()] = i;
	}
	printPath(ids, [&g](int id) {return g.node(id).sequence();}, v);
}

void printPath(const GfaGraph& g, const vg::Alignment& v)
{
	std::map<int, int> ids;
	for (int i = 0; i < g.nodes.size(); i++)
	{
		ids[g.nodes[i].id] = i;
	}
	printPath(ids, [&g](int id) {return g.nodes[id].sequence;}, v);
}

int main(int argc, char** argv)
{
	std::string graphfilename {argv[1]};
	if (graphfilename.substr(graphfilename.size()-3) == ".vg")
	{
		vg::Graph graph = CommonUtils::LoadVGGraph(argv[1]);
		{
			std::ifstream graphfile { argv[2], std::ios::in | std::ios::binary };
			std::function<void(vg::Alignment&)> lambda = [&graph](vg::Alignment& g) {
				std::cerr << g.name() << std::endl;
				printPath(graph, g);
			};
			stream::for_each(graphfile, lambda);
		}
	}
	else if (graphfilename.substr(graphfilename.size() - 4) == ".gfa")
	{
		GfaGraph graph = GfaGraph::LoadFromFile(argv[1]);
		{
			std::ifstream graphfile { argv[2], std::ios::in | std::ios::binary };
			std::function<void(vg::Alignment&)> lambda = [&graph](vg::Alignment& g) {
				std::cerr << g.name() << std::endl;
				printPath(graph, g);
			};
			stream::for_each(graphfile, lambda);
		}
	}

}