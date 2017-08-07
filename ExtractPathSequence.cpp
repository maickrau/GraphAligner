#include <algorithm>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include "vg.pb.h"
#include "stream.hpp"
#include "CommonUtils.h"

void printPath(const vg::Graph& g, const vg::Alignment& v)
{
	std::map<int, int> ids;
	for (int i = 0; i < g.node_size(); i++)
	{
		ids[g.node(i).id()] = i;
	}
	std::cout << ">" << v.name() << std::endl;
	for (int i = 0; i < v.path().mapping_size(); i++)
	{
		auto nodeid = v.path().mapping(i).position().node_id();
		auto sequence = g.node(ids[nodeid]).sequence();
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

int main(int argc, char** argv)
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