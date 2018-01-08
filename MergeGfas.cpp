#include <fstream>
#include <iostream>
#include <queue>
#include "GfaGraph.h"
#include "CommonUtils.h"

int main(int argc, char** argv)
{
	if (argc < 3) return -1;
	std::string outfile {argv[1]};
	auto graph = GfaGraph::LoadFromFile(argv[2]);
	for (int i = 3; i < argc; i++)
	{
		auto nextgraph = GfaGraph::LoadFromFile(argv[i]);
		graph.AddSubgraph(nextgraph);
	}
	graph.SaveToFile(outfile);
}