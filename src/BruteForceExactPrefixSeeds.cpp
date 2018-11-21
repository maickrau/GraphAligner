#include <cmath>
#include <fstream>
#include "vg.pb.h"
#include "stream.hpp"
#include "GfaGraph.h"
#include "fastqloader.h"

void bruteForceAddPaths(const GfaGraph& graph, std::vector<std::vector<std::pair<int, size_t>>>& result, int node, size_t offset, size_t label, size_t length, size_t k)
{
	label <<= 2;
	switch(graph.nodes.at(node)[offset])
	{
		case 'A':
		case 'a':
			label += 0;
			break;
		case 'C':
		case 'c':
			label += 1;
			break;
		case 'G':
		case 'g':
			label += 2;
			break;
		case 'T':
		case 't':
			label += 3;
			break;
		default:
			assert(false);
	}
	if (length == k - 1)
	{
		assert(label < result.size());
		result[label].emplace_back(node, offset);
		return;
	}
	if (offset < graph.nodes.at(node).size() - 1)
	{
		bruteForceAddPaths(graph, result, node, offset+1, label, length+1, k);
	}
	else
	{
		if (graph.edges.count(NodePos{node, true}) == 1)
		{
			for (auto edge : graph.edges.at(NodePos{node, true}))
			{
				bruteForceAddPaths(graph, result, edge.id, graph.edgeOverlap, label, length+1, k);
			}
		}
	}
}

std::vector<std::vector<std::pair<int, size_t>>> buildBruteForcePathIndex(const GfaGraph& graph, const int k)
{
	std::vector<std::vector<std::pair<int, size_t>>> result;
	result.resize(pow(4, k));
	for (auto node : graph.nodes)
	{
		for (size_t i = 0; i < node.second.size(); i++)
		{
			bruteForceAddPaths(graph, result, node.first, i, 0, 0, k);
		}
	}
	return result;
}

int main(int argc, char** argv)
{
	std::string graphFile { argv[1] };
	std::string readFile { argv[2] };
	int k = std::stoi(argv[3]);
	std::string outputSeedFile { argv[4] };

	std::cerr << "load graph" << std::endl;
	auto graph = GfaGraph::LoadFromFile(graphFile);
	std::cerr << "build index" << std::endl;
	auto index = buildBruteForcePathIndex(graph, k);

	std::cerr << "load reads" << std::endl;
	auto reads = loadFastqFromFile(readFile);
	std::vector<std::pair<std::string, size_t>> readLabels;
	readLabels.reserve(reads.size());
	std::vector<vg::Alignment> seeds;
	size_t numSeeds = 0;
	std::cerr << "count seeds" << std::endl;
	for (auto read : reads)
	{
		if (read.sequence.size() < k) continue;
		size_t label = 0;
		for (int i = 0; i < k; i++)
		{
			label <<= 2;
			switch(read.sequence[i])
			{
				case 'A':
				case 'a':
					label += 0;
					break;
				case 'C':
				case 'c':
					label += 1;
					break;
				case 'T':
				case 't':
					label += 2;
					break;
				case 'G':
				case 'g':
					label += 3;
					break;
				default:
					break;
			}
		}
		readLabels.emplace_back(read.seq_id, label);
		numSeeds += index[label].size();
	}
	std::cerr << numSeeds << " seeds" << std::endl;
	seeds.reserve(numSeeds);
	std::cerr << "get seeds" << std::endl;
	for (auto pair : readLabels)
	{
		for (auto pos : index[pair.second])
		{
			vg::Alignment seed;
			seed.set_name(pair.first);
			auto mapping = seed.mutable_path()->add_mapping();
			auto edit = mapping->add_edit();
			edit->set_from_length(k);
			edit->set_to_length(k);
			mapping->mutable_position()->set_node_id(pos.first);
			mapping->mutable_position()->set_offset(pos.second);
			seeds.push_back(seed);
			seed.set_query_position(k-1);
		}
	}

	std::cerr << "write seeds" << std::endl;
	std::ofstream outFile { outputSeedFile, std::ios::binary };
	stream::write_buffered(outFile, seeds, 0);
}