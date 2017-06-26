#include <chrono>
#include <iostream>
#include <random>
#include <fstream>
#include "vg.pb.h"
#include "stream.hpp"

std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0,1.0);

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

std::string introduceErrors(std::string real, double substitutionErrorRate, double insertionErrorRate, double deletionErrorRate)
{
	std::string result;
	for (size_t i = 0; i < real.size(); i++)
	{
		if (distribution(generator) < deletionErrorRate)
		{
		}
		else
		{
			if (distribution(generator) < substitutionErrorRate)
			{
				result += "ATCG"[rand() % 4];
			}
			else
			{
				result += real[i];
			}
		}
		if (distribution(generator) < insertionErrorRate / 10.0)
		{
			int length = rand() % 20;
			for (int j = 0; j < length; j++)
			{
				result += "ATCG"[rand() % 4];
			}
		}
	}
	return result;
}

bool is_file_exist(const char *fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}

std::string reverseComplement(std::string str)
{
	std::string result;
	for (int i = str.size()-1; i >= 0; i--)
	{
		switch (str[i])
		{
			case 'A':
			case 'a':
			result += 'T';
			break;
			case 'C':
			case 'c':
			result += 'G';
			break;
			case 'T':
			case 't':
			result += 'A';
			break;
			case 'G':
			case 'g':
			result += 'C';
			break;
			case 'N':
			case 'n':
			result += 'N';
			break;
		}
	}
	return result;
}

std::tuple<vg::Alignment, std::string, vg::Alignment> simulateOneRead(const vg::Graph& g, int length, double substitutionErrorRate, double insertionErrorRate, double deletionErrorRate, const std::map<size_t, std::vector<std::pair<size_t, bool>>>& outEdgesRight, const std::map<size_t, std::vector<std::pair<size_t, bool>>>& outEdgesLeft)
{
	bool reverse = false;
	if (distribution(generator) < 0.5) reverse = true;

	std::vector<std::pair<int, bool>> realNodes;

	int currentNode = rand() % g.node_size();
	int startNode = g.node(currentNode).id();
	int startPos = rand() % g.node(currentNode).sequence().size();
	std::string realsequence;
	if (reverse)
	{
		realsequence = reverseComplement(g.node(currentNode).sequence().substr(0, startPos));	
	}
	else
	{
		realsequence = g.node(currentNode).sequence().substr(startPos);	
	}
	while (realsequence.size() < length)
	{
		if (currentNode == 0) return simulateOneRead(g, length, substitutionErrorRate, insertionErrorRate, deletionErrorRate, outEdgesRight, outEdgesLeft);
		realNodes.emplace_back(g.node(currentNode).id(), reverse);
		if (reverse)
		{
			if (outEdgesLeft.count(currentNode) == 0) return simulateOneRead(g, length, substitutionErrorRate, insertionErrorRate, deletionErrorRate, outEdgesRight, outEdgesLeft);
			if (outEdgesLeft.at(currentNode).size() == 0) return simulateOneRead(g, length, substitutionErrorRate, insertionErrorRate, deletionErrorRate, outEdgesRight, outEdgesLeft);
			auto pickedIndex = rand() % outEdgesLeft.at(currentNode).size();
			reverse = outEdgesLeft.at(currentNode)[pickedIndex].second;
			currentNode = outEdgesLeft.at(currentNode)[pickedIndex].first;
		}
		else
		{
			if (outEdgesRight.count(currentNode) == 0) return simulateOneRead(g, length, substitutionErrorRate, insertionErrorRate, deletionErrorRate, outEdgesRight, outEdgesLeft);
			if (outEdgesRight.at(currentNode).size() == 0) return simulateOneRead(g, length, substitutionErrorRate, insertionErrorRate, deletionErrorRate, outEdgesRight, outEdgesLeft);
			auto pickedIndex = rand() % outEdgesRight.at(currentNode).size();
			reverse = outEdgesRight.at(currentNode)[pickedIndex].second;
			currentNode = outEdgesRight.at(currentNode)[pickedIndex].first;
		}
		if (reverse)
		{
			realsequence += reverseComplement(g.node(currentNode).sequence());
		}
		else
		{
			realsequence += g.node(currentNode).sequence();
		}
	}
	realNodes.emplace_back(g.node(currentNode).id(), reverse);
	realsequence = realsequence.substr(0, length);
	auto errorSequence = introduceErrors(realsequence, substitutionErrorRate, insertionErrorRate, deletionErrorRate);

	vg::Alignment result;
	result.set_name("read_" + std::to_string(rand()));
	result.set_sequence(realsequence);
	vg::Path* path = new vg::Path;
	result.set_allocated_path(path);
	for (int i = 0; i < realNodes.size(); i++)
	{
		auto mapping = path->add_mapping();
		auto position = new vg::Position;
		mapping->set_allocated_position(position);
		position->set_node_id(realNodes[i].first);
		position->set_is_reverse(realNodes[i].second);
		if (i == 0) position->set_offset(startPos);
	}

	vg::Alignment seed;
	seed.set_query_position(1);
	seed.set_name(result.name());
	vg::Path* seedpath = new vg::Path;
	seed.set_allocated_path(seedpath);
	auto seedmapping = seedpath->add_mapping();
	auto seedposition = new vg::Position;
	seedmapping->set_allocated_position(seedposition);
	seedposition->set_node_id(startNode);

	return std::make_tuple(result, errorSequence, seed);
}

int main(int argc, char** argv)
{
	generator.seed(std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1));

	vg::Graph graph;
	{
		if (is_file_exist(argv[1])){
			std::cout << "load graph from " << argv[1] << std::endl;
			std::ifstream graphfile { argv[1], std::ios::in | std::ios::binary };
			std::vector<vg::Graph> parts;
			std::function<void(vg::Graph&)> lambda = [&parts](vg::Graph& g) {
				parts.push_back(g);
			};
			stream::for_each(graphfile, lambda);
			graph = mergeGraphs(parts);
		}
		else{
			std::cout << "No graph file exists" << std::endl;
			std::exit(0);
		}
	}

	int numReads = std::stoi(argv[4]);
	int length = std::stoi(argv[5]);
	double substitution = std::stod(argv[6]);
	double insertions = std::stod(argv[7]);
	double deletions = std::stod(argv[9]);

	std::map<int, size_t> ids;
	for (int i = 0; i < graph.node_size(); i++)
	{
		ids[graph.node(i).id()] = i;
	}
	std::map<size_t, std::vector<std::pair<size_t, bool>>> outEdgesRight;
	std::map<size_t, std::vector<std::pair<size_t, bool>>> outEdgesLeft;
	for (int i = 0; i < graph.edge_size(); i++)
	{
		if (graph.edge(i).from_start())
		{
			bool direction = graph.edge(i).to_end();
			outEdgesLeft[ids[graph.edge(i).from()]].emplace_back(ids[graph.edge(i).to()], direction);
		}
		else
		{
			bool direction = graph.edge(i).to_end();
			outEdgesRight[ids[graph.edge(i).from()]].emplace_back(ids[graph.edge(i).to()], direction);
		}
		if (graph.edge(i).to_end())
		{
			bool direction = graph.edge(i).from_start();
			outEdgesRight[ids[graph.edge(i).to()]].emplace_back(ids[graph.edge(i).from()], !direction);
		}
		else
		{
			bool direction = graph.edge(i).from_start();
			outEdgesLeft[ids[graph.edge(i).to()]].emplace_back(ids[graph.edge(i).from()], !direction);
		}
	}

	std::vector<std::tuple<vg::Alignment, std::string, vg::Alignment>> reads;
	std::vector<vg::Alignment> truth;
	std::vector<vg::Alignment> seeds;
	for (int i = 0; i < numReads; i++)
	{
		reads.push_back(simulateOneRead(graph, length, substitution, insertions, deletions, outEdgesRight, outEdgesLeft));
		truth.emplace_back(std::get<0>(reads[i]));
		seeds.emplace_back(std::get<2>(reads[i]));
	}

	std::ofstream alignmentOut { argv[2], std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, truth, 0);

	std::ofstream seedsOut { argv[8], std::ios::out | std::ios::binary };
	stream::write_buffered(seedsOut, seeds, 0);

	std::ofstream fastqOut {argv[3]};
	for (int i = 0; i < numReads; i++)
	{
		fastqOut << "@" << std::get<0>(reads[i]).name() << "\n";
		fastqOut << std::get<1>(reads[i]) << "\n";
		fastqOut << "+" << "\n";
		fastqOut << std::string(std::get<1>(reads[i]).size(), '!') << "\n";
	}

}