#include <chrono>
#include <iostream>
#include <random>
#include <fstream>
#include "CommonUtils.h"
#include "vg.pb.h"
#include "stream.hpp"
#include "GfaGraph.h"

std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0,1.0);

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

bool is_file_exist(std::string fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}

std::tuple<vg::Alignment, std::string, vg::Alignment> simulateOneRead(const std::vector<std::string>& nodeSequences, const std::vector<int>& nodeIds, int overlap, int length, double substitutionErrorRate, double insertionErrorRate, double deletionErrorRate, const std::map<size_t, std::vector<std::pair<size_t, bool>>>& outEdgesRight, const std::map<size_t, std::vector<std::pair<size_t, bool>>>& outEdgesLeft)
{
	bool reverse = false;
	if (distribution(generator) < 0.5) reverse = true;

	std::vector<std::pair<int, bool>> realNodes;
	std::vector<size_t> nodelens;

	int currentNode = rand() % nodeSequences.size();
	int startNode = nodeIds[currentNode];
	assert(nodeSequences[currentNode].size() > overlap);
	int startPos = rand() % (nodeSequences[currentNode].size() - overlap);
	bool startReverse = reverse;
	std::string realsequence;
	if (reverse)
	{
		realsequence = CommonUtils::ReverseComplement(nodeSequences[currentNode]).substr(startPos);
	}
	else
	{
		realsequence = nodeSequences[currentNode].substr(startPos);
	}
	assert(realsequence.size() > overlap);
	realsequence.erase(realsequence.end()-overlap, realsequence.end());
	while (realsequence.size() < length)
	{
		if (currentNode == 0) return simulateOneRead(nodeSequences, nodeIds, overlap, length, substitutionErrorRate, insertionErrorRate, deletionErrorRate, outEdgesRight, outEdgesLeft);
		realNodes.emplace_back(nodeIds[currentNode], reverse);
		if (nodelens.size() == 0)
		{
			nodelens.emplace_back(nodeSequences[currentNode].size() - overlap - startPos);
		}
		else
		{
			nodelens.emplace_back(nodeSequences[currentNode].size() - overlap);
		}
		if (reverse)
		{
			if (outEdgesLeft.count(currentNode) == 0) return simulateOneRead(nodeSequences, nodeIds, overlap, length, substitutionErrorRate, insertionErrorRate, deletionErrorRate, outEdgesRight, outEdgesLeft);
			if (outEdgesLeft.at(currentNode).size() == 0) return simulateOneRead(nodeSequences, nodeIds, overlap, length, substitutionErrorRate, insertionErrorRate, deletionErrorRate, outEdgesRight, outEdgesLeft);
			auto pickedIndex = rand() % outEdgesLeft.at(currentNode).size();
			reverse = !outEdgesLeft.at(currentNode)[pickedIndex].second;
			currentNode = outEdgesLeft.at(currentNode)[pickedIndex].first;
		}
		else
		{
			if (outEdgesRight.count(currentNode) == 0) return simulateOneRead(nodeSequences, nodeIds, overlap, length, substitutionErrorRate, insertionErrorRate, deletionErrorRate, outEdgesRight, outEdgesLeft);
			if (outEdgesRight.at(currentNode).size() == 0) return simulateOneRead(nodeSequences, nodeIds, overlap, length, substitutionErrorRate, insertionErrorRate, deletionErrorRate, outEdgesRight, outEdgesLeft);
			auto pickedIndex = rand() % outEdgesRight.at(currentNode).size();
			reverse = !outEdgesRight.at(currentNode)[pickedIndex].second;
			currentNode = outEdgesRight.at(currentNode)[pickedIndex].first;
		}
		if (reverse)
		{
			realsequence += CommonUtils::ReverseComplement(nodeSequences[currentNode]);
		}
		else
		{
			realsequence += nodeSequences[currentNode];
		}
		assert(realsequence.size() > overlap);
		realsequence.erase(realsequence.end()-overlap, realsequence.end());
	}
	realNodes.emplace_back(nodeIds[currentNode], reverse);
	nodelens.emplace_back(nodeSequences[currentNode].size() - overlap);
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
		auto edit = mapping->add_edit();
		edit->set_from_length(nodelens[i]);
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
	seedposition->set_is_reverse(startReverse);

	return std::make_tuple(result, errorSequence, seed);
}

int main(int argc, char** argv)
{
	std::string graphFile {argv[1]};
	std::string alignmentOutFile {argv[2]};
	std::string fastqOutFile {argv[3]};
	int numReads = std::stoi(argv[4]);
	int length = std::stoi(argv[5]);
	double substitution = std::stod(argv[6]);
	double insertions = std::stod(argv[7]);
	std::string seedsOutFile {argv[8]};
	double deletions = std::stod(argv[9]);

	// generator.seed(std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1));
	generator.seed(0);
	srand(0);

	if (is_file_exist(graphFile)){
		std::cout << "load graph from " << graphFile << std::endl;
	}
	else{
		std::cout << "No graph file exists" << std::endl;
		std::exit(0);
	}
	std::vector<std::string> nodeSequences;
	std::vector<int> nodeIds;
	int overlap = 0;
	std::map<size_t, std::vector<std::pair<size_t, bool>>> outEdgesRight;
	std::map<size_t, std::vector<std::pair<size_t, bool>>> outEdgesLeft;

	if (graphFile.substr(graphFile.size()-3) == ".vg")
	{
		vg::Graph graph = CommonUtils::LoadVGGraph(graphFile);
		std::map<int, size_t> ids;
		for (int i = 0; i < graph.node_size(); i++)
		{
			nodeSequences.push_back(graph.node(i).sequence());
			nodeIds.push_back(graph.node(i).id());
			ids[graph.node(i).id()] = i;
		}
		for (int i = 0; i < graph.edge_size(); i++)
		{
			if (graph.edge(i).from_start())
			{
				bool direction = !graph.edge(i).to_end();
				outEdgesLeft[ids[graph.edge(i).from()]].emplace_back(ids[graph.edge(i).to()], direction);
			}
			else
			{
				bool direction = !graph.edge(i).to_end();
				outEdgesRight[ids[graph.edge(i).from()]].emplace_back(ids[graph.edge(i).to()], direction);
			}
			if (graph.edge(i).to_end())
			{
				bool direction = !graph.edge(i).from_start();
				outEdgesRight[ids[graph.edge(i).to()]].emplace_back(ids[graph.edge(i).from()], !direction);
			}
			else
			{
				bool direction = !graph.edge(i).from_start();
				outEdgesLeft[ids[graph.edge(i).to()]].emplace_back(ids[graph.edge(i).from()], !direction);
			}
		}
	}
	else
	{
		GfaGraph graph = GfaGraph::LoadFromFile(graphFile);
		std::map<int, size_t> ids;
		size_t nodenum = 0;
		overlap = graph.edgeOverlap;
		for (const auto& pair : graph.nodes)
		{
			nodeSequences.push_back(pair.second);
			nodeIds.push_back(pair.first);
			ids[pair.first] = nodenum;
			nodenum++;
		}
		for (auto edge : graph.edges)
		{
			auto from = edge.first;
			for (auto to : edge.second)
			{
				if (from.end)
				{
					bool direction = to.end;
					outEdgesRight[ids[from.id]].emplace_back(ids[to.id], direction);
				}
				else
				{
					bool direction = to.end;
					outEdgesLeft[ids[from.id]].emplace_back(ids[to.id], direction);
				}
				if (to.end)
				{
					bool direction = from.end;
					outEdgesLeft[ids[to.id]].emplace_back(ids[from.id], !direction);
				}
				else
				{
					bool direction = from.end;
					outEdgesRight[ids[to.id]].emplace_back(ids[from.id], !direction);
				}
			}
		}
	}

	std::vector<std::tuple<vg::Alignment, std::string, vg::Alignment>> reads;
	std::vector<vg::Alignment> truth;
	std::vector<vg::Alignment> seeds;
	for (int i = 0; i < numReads; i++)
	{
		reads.push_back(simulateOneRead(nodeSequences, nodeIds, overlap, length, substitution, insertions, deletions, outEdgesRight, outEdgesLeft));
		truth.emplace_back(std::get<0>(reads[i]));
		seeds.emplace_back(std::get<2>(reads[i]));
	}

	std::ofstream alignmentOut { alignmentOutFile, std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, truth, 0);

	std::ofstream seedsOut { seedsOutFile, std::ios::out | std::ios::binary };
	stream::write_buffered(seedsOut, seeds, 0);

	std::ofstream fastqOut {fastqOutFile};
	for (int i = 0; i < numReads; i++)
	{
		fastqOut << "@" << std::get<0>(reads[i]).name() << "\n";
		fastqOut << std::get<1>(reads[i]) << "\n";
		fastqOut << "+" << "\n";
		fastqOut << std::string(std::get<1>(reads[i]).size(), '!') << "\n";
	}

}