#include <fstream>
#include <algorithm>
#include <set>
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

double idendityPercent(std::tuple<int, int, int> result)
{
	return (double)std::get<0>(result) / (double)(std::get<0>(result) + std::get<1>(result) + std::get<2>(result));
}

std::tuple<int, int, int> alignmentIdentity(vg::Alignment real, vg::Alignment predicted, const std::map<int, int>& nodeSizes)
{
	std::set<int> leftNodes;
	std::set<int> rightNodes;
	for (int i = 0; i < real.path().mapping_size(); i++)
	{
		leftNodes.insert(real.path().mapping(i).position().node_id());
	}
	for (int i = 0; i < predicted.path().mapping_size(); i++)
	{
		rightNodes.insert(predicted.path().mapping(i).position().node_id());
	}
	std::vector<int> intersection;
	std::set_intersection(leftNodes.begin(), leftNodes.end(), rightNodes.begin(), rightNodes.end(), std::insert_iterator<std::vector<int>>(intersection, intersection.end()));
	int commonBP = 0;
	for (size_t i = 0; i < intersection.size(); i++)
	{
		commonBP += nodeSizes.at(intersection[i]);
	}
	int falseNegativeBP = -commonBP;
	for (int i = 0; i < real.path().mapping_size(); i++)
	{
		falseNegativeBP += nodeSizes.at(real.path().mapping(i).position().node_id());
	}
	int falsePositiveBP = -commonBP;
	for (int i = 0; i < predicted.path().mapping_size(); i++)
	{
		falsePositiveBP += nodeSizes.at(predicted.path().mapping(i).position().node_id());
	}
	auto result = std::make_tuple(commonBP, falseNegativeBP, falsePositiveBP);
	std::cout << real.name() << ": " << commonBP << "bp common, " << falseNegativeBP << "bp false negative, " << falsePositiveBP << "bp false positive (" << idendityPercent(result) << ") " << predicted.score() << " mismatches, read length " << predicted.sequence().size() << " (" << ((double)predicted.score() / (double)predicted.sequence().size()) << ")" << std::endl;
	return result;
}

int main(int argc, char** argv)
{
	std::map<std::string, vg::Alignment> real;
	std::map<std::string, vg::Alignment> predicted;
	std::map<int, int> nodeSizes;
	{
		std::ifstream graphfile { argv[3], std::ios::in | std::ios::binary };
		std::vector<vg::Graph> parts;
		std::function<void(vg::Graph&)> lambda = [&parts](vg::Graph& g) {
			parts.push_back(g);
		};
		stream::for_each(graphfile, lambda);
		auto graph = mergeGraphs(parts);
		for (int i = 0; i < graph.node_size(); i++)
		{
			nodeSizes[graph.node(i).id()] = graph.node(i).sequence().size();
		}
	}

	{
		std::ifstream truthfile { argv[1], std::ios::in | std::ios::binary };
		std::function<void(vg::Alignment&)> lambda = [&real](vg::Alignment& g) {
			real[g.name()] = g;
		};
		stream::for_each(truthfile, lambda);
	}

	{
		std::ifstream predictedfile { argv[2], std::ios::in | std::ios::binary };
		std::function<void(vg::Alignment&)> lambda2 = [&predicted](vg::Alignment& g) {
			predicted[g.name()] = g;
		};
		stream::for_each(predictedfile, lambda2);
	}

	int goodMatches = 0;
	int badMatches = 0;
	for (auto x : real)
	{
		if (predicted.count(x.first) == 0)
		{
			badMatches++;
			continue;
		}
		auto match = alignmentIdentity(x.second, predicted[x.first], nodeSizes);
		if (idendityPercent(match) < 0.7)
		{
			badMatches++;
		}
		else
		{
			goodMatches++;
		}
	}
	for (auto x : predicted)
	{
		if (real.count(x.first) == 0) badMatches++;
	}
	std::cout << "good matches: " << goodMatches << std::endl;
	std::cout << "bad matches: " << badMatches << std::endl;
}