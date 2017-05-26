#include <fstream>
#include <algorithm>
#include <set>
#include "vg.pb.h"
#include "stream.hpp"

double idendityPercent(std::tuple<int, int, int> result)
{
	return (double)std::get<0>(result) / (double)(std::get<0>(result) + std::get<1>(result) + std::get<2>(result));
}

std::tuple<int, int, int> alignmentIdentity(vg::Alignment real, vg::Alignment predicted)
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
	int common = intersection.size();
	int falseNegative = real.path().mapping_size() - intersection.size();
	int falsePositive = predicted.path().mapping_size() - intersection.size();
	auto result = std::make_tuple(common, falseNegative, falsePositive);
	std::cout << real.name() << ": " << common << " common, " << falseNegative << " false negative, " << falsePositive << " false positive (" << idendityPercent(result) << ") " << predicted.score() << std::endl;
	return result;
}

int main(int argc, char** argv)
{
	std::ifstream graphfile { argv[1], std::ios::in | std::ios::binary };
	std::map<std::string, vg::Alignment> real;
	std::map<std::string, vg::Alignment> predicted;
	std::function<void(vg::Alignment&)> lambda = [&real](vg::Alignment& g) {
		real[g.name()] = g;
	};
	stream::for_each(graphfile, lambda);

	std::ifstream graphfile2 { argv[2], std::ios::in | std::ios::binary };
	std::function<void(vg::Alignment&)> lambda2 = [&predicted](vg::Alignment& g) {
		predicted[g.name()] = g;
	};
	stream::for_each(graphfile2, lambda2);

	int goodMatches = 0;
	int badMatches = 0;
	for (auto x : real)
	{
		if (predicted.count(x.first) == 0)
		{
			badMatches++;
			continue;
		}
		auto match = alignmentIdentity(x.second, predicted[x.first]);
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