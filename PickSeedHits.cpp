#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include "vg.pb.h"
#include "stream.hpp"

bool seedHitExists(const vg::Alignment& seedHit, const std::vector<vg::Alignment>& existing)
{
	return std::any_of(existing.begin(), existing.end(), [&seedHit](auto& alignment) {
		return alignment.path().mapping(0).position().node_id() == seedHit.path().mapping(0).position().node_id() && alignment.query_position() == seedHit.query_position();
	});
}

int main(int argc, char** argv)
{
	std::map<std::string, std::vector<vg::Alignment>> alignments;
	int maxseeds = std::stoi(argv[2]);
	for (int i = 3; i < argc; i++)
	{
		std::ifstream seedfile { argv[i], std::ios::in | std::ios::binary };
		std::function<void(vg::Alignment&)> alignmentLambda = [&alignments, maxseeds](vg::Alignment& a) {
			if (a.path().mapping(0).position().node_id() <= 1) return;
			if (seedHitExists(a, alignments[a.name()])) return;
			if (alignments[a.name()].size() < maxseeds) alignments[a.name()].push_back(a);
		};
		stream::for_each(seedfile, alignmentLambda);
	}
	std::vector<vg::Alignment> writeAlignments;
	for (auto pair : alignments)
	{
		writeAlignments.insert(writeAlignments.end(), pair.second.begin(), pair.second.end());
	}
	std::ofstream alignmentOut { argv[1], std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, writeAlignments, 0);
}