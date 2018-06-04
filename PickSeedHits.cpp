#include <queue>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include "vg.pb.h"
#include "stream.hpp"
#include "CommonUtils.h"

bool seedHitExists(const vg::Alignment& seedHit, const std::vector<vg::Alignment>& existing)
{
	return std::any_of(existing.begin(), existing.end(), [&seedHit](auto& alignment) {
		return alignment.path().mapping(0).position().node_id() == seedHit.path().mapping(0).position().node_id() && alignment.query_position() == seedHit.query_position();
	});
}

int seedLength(const vg::Alignment& seedHit)
{
	if (seedHit.path().mapping_size() == 0) return 0;
	if (seedHit.path().mapping(0).edit_size() == 0) return 0;
	return seedHit.path().mapping(0).edit(0).from_length();
}

class AlignmentLengthCompare
{
public:
	bool operator()(const vg::Alignment& left, const vg::Alignment& right) const
	{
		return seedLength(left) > seedLength(right);
	}
};

int main(int argc, char** argv)
{
	std::string outputFileName { argv[1] };
	int maxSeeds = std::stoi(argv[2]);
	std::unordered_map<std::string, std::priority_queue<vg::Alignment, std::vector<vg::Alignment>, AlignmentLengthCompare>> alignments;
	size_t numElems = 0;
	for (int i = 3; i < argc; i++)
	{
		std::ifstream graphfile { argv[i], std::ios::in | std::ios::binary };
		std::function<void(vg::Alignment&)> lambda = [&alignments, maxSeeds, &numElems](vg::Alignment& g) {
			if (alignments[g.name()].size() < maxSeeds)
			{
				alignments[g.name()].emplace(g);
				numElems++;
				return;
			}
			if (AlignmentLengthCompare{}(g, alignments[g.name()].top()))
			{
				alignments[g.name()].pop();
				alignments[g.name()].emplace(g);
				return;
			}
		};
		stream::for_each(graphfile, lambda);
	}
	std::vector<vg::Alignment> writeAlignments;
	writeAlignments.reserve(numElems);
	for (auto& pair : alignments)
	{
		std::vector<vg::Alignment> insertAlns;
		insertAlns.reserve(pair.second.size());
		while (pair.second.size() > 0)
		{
			insertAlns.push_back(pair.second.top());
			pair.second.pop();
		}
		std::reverse(insertAlns.begin(), insertAlns.end());
		writeAlignments.insert(writeAlignments.end(), insertAlns.begin(), insertAlns.end());
	}
	assert(writeAlignments.size() == numElems);
	std::ofstream alignmentOut { outputFileName, std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, writeAlignments, 0);
}