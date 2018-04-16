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

int main(int argc, char** argv)
{
	std::string outputFileName { argv[1] };
	int maxSeeds = std::stoi(argv[2]);
	std::unordered_map<std::string, std::vector<vg::Alignment>> alignments;
	for (int i = 3; i < argc; i++)
	{
		auto fileAlns = CommonUtils::LoadVGAlignments(argv[i]);
		for (auto a : fileAlns)
		{
			alignments[a.name()].push_back(a);
		}
	}
	std::vector<vg::Alignment> writeAlignments;
	for (auto& pair : alignments)
	{
		std::stable_sort(pair.second.begin(), pair.second.end(), [](const vg::Alignment& left, const vg::Alignment& right) { return seedLength(left) > seedLength(right); });
		std::vector<vg::Alignment> result;
		for (int i = 0; result.size() < maxSeeds && i < pair.second.size(); i++)
		{
			if (!seedHitExists(pair.second[i], result)) result.push_back(pair.second[i]);
		}
		writeAlignments.insert(writeAlignments.end(), result.begin(), result.end());
	}
	std::ofstream alignmentOut { outputFileName, std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, writeAlignments, 0);
}