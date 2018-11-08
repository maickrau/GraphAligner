#include "CommonUtils.h"
#include "MummerSeeder.h"

char lowercase(char c)
{
	switch(c)
	{
		case 'a':
		case 'A':
			return 'a';
		case 'c':
		case 'C':
			return 'c';
		case 'g':
		case 'G':
			return 'g';
		case 't':
		case 'T':
			return 't';
	}
	return std::numeric_limits<char>::max();
}

MummerSeeder::MummerSeeder(const GfaGraph& graph, size_t minL)
{
	constructTree(graph, minL);
	minLen = minL;
}

void MummerSeeder::constructTree(const GfaGraph& graph, size_t minLen)
{
	std::string seq;
	for (auto node : graph.nodes)
	{
		nodePositions.push_back(seq.size());
		nodeIDs.push_back(node.first);
		nodeReverse.push_back(false);
		seq += node.second;
		seq += '$';
		nodePositions.push_back(seq.size());
		nodeIDs.push_back(node.first);
		nodeReverse.push_back(true);
		seq += CommonUtils::ReverseComplement(node.second);
		seq += '$';
	}
	nodePositions.push_back(seq.size());
	for (size_t i = 0; i < seq.size(); i++)
	{
		seq[i] = lowercase(seq[i]);
	}
	matcher = std::make_unique<mummer::mummer::sparseSA>(mummer::mummer::sparseSA::create_auto(seq.c_str(), seq.size(), 1, true));
	matcher->construct();
}

size_t MummerSeeder::getNodeIndex(size_t indexPos) const
{
	auto next = std::upper_bound(nodePositions.begin(), nodePositions.end(), indexPos);
	assert(next != nodePositions.begin());
	size_t index = (next - nodePositions.begin()) - 1;
	assert(index < nodePositions.size()-1);
	return index;
}

std::vector<SeedHit> MummerSeeder::getMumSeeds(std::string sequence, size_t maxCount) const
{
	for (size_t i = 0; i < sequence.size(); i++)
	{
		sequence[i] = lowercase(sequence[i]);
	}
	assert(matcher != nullptr);
	std::vector<mummer::mummer::match_t> MAMs;
	matcher->MAM(sequence, minLen, false, MAMs);
	std::cerr << MAMs.size() << " seeds" << std::endl;
	std::sort(MAMs.begin(), MAMs.end(), [](const mummer::mummer::match_t& left, const mummer::mummer::match_t& right) { return left.len > right.len; });
	if (MAMs.size() > maxCount)
	{
		MAMs.erase(MAMs.begin() + maxCount, MAMs.end());
	}
	std::cerr << MAMs.size() << " seeds" << std::endl;
	return matchesToSeeds(MAMs);
}

std::vector<SeedHit> MummerSeeder::matchesToSeeds(const std::vector<mummer::mummer::match_t>& matches) const
{
	std::vector<SeedHit> result;
	result.reserve(matches.size());
	for (auto match : matches)
	{
		auto index = getNodeIndex(match.ref);
		int nodeID = nodeIDs[index];
		size_t nodeOffset = match.ref - nodePositions[index];
		size_t seqPos = match.query;
		size_t matchLen = match.len;
		bool reverse = nodeReverse[index];
		result.emplace_back(nodeID, nodeOffset, seqPos, matchLen, reverse);
	}
	return result;
}
