#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
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
		case 'u':
		case 'U':
		case 't':
		case 'T':
			return 't';
		default:
		case '`':
			return '`';
	}
	assert(false);
	return std::numeric_limits<char>::max();
}

bool fileExists(const std::string& fileName)
{
	std::ifstream file { fileName };
	return file.good();
}

MummerSeeder::MummerSeeder(const GfaGraph& graph, const std::string& cachePrefix)
{
	if (cachePrefix.size() > 0 && fileExists(cachePrefix + ".aux"))
	{
		loadFrom(cachePrefix);
	}
	else
	{
		initTree(graph);
		if (cachePrefix.size() > 0) saveTo(cachePrefix);
	}
}

MummerSeeder::MummerSeeder(const vg::Graph& graph, const std::string& cachePrefix)
{
	if (cachePrefix.size() > 0 && fileExists(cachePrefix + ".aux"))
	{
		loadFrom(cachePrefix);
	}
	else
	{
		initTree(graph);
		if (cachePrefix.size() > 0) saveTo(cachePrefix);
	}
}

void MummerSeeder::initTree(const GfaGraph& graph)
{
	for (auto node : graph.nodes)
	{
		nodePositions.push_back(seq.size());
		nodeIDs.push_back(node.first);
		seq += node.second;
		seq += '`';
	}
	seq.pop_back();
	nodePositions.push_back(seq.size());
	for (size_t i = 0; i < seq.size(); i++)
	{
		seq[i] = lowercase(seq[i]);
	}
	seq.shrink_to_fit();
	matcher = std::make_unique<mummer::mummer::sparseSA>(mummer::mummer::sparseSA::create_auto(seq.c_str(), seq.size(), 0, true));
}

void MummerSeeder::initTree(const vg::Graph& graph)
{
	for (int i = 0; i < graph.node_size(); i++)
	{
		nodePositions.push_back(seq.size());
		nodeIDs.push_back(graph.node(i).id());
		seq += graph.node(i).sequence();
		seq += '`';
	}
	seq.pop_back();
	nodePositions.push_back(seq.size());
	for (size_t i = 0; i < seq.size(); i++)
	{
		seq[i] = lowercase(seq[i]);
	}
	seq.shrink_to_fit();
	matcher = std::make_unique<mummer::mummer::sparseSA>(mummer::mummer::sparseSA::create_auto(seq.c_str(), seq.size(), 0, true));
}

size_t MummerSeeder::getNodeIndex(size_t indexPos) const
{
	auto next = std::upper_bound(nodePositions.begin(), nodePositions.end(), indexPos);
	assert(next != nodePositions.begin());
	size_t index = (next - nodePositions.begin()) - 1;
	assert(index < nodePositions.size()-1);
	return index;
}

void MummerSeeder::saveTo(const std::string& prefix) const
{
	std::ofstream file { prefix + ".aux", std::ios::binary };
	{
		boost::archive::text_oarchive oa(file);
		oa << seq;
		oa << nodePositions;
		oa << nodeIDs;
	}
	matcher->save(prefix + "_index");
}

void MummerSeeder::loadFrom(const std::string& prefix)
{
	std::ifstream file { prefix + ".aux", std::ios::binary };
	{
		boost::archive::text_iarchive ia(file);
		ia >> seq;
		ia >> nodePositions;
		ia >> nodeIDs;
	}
	// same params that create_auto with minlen=0 passes
	matcher = std::make_unique<mummer::mummer::sparseSA>(seq, false, 1, true, false, false, 1, 0, true);
	matcher->load(prefix + "_index");
}

std::vector<SeedHit> MummerSeeder::getMumSeeds(std::string sequence, size_t maxCount, size_t minLen) const
{
	for (size_t i = 0; i < sequence.size(); i++)
	{
		sequence[i] = lowercase(sequence[i]);
	}
	assert(matcher != nullptr);
	std::vector<mummer::mummer::match_t> MAMs;
	matcher->MAM(sequence, minLen, false, MAMs);
	revcompInPlace(sequence);
	std::vector<mummer::mummer::match_t> bwMAMs;
	matcher->MAM(sequence, minLen, false, bwMAMs);
	auto seeds = matchesToSeeds(sequence.size(), MAMs, bwMAMs);
	std::sort(seeds.begin(), seeds.end(), [](const SeedHit& left, const SeedHit& right) { return left.matchLen > right.matchLen; });
	if (seeds.size() > maxCount)
	{
		seeds.erase(seeds.begin() + maxCount, seeds.end());
	}
	return seeds;
}

std::vector<SeedHit> MummerSeeder::getMemSeeds(std::string sequence, size_t maxCount, size_t minLen) const
{
	for (size_t i = 0; i < sequence.size(); i++)
	{
		sequence[i] = lowercase(sequence[i]);
	}
	assert(matcher != nullptr);
	std::vector<mummer::mummer::match_t> MEMs;
	matcher->MEM(sequence, minLen, false, MEMs);
	revcompInPlace(sequence);
	std::vector<mummer::mummer::match_t> bwMEMs;
	matcher->MEM(sequence, minLen, false, bwMEMs);
	auto seeds = matchesToSeeds(sequence.size(), MEMs, bwMEMs);
	std::sort(seeds.begin(), seeds.end(), [](const SeedHit& left, const SeedHit& right) { return left.matchLen > right.matchLen; });
	if (seeds.size() > maxCount)
	{
		seeds.erase(seeds.begin() + maxCount, seeds.end());
	}
	return seeds;
}

std::vector<SeedHit> MummerSeeder::matchesToSeeds(size_t seqLen, const std::vector<mummer::mummer::match_t>& fwmatches, const std::vector<mummer::mummer::match_t>& bwmatches) const
{
	std::vector<SeedHit> result;
	result.reserve(fwmatches.size() + bwmatches.size());
	for (auto match : fwmatches)
	{
		auto index = getNodeIndex(match.ref);
		int nodeID = nodeIDs[index];
		size_t nodeOffset = match.ref - nodePositions[index];
		size_t seqPos = match.query;
		size_t matchLen = match.len;
		result.emplace_back(nodeID, nodeOffset, seqPos, matchLen, false);
	}
	for (auto match : bwmatches)
	{
		auto index = getNodeIndex(match.ref);
		int nodeID = nodeIDs[index];
		size_t nodeOffset = match.ref - nodePositions[index];
		size_t seqPos = match.query;
		size_t matchLen = match.len;
		assert(match.len > 0);
		assert(nodeOffset + matchLen <= nodeLength(index));
		assert(seqPos + matchLen <= seqLen);
		nodeOffset = nodeLength(index) - nodeOffset - matchLen;
		seqPos = seqLen - seqPos - matchLen;
		assert(nodeOffset < nodeLength(index));
		assert(seqPos < seqLen);
		result.emplace_back(nodeID, nodeOffset, seqPos, matchLen, true);
	}
	return result;
}

size_t MummerSeeder::nodeLength(size_t indexPos) const
{
	//-1 for separator
	return nodePositions[indexPos+1] - nodePositions[indexPos] - 1;
}

void MummerSeeder::revcompInPlace(std::string& seq) const
{
	std::reverse(seq.begin(), seq.end());
	for (size_t i = 0; i < seq.size(); i++)
	{
		switch(seq[i])
		{
			case 'a':
				seq[i] = 't';
				break;
			case 'u':
			case 't':
				seq[i] = 'a';
				break;
			case 'c':
				seq[i] = 'g';
				break;
			case 'g':
				seq[i] = 'c';
				break;
			default:
				assert(false);
		}
	}
}
