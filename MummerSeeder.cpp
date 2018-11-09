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
		case 't':
		case 'T':
			return 't';
	}
	return std::numeric_limits<char>::max();
}

bool fileExists(const std::string& fileName)
{
	std::ifstream file { fileName };
	return file.good();
}

MummerSeeder::MummerSeeder(const GfaGraph& graph, size_t minL, const std::string& cachePrefix)
{
	if (cachePrefix.size() > 0 && fileExists(cachePrefix + ".aux"))
	{
		loadFrom(cachePrefix);
	}
	else
	{
		initTree(graph, minL);
		minLen = minL;
		if (cachePrefix.size() > 0) saveTo(cachePrefix);
	}
}

MummerSeeder::MummerSeeder(const vg::Graph& graph, size_t minL, const std::string& cachePrefix)
{
	if (cachePrefix.size() > 0 && fileExists(cachePrefix + ".aux"))
	{
		loadFrom(cachePrefix);
	}
	else
	{
		initTree(graph, minL);
		minLen = minL;
		if (cachePrefix.size() > 0) saveTo(cachePrefix);
	}
}

void MummerSeeder::initTree(const GfaGraph& graph, size_t minLen)
{
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
	seq.shrink_to_fit();
	matcher = std::make_unique<mummer::mummer::sparseSA>(mummer::mummer::sparseSA::create_auto(seq.c_str(), seq.size(), 1, true));
}

void MummerSeeder::initTree(const vg::Graph& graph, size_t minLen)
{
	for (int i = 0; i < graph.node_size(); i++)
	{
		nodePositions.push_back(seq.size());
		nodeIDs.push_back(graph.node(i).id());
		nodeReverse.push_back(false);
		seq += graph.node(i).sequence();
		seq += '$';
		nodePositions.push_back(seq.size());
		nodeIDs.push_back(graph.node(i).id());
		nodeReverse.push_back(true);
		seq += CommonUtils::ReverseComplement(graph.node(i).sequence());
		seq += '$';
	}
	nodePositions.push_back(seq.size());
	for (size_t i = 0; i < seq.size(); i++)
	{
		seq[i] = lowercase(seq[i]);
	}
	seq.shrink_to_fit();
	matcher = std::make_unique<mummer::mummer::sparseSA>(mummer::mummer::sparseSA::create_auto(seq.c_str(), seq.size(), 1, true));
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
		oa << minLen;
		oa << seq;
		oa << nodePositions;
		oa << nodeIDs;
		oa << nodeReverse;
	}
	matcher->save(prefix + "_index");
}

void MummerSeeder::loadFrom(const std::string& prefix)
{
	std::ifstream file { prefix + ".aux", std::ios::binary };
	{
		boost::archive::text_iarchive ia(file);
		ia >> minLen;
		ia >> seq;
		ia >> nodePositions;
		ia >> nodeIDs;
		ia >> nodeReverse;
	}
	size_t kmer = minLen;
	if (kmer > 10) kmer = 10;
	// same params that create_auto passes
	matcher = std::make_unique<mummer::mummer::sparseSA>(seq, false, 1, true, false, kmer>0, 1, kmer, true);
	matcher->load(prefix + "_index");
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
	std::sort(MAMs.begin(), MAMs.end(), [](const mummer::mummer::match_t& left, const mummer::mummer::match_t& right) { return left.len > right.len; });
	if (MAMs.size() > maxCount)
	{
		MAMs.erase(MAMs.begin() + maxCount, MAMs.end());
	}
	return matchesToSeeds(MAMs);
}

std::vector<SeedHit> MummerSeeder::getMemSeeds(std::string sequence, size_t maxCount) const
{
	for (size_t i = 0; i < sequence.size(); i++)
	{
		sequence[i] = lowercase(sequence[i]);
	}
	assert(matcher != nullptr);
	std::vector<mummer::mummer::match_t> MEMs;
	matcher->MEM(sequence, minLen, false, MEMs);
	std::sort(MEMs.begin(), MEMs.end(), [](const mummer::mummer::match_t& left, const mummer::mummer::match_t& right) { return left.len > right.len; });
	if (MEMs.size() > maxCount)
	{
		MEMs.erase(MEMs.begin() + maxCount, MEMs.end());
	}
	return matchesToSeeds(MEMs);
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
