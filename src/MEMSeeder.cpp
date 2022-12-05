#include <fstream>
#include "CommonUtils.h"
#include "MEMSeeder.h"
#include "Serialize.h"

bool fileExists(const std::string& fileName)
{
	std::ifstream file { fileName };
	return file.good();
}

MEMSeeder::MEMSeeder(const GfaGraph& graph, const std::string& cachePrefix, const double uniqueBonusFactor, const bool lowMemoryMEMIndexConstruction, const bool useWaveletTree, const size_t windowSize) :
	uniqueBonusFactor(uniqueBonusFactor),
	windowSize(windowSize)
{
	if (cachePrefix.size() > 0 && fileExists(cachePrefix + ".index"))
	{
		loadFrom(cachePrefix);
		for (size_t i = 0; i < nodeIDs.size(); i++)
		{
			assert(nodePositions[i+1] - nodePositions[i] == graph.nodes.at(nodeIDs[i]).size()+1);
		}
	}
	else
	{
		initTree(graph, lowMemoryMEMIndexConstruction, useWaveletTree);
		if (cachePrefix.size() > 0) saveTo(cachePrefix);
	}
}

MEMSeeder::MEMSeeder(const vg::Graph& graph, const std::string& cachePrefix, const double uniqueBonusFactor, const bool lowMemoryMEMIndexConstruction, const bool useWaveletTree, const size_t windowSize) :
	uniqueBonusFactor(uniqueBonusFactor),
	windowSize(windowSize)
{
	if (cachePrefix.size() > 0 && fileExists(cachePrefix + ".index"))
	{
		loadFrom(cachePrefix);
	}
	else
	{
		initTree(graph, lowMemoryMEMIndexConstruction, useWaveletTree);
		if (cachePrefix.size() > 0) saveTo(cachePrefix);
	}
}

void MEMSeeder::saveTo(const std::string& prefix) const
{
	std::ofstream file { prefix + ".index", std::ios::binary };
	index.save(file);
	serialize(file, nodePositions);
	serialize(file, nodeIDs);
}

void MEMSeeder::loadFrom(const std::string& prefix)
{
	std::ifstream file { prefix + ".index", std::ios::binary };
	index.load(file);
	deserialize(file, nodePositions);
	deserialize(file, nodeIDs);
	prefixIndex = MEMfinder::buildPrefixIndex(index, 10);
}

void MEMSeeder::initTree(const GfaGraph& graph, const bool lowMemoryMEMIndexConstruction, const bool useWaveletTree)
{
	std::string seq;
	for (size_t i = 0; i < graph.nodes.size(); i++)
	{
		nodePositions.push_back(seq.size());
		nodeIDs.push_back(i);
		seq += graph.nodes[i].toString();
		seq += '`';
	}
	for (size_t i = 0; i < seq.size(); i++)
	{
		switch(seq[i])
		{
			case 'a':
			case 'A':
				seq[i] = 2;
				break;
			case 'c':
			case 'C':
				seq[i] = 3;
				break;
			case 'g':
			case 'G':
				seq[i] = 4;
				break;
			case 't':
			case 'T':
				seq[i] = 5;
				break;
			default:
				seq[i] = 1;
				break;
		}
	}
	nodePositions.push_back(seq.size());
	seq.push_back(0);
	if (lowMemoryMEMIndexConstruction)
	{
		index.initializeLowMemory(std::move(seq), 16, useWaveletTree);
	}
	else
	{
		index.initialize(std::move(seq), 16, useWaveletTree);
	}
	prefixIndex = MEMfinder::buildPrefixIndex(index, 10);
}

void MEMSeeder::initTree(const vg::Graph& graph, const bool lowMemoryMEMIndexConstruction, const bool useWaveletTree)
{
	std::string seq;
	for (int i = 0; i < graph.node_size(); i++)
	{
		nodePositions.push_back(seq.size());
		nodeIDs.push_back(graph.node(i).id());
		seq += graph.node(i).sequence();
		seq += '`';
	}
	for (size_t i = 0; i < seq.size(); i++)
	{
		switch(seq[i])
		{
			case 'a':
			case 'A':
				seq[i] = 2;
				break;
			case 'c':
			case 'C':
				seq[i] = 3;
				break;
			case 'g':
			case 'G':
				seq[i] = 4;
				break;
			case 't':
			case 'T':
				seq[i] = 5;
				break;
			default:
				seq[i] = 1;
				break;
		}
	}
	nodePositions.push_back(seq.size());
	seq.push_back(0);
	if (lowMemoryMEMIndexConstruction)
	{
		index.initializeLowMemory(std::move(seq), 16, useWaveletTree);
	}
	else
	{
		index.initialize(std::move(seq), 16, useWaveletTree);
	}
	prefixIndex = MEMfinder::buildPrefixIndex(index, 10);
}

size_t MEMSeeder::getNodeIndex(size_t indexPos) const
{
	assert(indexPos < nodePositions.back());
	auto next = std::upper_bound(nodePositions.begin(), nodePositions.end(), indexPos);
	assert(next != nodePositions.begin());
	size_t index = (next - nodePositions.begin()) - 1;
	assert(index < nodePositions.size()-1);
	assert(nodePositions[index+1] > indexPos);
	assert(nodePositions[index] <= indexPos);
	return index;
}

std::vector<SeedHit> MEMSeeder::getMumSeeds(const std::string& sequence, size_t maxCount, size_t minLen) const
{
	assert(index.initialized());
	std::vector<SeedHit> result;
	auto matches = MEMfinder::getBestFwBwMUMs(index, sequence, minLen, maxCount);
	assert(matches.size() <= maxCount);
	result.reserve(matches.size());
	for (size_t i = 0; i < matches.size(); i++)
	{
		result.push_back(matchToSeed(matches[i]));
	}
	assert(result.size() == matches.size());
	assert(result.size() <= maxCount);
	return result;
}

std::vector<SeedHit> MEMSeeder::getMemSeeds(const std::string& sequence, size_t maxCount, size_t minLen) const
{
	assert(index.initialized());
	std::vector<SeedHit> result;
	std::vector<MEMfinder::Match> matches;
	size_t window = windowSize;
	if (window == 0) window = sequence.size();
	if (minLen < 10)
	{
		matches = MEMfinder::getBestFwBwMEMs(index, sequence, minLen, maxCount, uniqueBonusFactor, window);
	}
	else
	{
		matches = MEMfinder::getBestFwBwMEMs(index, sequence, minLen, maxCount, uniqueBonusFactor, prefixIndex, 10, window);
	}
	assert(matches.size() <= maxCount);
	result.reserve(matches.size());
	for (size_t i = 0; i < matches.size(); i++)
	{
		result.push_back(matchToSeed(matches[i]));
	}
	assert(result.size() == matches.size());
	assert(result.size() <= maxCount);
	return result;
}

SeedHit MEMSeeder::matchToSeed(MEMfinder::Match match) const
{
	assert(match.length > 0);
	assert(match.refPos + match.length <= nodePositions.back());
	auto index = getNodeIndex(match.refPos);
	assert(getNodeIndex(match.refPos + match.length) == index);
	int nodeID = nodeIDs[index];
	size_t nodeOffset = match.refPos - nodePositions[index];
	if (!match.fw)
	{
		nodeOffset = nodeLength(index) - (nodeOffset + match.length);
	}
	size_t seqPos = match.queryPos;
	size_t matchLen = match.length;
	return SeedHit { nodeID, nodeOffset, seqPos, matchLen, matchLen, !match.fw };
}

size_t MEMSeeder::nodeLength(size_t indexPos) const
{
	//-1 for separator
	return nodePositions[indexPos+1] - nodePositions[indexPos] - 1;
}
