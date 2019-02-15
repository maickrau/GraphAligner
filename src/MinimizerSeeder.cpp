#include <cmath>
#include "CommonUtils.h"
#include "MinimizerSeeder.h"

size_t charToInt(char c)
{
	switch(c)
	{
		case 'a':
		case 'A':
			return 0;
		case 'c':
		case 'C':
			return 1;
		case 'g':
		case 'G':
			return 2;
		case 't':
		case 'T':
			return 3;
	}
	assert(false);
	return 0;
}

std::vector<bool> getValidChars()
{
	std::vector<bool> result;
	result.resize(256, false);
	result['a'] = true;
	result['A'] = true;
	result['c'] = true;
	result['C'] = true;
	result['g'] = true;
	result['G'] = true;
	result['t'] = true;
	result['T'] = true;
	return result;
}

std::vector<bool> validChar = getValidChars();

template <typename OrderF, typename CallbackF>
void iterateMinimizers(const std::string& str, OrderF kmerOrder, size_t minimizerLength, size_t windowSize, CallbackF callback)
{
	assert(minimizerLength * 2 <= sizeof(size_t) * 8);
	assert(minimizerLength <= windowSize);
	if (str.size() < minimizerLength) return;
	size_t mask = ~(0xFFFFFFFFFFFFFFFF << (minimizerLength * 2));
	assert(mask == pow(4, minimizerLength)-1);
	size_t offset = 0;
	std::vector<size_t> window;
	window.resize(windowSize - minimizerLength + 1);
start:
	while (offset < str.size() && !validChar[str[offset]]) offset++;
	size_t kmer = 0;
	size_t minOrder = kmerOrder(kmer);
	size_t minPos = offset;
	size_t minKmer = kmer;
	size_t minWindowPos = 0;
	for (size_t i = 0; i < minimizerLength && offset+i < str.size(); i++)
	{
		if (!validChar[str[offset+i]])
		{
			offset = offset+i;
			goto start;
		}
		kmer <<= 2;
		kmer |= charToInt(str[offset+i]);
	}
	window[minimizerLength % window.size()] = kmer;
	for (size_t i = minimizerLength; i < windowSize && offset+i < str.size(); i++)
	{
		if (!validChar[str[offset+i]])
		{
			offset = offset+i;
			goto start;
		}
		kmer <<= 2;
		kmer &= mask;
		kmer |= charToInt(str[offset+i]);
		window[i % window.size()] = kmer;
		if (kmerOrder(kmer) < minOrder)
		{
			minOrder = kmerOrder(kmer);
			minPos = offset + i - minimizerLength + 1;
			minKmer = kmer;
			minWindowPos = i % window.size();
		}
	}
	if (minOrder != std::numeric_limits<size_t>::max()) callback(minPos, minKmer);
	for (size_t i = windowSize; offset+i < str.size(); i++)
	{
		if (!validChar[str[offset+i]])
		{
			offset = offset+i;
			goto start;
		}
		kmer <<= 2;
		kmer &= mask;
		kmer |= charToInt(str[offset+i]);
		window[i % window.size()] = kmer;
		if (minWindowPos == i % window.size())
		{
			minOrder = kmerOrder(window[0]);
			minWindowPos = 0;
			minPos = offset + i - (i % window.size());
			minKmer = window[0];
			for (size_t j = 1; j < window.size(); j++)
			{
				if (kmerOrder(window[j]) < minOrder)
				{
					minOrder = kmerOrder(window[j]);
					minWindowPos = j;
					minPos = offset + i - ((i - j) % window.size());
					minKmer = window[j];
				}
			}
			if (minOrder != std::numeric_limits<size_t>::max())
			{
				for (size_t j = 0; j < window.size(); j++)
				{
					if (kmerOrder(window[j]) == minOrder)
					{
						callback(offset + i - ((i - j) % window.size()), minKmer);
					}
				}
			}
		}
		else if (kmerOrder(kmer) <= minOrder)
		{
			minOrder = kmerOrder(kmer);
			minPos = offset + i - minimizerLength + 1;
			minKmer = kmer;
			minWindowPos = i % window.size();
			if (minOrder != std::numeric_limits<size_t>::max()) callback(minPos, minKmer);
		}
	}
}

MinimizerSeeder::MinimizerSeeder(const GfaGraph& graph, size_t minimizerLength, size_t windowSize) :
minimizerIndex(),
kmerOrder(),
minimizerLength(minimizerLength),
windowSize(windowSize),
maxCount(0),
rd(),
gen(rd()),
dis(0, std::numeric_limits<size_t>::max()-1)
{
	kmerOrder.set_empty_key(std::numeric_limits<size_t>::max());
	minimizerIndex.set_empty_key(std::numeric_limits<size_t>::max());
	assert(minimizerLength * 2 <= sizeof(size_t) * 8);
	assert(minimizerLength <= windowSize);
	for (auto node : graph.nodes)
	{
		addMinimizers(node.second, node.first*2);
		addMinimizers(CommonUtils::ReverseComplement(node.second), node.first*2+1);
	}
	finalizeOrder();
	initMaxCount();
}

MinimizerSeeder::MinimizerSeeder(const vg::Graph& graph, size_t minimizerLength, size_t windowSize) :
minimizerIndex(),
kmerOrder(),
minimizerLength(minimizerLength),
windowSize(windowSize),
maxCount(0),
rd(),
gen(rd()),
dis(0, std::numeric_limits<size_t>::max()-1)
{
	kmerOrder.set_empty_key(std::numeric_limits<size_t>::max());
	minimizerIndex.set_empty_key(std::numeric_limits<size_t>::max());
	assert(minimizerLength * 2 <= sizeof(size_t) * 8);
	assert(minimizerLength <= windowSize);
	for (int i = 0; i < graph.node_size(); i++)
	{
		addMinimizers(graph.node(i).sequence(), graph.node(i).id()*2);
		addMinimizers(CommonUtils::ReverseComplement(graph.node(i).sequence()), graph.node(i).id()*2+1);
	}
	finalizeOrder();
	initMaxCount();
}

std::vector<SeedHit> MinimizerSeeder::getSeeds(const std::string& sequence, size_t maxCount) const
{
	std::vector<std::tuple<size_t, size_t, size_t>> matchIndices;
	iterateMinimizers(sequence, [this](size_t kmer){ return getOrder(kmer); }, minimizerLength, windowSize, [this, &matchIndices](size_t pos, size_t kmer)
	{
		auto found = minimizerIndex.find(kmer);
		assert(found != minimizerIndex.end());
		matchIndices.emplace_back(pos, kmer, found->second.size());
	});
	//prefer less common minimizers
	std::sort(matchIndices.begin(), matchIndices.end(), [this](const std::tuple<size_t, size_t, size_t>& left, const std::tuple<size_t, size_t, size_t>& right)
	{
		return std::get<2>(left) < std::get<2>(right);
	});
	std::vector<SeedHit> result;
	for (auto match : matchIndices)
	{
		auto found = minimizerIndex.find(std::get<1>(match));
		assert(found != minimizerIndex.end());
		for (auto index : found->second)
		{
			if (result.size() >= maxCount) break;
			result.push_back(matchToSeedHit(index.first, index.second, std::get<0>(match), std::get<2>(match)));
		}
		if (result.size() >= maxCount) break;
	}
	return result;
}

SeedHit MinimizerSeeder::matchToSeedHit(int nodeId, size_t nodeOffset, size_t seqPos, int count) const
{
	SeedHit result { nodeId/2, nodeOffset, seqPos, maxCount - count, nodeId % 2 };
	return result;
}

void MinimizerSeeder::addMinimizers(const std::string& str, int nodeId)
{
	iterateMinimizers(str, [this](size_t kmer){ return getOrAddOrder(kmer); }, minimizerLength, windowSize, [this, nodeId](size_t pos, size_t kmer)
	{
		minimizerIndex[kmer].emplace_back(nodeId, pos);
	});
}

size_t MinimizerSeeder::getOrAddOrder(size_t kmer)
{
	auto found = kmerOrder.find(kmer);
	if (found != kmerOrder.end()) return found->second;
	size_t result = dis(gen);
	kmerOrder[kmer] = result;
	return result;
}

size_t MinimizerSeeder::getOrder(size_t kmer) const
{
	auto found = kmerOrder.find(kmer);
	if (found != kmerOrder.end()) return found->second;
	return std::numeric_limits<size_t>::max();
}

void MinimizerSeeder::initMaxCount()
{
	maxCount = 0;
	for (auto pair : minimizerIndex)
	{
		maxCount = std::max(maxCount, pair.second.size());
	}
	maxCount += 1;
}

void MinimizerSeeder::finalizeOrder()
{
	decltype(kmerOrder) newOrder;
	newOrder.set_empty_key(std::numeric_limits<size_t>::max());
	for (auto pair : kmerOrder)
	{
		auto found = minimizerIndex.find(pair.first);
		if (found != minimizerIndex.end() && found->second.size() > 0)
		{
			newOrder[pair.first] = pair.second;
		}
	}
	std::swap(kmerOrder, newOrder);
}
