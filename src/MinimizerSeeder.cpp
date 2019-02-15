#include <random>
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

template <typename F>
void iterateMinimizers(const std::string& str, const std::vector<size_t>& kmerOrder, size_t minimizerLength, size_t windowSize, F f)
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
	size_t minOrder = kmerOrder[kmer];
	size_t minPos = offset;
	size_t minKmer = kmer;
	size_t minWindowPos = 0;
	for (size_t i = 0; i < minimizerLength && offset+i < str.size(); i++)
	{
		if (!validChar[str[offset+i]]) goto start;
		kmer <<= 2;
		kmer |= charToInt(str[offset+i]);
	}
	window[minimizerLength % window.size()] = kmer;
	for (size_t i = minimizerLength; i < windowSize && offset+i < str.size(); i++)
	{
		if (!validChar[str[offset+i]]) goto start;
		kmer <<= 2;
		kmer &= mask;
		kmer |= charToInt(str[offset+i]);
		window[i % window.size()] = kmer;
		if (kmerOrder[kmer] < minOrder)
		{
			minOrder = kmerOrder[kmer];
			minPos = offset + i - minimizerLength + 1;
			minKmer = kmer;
			minWindowPos = i % window.size();
		}
	}
	f(minPos, minKmer);
	for (size_t i = windowSize; offset+i < str.size(); i++)
	{
		if (!validChar[str[offset+i]]) goto start;
		kmer <<= 2;
		kmer &= mask;
		kmer |= charToInt(str[offset+i]);
		window[i % window.size()] = kmer;
		if (minWindowPos == i % window.size())
		{
			minOrder = kmerOrder[window[0]];
			minWindowPos = 0;
			minPos = offset + i - (i % window.size());
			minKmer = window[0];
			for (size_t j = 1; j < window.size(); j++)
			{
				if (kmerOrder[window[j]] < minOrder)
				{
					minOrder = kmerOrder[window[j]];
					minWindowPos = j;
					minPos = offset + i - ((i - j) % window.size());
					minKmer = window[j];
				}
			}
			for (size_t j = 0; j < window.size(); j++)
			{
				if (kmerOrder[window[j]] == minOrder)
				{
					assert(window[j] == minKmer);
					f(offset + i - ((i - j) % window.size()), minKmer);
				}
			}
		}
		else if (kmerOrder[kmer] <= minOrder)
		{
			minOrder = kmerOrder[kmer];
			minPos = offset + i - minimizerLength + 1;
			minKmer = kmer;
			minWindowPos = i % window.size();
			f(minPos, minKmer);
		}
	}
}

MinimizerSeeder::MinimizerSeeder(const GfaGraph& graph, size_t minimizerLength, size_t windowSize) :
minimizerIndex(),
kmerOrder(),
minimizerLength(minimizerLength),
windowSize(windowSize),
maxCount(0)
{
	initOrder();
	initEmptyIndex();
	for (auto node : graph.nodes)
	{
		addMinimizers(node.second, node.first*2);
		addMinimizers(CommonUtils::ReverseComplement(node.second), node.first*2+1);
	}
	initMaxCount();
}

MinimizerSeeder::MinimizerSeeder(const vg::Graph& graph, size_t minimizerLength, size_t windowSize) :
minimizerIndex(),
kmerOrder(),
minimizerLength(minimizerLength),
windowSize(windowSize),
maxCount()
{
	initOrder();
	initEmptyIndex();
	for (int i = 0; i < graph.node_size(); i++)
	{
		addMinimizers(graph.node(i).sequence(), graph.node(i).id()*2);
		addMinimizers(CommonUtils::ReverseComplement(graph.node(i).sequence()), graph.node(i).id()*2+1);
	}
	initMaxCount();
}

std::vector<SeedHit> MinimizerSeeder::getSeeds(const std::string& sequence, size_t maxCount) const
{
	std::vector<std::pair<size_t, size_t>> matchIndices;
	iterateMinimizers(sequence, kmerOrder, minimizerLength, windowSize, [this, &matchIndices](size_t pos, size_t kmer)
	{
		if (minimizerIndex[kmer].size() > 0) matchIndices.emplace_back(pos, kmer);
	});
	//prefer less common minimizers
	std::sort(matchIndices.begin(), matchIndices.end(), [this](const std::pair<size_t, size_t>& left, const std::pair<size_t, size_t>& right)
	{
		return minimizerIndex[left.second].size() < minimizerIndex[right.second].size();
	});
	std::vector<SeedHit> result;
	for (auto match : matchIndices)
	{
		assert(minimizerIndex[match.second].size() >= 1);
		int count = minimizerIndex[match.second].size();
		for (auto index : minimizerIndex[match.second])
		{
			if (result.size() >= maxCount) break;
			result.push_back(matchToSeedHit(index.first, index.second, match.first, count));
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
	iterateMinimizers(str, kmerOrder, minimizerLength, windowSize, [this, nodeId](size_t pos, size_t kmer)
	{
		minimizerIndex[kmer].emplace_back(nodeId, pos);
	});
}

void MinimizerSeeder::initOrder()
{
	kmerOrder.resize(pow(4, minimizerLength));
	for (size_t i = 0; i < kmerOrder.size(); i++)
	{
		kmerOrder[i] = i;
	}
	std::random_device rd;
	std::mt19937 gen(rd());
	for (size_t i = kmerOrder.size()-1; i > 0; i--)
	{
		std::uniform_int_distribution<> dis(0, i);
		std::swap(kmerOrder[i], kmerOrder[dis(gen)]);
	}
}

void MinimizerSeeder::initMaxCount()
{
	maxCount = 0;
	for (size_t i = 0; i < minimizerIndex.size(); i++)
	{
		maxCount = std::max(maxCount, minimizerIndex[i].size());
	}
	maxCount += 1;
}

void MinimizerSeeder::initEmptyIndex()
{
	minimizerIndex.resize(pow(4, minimizerLength));
}
