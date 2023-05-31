#include <map>
#include "DiploidHeuristic.h"

uint64_t hash(uint64_t key);

void switchBasePairsToNumbers(std::string& str)
{
	for (size_t i = 0; i < str.size(); i++)
	{
		switch(str[i])
		{
			case 'A':
				str[i] = 0;
				break;
			case 'C':
				str[i] = 1;
				break;
			case 'G':
				str[i] = 2;
				break;
			case 'T':
				str[i] = 3;
				break;
		}
	}
}

template <typename F>
void iterateSplittedSeqs(const std::string& str, F callback)
{
	size_t lastValid = 0;
	for (size_t i = 0; i < str.size(); i++)
	{
		if (str[i] == 'A' || str[i] == 'C' || str[i] == 'G' || str[i] == 'T') continue;
		if (i == lastValid)
		{
			lastValid = i+1;
			continue;
		}
		std::string sub = str.substr(lastValid, i-lastValid);
		switchBasePairsToNumbers(sub);
		callback(sub);
		lastValid = i+1;
	}
	if (lastValid < str.size())
	{
		std::string sub = str.substr(lastValid);
		switchBasePairsToNumbers(sub);
		callback(sub);
	}
}

template <typename F>
void iterateSmallSyncmers(const size_t k, const size_t s, const std::string& sequence, F callback)
{
	if (sequence.size() < k) return;
	uint64_t kmer = 0;
	uint64_t smer = 0;
	const uint64_t smerMask = ~(0xFFFFFFFFFFFFFFFF << (s * 2));
	const uint64_t kmerMask = ~(0xFFFFFFFFFFFFFFFF << (k * 2));
	std::vector<std::pair<size_t, size_t>> window;
	for (size_t i = 0; i < s; i++)
	{
		assert(sequence[i] <= 3);
		smer <<= 2;
		smer += sequence[i];
	}
	window.emplace_back(0, hash(smer));
	for (size_t i = s; i < k; i++)
	{
		smer <<= 2;
		smer += sequence[i];
		smer &= smerMask;
		size_t h = hash(smer);
		while (window.size() >= 1 && window.back().second > h) window.pop_back();
		window.emplace_back(i-s+1, h);
	}
	for (size_t i = 0; i < k; i++)
	{
		assert(sequence[i] <= 3);
		kmer <<= 2;
		kmer += sequence[i];
	}
	assert(window.size() >= 1);
	if (window[0].first == 0 || window[0].first == k-s) callback(0, kmer);
	for (size_t i = k; i < sequence.size(); i++)
	{
		smer <<= 2;
		smer += sequence[i];
		smer &= smerMask;
		kmer <<= 2;
		kmer += sequence[i];
		kmer &= kmerMask;
		while (window.size() >= 1 && window[0].first <= i-k) window.erase(window.begin());
		size_t h = hash(smer);
		while (window.size() >= 1 && window.back().second > h) window.pop_back();
		window.emplace_back(i-s+1, h);
		assert(window.size() >= 1);
		if (window[0].first == i-k+1 || window[0].first == i-s+1) callback(i-k+1, kmer);
	}
}

void countKmers(phmap::flat_hash_map<uint64_t, uint8_t>& kmerCount, const std::string& seq, const size_t maxCount, const size_t k)
{
	iterateSmallSyncmers(k, 5, seq, [&kmerCount, maxCount](size_t pos, uint64_t kmer)
	{
		if (kmerCount[kmer] < maxCount) kmerCount[kmer] += 1;
	});
}

void getKmerPositions(const phmap::flat_hash_map<uint64_t, uint8_t>& kmerCount, phmap::flat_hash_map<uint64_t, std::pair<size_t, size_t>>& kmerPosition, const std::string& seq, const size_t node, const size_t k)
{
	iterateSmallSyncmers(k, 5, seq, [&kmerCount, &kmerPosition, node](size_t pos, uint64_t kmer)
	{
		if (kmerCount.at(kmer) > 2) return;
		auto found = kmerPosition.find(kmer);
		if (found == kmerPosition.end())
		{
			kmerPosition[kmer] = std::make_pair(node, std::numeric_limits<size_t>::max());
		}
		else
		{
			assert(kmerPosition[kmer].second == std::numeric_limits<size_t>::max());
			kmerPosition[kmer].second = node;
		}
	});
}

void DiploidHeuristicSplitter::getHomologyPairs(const phmap::flat_hash_map<uint64_t, uint8_t>& kmerCount, const phmap::flat_hash_map<uint64_t, std::pair<size_t, size_t>>& kmerPosition, const AlignmentGraph& graph)
{
	std::vector<uint64_t> currentString;
	std::vector<std::vector<uint64_t>> validStrings;
	for (size_t i = 0; i < graph.BigraphNodeCount(); i++)
	{
		std::string seq = graph.BigraphNodeSeq(i);
		iterateSplittedSeqs(seq, [this, &currentString, &validStrings, &kmerCount, &kmerPosition, &graph](std::string str)
		{
			iterateSmallSyncmers(k, 5, str, [this, &currentString, &validStrings, &kmerCount, &kmerPosition, &graph](size_t pos, uint64_t kmer)
			{
				if (kmerCount.at(kmer) > 2)
				{
					currentString.clear();
					return;
				}
				if (kmerCount.at(kmer) == 1)
				{
					if (currentString.size() >= 1)
					{
						currentString.push_back(kmer);
						return;
					}
					currentString.clear();
					return;
				}
				if (kmerCount.at(kmer) == 2)
				{
					if (kmerPosition.at(kmer).first == kmerPosition.at(kmer).second)
					{
						currentString.clear();
						return;
					}
					if (graph.BigraphNodeName(kmerPosition.at(kmer).first) == graph.BigraphNodeName(kmerPosition.at(kmer).second))
					{
						currentString.clear();
						return;
					}
					if (currentString.size() == 0)
					{
						currentString.push_back(kmer);
						return;
					}
					if (currentString.size() == 1)
					{
						currentString[0] = kmer;
						return;
					}
					assert(currentString.size() >= 2);
					currentString.push_back(kmer);
					if (kmerPosition.at(currentString[0]) != kmerPosition.at(currentString.back()))
					{
						currentString.clear();
						currentString.push_back(kmer);
						return;
					}
					validStrings.push_back(currentString);
					currentString.clear();
				}
			});
		});
	}
	std::map<std::pair<size_t, size_t>, std::pair<size_t, size_t>> potentialPairs;
	for (size_t i = 0; i < validStrings.size(); i++)
	{
		assert(validStrings[i].size() >= 3);
		std::pair<size_t, size_t> key { validStrings[i][0], validStrings[i].back() };
		assert(kmerPosition.at(validStrings[i][0]) == kmerPosition.at(validStrings[i].back()));
		if (potentialPairs.count(key) == 0)
		{
			potentialPairs[key] = std::make_pair(i, std::numeric_limits<size_t>::max());
		}
		else if (potentialPairs.at(key).second == std::numeric_limits<size_t>::max())
		{
			potentialPairs[key].second = i;
		}
		else
		{
			potentialPairs[key] = std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
		}
	}
	phmap::flat_hash_set<std::pair<size_t, size_t>> pairs;
	for (auto pair : potentialPairs)
	{
		if (pair.second.second == std::numeric_limits<size_t>::max()) continue;
		size_t firstNode = kmerPosition.at(validStrings[pair.second.first][1]).first;
		size_t secondNode = kmerPosition.at(validStrings[pair.second.second][1]).first;
		assert(secondNode > firstNode);
		for (size_t i = 1; i < validStrings[pair.second.first].size()-1; i++)
		{
			assert(kmerPosition.at(validStrings[pair.second.first][i]).first == firstNode);
			kmerForbidsNode[validStrings[pair.second.first][i]] = secondNode;
		}
		for (size_t i = 1; i < validStrings[pair.second.second].size()-1; i++)
		{
			assert(kmerPosition.at(validStrings[pair.second.second][i]).first == secondNode);
			kmerForbidsNode[validStrings[pair.second.second][i]] = firstNode;
		}
		pairs.emplace(firstNode, secondNode);
	}
	for (auto pair : pairs)
	{
		homologyPairs.emplace_back(pair);
	}
}

void DiploidHeuristicSplitter::initializePairs(const AlignmentGraph& graph, size_t k)
{
	this->k = k;
	phmap::flat_hash_map<uint64_t, uint8_t> kmerCount;
	for (size_t i = 0; i < graph.BigraphNodeCount(); i++)
	{
		std::string seq = graph.BigraphNodeSeq(i);
		iterateSplittedSeqs(seq, [&kmerCount, k](std::string str)
		{
			countKmers(kmerCount, str, 3, k);
		});
	}
	size_t countOnes = 0;
	size_t countTwos = 0;
	for (auto pair : kmerCount)
	{
		if (pair.second == 1) countOnes += 1;
		if (pair.second == 2) countTwos += 1;
	}
	phmap::flat_hash_map<uint64_t, std::pair<size_t, size_t>> kmerPosition;
	for (size_t i = 0; i < graph.BigraphNodeCount(); i++)
	{
		std::string seq = graph.BigraphNodeSeq(i);
		iterateSplittedSeqs(seq, [&kmerCount, &kmerPosition, i, k](std::string str)
		{
			getKmerPositions(kmerCount, kmerPosition, str, i, k);
		});
	}
	getHomologyPairs(kmerCount, kmerPosition, graph);
}

std::vector<size_t> DiploidHeuristicSplitter::getForbiddenNodes(std::string sequence) const
{
	if (homologyPairs.size() == 0) return std::vector<size_t> {};
	phmap::flat_hash_map<size_t, size_t> forbidCount;
	iterateSplittedSeqs(sequence, [this, &forbidCount](std::string str)
	{
		iterateSmallSyncmers(k, 5, str, [this, &forbidCount](size_t pos, uint64_t kmer)
		{
			if (kmerForbidsNode.count(kmer) == 0) return;
			forbidCount[kmerForbidsNode.at(kmer)] += 1;
		});
	});
	std::vector<size_t> result;
	for (auto pair : homologyPairs)
	{
		if (forbidCount[pair.first] <= 3 && forbidCount[pair.second] <= 3) continue;
		if (forbidCount[pair.first] >= 1 && forbidCount[pair.second] >= 1 && forbidCount[pair.first] < forbidCount[pair.second]*10 && forbidCount[pair.second] < forbidCount[pair.first]*10) continue;
		if (forbidCount[pair.first] > forbidCount[pair.second]) result.push_back(pair.first);
		if (forbidCount[pair.second] > forbidCount[pair.first]) result.push_back(pair.second);
	}
	return result;
}
