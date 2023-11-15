#include <map>
#include "CommonUtils.h"
#include "DiploidHeuristic.h"
#include "Serialize.h"

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
	__uint128_t kmer = 0;
	__uint128_t smer = 0;
	const __uint128_t smerMask = ((__uint128_t)1 << (__uint128_t)(s*2)) - (__uint128_t)1;
	const __uint128_t kmerMask = ((__uint128_t)1 << (__uint128_t)(k*2)) - (__uint128_t)1;
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

void countKmers(phmap::flat_hash_map<__uint128_t, uint8_t>& kmerCount, const std::string& seq, const size_t maxCount, const size_t k)
{
	iterateSmallSyncmers(k, 5, seq, [&kmerCount, maxCount](size_t pos, __uint128_t kmer)
	{
		if (kmerCount[kmer] < maxCount) kmerCount[kmer] += 1;
	});
}

void getKmerPositions(const phmap::flat_hash_map<__uint128_t, uint8_t>& kmerCount, phmap::flat_hash_map<__uint128_t, std::pair<size_t, size_t>>& kmerPosition, const std::string& seq, const size_t node, const size_t k)
{
	iterateSmallSyncmers(k, 5, seq, [&kmerCount, &kmerPosition, node](size_t pos, __uint128_t kmer)
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

size_t DiploidHeuristicSplitterOneK::getk() const
{
	return k;
}

void DiploidHeuristicSplitterOneK::getHomologyPairs(const phmap::flat_hash_map<__uint128_t, uint8_t>& kmerCount, const phmap::flat_hash_map<__uint128_t, std::pair<size_t, size_t>>& kmerPosition, const AlignmentGraph& graph)
{
	std::vector<__uint128_t> currentString;
	std::vector<std::vector<__uint128_t>> validStrings;
	for (size_t i = 0; i < graph.BigraphNodeCount(); i++)
	{
		std::string seq = graph.BigraphNodeSeq(i);
		iterateSplittedSeqs(seq, [this, &currentString, &validStrings, &kmerCount, &kmerPosition, &graph](std::string str)
		{
			iterateSmallSyncmers(k, 5, str, [this, &currentString, &validStrings, &kmerCount, &kmerPosition, &graph](size_t pos, __uint128_t kmer)
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
	std::map<std::pair<__uint128_t, __uint128_t>, std::pair<size_t, size_t>> potentialPairs;
	for (size_t i = 0; i < validStrings.size(); i++)
	{
		assert(validStrings[i].size() >= 3);
		std::pair<__uint128_t, __uint128_t> key { validStrings[i][0], validStrings[i].back() };
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

void DiploidHeuristicSplitterOneK::initializePairs(const AlignmentGraph& graph, size_t k)
{
	this->k = k;
	phmap::flat_hash_map<__uint128_t, uint8_t> kmerCount;
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
	phmap::flat_hash_map<__uint128_t, std::pair<size_t, size_t>> kmerPosition;
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

phmap::flat_hash_set<size_t> DiploidHeuristicSplitterOneK::getForbiddenNodes(std::string sequence) const
{
	if (homologyPairs.size() == 0) return phmap::flat_hash_set<size_t> {};
	phmap::flat_hash_map<size_t, size_t> forbidCount;
	iterateSplittedSeqs(sequence, [this, &forbidCount](std::string str)
	{
		iterateSmallSyncmers(k, 5, str, [this, &forbidCount](size_t pos, __uint128_t kmer)
		{
			if (kmerForbidsNode.count(kmer) == 0) return;
			forbidCount[kmerForbidsNode.at(kmer)] += 1;
		});
	});
	phmap::flat_hash_set<size_t> result;
	for (auto pair : homologyPairs)
	{
		if (forbidCount[pair.first] <= 3 && forbidCount[pair.second] <= 3) continue;
		if (forbidCount[pair.first] >= 1 && forbidCount[pair.second] >= 1 && forbidCount[pair.first] < forbidCount[pair.second]*10 && forbidCount[pair.second] < forbidCount[pair.first]*10) continue;
		if (forbidCount[pair.first] > forbidCount[pair.second]) result.emplace(pair.first);
		if (forbidCount[pair.second] > forbidCount[pair.first]) result.emplace(pair.second);
	}
	return result;
}

void DiploidHeuristicSplitterOneK::write(std::ostream& file) const
{
	serialize(file, (uint64_t)k);
	serialize(file, (uint64_t)kmerForbidsNode.size());
	for (auto pair : kmerForbidsNode)
	{
		serialize(file, (__uint128_t)pair.first);
		serialize(file, (uint64_t)pair.second);
	}
	serialize(file, homologyPairs.size());
	for (size_t i = 0; i < (uint64_t)homologyPairs.size(); i++)
	{
		serialize(file, (uint64_t)homologyPairs[i].first);
		serialize(file, (uint64_t)homologyPairs[i].second);
	}
}

void DiploidHeuristicSplitterOneK::read(std::istream& file)
{
	assert(file.good());
	uint64_t numItems;
	deserialize(file, numItems);
	k = numItems;
	deserialize(file, numItems);
	for (size_t i = 0; i < numItems; i++)
	{
		__uint128_t key;
		uint64_t value;
		deserialize(file, key);
		deserialize(file, value);
		kmerForbidsNode[key] = value;
	}
	deserialize(file, numItems);
	homologyPairs.resize(numItems);
	for (size_t i = 0; i < numItems; i++)
	{
		uint64_t val;
		deserialize(file, val);
		homologyPairs[i].first = val;
		deserialize(file, val);
		homologyPairs[i].second = val;
	}
	assert(file.good());
}

void DiploidHeuristicSplitter::initializePairs(const AlignmentGraph& graph, const std::vector<size_t>& kValues)
{
	splitters.resize(kValues.size());
	for (size_t i = 0; i < kValues.size(); i++)
	{
		splitters[i].initializePairs(graph, kValues[i]);
	}
}

std::vector<size_t> DiploidHeuristicSplitter::getForbiddenNodes(std::string sequence) const
{
	phmap::flat_hash_set<size_t> result;
	for (size_t i = 0; i < splitters.size(); i++)
	{
		phmap::flat_hash_set<size_t> forbidden = splitters[i].getForbiddenNodes(sequence);
		result.insert(forbidden.begin(), forbidden.end());
	}
	std::vector<size_t> resultVec { result.begin(), result.end() };
	return resultVec;
}

void DiploidHeuristicSplitter::write(std::string filename) const
{
	std::ofstream file { filename, std::ios::binary };
	uint64_t splitterCount = splitters.size();
	serialize(file, (uint64_t)splitterCount);
	for (size_t i = 0; i < splitters.size(); i++)
	{
		splitters[i].write(file);
	}
}

void DiploidHeuristicSplitter::read(std::string filename)
{
	std::ifstream file { filename, std::ios::binary };
	uint64_t splitterCount;
	deserialize(file, splitterCount);
	splitters.resize(splitterCount);
	for (size_t i = 0; i < splitters.size(); i++)
	{
		splitters[i].read(file);
	}
}

std::vector<size_t> DiploidHeuristicSplitter::getKValues() const
{
	std::vector<size_t> result;
	for (size_t i = 0; i < splitters.size(); i++)
	{
		result.push_back(splitters[i].getk());
	}
	return result;
}
