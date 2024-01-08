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
	std::vector<std::pair<__uint128_t, size_t>> currentString;
	std::vector<std::pair<std::vector<std::pair<__uint128_t, size_t>>, size_t>> validStrings;
	for (size_t i = 0; i < graph.BigraphNodeCount(); i++)
	{
		std::string seq = graph.BigraphNodeSeq(i);
		nodeLengths[i] = seq.size();
		iterateSplittedSeqs(seq, [this, &currentString, i, &validStrings, &kmerCount, &kmerPosition, &graph](std::string str)
		{
			iterateSmallSyncmers(k, 5, str, [this, &currentString, i, &validStrings, &kmerCount, &kmerPosition, &graph](size_t pos, __uint128_t kmer)
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
						currentString.emplace_back(kmer, pos);
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
						currentString.emplace_back(kmer, pos);
						return;
					}
					if (currentString.size() == 1)
					{
						currentString[0] = std::make_pair(kmer, pos);
						return;
					}
					assert(currentString.size() >= 2);
					currentString.emplace_back(kmer, pos);
					if (kmerPosition.at(currentString[0].first) != kmerPosition.at(currentString.back().first))
					{
						currentString.clear();
						currentString.emplace_back(kmer, pos);
						return;
					}
					validStrings.emplace_back(currentString, i);
					currentString.clear();
				}
			});
		});
	}
	std::map<std::pair<__uint128_t, __uint128_t>, std::pair<size_t, size_t>> potentialPairs;
	for (size_t i = 0; i < validStrings.size(); i++)
	{
		assert(validStrings[i].first.size() >= 3);
		std::pair<__uint128_t, __uint128_t> key { validStrings[i].first[0].first, validStrings[i].first.back().first };
		assert(kmerPosition.at(validStrings[i].first[0].first) == kmerPosition.at(validStrings[i].first.back().first));
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
	phmap::flat_hash_map<size_t, phmap::flat_hash_set<size_t>> pairs;
	for (auto pair : potentialPairs)
	{
		if (pair.second.second == std::numeric_limits<size_t>::max()) continue;
		size_t firstNode = validStrings[pair.second.first].second;
		assert(kmerPosition.at(validStrings[pair.second.first].first[1].first).first == firstNode);
		size_t secondNode = validStrings[pair.second.second].second;
		assert(kmerPosition.at(validStrings[pair.second.second].first[1].first).first == secondNode);
		assert(secondNode > firstNode);
		for (size_t i = 1; i < validStrings[pair.second.first].first.size()-1; i++)
		{
			assert(kmerPosition.at(validStrings[pair.second.first].first[i].first).first == firstNode);
			kmerImpliesNode[validStrings[pair.second.first].first[i].first] = std::make_pair(firstNode, validStrings[pair.second.first].first[i].second);
		}
		for (size_t i = 1; i < validStrings[pair.second.second].first.size()-1; i++)
		{
			assert(kmerPosition.at(validStrings[pair.second.second].first[i].first).first == secondNode);
			kmerImpliesNode[validStrings[pair.second.second].first[i].first] = std::make_pair(secondNode, validStrings[pair.second.second].first[i].second);
		}
		pairs[firstNode].emplace(secondNode);
		pairs[secondNode].emplace(firstNode);
	}
	for (auto pair : pairs)
	{
		conflictPairs[pair.first].insert(conflictPairs[pair.first].end(), pair.second.begin(), pair.second.end());
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

phmap::flat_hash_set<std::tuple<size_t, int, int>> DiploidHeuristicSplitterOneK::getForbiddenNodes(std::string sequence) const
{
	if (conflictPairs.size() == 0) return phmap::flat_hash_set<std::tuple<size_t, int, int>> {};
	phmap::flat_hash_map<size_t, std::vector<int>> forbidPositions;
	iterateSplittedSeqs(sequence, [this, &forbidPositions](std::string str)
	{
		iterateSmallSyncmers(k, 5, str, [this, &forbidPositions](size_t pos, __uint128_t kmer)
		{
			if (kmerImpliesNode.count(kmer) == 0) return;
			std::pair<size_t, size_t> val = kmerImpliesNode.at(kmer);
			forbidPositions[val.first].emplace_back((int)pos - (int)val.second);
		});
	});
	phmap::flat_hash_map<size_t, std::vector<std::pair<int, int>>> solidPositions;
	for (auto& pair : forbidPositions)
	{
		std::sort(pair.second.begin(), pair.second.end());
		size_t clusterStart = 0;
		for (size_t i = 1; i < pair.second.size(); i++)
		{
			if (pair.second[i] > pair.second[clusterStart] + 100)
			{
				if (i - clusterStart >= 3)
				{
					int pos = (pair.second[clusterStart] + pair.second[i-1])/2;
					solidPositions[pair.first].emplace_back(pos, pos + nodeLengths.at(pair.first));
				}
				clusterStart = i;
			}
		}
		if (pair.second.size() - clusterStart >= 3)
		{
			int pos = (pair.second[clusterStart] + pair.second.back())/2;
			solidPositions[pair.first].emplace_back(pos, pos + nodeLengths.at(pair.first));
		}
	}
	phmap::flat_hash_map<size_t, std::vector<std::pair<int, int>>> forbiddenSpans;
	for (const auto& pair : solidPositions)
	{
		if (conflictPairs.count(pair.first) == 0) continue;
		for (const auto pos : pair.second)
		{
			for (const size_t otherNode : conflictPairs.at(pair.first))
			{
				bool canBlock = true;
				if (solidPositions.count(otherNode) == 1)
				{
					for (const auto& pos2 : solidPositions.at(otherNode))
					{
						if (pos2.first < pos.second && pos2.second > pos.first)
						{
							canBlock = false;
							break;
						}
					}
				}
				if (!canBlock) continue;
				forbiddenSpans[otherNode].emplace_back(pos.first, pos.second);
			}
		}
	}
	phmap::flat_hash_set<std::tuple<size_t, int, int>> result;
	for (auto& pair : forbiddenSpans)
	{
		std::sort(pair.second.begin(), pair.second.end(), [](auto left, auto right) { return left.first < right.first; });
		int currentSpanStart = 0;
		int currentSpanEnd = 0;
		for (size_t i = 1; i < pair.second.size(); i++)
		{
			if (pair.second[i].first > currentSpanEnd)
			{
				if (currentSpanEnd > currentSpanStart) result.emplace(pair.first, currentSpanStart, currentSpanEnd);
				currentSpanStart = pair.second[i].first;
				currentSpanEnd = pair.second[i].second;
			}
			else
			{
				currentSpanEnd = std::max(currentSpanEnd, pair.second[i].second);
			}
		}
		if (currentSpanEnd > currentSpanStart) result.emplace(pair.first, currentSpanStart, currentSpanEnd);
	}
	return result;
}

void DiploidHeuristicSplitterOneK::write(std::ostream& file) const
{
	serialize(file, (uint64_t)k);
	serialize(file, (uint64_t)kmerImpliesNode.size());
	for (auto pair : kmerImpliesNode)
	{
		serialize(file, (__uint128_t)pair.first);
		serialize(file, (uint64_t)pair.second.first);
		serialize(file, (uint64_t)pair.second.second);
	}
	serialize(file, (uint64_t)conflictPairs.size());
	for (const auto& pair : conflictPairs)
	{
		serialize(file, (uint64_t)pair.first);
		serialize(file, (uint64_t)pair.second.size());
		for (auto node : pair.second)
		{
			serialize(file, (uint64_t)node);
		}
	}
	serialize(file, (uint64_t)nodeLengths.size());
	for (auto pair : nodeLengths)
	{
		serialize(file, (uint64_t)pair.first);
		serialize(file, (uint64_t)pair.second);
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
		uint64_t value1, value2;
		deserialize(file, key);
		deserialize(file, value1);
		deserialize(file, value2);
		kmerImpliesNode[key] = std::make_pair(value1, value2);
	}
	deserialize(file, numItems);
	for (size_t i = 0; i < numItems; i++)
	{
		uint64_t val;
		deserialize(file, val);
		uint64_t count;
		deserialize(file, count);
		for (size_t j = 0; j < count; j++)
		{
			uint64_t othernode;
			deserialize(file, othernode);
			conflictPairs[val].emplace_back(othernode);
		}
	}
	deserialize(file, numItems);
	for (size_t i = 0; i < numItems; i++)
	{
		uint64_t key, value;
		deserialize(file, key);
		deserialize(file, value);
		nodeLengths[key] = value;
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

std::vector<std::tuple<size_t, int, int>> DiploidHeuristicSplitter::getForbiddenNodes(std::string sequence) const
{
	phmap::flat_hash_set<std::tuple<size_t, int, int>> result;
	for (size_t i = 0; i < splitters.size(); i++)
	{
		phmap::flat_hash_set<std::tuple<size_t, int, int>> forbidden = splitters[i].getForbiddenNodes(sequence);
		result.insert(forbidden.begin(), forbidden.end());
	}
	std::vector<std::tuple<size_t, int, int>> resultVec { result.begin(), result.end() };
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
