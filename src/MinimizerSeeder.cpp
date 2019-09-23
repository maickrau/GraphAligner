#include <thread>
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

// https://naml.us/post/inverse-of-a-hash-function/
uint64_t hash(uint64_t key) {
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}

std::vector<bool> validChar = getValidChars();

template <typename CallbackF>
void iterateMinimizers(const std::string& str, size_t minimizerLength, size_t windowSize, CallbackF callback)
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
	if (offset + minimizerLength > str.size()) return;
	size_t kmer = 0;
	for (size_t i = 0; i < minimizerLength; i++)
	{
		assert(offset+i < str.size());
		if (!validChar[str[offset+i]])
		{
			offset = offset+i;
			goto start;
		}
		kmer <<= 2;
		kmer |= charToInt(str[offset+i]);
	}
	size_t minOrder = hash(kmer);
	size_t minPos = offset + minimizerLength - 1;
	size_t minKmer = kmer;
	size_t minWindowPos = minimizerLength % window.size();
	window[(minimizerLength-1) % window.size()] = kmer;
	for (size_t i = minimizerLength; i < windowSize && offset+i < str.size(); i++)
	{
		if (!validChar[str[offset+i]])
		{
			if (minOrder != std::numeric_limits<size_t>::max())
			{
				for (size_t j = 0; j < window.size(); j++)
				{
					size_t seqPos = ((j + window.size() - minimizerLength) % window.size()) + offset + minimizerLength;
					if (seqPos >= str.size()) continue;
					if (hash(window[j]) == minOrder)
					{
						callback(seqPos, window[j]);
					}
				}
			}
			offset = offset+i;
			goto start;
		}
		kmer <<= 2;
		kmer &= mask;
		kmer |= charToInt(str[offset+i]);
		window[i % window.size()] = kmer;
		if (hash(kmer) < minOrder)
		{
			minOrder = hash(kmer);
			minPos = offset + i;
			minKmer = kmer;
			minWindowPos = i % window.size();
		}
	}
	if (minOrder != std::numeric_limits<size_t>::max())
	{
		for (size_t j = 0; j < window.size(); j++)
		{
			size_t seqPos = ((j + window.size() - minimizerLength) % window.size()) + offset + minimizerLength;
			if (seqPos >= str.size()) continue;
			if (hash(window[j]) == minOrder)
			{
				callback(seqPos, window[j]);
			}
		}
	}
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
			minOrder = hash(window[0]);
			minWindowPos = 0;
			minPos = offset + i - (i % window.size());
			minKmer = window[0];
			for (size_t j = 1; j < window.size(); j++)
			{
				if (hash(window[j]) < minOrder)
				{
					minOrder = hash(window[j]);
					minWindowPos = j;
					minPos = offset + i - ((i - j) % window.size());
					minKmer = window[j];
				}
			}
			if (minOrder != std::numeric_limits<size_t>::max())
			{
				for (size_t j = 0; j < window.size(); j++)
				{
					if (hash(window[j]) == minOrder)
					{
						callback(offset + i - ((i - j) % window.size()), window[j]);
					}
				}
			}
		}
		else if (hash(kmer) <= minOrder)
		{
			minOrder = hash(kmer);
			minPos = offset + i;
			minKmer = kmer;
			minWindowPos = i % window.size();
			if (minOrder != std::numeric_limits<size_t>::max()) callback(minPos, minKmer);
		}
	}
}

MinimizerSeeder::MinimizerSeeder(const GfaGraph& graph, size_t minimizerLength, size_t windowSize, size_t numThreads) :
locator(),
starts(),
positions(),
minimizerLength(minimizerLength),
windowSize(windowSize),
maxCount(0),
rd(),
gen((size_t)rd() ^ (size_t)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()),
dis(0, std::numeric_limits<size_t>::max()-1)
{
	assert(minimizerLength * 2 <= sizeof(size_t) * 8);
	assert(minimizerLength <= windowSize);
	initMinimizers(graph, numThreads);
	initMaxCount();
}

MinimizerSeeder::MinimizerSeeder(const vg::Graph& graph, size_t minimizerLength, size_t windowSize, size_t numThreads) :
locator(),
starts(),
positions(),
minimizerLength(minimizerLength),
windowSize(windowSize),
maxCount(0),
rd(),
gen((size_t)rd() ^ (size_t)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()),
dis(0, std::numeric_limits<size_t>::max()-1)
{
	assert(minimizerLength * 2 <= sizeof(size_t) * 8);
	assert(minimizerLength <= windowSize);
	initMinimizers(graph, numThreads);
	initMaxCount();
}

void MinimizerSeeder::initMinimizers(const GfaGraph& graph, size_t numThreads)
{
	auto nodeIter = graph.nodes.begin();
	std::mutex nodeMutex;
	std::vector<std::thread> threads;
	std::vector<std::vector<std::tuple<size_t, int, size_t>>> resultPerThread;
	resultPerThread.resize(numThreads);

	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([this, &graph, &resultPerThread, i, &nodeMutex, &nodeIter](){
			while (true)
			{
				auto iter = graph.nodes.end();
				{
					std::lock_guard<std::mutex> guard { nodeMutex };
					iter = nodeIter;
					if (nodeIter != graph.nodes.end()) ++nodeIter;
				}
				if (iter == graph.nodes.end()) break;
				int nodeId = iter->first * 2;
				iterateMinimizers(iter->second, minimizerLength, windowSize, [&resultPerThread, i, nodeId](size_t pos, size_t kmer)
				{
					resultPerThread[i].emplace_back(kmer, nodeId, pos);
				});
				nodeId = iter->first * 2 + 1;
				iterateMinimizers(CommonUtils::ReverseComplement(iter->second), minimizerLength, windowSize, [&resultPerThread, i, nodeId](size_t pos, size_t kmer)
				{
					resultPerThread[i].emplace_back(kmer, nodeId, pos);
				});
			}
			std::sort(resultPerThread[i].begin(), resultPerThread[i].end(), [](const std::tuple<size_t, int, size_t>& left, const std::tuple<size_t, int, size_t>& right) { return std::get<0>(left) < std::get<0>(right); });
		});
	}

	size_t totalSize = 0;

	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
		totalSize += resultPerThread[i].size();
	}
	threads.clear();

	positions.reserve(totalSize);
	std::vector<size_t> threadIndex;
	threadIndex.resize(numThreads, 0);
	size_t current = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < numThreads; i++)
	{
		if (resultPerThread[i].size() == 0) continue;
		current = std::min(current, std::get<0>(resultPerThread[i][0]));
	}
	while (true)
	{
		size_t next = std::numeric_limits<size_t>::max();
		size_t numFinished = 0;
		starts.emplace_back(positions.size());
		locator[current] = starts.size()-1;
		for (size_t i = 0; i < numThreads; i++)
		{
			while (threadIndex[i] < resultPerThread[i].size() && std::get<0>(resultPerThread[i][threadIndex[i]]) == current)
			{
				positions.emplace_back(std::get<1>(resultPerThread[i][threadIndex[i]]), std::get<2>(resultPerThread[i][threadIndex[i]]));
				threadIndex[i] += 1;
			}
			if (threadIndex[i] < resultPerThread[i].size())
			{
				assert(std::get<0>(resultPerThread[i][threadIndex[i]]) > current);
				next = std::min(next, std::get<0>(resultPerThread[i][threadIndex[i]]));
			}
			else
			{
				numFinished += 1;
			}
		}
		if (numFinished == numThreads) break;
		assert(next > current);
		current = next;
	}
	starts.emplace_back(positions.size());
	assert(starts.size() == locator.size()+1);
	assert(positions.size() == totalSize);
	std::cerr << "locator: " << locator.size() << " " << locator.capacity() << std::endl;
	std::cerr << "starts: " << starts.size() << " " << starts.capacity() << std::endl;
	std::cerr << "positions: " << positions.size() << " " << positions.capacity() << std::endl;
}

void MinimizerSeeder::initMinimizers(const vg::Graph& graph, size_t numThreads)
{
	int nodeI = 0;
	std::mutex nodeMutex;
	std::vector<std::thread> threads;
	std::vector<std::vector<std::tuple<size_t, int, size_t>>> resultPerThread;
	resultPerThread.resize(numThreads);

	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([this, &graph, &resultPerThread, i, &nodeMutex, &nodeI](){
			while (true)
			{
				auto index = graph.node_size();
				{
					std::lock_guard<std::mutex> guard { nodeMutex };
					index = nodeI;
					if (nodeI < graph.node_size()) ++nodeI;
				}
				if (index == graph.node_size()) break;
				int nodeId = graph.node(index).id() * 2;
				iterateMinimizers(graph.node(index).sequence(), minimizerLength, windowSize, [&resultPerThread, i, nodeId](size_t pos, size_t kmer)
				{
					resultPerThread[i].emplace_back(kmer, nodeId, pos);
				});
				nodeId = graph.node(index).id() * 2 + 1;
				iterateMinimizers(CommonUtils::ReverseComplement(graph.node(index).sequence()), minimizerLength, windowSize, [&resultPerThread, i, nodeId](size_t pos, size_t kmer)
				{
					resultPerThread[i].emplace_back(kmer, nodeId, pos);
				});
			}
			std::sort(resultPerThread[i].begin(), resultPerThread[i].end(), [](const std::tuple<size_t, int, size_t>& left, const std::tuple<size_t, int, size_t>& right) { return std::get<0>(left) < std::get<0>(right); });
		});
	}

	size_t totalSize = 0;

	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
		totalSize += resultPerThread[i].size();
	}
	threads.clear();

	positions.reserve(totalSize);
	std::vector<size_t> threadIndex;
	threadIndex.resize(numThreads, 0);
	size_t current = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < numThreads; i++)
	{
		if (resultPerThread[i].size() == 0) continue;
		current = std::min(current, std::get<0>(resultPerThread[i][0]));
	}
	while (true)
	{
		size_t next = std::numeric_limits<size_t>::max();
		size_t numFinished = 0;
		starts.emplace_back(positions.size());
		locator[current] = starts.size()-1;
		for (size_t i = 0; i < numThreads; i++)
		{
			while (threadIndex[i] < resultPerThread[i].size() && std::get<0>(resultPerThread[i][threadIndex[i]]) == current)
			{
				positions.emplace_back(std::get<1>(resultPerThread[i][threadIndex[i]]), std::get<2>(resultPerThread[i][threadIndex[i]]));
				threadIndex[i] += 1;
			}
			if (threadIndex[i] < resultPerThread[i].size())
			{
				assert(std::get<0>(resultPerThread[i][threadIndex[i]]) > current);
				next = std::min(next, std::get<0>(resultPerThread[i][threadIndex[i]]));
			}
			else
			{
				numFinished += 1;
			}
		}
		if (numFinished == numThreads) break;
		assert(next > current);
		current = next;
	}
	starts.emplace_back(positions.size());
	assert(starts.size() == locator.size()+1);
	assert(positions.size() == totalSize);
}

std::vector<SeedHit> MinimizerSeeder::getSeeds(const std::string& sequence, size_t maxCount, size_t chunkSize) const
{
	size_t numChunks = (sequence.size() + chunkSize - 1) / chunkSize;
	size_t bpPerChunk = (sequence.size() + numChunks - 1) / numChunks;
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> matchIndices;
	matchIndices.resize(numChunks);
	iterateMinimizers(sequence, minimizerLength, windowSize, [this, &matchIndices, bpPerChunk](size_t pos, size_t kmer)
	{
		auto found = locator.find(kmer);
		if (found == locator.end()) return;
		size_t chunk = pos / bpPerChunk;
		assert(chunk < matchIndices.size());
		matchIndices[chunk].emplace_back(pos, kmer, starts[found->second+1] - starts[found->second]);
	});
	//prefer less common minimizers
	for (size_t i = 0; i < numChunks; i++)
	{
		std::sort(matchIndices[i].begin(), matchIndices[i].end(), [this](const std::tuple<size_t, size_t, size_t>& left, const std::tuple<size_t, size_t, size_t>& right)
		{
			return std::get<2>(left) < std::get<2>(right);
		});
	}
	std::vector<SeedHit> result;
	for (size_t i = 0; i < numChunks; i++)
	{
		size_t seedsHere = 0;
		for (auto match : matchIndices[i])
		{
			auto found = locator.find(std::get<1>(match));
			assert(found != locator.end());
			for (size_t i = starts[found->second]; i < starts[found->second+1]; i++)
			{
				if (seedsHere >= maxCount) break;
				result.push_back(matchToSeedHit(positions[i].first, positions[i].second, std::get<0>(match), std::get<2>(match)));
				seedsHere += 1;
			}
			if (seedsHere >= maxCount) break;
		}
	}
	return result;
}

SeedHit MinimizerSeeder::matchToSeedHit(int nodeId, size_t nodeOffset, size_t seqPos, int count) const
{
	SeedHit result { nodeId/2, nodeOffset, seqPos, maxCount - count, nodeId % 2 };
	return result;
}

void MinimizerSeeder::initMaxCount()
{
	maxCount = 0;
	for (size_t i = 0; i < starts.size()-1; i++)
	{
		maxCount = std::max(maxCount, starts[i+1] - starts[i]);
	}
	maxCount += 1;
}
