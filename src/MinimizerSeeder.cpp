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

MinimizerSeeder::MinimizerSeeder(const AlignmentGraph& graph, size_t minimizerLength, size_t windowSize, size_t numThreads) :
graph(graph),
locator(),
starts(),
positions(),
minimizerLength(minimizerLength),
windowSize(windowSize),
maxCount(0)
{
	assert(minimizerLength * 2 <= sizeof(size_t) * 8);
	assert(minimizerLength <= windowSize);
	initMinimizers(numThreads);
	initMaxCount();
}

void MinimizerSeeder::initMinimizers(size_t numThreads)
{
	size_t positionSize = log2(graph.nodeIDs.size()) + 1;
	positions.width(positionSize+6);
	auto nodeIter = graph.nodeLookup.begin();
	std::mutex nodeMutex;
	std::vector<std::thread> threads;
	std::vector<sdsl::int_vector<0>> resultPerThread;
	resultPerThread.resize(numThreads);
	for (size_t i = 0; i < numThreads; i++)
	{
		resultPerThread[i].width(2 * minimizerLength + positionSize + 1);
	}

	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([this, &resultPerThread, i, &nodeMutex, &nodeIter, positionSize](){
			size_t vecPos = 0;
			resultPerThread[i].resize(10);
			while (true)
			{
				auto iter = graph.nodeLookup.end();
				{
					std::lock_guard<std::mutex> guard { nodeMutex };
					iter = nodeIter;
					if (nodeIter != graph.nodeLookup.end()) ++nodeIter;
				}
				if (iter == graph.nodeLookup.end()) break;
				int nodeId = iter->first;
				int firstNode = iter->second[0];
				std::string sequence;
				sequence.resize(graph.originalNodeSize.at(nodeId));
				for (size_t i = 0; i < sequence.size(); i++)
				{
					size_t nodeidHere = graph.GetUnitigNode(nodeId, i);
					sequence[i] = graph.NodeSequences(nodeidHere, i - graph.nodeOffset[nodeidHere]);
				}
				iterateMinimizers(sequence, minimizerLength, windowSize, [this, &resultPerThread, &vecPos, positionSize, i, nodeId](size_t pos, size_t kmer)
				{
					size_t splitNode = graph.GetUnitigNode(nodeId, pos);
					size_t remainingOffset = pos - graph.nodeOffset[splitNode];
					assert(remainingOffset < 64);
					unsigned __int128 storeThis = kmer;
					storeThis <<= positionSize;
					storeThis += splitNode;
					storeThis <<= 6;
					storeThis += remainingOffset;
					if (vecPos == resultPerThread[i].size()) resultPerThread[i].resize(resultPerThread[i].size() * 2);
					resultPerThread[i][vecPos] = storeThis;
					vecPos += 1;
				});
			}
			std::sort(resultPerThread[i].begin(), resultPerThread[i].end(), [positionSize](unsigned __int128 left, unsigned __int128 right) { return left < right; });
		});
	}

	size_t totalSize = 0;

	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
		totalSize += resultPerThread[i].size();
	}
	threads.clear();

	positions.resize(totalSize);
	size_t posi = 0;
	std::vector<size_t> threadIndex;
	threadIndex.resize(numThreads, 0);
	size_t current = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < numThreads; i++)
	{
		if (resultPerThread[i].size() == 0) continue;
		current = std::min(current, (resultPerThread[i][0] >> (positionSize + 6)));
	}
	std::vector<size_t> locatorKeys;
	std::vector<size_t> locatorValues;
	size_t startPos = 0;
	while (true)
	{
		size_t next = std::numeric_limits<size_t>::max();
		size_t numFinished = 0;
		locatorKeys.emplace_back(current);
		locatorValues.emplace_back(startPos);
		startPos += 1;
		for (size_t i = 0; i < numThreads; i++)
		{
			while (threadIndex[i] < resultPerThread[i].size() && (resultPerThread[i][threadIndex[i]] >> (positionSize + 6)) == current)
			{
				unsigned __int128 mergedPos = resultPerThread[i][threadIndex[i]];
				unsigned __int128 insert = mergedPos & ~(-1 << (positionSize + 6));
				assert(insert < (unsigned __int128)(1 << (positionSize + 6)));
				positions[posi] = insert;
				posi += 1;
				threadIndex[i] += 1;
			}
			if (threadIndex[i] < resultPerThread[i].size())
			{
				assert((resultPerThread[i][threadIndex[i]] >> (positionSize + 6)) > current);
				next = std::min(next, (resultPerThread[i][threadIndex[i]] >> (positionSize + 6)));
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
	locator.build(locatorKeys, locatorValues, numThreads);
	starts.resize(positions.size());
	sdsl::util::set_to_value(starts, 0);
	for (auto val : locatorValues)
	{
		starts[val] = 1;
	}
	sdsl::util::init_support(startSelector, &starts);
	assert(posi == positions.size());
	assert(positions.size() == totalSize);
	std::cerr << "locator: " << locator.size() << std::endl;
	std::cerr << "starts: " << starts.size() << " " << starts.capacity() << std::endl;
	std::cerr << "positions: " << positions.size() << " " << positions.capacity() << std::endl;
}

std::vector<SeedHit> MinimizerSeeder::getSeeds(const std::string& sequence, size_t maxCount, size_t chunkSize) const
{
	size_t numChunks = (sequence.size() + chunkSize - 1) / chunkSize;
	size_t bpPerChunk = (sequence.size() + numChunks - 1) / numChunks;
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> matchIndices;
	matchIndices.resize(numChunks);
	iterateMinimizers(sequence, minimizerLength, windowSize, [this, &matchIndices, bpPerChunk](size_t pos, size_t kmer)
	{
		if (!locator.contains(kmer)) return;
		auto found = locator[kmer];
		size_t chunk = pos / bpPerChunk;
		assert(chunk < matchIndices.size());
		matchIndices[chunk].emplace_back(pos, kmer, getStart(found+1) - getStart(found));
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
			auto found = locator[std::get<1>(match)];
			size_t end = getStart(found+1);
			for (size_t i = getStart(found); i < end; i++)
			{
				if (seedsHere >= maxCount) break;
				volatile size_t mergepos = positions[i];
				result.push_back(matchToSeedHit((size_t)(mergepos / 64), mergepos % 64, std::get<0>(match), std::get<2>(match)));
				seedsHere += 1;
			}
			if (seedsHere >= maxCount) break;
		}
	}
	return result;
}

SeedHit MinimizerSeeder::matchToSeedHit(int nodeId, size_t nodeOffset, size_t seqPos, int count) const
{
	assert(nodeId < graph.nodeIDs.size());
	assert(nodeId < graph.nodeOffset.size());
	assert(nodeId < graph.reverse.size());
	SeedHit result { graph.nodeIDs[nodeId]/2, nodeOffset + graph.nodeOffset[nodeId], seqPos, maxCount - count, graph.reverse[nodeId] };
	return result;
}

void MinimizerSeeder::initMaxCount()
{
	maxCount = 0;
	for (size_t i = 0; i < locator.size()-1; i++)
	{
		maxCount = std::max(maxCount, getStart(i+1) - getStart(i));
	}
	maxCount += 1;
}

size_t MinimizerSeeder::getStart(size_t index) const
{
	return startSelector(index+1);
}