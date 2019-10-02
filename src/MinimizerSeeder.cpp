#include <thread>
#include <cmath>
#include <concurrentqueue.h>
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
buckets(),
minimizerLength(minimizerLength),
windowSize(windowSize),
maxCount(0)
{
	assert(minimizerLength * 2 <= sizeof(size_t) * 8);
	assert(minimizerLength <= windowSize);
	initMinimizers(numThreads);
	initMaxCount();
}

template <typename F>
void iteratePerThreads(const std::vector<sdsl::int_vector<0>>& resultPerThread, size_t positionSize, F callback)
{
	size_t numThreads = resultPerThread.size();
	size_t current = std::numeric_limits<size_t>::max();
	std::vector<size_t> threadIndex;
	threadIndex.resize(numThreads, 0);
	for (size_t i = 0; i < numThreads; i++)
	{
		if (resultPerThread[i].size() == 0) continue;
		current = std::min(current, (resultPerThread[i][0] >> (positionSize + 6)));
	}
	while (true)
	{
		size_t next = std::numeric_limits<size_t>::max();
		size_t numFinished = 0;
		for (size_t i = 0; i < numThreads; i++)
		{
			while (threadIndex[i] < resultPerThread[i].size() && (resultPerThread[i][threadIndex[i]] >> (positionSize + 6)) == current)
			{
				callback(resultPerThread[i][threadIndex[i]]);
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
}

void MinimizerSeeder::initMinimizers(size_t numThreads)
{
	size_t positionSize = log2(graph.nodeIDs.size()) + 1;
	assert(minimizerLength * 2 + positionSize + 6 < 128);
	auto nodeIter = graph.nodeLookup.begin();
	std::mutex nodeMutex;
	std::vector<std::thread> threads;
	std::vector<sdsl::int_vector<0>> resultPerThread;
	std::vector<moodycamel::ConcurrentQueue<unsigned __int128>> positionDistributor;
	resultPerThread.resize(numThreads);
	positionDistributor.resize(numThreads);
	buckets.resize(numThreads);
	std::atomic<size_t> threadsDone;
	threadsDone = 0;
	for (size_t i = 0; i < numThreads; i++)
	{
		resultPerThread[i].width(2 * minimizerLength + positionSize + 6);
		buckets[i].positions.width(positionSize + 6);
	}

	for (size_t thread = 0; thread < numThreads; thread++)
	{
		threads.emplace_back([this, &positionDistributor, &threadsDone, &resultPerThread, thread, numThreads, &nodeMutex, &nodeIter, positionSize](){
			size_t vecPos = 0;
			resultPerThread[thread].resize(10);
			unsigned __int128 readThis;
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
				for (size_t pos = 0; pos < sequence.size(); pos++)
				{
					size_t nodeidHere = graph.GetUnitigNode(nodeId, pos);
					sequence[pos] = graph.NodeSequences(nodeidHere, pos - graph.nodeOffset[nodeidHere]);
				}
				iterateMinimizers(sequence, minimizerLength, windowSize, [this, &positionDistributor, &resultPerThread, &vecPos, positionSize, thread, nodeId](size_t pos, size_t kmer)
				{
					size_t splitNode = graph.GetUnitigNode(nodeId, pos);
					size_t remainingOffset = pos - graph.nodeOffset[splitNode];
					assert(remainingOffset < 64);
					unsigned __int128 readThis;
					while (positionDistributor[thread].try_dequeue(readThis))
					{
						size_t kmer = readThis >> (positionSize + 6);
						assert(getBucket(kmer) == thread);
						if (vecPos == resultPerThread[thread].size()) resultPerThread[thread].resize(resultPerThread[thread].size() * 2);
						resultPerThread[thread][vecPos] = readThis;
						assert((unsigned __int128)resultPerThread[thread][vecPos] == readThis);
						assert(getBucket((unsigned __int128)resultPerThread[thread][vecPos] >> (positionSize + 6)) == thread);
						vecPos += 1;
					}
					unsigned __int128 storeThis = kmer;
					size_t bucket = getBucket(kmer);
					storeThis <<= positionSize;
					storeThis += splitNode;
					storeThis <<= 6;
					storeThis += remainingOffset;
					positionDistributor[bucket].enqueue(storeThis);
				});
			}
			threadsDone += 1;
			while (threadsDone < numThreads)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				while (positionDistributor[thread].try_dequeue(readThis))
				{
					size_t kmer = readThis >> (positionSize + 6);
					assert(getBucket(kmer) == thread);
					if (vecPos == resultPerThread[thread].size()) resultPerThread[thread].resize(resultPerThread[thread].size() * 2);
					resultPerThread[thread][vecPos] = readThis;
					assert((unsigned __int128)resultPerThread[thread][vecPos] == readThis);
					assert(getBucket((unsigned __int128)resultPerThread[thread][vecPos] >> (positionSize + 6)) == thread);
					vecPos += 1;
				}
			}
			while (positionDistributor[thread].try_dequeue(readThis))
			{
				size_t kmer = readThis >> (positionSize + 6);
				assert(getBucket(kmer) == thread);
				if (vecPos == resultPerThread[thread].size()) resultPerThread[thread].resize(resultPerThread[thread].size() * 2);
				resultPerThread[thread][vecPos] = readThis;
				assert((unsigned __int128)resultPerThread[thread][vecPos] == readThis);
				assert(getBucket((unsigned __int128)resultPerThread[thread][vecPos] >> (positionSize + 6)) == thread);
				vecPos += 1;
			}
			resultPerThread[thread].resize(vecPos);
			std::sort(resultPerThread[thread].begin(), resultPerThread[thread].end(), [positionSize](unsigned __int128 left, unsigned __int128 right) { return left < right; });
			size_t current = std::numeric_limits<size_t>::max();
			std::vector<size_t> locatorKeys;
			for (unsigned __int128 value : resultPerThread[thread])
			{
				size_t kmer = value >> (positionSize + 6);
				assert(getBucket(kmer) == thread);
				if (kmer == current)
				{
					continue;
				}
				current = kmer;
				locatorKeys.push_back(current);
			}
			buckets[thread].locator = new boomphf::mphf<size_t,KmerBucket::hasher_t>(locatorKeys.size(), locatorKeys, 1, 2, true, false);
			{
				decltype(locatorKeys) tmp;
				std::swap(tmp, locatorKeys);
			}
			sdsl::int_vector<0> counts;
			counts.width(log2(resultPerThread[thread].size()) + 1);
			counts.resize(buckets[thread].locator->nbKeys());
			buckets[thread].kmerCheck.width(minimizerLength * 2);
			buckets[thread].kmerCheck.resize(buckets[thread].locator->nbKeys());
			sdsl::util::set_to_value(counts, 0);
			for (unsigned __int128 value : resultPerThread[thread])
			{
				size_t kmer = value >> (positionSize + 6);
				size_t index = buckets[thread].locator->lookup(kmer);
				counts[index] += 1;
				buckets[thread].kmerCheck[index] = kmer;
			}
			buckets[thread].starts.resize(resultPerThread[thread].size());
			sdsl::util::set_to_value(buckets[thread].starts, 0);
			buckets[thread].starts[0] = 1;
			buckets[thread].starts[counts[0]] = 1;
			for (size_t i = 1; i < counts.size()-1; i++)
			{
				counts[i] += counts[i-1];
				buckets[thread].starts[counts[i]] = 1;
			}
			if (counts.size() >= 2) counts[counts.size()-1] += counts[counts.size()-2];
			assert(counts[counts.size()-1] == resultPerThread[thread].size());
			buckets[thread].positions.resize(resultPerThread[thread].size());
			for (auto value : resultPerThread[thread])
			{
				size_t kmer = value >> (positionSize + 6);
				size_t index = buckets[thread].locator->lookup(kmer);
				assert(counts[index] > 0);
				counts[index] -= 1;
				size_t pos = counts[index];
				unsigned __int128 insert = value & ~(-1 << (positionSize + 6));
				assert(insert < (unsigned __int128)(1 << (positionSize + 6)));
				buckets[thread].positions[pos] = insert;
			};

			sdsl::util::init_support(buckets[thread].startSelector, &buckets[thread].starts);
			assert(buckets[thread].positions.size() == resultPerThread[thread].size());
		});
	}

	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
	threads.clear();
}

std::vector<SeedHit> MinimizerSeeder::getSeeds(const std::string& sequence, size_t maxCount, size_t chunkSize) const
{
	size_t numChunks = (sequence.size() + chunkSize - 1) / chunkSize;
	size_t bpPerChunk = (sequence.size() + numChunks - 1) / numChunks;
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> matchIndices;
	matchIndices.resize(numChunks);
	iterateMinimizers(sequence, minimizerLength, windowSize, [this, &matchIndices, bpPerChunk](size_t pos, size_t kmer)
	{
		volatile size_t bucket = getBucket(kmer);
		assert(bucket < buckets.size());
		volatile size_t index = buckets[bucket].locator->lookup(kmer);
		if (index == ULLONG_MAX) return;
		assert(index < buckets[bucket].kmerCheck.size());
		if (buckets[bucket].kmerCheck[(size_t)index] != kmer) return;
		size_t chunk = pos / bpPerChunk;
		assert(chunk < matchIndices.size());
		matchIndices[chunk].emplace_back(pos, kmer, getStart(bucket, index+1) - getStart(bucket, index));
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
			size_t bucket = getBucket(std::get<1>(match));
			auto found = buckets[bucket].locator->lookup(std::get<1>(match));
			size_t end = getStart(bucket, found+1);
			for (size_t i = getStart(bucket, found); i < end; i++)
			{
				if (seedsHere >= maxCount) break;
				size_t mergepos = buckets[bucket].positions[i];
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
	for (size_t bucket = 0; bucket < buckets.size(); bucket++)
	{
		for (size_t i = 0; i < buckets[bucket].locator->nbKeys()-1; i++)
		{
			maxCount = std::max(maxCount, getStart(bucket, i+1) - getStart(bucket, i));
		}
	}
	maxCount += 1;
}

size_t MinimizerSeeder::getStart(size_t bucket, size_t index) const
{
	return buckets[bucket].startSelector(index+1);
}

size_t MinimizerSeeder::getBucket(size_t hash) const
{
	return hash % buckets.size();
}

MinimizerSeeder::KmerBucket::KmerBucket() :
	locator(nullptr)
{}

MinimizerSeeder::KmerBucket::~KmerBucket()
{
	if (locator != nullptr) delete locator;
}