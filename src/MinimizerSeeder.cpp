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

void MinimizerSeeder::initMinimizers(size_t numThreads)
{
	size_t positionSize = log2(graph.nodeIDs.size()) + 1;
	assert(positionSize + 6 < 64);
	assert(minimizerLength * 2 < 64);
	auto nodeIter = graph.nodeLookup.begin();
	std::mutex nodeMutex;
	std::vector<std::thread> threads;
	std::vector<sdsl::int_vector<0>> kmerPerBucket;
	std::vector<sdsl::int_vector<0>> positionPerBucket;
	std::vector<moodycamel::ConcurrentQueue<std::pair<uint64_t, uint64_t>>> positionDistributor;
	kmerPerBucket.resize(numThreads);
	positionPerBucket.resize(numThreads);
	positionDistributor.resize(numThreads);
	buckets.resize(numThreads);
	std::atomic<size_t> threadsDone;
	threadsDone = 0;
	for (size_t i = 0; i < numThreads; i++)
	{
		kmerPerBucket[i].width(minimizerLength * 2);
		positionPerBucket[i].width(positionSize + 6);
		buckets[i].positions.width(positionSize + 6);
	}

	for (size_t thread = 0; thread < numThreads; thread++)
	{
		threads.emplace_back([this, &positionDistributor, &threadsDone, &kmerPerBucket, &positionPerBucket, thread, numThreads, &nodeMutex, &nodeIter, positionSize](){
			size_t vecPos = 0;
			kmerPerBucket[thread].resize(10);
			positionPerBucket[thread].resize(10);
			std::pair<uint64_t, uint64_t> readThis;
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
				iterateMinimizers(sequence, minimizerLength, windowSize, [this, &positionDistributor, &kmerPerBucket, &positionPerBucket, &vecPos, positionSize, thread, nodeId](size_t pos, size_t kmer)
				{
					size_t splitNode = graph.GetUnitigNode(nodeId, pos);
					assert(splitNode < (size_t)1 << positionSize);
					size_t remainingOffset = pos - graph.nodeOffset[splitNode];
					assert(remainingOffset < 64);
					std::pair<uint64_t, uint64_t> readThis;
					while (positionDistributor[thread].try_dequeue(readThis))
					{
						assert(readThis.first < ((uint64_t)1) << (kmerPerBucket[thread].width()));
						assert(readThis.second < ((uint64_t)1) << (positionPerBucket[thread].width()));
						assert(getBucket(readThis.first) == thread);
						if (vecPos == kmerPerBucket[thread].size())
						{
							kmerPerBucket[thread].resize(kmerPerBucket[thread].size() * 2);
							positionPerBucket[thread].resize(kmerPerBucket[thread].size() * 2);
						}
						kmerPerBucket[thread][vecPos] = readThis.first;
						positionPerBucket[thread][vecPos] = readThis.second;
						vecPos += 1;
					}
					std::pair<uint64_t, uint64_t> storeThis;
					storeThis.first = kmer;
					size_t bucket = getBucket(kmer);
					storeThis.second = splitNode;
					storeThis.second <<= 6;
					storeThis.second += remainingOffset;
					positionDistributor[bucket].enqueue(storeThis);
				});
			}
			threadsDone += 1;
			while (threadsDone < numThreads)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				while (positionDistributor[thread].try_dequeue(readThis))
				{
					assert(readThis.first < ((uint64_t)1) << (kmerPerBucket[thread].width()));
					assert(readThis.second < ((uint64_t)1) << (positionPerBucket[thread].width()));
					assert(getBucket(readThis.first) == thread);
					if (vecPos == kmerPerBucket[thread].size())
					{
						kmerPerBucket[thread].resize(kmerPerBucket[thread].size() * 2);
						positionPerBucket[thread].resize(kmerPerBucket[thread].size() * 2);
					}
					kmerPerBucket[thread][vecPos] = readThis.first;
					positionPerBucket[thread][vecPos] = readThis.second;
					vecPos += 1;
				}
			}
			while (positionDistributor[thread].try_dequeue(readThis))
			{
				assert(readThis.first < ((uint64_t)1) << (kmerPerBucket[thread].width()));
				assert(readThis.second < ((uint64_t)1) << (positionPerBucket[thread].width()));
				assert(getBucket(readThis.first) == thread);
				if (vecPos == kmerPerBucket[thread].size())
				{
					kmerPerBucket[thread].resize(kmerPerBucket[thread].size() * 2);
					positionPerBucket[thread].resize(kmerPerBucket[thread].size() * 2);
				}
				kmerPerBucket[thread][vecPos] = readThis.first;
				positionPerBucket[thread][vecPos] = readThis.second;
				vecPos += 1;
			}
			kmerPerBucket[thread].resize(vecPos);
			positionPerBucket[thread].resize(vecPos);
			{
				std::vector<uint64_t> locatorKeys;
				{
					sdsl::int_vector<0> sortedKmers;
					sortedKmers = kmerPerBucket[thread];
					std::sort(sortedKmers.begin(), sortedKmers.end(), [positionSize](uint64_t left, uint64_t right) { return left < right; });
					size_t current = std::numeric_limits<size_t>::max();
					for (uint64_t kmer : sortedKmers)
					{
						assert(getBucket(kmer) == thread);
						if (kmer == current)
						{
							continue;
						}
						current = kmer;
						locatorKeys.push_back(current);
					}
				}
				buckets[thread].locator = new boomphf::mphf<uint64_t,KmerBucket::hasher_t>(locatorKeys.size(), locatorKeys, 1, 2, true, false);
			}
			sdsl::int_vector<0> counts;
			counts.width(log2(kmerPerBucket[thread].size()) + 1);
			counts.resize(buckets[thread].locator->nbKeys());
			buckets[thread].kmerCheck.width(minimizerLength * 2);
			buckets[thread].kmerCheck.resize(buckets[thread].locator->nbKeys());
			sdsl::util::set_to_value(counts, 0);
			for (size_t i = 0; i < kmerPerBucket[thread].size(); i++)
			{
				uint64_t kmer = kmerPerBucket[thread][i];
				size_t index = buckets[thread].locator->lookup(kmer);
				counts[index] += 1;
				buckets[thread].kmerCheck[index] = kmer;
			}
			buckets[thread].starts.resize(kmerPerBucket[thread].size()+1);
			sdsl::util::set_to_value(buckets[thread].starts, 0);
			buckets[thread].starts[0] = 1;
			buckets[thread].starts[counts[0]] = 1;
			for (size_t i = 1; i < counts.size(); i++)
			{
				counts[i] += counts[i-1];
				buckets[thread].starts[counts[i]] = 1;
			}
			assert(counts[counts.size()-1] == kmerPerBucket[thread].size());
			assert(buckets[thread].starts[kmerPerBucket[thread].size()]);
			buckets[thread].positions.resize(kmerPerBucket[thread].size());
			for (size_t i = 0; i < kmerPerBucket[thread].size(); i++)
			{
				size_t kmer = kmerPerBucket[thread][i];
				size_t index = buckets[thread].locator->lookup(kmer);
				assert(counts[index] > 0);
				counts[index] -= 1;
				size_t pos = counts[index];
				uint64_t insert = positionPerBucket[thread][i];
				buckets[thread].positions[pos] = insert;
			};

			sdsl::util::init_support(buckets[thread].startSelector, &buckets[thread].starts);
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
	std::vector<std::vector<std::tuple<size_t, size_t, size_t, size_t>>> matchIndices;
	matchIndices.resize(numChunks);
	size_t lastChainNodeId = 0;
	size_t lastChainNodeOffset = 0;
	size_t lastChainReadPos = 0;
	size_t lastChainCount = 0;
	size_t lastChainKmer = 0;
	size_t lastChainChunk = 0;
	iterateMinimizers(sequence, minimizerLength, windowSize, [this, &matchIndices, bpPerChunk, &lastChainNodeId, &lastChainNodeOffset, &lastChainReadPos, &lastChainCount, &lastChainKmer, &lastChainChunk](size_t pos, size_t kmer)
	{
		volatile size_t bucket = getBucket(kmer);
		assert(bucket < buckets.size());
		volatile size_t index = buckets[bucket].locator->lookup(kmer);
		if (index == ULLONG_MAX) return;
		assert(index < buckets[bucket].kmerCheck.size());
		if (buckets[bucket].kmerCheck[(size_t)index] != kmer) return;
		size_t chunk = pos / bpPerChunk;
		assert(chunk < matchIndices.size());
		bool canChain = false;
		size_t count = getStart(bucket, index+1) - getStart(bucket, index);
		if (count > 1)
		{
			if (lastChainCount > 0) matchIndices[lastChainChunk].emplace_back(lastChainReadPos, lastChainKmer, 1, lastChainCount);
			matchIndices[chunk].emplace_back(pos, kmer, count, 1);
			lastChainCount = 0;
			return;
		}
		size_t splitpos = buckets[bucket].positions[getStart(bucket, index)];
		size_t splitNodeId = pos >> 6;
		size_t splitOffset = pos & 63;
		size_t nodeOffset = graph.nodeOffset[splitNodeId] + splitOffset;
		size_t nodeId = graph.nodeIDs[splitNodeId];
		do
		{
			// unique
			if (count != 1) break;
			// right orientation in read
			if (pos <= lastChainReadPos) break;
			// same node
			if (nodeId != lastChainNodeId) break;
			// right orientation in graph
			if (nodeOffset <= lastChainNodeOffset) break;
			// if the kmer overlaps with the last one, the offsets have to match
			if ((nodeOffset - lastChainNodeOffset < minimizerLength || (pos - lastChainReadPos < minimizerLength)) && nodeOffset - lastChainNodeOffset != pos - lastChainReadPos) break;
			// if it doesn't overlap the offsets have to be within 50bp arbitrarily
			int readDiff = pos - lastChainReadPos;
			int nodeDiff = nodeOffset - lastChainNodeOffset;
			if ((nodeOffset - lastChainNodeOffset >= minimizerLength && (pos - lastChainReadPos >= minimizerLength)) && (readDiff - nodeDiff > 50 || readDiff - nodeDiff < -50)) break;
			canChain = true;
		} while (false);
		if (!canChain)
		{
			if (lastChainCount > 0) matchIndices[lastChainChunk].emplace_back(lastChainReadPos, lastChainKmer, 1, lastChainCount);
			lastChainCount = 0;
		}
		lastChainNodeId = nodeId;
		lastChainNodeOffset = nodeOffset;
		lastChainReadPos = pos;
		lastChainCount += 1;
		lastChainChunk = chunk;
		lastChainKmer = kmer;
	});
	if (lastChainCount > 0) matchIndices[lastChainChunk].emplace_back(lastChainReadPos, lastChainKmer, 1, lastChainCount);
	//prefer longer chains first, less common minimizers second
	for (size_t i = 0; i < numChunks; i++)
	{
		std::sort(matchIndices[i].begin(), matchIndices[i].end(), [this](const std::tuple<size_t, size_t, size_t, size_t>& left, const std::tuple<size_t, size_t, size_t, size_t>& right)
		{
			return std::get<3>(left) < std::get<3>(right) || (std::get<3>(left) == std::get<3>(right) && std::get<2>(left) < std::get<2>(right));
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
				size_t nodeId = mergepos >> 6;
				size_t offset = mergepos & 63;
				result.push_back(matchToSeedHit(nodeId, offset, std::get<0>(match), std::get<2>(match)));
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