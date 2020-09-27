#include <queue>
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
void iterateKmers(const std::string& str, size_t kmerLength, size_t windowSize, CallbackF callback)
{
	const size_t realWindow = windowSize - kmerLength + 1;
	assert(kmerLength * 2 <= sizeof(size_t) * 8);
	if (str.size() < kmerLength) return;
	const size_t mask = ~(0xFFFFFFFFFFFFFFFF << (kmerLength * 2));
	assert(mask == pow(4, kmerLength)-1);
	size_t offset = 0;
start:
	while (offset < str.size() && !validChar[str[offset]]) offset++;
	if (offset + kmerLength > str.size()) return;
	size_t kmer = 0;
	for (size_t i = 0; i < kmerLength; i++)
	{
		if (!validChar[str[offset+i]])
		{
			offset += i;
			goto start;
		}
		kmer <<= 2;
		kmer |= charToInt(str[offset+i]);
	}
	callback(offset + kmerLength-1, kmer);
	size_t lastKmer = kmer;
	size_t lastPos = offset + kmerLength-1;
	for (size_t i = kmerLength; offset+i < str.size(); i++)
	{
		if (!validChar[str[offset+i]])
		{
			offset += i;
			goto start;
		}
		kmer <<= 2;
		kmer &= mask;
		kmer |= charToInt(str[offset + i]);
		if (lastKmer != kmer || lastPos <= offset + i - realWindow)
		{
			callback(offset + i, kmer);
			lastKmer = kmer;
			lastPos = offset + i;
		}
	}
}

template <typename CallbackF>
void iterateMinimizersReal(const std::string& str, size_t minimizerLength, size_t windowSize, CallbackF callback)
{
	assert(minimizerLength * 2 <= sizeof(size_t) * 8);
	assert(minimizerLength <= windowSize);
	if (str.size() < minimizerLength) return;
	const size_t realWindow = windowSize - minimizerLength + 1;
	const size_t mask = ~(0xFFFFFFFFFFFFFFFF << (minimizerLength * 2));
	assert(mask == pow(4, minimizerLength)-1);
	size_t offset = 0;
	std::deque<std::tuple<size_t, size_t, size_t>> window;
start:
	while (offset < str.size() && !validChar[str[offset]]) offset++;
	if (offset + windowSize > str.size()) return;
	size_t kmer = 0;
	for (size_t i = 0; i < minimizerLength; i++)
	{
		if (!validChar[str[offset+i]])
		{
			offset += i;
			goto start;
		}
		kmer <<= 2;
		kmer |= charToInt(str[offset+i]);
	}
	window.clear();
	window.emplace_back(offset+minimizerLength-1, kmer, hash(kmer));
	for (size_t i = minimizerLength; i < minimizerLength + realWindow; i++)
	{
		if (!validChar[str[offset + i]])
		{
			offset += i;
			goto start;
		}
		kmer <<= 2;
		kmer &= mask;
		kmer |= charToInt(str[offset + i]);
		auto hashed = hash(kmer);
		while (!window.empty() && std::get<2>(window.back()) > hashed) window.pop_back();
		window.emplace_back(offset+i, kmer, hashed);
	}
	auto iter = window.begin();
	while (iter != window.end() && std::get<2>(*iter) == std::get<2>(window.front()))
	{
		callback(std::get<0>(*iter), std::get<1>(*iter));
		++iter;
	}
	for (size_t i = minimizerLength + realWindow; offset+i < str.size(); i++)
	{
		if (!validChar[str[offset+i]])
		{
			offset += i;
			goto start;
		}
		kmer <<= 2;
		kmer &= mask;
		kmer |= charToInt(str[offset + i]);
		auto hashed = hash(kmer);
		size_t oldMinimum = std::get<2>(window.front());
		bool frontPopped = false;
		while (!window.empty() && std::get<0>(window.front()) <= offset + i - realWindow)
		{
			frontPopped = true;
			window.pop_front();
		}
		if (frontPopped)
		{
			while (window.size() >= 2 && std::get<2>(window.front()) == std::get<2>(*(window.begin()+1))) window.pop_front();
		}
		while (!window.empty() && std::get<2>(window.back()) > hashed) window.pop_back();
		window.emplace_back(offset+i, kmer, hashed);
		if (std::get<2>(window.front()) != oldMinimum)
		{
			auto iter = window.begin();
			while (iter != window.end() && std::get<2>(*iter) == std::get<2>(window.front()))
			{
				callback(std::get<0>(*iter), std::get<1>(*iter));
				++iter;
			}
		}
		else if (std::get<2>(window.back()) == std::get<2>(window.front()))
		{
			callback(std::get<0>(window.back()), std::get<1>(window.back()));
		}
	}
}

#ifndef EXTRACORRECTNESSASSERTIONS

template <typename CallbackF>
void iterateMinimizers(const std::string& str, size_t minimizerLength, size_t windowSize, CallbackF callback)
{
	iterateMinimizersReal(str, minimizerLength, windowSize, callback);
}

#else

template <typename CallbackF>
void iterateMinimizersSimple(const std::string& str, size_t minimizerLength, size_t windowSize, CallbackF callback)
{
	if (str.size() <= windowSize)
	{
		size_t windowMinimum = std::numeric_limits<size_t>::max();
		for (size_t i = minimizerLength-1; i < str.size(); i++)
		{
			size_t kmer = 0;
			for (size_t j = i + 1 - minimizerLength; j <= i; j++)
			{
				kmer <<= 2;
				kmer |= charToInt(str[j]);
			}
			windowMinimum = std::min(windowMinimum, hash(kmer));
		}
		for (size_t i = minimizerLength-1; i < str.size(); i++)
		{
			size_t kmer = 0;
			for (size_t j = i + 1 - minimizerLength; j <= i; j++)
			{
				kmer <<= 2;
				kmer |= charToInt(str[j]);
			}
			if (hash(kmer) == windowMinimum) callback(i, kmer);
		}
		return;
	}
	std::vector<bool> isMinimizer;
	isMinimizer.resize(str.size(), false);
	std::vector<size_t> kmer;
	std::vector<size_t> kmerHash;
	kmer.resize(str.size(), 0);
	kmerHash.resize(str.size(), std::numeric_limits<size_t>::max());
	for (size_t i = minimizerLength-1; i < str.size(); i++)
	{
		kmer[i] = 0;
		for (size_t j = i + 1 - minimizerLength; j <= i; j++)
		{
			kmer[i] <<= 2;
			kmer[i] |= charToInt(str[j]);
		}
		kmerHash[i] = hash(kmer[i]);
	}
	for (size_t i = windowSize-1; i < str.size(); i++)
	{
		size_t windowMinimum = kmerHash[i];
		for (size_t j = i + minimizerLength - windowSize; j <= i; j++)
		{
			windowMinimum = std::min(windowMinimum, kmerHash[j]);
		}
		for (size_t j = i + minimizerLength - windowSize; j <= i; j++)
		{
			if (kmerHash[j] == windowMinimum)
			{
				isMinimizer[j] = true;
			}
		}
	}
	for (size_t i = 0; i < str.size(); i++)
	{
		if (!isMinimizer[i]) continue;
		callback(i, kmer[i]);
	}
}

template <typename CallbackF>
void iterateMinimizers(const std::string& str, size_t minimizerLength, size_t windowSize, CallbackF callback)
{
	std::vector<std::pair<size_t, size_t>> simpleMinimizers;
	std::vector<std::pair<size_t, size_t>> otherMinimizers;
	iterateMinimizersSimple(str, minimizerLength, windowSize, [&simpleMinimizers](size_t pos, size_t kmer) { simpleMinimizers.emplace_back(pos, kmer); });
	iterateMinimizersReal(str, minimizerLength, windowSize, [&otherMinimizers](size_t pos, size_t kmer) { otherMinimizers.emplace_back(pos, kmer); });
	std::sort(otherMinimizers.begin(), otherMinimizers.end(), [](std::pair<size_t, size_t> left, std::pair<size_t, size_t> right) { return left.first < right.first; });
	std::sort(simpleMinimizers.begin(), simpleMinimizers.end(), [](std::pair<size_t, size_t> left, std::pair<size_t, size_t> right) { return left.first < right.first; });
	// assert(simpleMinimizers.size() == otherMinimizers.size());
	for (size_t i = 0; i < otherMinimizers.size(); i++)
	{
		// assert(simpleMinimizers[i] == otherMinimizers[i]);
		callback(otherMinimizers[i].first, otherMinimizers[i].second);
	}
}

#endif

MinimizerSeeder::MinimizerSeeder(const AlignmentGraph& graph, size_t minimizerLength, size_t windowSize, size_t numThreads, double keepLeastFrequentFraction) :
graph(graph),
buckets(),
minimizerLength(minimizerLength),
windowSize(windowSize),
maxCount(0)
{
	assert(minimizerLength * 2 <= sizeof(size_t) * 8);
	assert(minimizerLength <= windowSize);
	initMinimizers(numThreads);
	initMaxCount(keepLeastFrequentFraction);
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

	std::unordered_map<size_t, size_t> nodeMinimizerStart;
	for (size_t i = 0; i < graph.NodeSize(); i++)
	{
		nodeMinimizerStart[graph.nodeIDs[i]] = std::max(nodeMinimizerStart[graph.nodeIDs[i]], (size_t)0);
		bool skipStart = false;
		for (auto n : graph.inNeighbors[i])
		{
			if (graph.nodeIDs[n] != graph.nodeIDs[i])
			{
				skipStart = true;
				break;
			}
		}
		if (skipStart)
		{
			nodeMinimizerStart[graph.nodeIDs[i]] = std::max(nodeMinimizerStart[graph.nodeIDs[i]], graph.nodeOffset[i]);
		}
	}

	for (size_t thread = 0; thread < numThreads; thread++)
	{
		threads.emplace_back([this, &nodeMinimizerStart, &positionDistributor, &threadsDone, &kmerPerBucket, &positionPerBucket, thread, numThreads, &nodeMutex, &nodeIter, positionSize](){
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
				iterateMinimizers(sequence, minimizerLength, windowSize, [this, &nodeMinimizerStart, &positionDistributor, &kmerPerBucket, &positionPerBucket, &vecPos, positionSize, thread, nodeId](size_t pos, size_t kmer)
				{
					if (pos < nodeMinimizerStart.at(nodeId)) return;
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
			buckets[thread].startPos.width(log2(kmerPerBucket[thread].size())+1);
			buckets[thread].startPos.resize(buckets[thread].locator->nbKeys() + 1);
			buckets[thread].kmerCheck.width(minimizerLength * 2);
			buckets[thread].kmerCheck.resize(buckets[thread].locator->nbKeys());
			sdsl::util::set_to_value(buckets[thread].startPos, 0);
			for (size_t i = 0; i < kmerPerBucket[thread].size(); i++)
			{
				uint64_t kmer = kmerPerBucket[thread][i];
				size_t index = buckets[thread].locator->lookup(kmer);
				buckets[thread].startPos[index] += 1;
				buckets[thread].kmerCheck[index] = kmer;
			}
			if (buckets[thread].startPos.size() > 0)
			{
				for (size_t i = 1; i < buckets[thread].startPos.size(); i++)
				{
					buckets[thread].startPos[i] += buckets[thread].startPos[i-1];
				}
				assert(buckets[thread].startPos[buckets[thread].startPos.size()-1] == kmerPerBucket[thread].size());
				buckets[thread].positions.resize(kmerPerBucket[thread].size());
				for (size_t i = 0; i < kmerPerBucket[thread].size(); i++)
				{
					size_t kmer = kmerPerBucket[thread][i];
					size_t index = buckets[thread].locator->lookup(kmer);
					assert(buckets[thread].startPos[index] > 0);
					buckets[thread].startPos[index] -= 1;
					size_t pos = buckets[thread].startPos[index];
					uint64_t insert = positionPerBucket[thread][i];
					buckets[thread].positions[pos] = insert;
				};
			}
		});
	}

	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
	threads.clear();
}

void MinimizerSeeder::addMinimizers(std::vector<SeedHit>& result, std::vector<std::tuple<size_t, size_t, size_t, size_t>>& matchIndices, size_t maxCount) const
{
	//prefer less common minimizers
	std::sort(matchIndices.begin(), matchIndices.end(), [this](const std::tuple<size_t, size_t, size_t, size_t>& left, const std::tuple<size_t, size_t, size_t, size_t>& right)
	{
		return std::get<3>(left) < std::get<3>(right);
	});
	size_t seedsHere = 0;
	size_t allowedCount = 0;
	for (auto match : matchIndices)
	{
		size_t bucket = std::get<1>(match);
		size_t start = std::get<2>(match);
		size_t end = start + std::get<3>(match);
		assert(end - start >= allowedCount);
		if (seedsHere >= maxCount && end - start > allowedCount) break;
		allowedCount = end - start;
		for (size_t i = start; i < end; i++)
		{
			size_t mergepos = buckets[bucket].positions[i];
			size_t nodeId = mergepos >> 6;
			size_t offset = mergepos & 63;
			result.push_back(matchToSeedHit(nodeId, offset, std::get<0>(match), std::get<3>(match)));
		}
		seedsHere += end - start;
	}
}

std::vector<SeedHit> MinimizerSeeder::getSeeds(const std::string& sequence, double density) const
{
	std::vector<std::tuple<size_t, size_t, size_t, size_t>> matchIndices;
	iterateKmers(sequence, minimizerLength, windowSize, [this, &matchIndices](size_t pos, size_t kmer)
	{
		size_t bucket = getBucket(kmer);
		assert(bucket < buckets.size());
		size_t index = buckets[bucket].locator->lookup(kmer);
		if (index == ULLONG_MAX) return;
		assert(index < buckets[bucket].kmerCheck.size());
		if (buckets[bucket].kmerCheck[(size_t)index] != kmer) return;
		size_t start = getStart(bucket, index);
		size_t end = getStart(bucket, index+1);
		size_t count = end - start;
		if (count >= maxCount) return;
		matchIndices.emplace_back(pos, bucket, start, count);
	});
	std::vector<SeedHit> result;
	size_t maxHits = sequence.size() * density;
	if (density == -1) maxHits = std::numeric_limits<size_t>::max();
	addMinimizers(result, matchIndices, maxHits);
	return result;
}

SeedHit MinimizerSeeder::matchToSeedHit(int nodeId, size_t nodeOffset, size_t seqPos, int count) const
{
	assert(nodeId < graph.nodeIDs.size());
	assert(nodeId < graph.nodeOffset.size());
	assert(nodeId < graph.reverse.size());
	SeedHit result { graph.nodeIDs[nodeId]/2, nodeOffset + graph.nodeOffset[nodeId], seqPos, minimizerLength, maxCount - count, graph.reverse[nodeId] };
	result.alignmentGraphNodeId = nodeId;
	result.alignmentGraphNodeOffset = nodeOffset;
	return result;
}

void MinimizerSeeder::initMaxCount(double keepLeastFrequentFraction)
{
	maxCount = 0;
	std::vector<size_t> counts;
	for (size_t bucket = 0; bucket < buckets.size(); bucket++)
	{
		if (buckets[bucket].locator->nbKeys() == 0) continue;
		for (size_t i = 0; i < buckets[bucket].locator->nbKeys()-1; i++)
		{
			counts.push_back(getStart(bucket, i+1) - getStart(bucket, i));
		}
	}
	std::sort(counts.begin(), counts.end());
	if (counts.size() == 0) return;
	size_t index = counts.size() * keepLeastFrequentFraction;
	if (index == counts.size()) index = counts.size()-1;
	maxCount = counts[index];
	maxCount += 1;
}

size_t MinimizerSeeder::getStart(size_t bucket, size_t index) const
{
	return buckets[bucket].startPos[index];
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

bool MinimizerSeeder::canSeed() const
{
	return maxCount > 0;
}