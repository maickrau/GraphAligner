#ifndef MinimizerSeeder_h
#define MinimizerSeeder_h

#include <random>
#include <vector>
#include <string>
#include <sdsl/int_vector.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <ParallelBB.h>
#include "AlignmentGraph.h"
#include "GraphAlignerWrapper.h"
#include "BooPHF.h"

class MinimizerSeeder
{
	struct KmerBucket
	{
		KmerBucket();
		KmerBucket(const KmerBucket& other) = delete;
		KmerBucket(KmerBucket&& other) = default;
		~KmerBucket();
		KmerBucket& operator=(const KmerBucket& other) = delete;
		KmerBucket& operator=(KmerBucket&& other) = default;
		typedef boomphf::SingleHashFunctor<size_t> hasher_t;
		typedef boomphf::mphf<size_t, hasher_t> boophf_t;
		boophf_t* locator;
		sdsl::int_vector<0> kmerCheck;
		sdsl::bit_vector starts;
		sdsl::bit_vector::select_1_type startSelector;
		sdsl::int_vector<0> positions;
	};
public:
	MinimizerSeeder(const AlignmentGraph& graph, size_t minimizerLength, size_t windowSize, size_t numThreads);
	std::vector<SeedHit> getSeeds(const std::string& sequence, size_t maxCount, size_t chunkSize) const;
private:
	size_t getStart(size_t bucket, size_t index) const;
	size_t getBucket(size_t hash) const;
	SeedHit matchToSeedHit(int nodeId, size_t nodeOffset, size_t seqPos, int count) const;
	void initMinimizers(size_t numThreads);
	void initMaxCount();
	const AlignmentGraph& graph;
	std::vector<KmerBucket> buckets;
	size_t minimizerLength;
	size_t windowSize;
	size_t maxCount;
};

#endif
