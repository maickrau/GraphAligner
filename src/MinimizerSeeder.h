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
		typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
		typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;
		boophf_t* locator;
		sdsl::int_vector<0> kmerCheck;
		sdsl::int_vector<0> startPos;
		sdsl::int_vector<0> positions;
	};
public:
	MinimizerSeeder(const AlignmentGraph& graph, size_t minimizerLength, size_t windowSize, size_t numThreads, double keepLeastFrequentFraction);
	std::vector<SeedHit> getSeeds(const std::string& sequence, double density) const;
	bool canSeed() const;
private:
	void addMinimizers(std::vector<SeedHit>& result, std::vector<std::tuple<size_t, size_t, size_t, size_t>>& matchIndices, size_t maxCount) const;
	size_t getStart(size_t bucket, size_t index) const;
	size_t getBucket(size_t hash) const;
	SeedHit matchToSeedHit(int nodeId, size_t nodeOffset, size_t seqPos, int count) const;
	void initMinimizers(size_t numThreads);
	void initMaxCount(double keepLeastFrequentFraction);
	const AlignmentGraph& graph;
	std::vector<KmerBucket> buckets;
	size_t minimizerLength;
	size_t windowSize;
	size_t maxCount;
};

#endif
