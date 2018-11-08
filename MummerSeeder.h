#ifndef MummerSeeder_h
#define MummerSeeder_h

#include <vector>
#include <string>
#include <mummer/sparseSA.hpp>
#include <mummer/fasta.hpp>
#include "GfaGraph.h"
#include "GraphAlignerWrapper.h"

class MummerSeeder
{
public:
	MummerSeeder(const GfaGraph& graph, size_t minLen);
	std::vector<SeedHit> getMumSeeds(std::string sequence, size_t maxCount) const;
private:
	std::vector<SeedHit> matchesToSeeds(const std::vector<mummer::mummer::match_t>& matches) const;
	size_t getNodeIndex(size_t indexPos) const;
	void constructTree(const GfaGraph& graph, size_t minLen);
	size_t minLen;
	std::unique_ptr<mummer::mummer::sparseSA> matcher;
	std::vector<size_t> nodePositions;
	std::vector<int> nodeIDs;
	std::vector<bool> nodeReverse;
};

#endif
