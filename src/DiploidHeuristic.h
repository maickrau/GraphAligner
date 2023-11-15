#ifndef DiploidHeuristic_h
#define DiploidHeuristic_h

#include <vector>
#include <tuple>
#include <cstddef>
#include <phmap.h>
#include <fstream>
#include "AlignmentGraph.h"

class DiploidHeuristicSplitterOneK
{
public:
	void initializePairs(const AlignmentGraph& graph, size_t k);
	phmap::flat_hash_set<size_t> getForbiddenNodes(std::string sequence) const;
	void write(std::ostream& file) const;
	void read(std::istream& file);
	size_t getk() const;
private:
	void getHomologyPairs(const phmap::flat_hash_map<__uint128_t, uint8_t>& kmerCounts, const phmap::flat_hash_map<__uint128_t, std::pair<size_t, size_t>>& kmerPositions, const AlignmentGraph& seq);
	size_t k;
	phmap::flat_hash_map<__uint128_t, size_t> kmerForbidsNode;
	std::vector<std::pair<size_t, size_t>> homologyPairs;
};

class DiploidHeuristicSplitter
{
public:
	void initializePairs(const AlignmentGraph& graph, const std::vector<size_t>& kValues);
	std::vector<size_t> getForbiddenNodes(std::string sequence) const;
	void write(std::string file) const;
	void read(std::string file);
	std::vector<size_t> getKValues() const;
private:
	std::vector<DiploidHeuristicSplitterOneK> splitters;
};

#endif
