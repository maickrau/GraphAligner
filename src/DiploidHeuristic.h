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
	phmap::flat_hash_set<std::tuple<size_t, int, int>> getForbiddenNodes(std::string sequence) const;
	void write(std::ostream& file) const;
	void read(std::istream& file);
	size_t getk() const;
private:
	void getHomologyPairs(const phmap::flat_hash_map<__uint128_t, uint8_t>& kmerCounts, const phmap::flat_hash_map<__uint128_t, std::pair<size_t, size_t>>& kmerPositions, const AlignmentGraph& seq);
	size_t k;
	phmap::flat_hash_map<__uint128_t, std::pair<size_t, size_t>> kmerImpliesNode;
	phmap::flat_hash_map<size_t, std::vector<size_t>> conflictPairs;
	phmap::flat_hash_map<size_t, size_t> nodeLengths;
};

class DiploidHeuristicSplitter
{
public:
	void initializePairs(const AlignmentGraph& graph, const std::vector<size_t>& kValues);
	std::vector<std::tuple<size_t, int, int>> getForbiddenNodes(std::string sequence) const;
	void write(std::string file) const;
	void read(std::string file);
	std::vector<size_t> getKValues() const;
private:
	std::vector<DiploidHeuristicSplitterOneK> splitters;
};

#endif
