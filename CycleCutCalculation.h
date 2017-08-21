#ifndef CycleCutCalculator_h
#define CycleCutCalculator_h

#include <vector>
#include <unordered_set>
#include <set>
#include <functional>
#include "AlignmentGraph.h"

class CycleCutCalculation
{
public:
	CycleCutCalculation(const AlignmentGraph& graph);
	AlignmentGraph::CycleCut GetCycleCutByIndex(size_t startNode, int wordSize) const;
	AlignmentGraph::CycleCut GetCycleCutBySupersequence(size_t startNode, int wordSize) const;
private:
	std::vector<size_t> getSupersequenceIndexingAndPredecessors(size_t cycleStart, int sizeLeft, std::vector<std::set<size_t>>& predecessors) const;
	std::vector<size_t> getCycleCuttersSupersequence(size_t cycleStart, int sizeLeft) const;
	void getPredecessorsFromSupersequence(size_t cycleStart, int sizeLeft, const std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const;
	void iterateOverCycleCuttingTreeRec(size_t cycleStart, size_t node, int sizeLeft, std::vector<size_t>& currentStack, std::function<void(const std::vector<size_t>&)> function) const;
	void iterateOverCycleCuttingTree(size_t cycleStart, int sizeLeft, std::function<void(const std::vector<size_t>&)> function) const;
	void iterateOverCycleCuttingCyclesRec(size_t cycleStart, size_t node, std::vector<size_t>& currentStack, std::unordered_set<size_t>& cycleful, std::unordered_set<size_t>& uncycleful, int sizeLeft, std::function<void(const std::vector<size_t>&)> cycleFunction) const;
	void iterateOverCycleCuttingCycles(size_t cycleStart, int sizeLeft, std::function<void(const std::vector<size_t>&)> cycleFunction, std::function<void(const std::set<size_t>&)> uncycleFunction) const;
	void getCycleCuttersByIndex(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const;
	void getCycleCuttersBySupersequence(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const;
	const AlignmentGraph& graph;
};

#endif
