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
	AlignmentGraph::CycleCut GetCycleCutSimplest(size_t startNode, int wordSize) const;
	AlignmentGraph::CycleCut GetCycleCutTooBig(size_t startNode, int wordSize) const;
	AlignmentGraph::CycleCut GetCycleCut(size_t startNode, int wordSize) const;
private:
	std::pair<std::unordered_set<size_t>, std::unordered_set<size_t>> splitCyclicAndNoncyclic(size_t cycleStart, int sizeLeft) const;
	std::unordered_set<size_t> getReachable(size_t cycleStart, size_t sizeLeft) const;
	void splitCyclicAndNoncyclicRec(std::vector<size_t>& stack, size_t currentNode, const std::unordered_set<size_t>& reachable, std::unordered_set<size_t>& visited, std::unordered_set<size_t>& cyclic) const;
	std::vector<size_t> getCycleCuttersSimplest(size_t cycleStart, int sizeLeft, std::vector<std::set<size_t>>& predecessors) const;
	std::vector<size_t> getCycleCuttersOrder(size_t cycleStart, int sizeLeft, std::vector<std::set<size_t>>& predecessors) const;
	void filterUnnecessaryCharacters(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors) const;
	std::vector<size_t> getCycleCuttersSupersequence(size_t cycleStart, int sizeLeft, const std::vector<std::vector<size_t>>& paths) const;
	void getPredecessorsFromSupersequenceOverEdgeCoveringPaths(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut, const std::vector<std::vector<size_t>>& paths) const;
	void iterateOverEdgeCoveringPaths(size_t cycleStart, int sizeLeft, std::function<void(const std::vector<size_t>&)> function) const;
	void getCycleCuttersSimplestff(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const;
	void getCycleCuttersTooBig(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const;
	void getCycleCutters(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut) const;
	const AlignmentGraph& graph;
};

#endif
