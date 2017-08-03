#ifndef TopologicalSort_H
#define TopologicalSort_H

#include <vector>
#include "vg.pb.h"
#include "BigraphToDigraph.h"

std::vector<size_t> topologicalSort(const vg::Graph& graph);
std::vector<size_t> topologicalSort(const DirectedGraph& graph);
std::vector<size_t> topologicalSort(const std::vector<std::vector<size_t>>& graph);

#endif