#ifndef SubgraphFromSeed_H
#define SubgraphFromSeed_H

#include "vg.pb.h"
vg::Graph ExtractSubgraph(const vg::Graph& graph, vg::Alignment seed, int length);

#endif