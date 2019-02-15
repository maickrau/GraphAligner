#ifndef Aligner_h
#define Aligner_h

#include <string>
#include <vector>
#include "AlignmentGraph.h"
#include "vg.pb.h"

struct AlignerParams
{
	std::string graphFile;
	std::vector<std::string> fastqFiles;
	size_t numThreads;
	size_t initialBandwidth;
	size_t rampBandwidth;
	int dynamicRowStart;
	size_t maxCellsPerSlice;
	std::vector<std::string> seedFiles;
	std::string outputAlignmentFile;
	bool verboseMode;
	bool tryAllSeeds;
	bool highMemory;
	size_t mxmLength;
	size_t mumCount;
	size_t memCount;
	bool outputAllAlns;
	std::string seederCachePrefix;
	bool forceGlobal;
	bool outputJSON;
	bool preciseClipping;
	size_t minimizerCount;
	size_t minimizerLength;
	size_t minimizerWindowSize;
};

void alignReads(AlignerParams params);
void replaceDigraphNodeIdsWithOriginalNodeIds(vg::Alignment& alignment, const AlignmentGraph& graph);

#endif
