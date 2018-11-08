#ifndef Aligner_h
#define Aligner_h

#include <string>
#include <vector>
#include "vg.pb.h"

struct AlignerParams
{
	std::string graphFile;
	std::string fastqFile;
	int numThreads;
	int initialBandwidth;
	int rampBandwidth;
	int dynamicRowStart;
	size_t maxCellsPerSlice;
	std::string seedFile;
	std::string outputAlignmentFile;
	bool quietMode;
	bool sloppyOptimizations;
	bool highMemory;
	bool useSubgraph;
	size_t mxmLength;
	size_t mumCount;
	size_t memCount;
	bool outputAllAlns;
	std::string seederCachePrefix;
};

void alignReads(AlignerParams params);

#endif
