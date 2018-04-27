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
};

void alignReads(AlignerParams params);
void wabiExperiments(AlignerParams params);

#endif
