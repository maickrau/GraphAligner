#ifndef Aligner_h
#define Aligner_h

#include <string>
#include <vector>

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
	bool linear;
};

void alignReads(AlignerParams params);
void wabiExperiments(AlignerParams params);

#endif
