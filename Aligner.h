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
	std::string alignmentFile;
	std::string auggraphFile;
	int dynamicRowStart;
	std::string seedFile;
};

void alignReads(AlignerParams params);

#endif
