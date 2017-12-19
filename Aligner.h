#ifndef Aligner_h
#define Aligner_h

#include <string>
#include <vector>
#include "vg.pb.h"

void alignReads(std::string graphFile, std::string fastqFile, int numThreads, int initialBandwidth, int rampBandwidth, std::string alignmentFile, std::string auggraphFile, int dynamicRowStart, std::string seedFile, bool sqrtSpace, bool stats);

#endif
