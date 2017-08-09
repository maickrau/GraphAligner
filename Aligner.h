#ifndef Aligner_h
#define Aligner_h

#include <string>
#include <vector>
#include "vg.pb.h"

void alignReads(std::string graphFile, std::string fastqFile, int numThreads, int dynamicWidth, std::string alignmentFile, std::string auggraphFile, int dynamicRowStart, std::string seedFile, int startBandwidth, std::string mfvsFilename, std::string orderFilename, std::string cycleCutFilename);

#endif
