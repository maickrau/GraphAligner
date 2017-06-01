#ifndef Aligner_h
#define Aligner_h

#include <string>
#include <vector>
#include "vg.pb.h"

void alignReads(std::string graphFile, std::string fastqFile, std::string seedFile, int numThreads, int bandwidth, std::string alignmentFile, std::string auggraphFile);

#endif
