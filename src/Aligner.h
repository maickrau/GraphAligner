#ifndef Aligner_h
#define Aligner_h

#include <string>
#include <vector>
#include "AlignmentGraph.h"
#include "vg.pb.h"
#include "AlignmentSelection.h"

struct AlignerParams
{
	std::string graphFile;
	std::vector<std::string> fastqFiles;
	size_t numThreads;
	size_t initialBandwidth;
	bool dynamicRowStart;
	size_t maxCellsPerSlice;
	std::vector<std::string> seedFiles;
	std::string outputGAMFile;
	std::string outputJSONFile;
	std::string outputGAFFile;
	std::string outputCorrectedFile;
	std::string outputCorrectedClippedFile;
	bool verboseMode;
	bool tryAllSeeds;
	size_t mxmLength;
	size_t mumCount;
	size_t memCount;
	std::string seederCachePrefix;
	AlignmentSelection::SelectionMethod alignmentSelectionMethod;
	double selectionECutoff;
	bool compressCorrected;
	bool compressClipped;
	size_t minimizerLength;
	size_t minimizerWindowSize;
	double minimizerSeedDensity;
	size_t seedClusterMinSize;
	double minimizerDiscardMostNumerousFraction;
	double seedExtendDensity;
	double preciseClippingIdentityCutoff;
	int Xdropcutoff;
	size_t DPRestartStride;
	bool multiseedDP;
	double multimapScoreFraction;
	bool cigarMatchMismatchMerge;
	double minAlignmentScore;
};

void alignReads(AlignerParams params);
void replaceDigraphNodeIdsWithOriginalNodeIds(vg::Alignment& alignment, const AlignmentGraph& graph);

#endif
