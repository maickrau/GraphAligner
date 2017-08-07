#ifndef CommonUtils_h
#define CommonUtils_h

#include <string>
#include <vector>
#include "vg.pb.h"

namespace CommonUtils
{
	vg::Graph LoadVGGraph(std::string filename);
	std::string ReverseComplement(std::string original);
}

#endif
