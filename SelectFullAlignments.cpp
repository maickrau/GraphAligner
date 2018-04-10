#include <fstream>
#include <string>
#include <unordered_map>
#include "CommonUtils.h"
#include "fastqloader.h"
#include "stream.hpp"

int main(int argc, char** argv)
{
	std::string alnfile { argv[1] };
	std::string readfile { argv[2] };
	std::string outfile { argv[3] };

	auto alns = CommonUtils::LoadVGAlignments(alnfile);
	auto reads = loadFastqFromFile(readfile);

	std::unordered_map<std::string, int> readLens;
	for (auto read : reads)
	{
		readLens[read.seq_id] = read.sequence.size();
	}

	std::vector<vg::Alignment> result;
	for (auto aln : alns)
	{
		if (aln.sequence().size() == readLens[aln.name()] - 1) result.push_back(aln);
	}

	std::ofstream resultFile { outfile, std::ios::out | std::ios::binary };
	stream::write_buffered(resultFile, result, 0);
}
