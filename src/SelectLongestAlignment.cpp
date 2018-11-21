#include <fstream>
#include <string>
#include <unordered_map>
#include "CommonUtils.h"
#include "fastqloader.h"
#include "stream.hpp"

int main(int argc, char** argv)
{
	std::string alnfile { argv[1] };
	std::string outfile { argv[2] };

	auto alns = CommonUtils::LoadVGAlignments(alnfile);

	std::unordered_map<std::string, vg::Alignment> result;
	for (auto aln : alns)
	{
		if (result.count(aln.name()) == 0) result[aln.name()] = aln;
		else if (aln.sequence().size() > result[aln.name()].sequence().size()) result[aln.name()] = aln;
		else if (aln.sequence().size() == result[aln.name()].sequence().size() && aln.score() < result[aln.name()].score()) result[aln.name()] = aln;
	}

	std::vector<vg::Alignment> writeAlns;
	writeAlns.reserve(result.size());
	for (auto pair : result)
	{
		writeAlns.push_back(pair.second);
	}

	std::ofstream resultFile { outfile, std::ios::out | std::ios::binary };
	stream::write_buffered(resultFile, writeAlns, 0);
}
