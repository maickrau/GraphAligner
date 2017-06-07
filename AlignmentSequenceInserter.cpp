#include <fstream>
#include "stream.hpp"
#include "vg.pb.h"
#include "fastqloader.h"

int main(int argc, char** argv)
{
	auto reads = loadFastqFromFile(argv[1]);
	std::vector<vg::Alignment> alignments;
	std::ifstream seedfile { argv[2], std::ios::in | std::ios::binary };
	std::function<void(vg::Alignment&)> alignmentLambda = [&alignments](vg::Alignment& a) {
		alignments.push_back(a);
	};
	stream::for_each(seedfile, alignmentLambda);
	std::map<std::string, FastQ> readsmap;
	for (size_t i = 0; i < reads.size(); i++)
	{
		readsmap[reads[i].seq_id] = reads[i];
	}
	for (size_t i = 0; i < alignments.size(); i++)
	{
		alignments[i].set_sequence(readsmap[alignments[i].name()].sequence);
	}

	std::ofstream alignmentOut { argv[3], std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, alignments, 0);
}