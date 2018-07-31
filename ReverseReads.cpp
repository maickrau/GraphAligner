#include <fstream>
#include "fastqloader.h"

int main(int argc, char** argv)
{
	auto reads = loadFastqFromFile(argv[1]);
	std::ofstream output {argv[2]};
	for (size_t i = 0; i < reads.size(); i++)
	{
		auto reverse = reads[i].reverseComplement();
		output << ">" << reverse.seq_id << "_Reverse" << "\n";
		output << reverse.sequence << "\n";
	}
}