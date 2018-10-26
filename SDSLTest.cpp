#include "STSeeder.h"
#include "GfaGraph.h"
#include "fastqloader.h"

int main(int argc, char** argv)
{
	std::string graphFile { argv[1] };
	std::string readFile { argv[2] };

	auto graph = GfaGraph::LoadFromFile(argv[1]);
	auto reads = loadFastqFromFile(readFile);
	STSeeder seeder { graph };
	for (auto read : reads)
	{
		seeder.getSeeds(read.sequence, 1);
	}
}