#include <chrono>
#include "GfaGraph.h"
#include "fastqloader.h"
#include "MummerSeeder.h"

int main(int argc, char** argv)
{
	std::string graphFile { argv[1] };
	std::string readFile { argv[2] };
	size_t minlen = std::stoi(argv[3]);

	auto graph = GfaGraph::LoadFromFile(graphFile);
	auto reads = loadFastqFromFile(readFile);

	std::cerr << "build seeder" << std::endl;
	auto timeStart = std::chrono::system_clock::now();
	MummerSeeder seeder { graph, minlen };
	auto timeEnd = std::chrono::system_clock::now();
	size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
	std::cerr << "seeder built, " << time << "ms" << std::endl;
	timeStart = std::chrono::system_clock::now();
	auto seeds = seeder.getMumSeeds(reads[0].sequence, std::numeric_limits<size_t>::max());
	timeEnd = std::chrono::system_clock::now();
	time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
	std::cerr << "seeds gotten, " << time << "ms" << std::endl;
	for (auto seed : seeds)
	{
		std::cout << seed.nodeID << (seed.reverse ? "-" : "+") << "(" << seed.nodeOffset << ")," << seed.seqPos << "," << seed.matchLen << std::endl;
	}

}
