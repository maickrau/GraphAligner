#include <iostream>
#include <unistd.h>
#include <fstream>
#include <limits>
#include "Aligner.h"
#include "stream.hpp"
#include "ThreadReadAssertion.h"
#include "ByteStuff.h"

int main(int argc, char** argv)
{
#ifndef NOBUILTINPOPCOUNT
	if (__builtin_cpu_supports("popcnt") == 0)
	{
		std::cerr << "CPU does not support builtin popcount operation" << std::endl;
		std::cerr << "recompile with -DNOBUILTINPOPCOUNT" << std::endl;
		std::abort();
	}
#endif

    struct sigaction act;
    act.sa_handler = ThreadReadAssertion::signal;
    sigemptyset(&act.sa_mask);
    act.sa_flags = 0;
    sigaction(SIGSEGV, &act, 0);

    AlignerParams params;
	params.graphFile = "";
	params.fastqFile = "";
	params.seedFile = "";
	params.numThreads = 0;
	params.initialBandwidth = 0;
	params.rampBandwidth = 0;
	params.dynamicRowStart = 64;
	params.maxCellsPerSlice = std::numeric_limits<decltype(params.maxCellsPerSlice)>::max();
	bool initialFullBand = false;
	params.quietMode = false;
	params.sloppyOptimizations = false;
	int c;

	while ((c = getopt(argc, argv, "g:f:")) != -1)
	{
		switch(c)
		{
			case 'g':
				params.graphFile = std::string(optarg);
				break;
			case 'f':
				params.fastqFile = std::string(optarg);
				break;
		}
	}

	ByteStuff::precalculateByteStuff();

	wabiExperiments(params);
	// alignReads(params);

	return 0;
}
