#include <iostream>
#include <unistd.h>
#include <fstream>
#include "Aligner.h"
#include "stream.hpp"
#include "ThreadReadAssertion.h"

int main(int argc, char** argv)
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;

    struct sigaction act;
    act.sa_handler = ThreadReadAssertion::signal;
    sigemptyset(&act.sa_mask);
    act.sa_flags = 0;
    sigaction(SIGSEGV, &act, 0);

	std::string graphFile = "";
	std::string fastqFile = "";
	std::string alignmentFile = "";
	std::string auggraphFile = "";
	std::string seedFile = "";
	int numThreads = 0;
	int initialBandwidth = 0;
	int rampBandwidth = 0;
	int dynamicRowStart = 64;
	int c;
	bool initialFullBand = false;

	while ((c = getopt(argc, argv, "g:f:a:t:B:A:is:d:MSb:")) != -1)
	{
		switch(c)
		{
			case 'g':
				graphFile = std::string(optarg);
				break;
			case 'f':
				fastqFile = std::string(optarg);
				break;
			case 'a':
				alignmentFile = std::string(optarg);
				break;
			case 't':
				numThreads = std::stoi(optarg);
				break;
			case 'b':
				initialBandwidth = std::stoi(optarg);
				break;
			case 'B':
				rampBandwidth = std::stoi(optarg);
				break;
			case 'A':
				auggraphFile = std::string(optarg);
				break;
			case 'i':
				initialFullBand = true;
				break;
			case 's':
				seedFile = std::string(optarg);
				break;
			case 'd':
				dynamicRowStart = std::stoi(optarg);
				break;
		}
	}

	if (dynamicRowStart % 64 != 0)
	{
		std::cerr << "dynamic row start has to be a multiple of 64" << std::endl;
		std::exit(0);
	}

	if (numThreads < 1)
	{
		std::cerr << "number of threads must be >= 1" << std::endl;
		std::exit(0);
	}

	if (initialBandwidth < 2)
	{
		std::cerr << "bandwidth must be >= 2" << std::endl;
		std::exit(0);
	}

	if (rampBandwidth != 0 && rampBandwidth <= initialBandwidth)
	{
		std::cerr << "backup bandwidth must be higher than initial bandwidth" << std::endl;
		std::exit(0);
	}

	if (!initialFullBand && seedFile == "")
	{
		std::cerr << "either initial full band or seed file must be set" << std::endl;
		std::exit(0);
	}

	alignReads(graphFile, fastqFile, numThreads, initialBandwidth, rampBandwidth, alignmentFile, auggraphFile, dynamicRowStart, seedFile);

	return 0;
}
