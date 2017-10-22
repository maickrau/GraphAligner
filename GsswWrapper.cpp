#include <iostream>
#include <unistd.h>
#include <fstream>
#include "Aligner.h"
#include "stream.hpp"

int main(int argc, char** argv)
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;
	std::string graphFile = "";
	std::string fastqFile = "";
	std::string alignmentFile = "";
	std::string auggraphFile = "";
	std::string seedFile = "";
	int startBandwidth = 0;
	int numThreads = 0;
	int dynamicWidth = 0;
	int dynamicRowStart = 64;
	int c;
	bool initialFullBand = false;
	bool sqrtSpace = false;

	while ((c = getopt(argc, argv, "g:f:a:t:B:A:b:is:d:M")) != -1)
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
			case 'B':
				dynamicWidth = std::stoi(optarg);
				break;
			case 'A':
				auggraphFile = std::string(optarg);
				break;
			case 'b':
				startBandwidth = std::stoi(optarg);
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
			case 'M':
				sqrtSpace = true;
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

	if (dynamicWidth < 2)
	{
		std::cerr << "dynamic bandwidth must be >= 2" << std::endl;
		std::exit(0);
	}

	if (!initialFullBand && (seedFile == "" || startBandwidth == 0))
	{
		std::cerr << "either initial full band or seed file and start bandwidth must be set" << std::endl;
		std::exit(0);
	}

	alignReads(graphFile, fastqFile, numThreads, dynamicWidth, alignmentFile, auggraphFile, dynamicRowStart, seedFile, startBandwidth, sqrtSpace);

	return 0;
}
