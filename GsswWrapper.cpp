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
	std::string seedFile = "";
	std::string alignmentFile = "";
	int numThreads = 0;
	int startBandwidth = 0;
	int dynamicWidth = 0;
	int c;

	while ((c = getopt(argc, argv, "g:f:s:a:t:b:B:")) != -1)
	{
		switch(c)
		{
			case 'g':
				graphFile = std::string(optarg);
				break;
			case 'f':
				fastqFile = std::string(optarg);
				break;
			case 's':
				seedFile = std::string(optarg);
				break;
			case 'a':
				alignmentFile = std::string(optarg);
				break;
			case 't':
				numThreads = std::stoi(optarg);
				break;
			case 'b':
				startBandwidth = std::stoi(optarg);
				break;
			case 'B':
				dynamicWidth = std::stoi(optarg);
				break;
		}
	}

	if (numThreads < 1)
	{
		std::cerr << "number of threads must be >= 1" << std::endl;
		std::exit(0);
	}

	if (startBandwidth < 2)
	{
		std::cerr << "starting bandwidth must be >= 2" << std::endl;
		std::exit(0);
	}

	if (dynamicWidth < 2)
	{
		std::cerr << "dynamic bandwidth must be >= 2" << std::endl;
		std::exit(0);
	}

	alignReads(graphFile, fastqFile, seedFile, numThreads, startBandwidth, dynamicWidth, alignmentFile);

	return 0;
}
