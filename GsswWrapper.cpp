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
	int bandwidth = 0;
	int c;

	while ((c = getopt(argc, argv, "g:f:s:a:t:b:")) != -1)
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
				bandwidth = std::stoi(optarg);
				break;
		}
	}

	if (numThreads < 1)
	{
		std::cerr << "number of threads must be >= 1" << std::endl;
		std::exit(0);
	}

	if (bandwidth < 2)
	{
		std::cerr << "bandwidth must be >= 2" << std::endl;
		std::exit(0);
	}

	alignReads(graphFile, fastqFile, seedFile, numThreads, bandwidth, alignmentFile);

	return 0;
}
