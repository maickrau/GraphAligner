#include <fstream>
#include "fastqloader.h"

std::vector<FastQ> loadFastqFromFile(std::string filename)
{
	std::ifstream file {filename};
	std::vector<FastQ> result;
	do
	{
		std::string line;
		std::getline(file, line);
		if (line[0] != '@') continue;
		std::getline(file, line);
		line = line.substr(0, 3);
		if (line.back() == '\r') line.pop_back();
		FastQ newread;
		newread.sequence = line;
		std::getline(file, line);
		std::getline(file, line);
		line = line.substr(0, 3);
		if (line.back() == '\r') line.pop_back();
		newread.quality = line;
		result.push_back(newread);
	} while (file.good());
	return result;
}