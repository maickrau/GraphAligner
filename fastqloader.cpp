#include <algorithm>
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
		FastQ newread;
		if (line.back() == '\r') line.pop_back();
		newread.seq_id = line.substr(1);
		std::getline(file, line);
		if (line.back() == '\r') line.pop_back();
		newread.sequence = line;
		std::getline(file, line);
		std::getline(file, line);
		if (line.back() == '\r') line.pop_back();
		newread.quality = line;
		result.push_back(newread);
	} while (file.good());
	return result;
}

FastQ FastQ::reverseComplement() const
{
	FastQ result;
	for (int i = sequence.size()-1; i >= 0; i--)
	{
		switch (sequence[i])
		{
			case 'A':
			case 'a':
				result.sequence += 'T';
				break;
			case 'C':
			case 'c':
				result.sequence += 'G';
				break;
			case 'T':
			case 't':
				result.sequence += 'A';
				break;
			case 'G':
			case 'g':
				result.sequence += 'C';
				break;
			case 'N':
			case 'n':
				result.sequence += 'N';
				break;
		}
	}
	result.quality = quality;
	std::reverse(result.quality.begin(), result.quality.end());
	return result;
}
