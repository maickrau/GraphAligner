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

std::string FastQ::reverseComplement(std::string str)
{
	std::string result;
	for (int i = str.size()-1; i >= 0; i--)
	{
		switch (str[i])
		{
			case 'A':
			case 'a':
				result += 'T';
				break;
			case 'C':
			case 'c':
				result += 'G';
				break;
			case 'T':
			case 't':
				result += 'A';
				break;
			case 'G':
			case 'g':
				result += 'C';
				break;
			case 'N':
			case 'n':
				result += 'N';
				break;
		}
	}
	return result;
}

FastQ FastQ::reverseComplement() const
{
	FastQ result;
	result.sequence = reverseComplement(sequence);
	result.seq_id = seq_id;
	result.quality = quality;
	std::reverse(result.quality.begin(), result.quality.end());
	return result;
}
