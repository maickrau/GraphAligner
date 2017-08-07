#include <algorithm>
#include <fstream>
#include "fastqloader.h"
#include "CommonUtils.h"

std::vector<FastQ> loadFastqFastqFromFile(std::string filename)
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

std::vector<FastQ> loadFastqFastaFromFile(std::string filename)
{
	std::ifstream file {filename};
	std::vector<FastQ> result;
	std::string line;
	std::getline(file, line);
	do
	{
		if (line[0] != '>')
		{
			std::getline(file, line);
			continue;
		}
		FastQ newread;
		if (line.back() == '\r') line.pop_back();
		newread.seq_id = line.substr(1);
		newread.sequence = "";
		do
		{
			std::getline(file, line);
			if (line[0] == '>') break;
			if (line.back() == '\r') line.pop_back();
			newread.sequence += line;
		} while (file.good());
		for (size_t i = 0; i < newread.sequence.size(); i++)
		{
			newread.quality += '!';
		}
		result.push_back(newread);
	} while (file.good());
	return result;
}

std::vector<FastQ> loadFastqFromFile(std::string filename)
{
	if (filename.substr(filename.size()-6) == ".fastq") return loadFastqFastqFromFile(filename);
	if (filename.substr(filename.size()-3) == ".fq") return loadFastqFastqFromFile(filename);
	if (filename.substr(filename.size()-6) == ".fasta") return loadFastqFastaFromFile(filename);
	if (filename.substr(filename.size()-3) == ".fa") return loadFastqFastaFromFile(filename);
	return std::vector<FastQ>{};
}

FastQ FastQ::reverseComplement() const
{
	FastQ result;
	result.sequence = CommonUtils::ReverseComplement(sequence);
	result.seq_id = seq_id;
	result.quality = quality;
	std::reverse(result.quality.begin(), result.quality.end());
	return result;
}
