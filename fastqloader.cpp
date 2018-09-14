#include <algorithm>
#include <fstream>
#include <zstr.hpp> //https://github.com/mateidavid/zstr
#include "fastqloader.h"
#include "CommonUtils.h"

std::vector<FastQ> loadFastqFastqFromStream(std::istream& file, bool includeQuality)
{
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
		if (includeQuality) newread.quality = line;
		result.push_back(newread);
	} while (file.good());
	return result;
}

std::vector<FastQ> loadFastqFastaFromStream(std::istream& file, bool includeQuality)
{
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
		if (includeQuality)
		{
			for (size_t i = 0; i < newread.sequence.size(); i++)
			{
				newread.quality += '!';
			}
		}
		result.push_back(newread);
	} while (file.good());
	return result;
}

std::vector<FastQ> loadFastqFastqFromFile(std::string filename, bool includeQuality)
{
	std::ifstream file {filename};
	return loadFastqFastqFromStream(file, includeQuality);
}

std::vector<FastQ> loadFastqFastaFromFile(std::string filename, bool includeQuality)
{
	std::ifstream file {filename};
	return loadFastqFastaFromStream(file, includeQuality);
}

std::vector<FastQ> loadFastqFastqFromGzippedFile(std::string filename, bool includeQuality)
{
	zstr::ifstream file { filename };
	return loadFastqFastqFromStream(file, includeQuality);
}

std::vector<FastQ> loadFastqFastaFromGzippedFile(std::string filename, bool includeQuality)
{
	zstr::ifstream file { filename };
	return loadFastqFastaFromStream(file, includeQuality);
}

std::vector<FastQ> loadFastqFromFile(std::string filename, bool includeQuality)
{
	bool gzipped = false;
	std::string originalFilename = filename;
	if (filename.substr(filename.size()-3) == ".gz")
	{
		gzipped = true;
		filename = filename.substr(0, filename.size()-3);
	}
	bool fastq = false;
	bool fasta = false;
	if (filename.substr(filename.size()-6) == ".fastq") fastq = true;
	if (filename.substr(filename.size()-3) == ".fq") fastq = true;
	if (filename.substr(filename.size()-6) == ".fasta") fasta = true;
	if (filename.substr(filename.size()-3) == ".fa") fasta = true;
	if (fasta)
	{
		if (gzipped)
		{
			return loadFastqFastaFromGzippedFile(originalFilename, includeQuality);
		}
		else
		{
			return loadFastqFastaFromFile(originalFilename, includeQuality);
		}
	}
	if (fastq)
	{
		if (gzipped)
		{
			return loadFastqFastqFromGzippedFile(originalFilename, includeQuality);
		}
		else
		{
			return loadFastqFastqFromFile(originalFilename, includeQuality);
		}
	}
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
