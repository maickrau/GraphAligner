#ifndef FastqLoader_H
#define FastqLoader_H

#include <exception>
#include <string>
#include <vector>
#include <zstr.hpp> //https://github.com/mateidavid/zstr
#include <api/BamReader.h>

class FastQ {
public:
	template <typename F>
	static void streamFastqFastqFromStream(std::istream& file, bool includeQuality, F f)
	{
		do
		{
			std::string line;
			std::getline(file, line);
			if (!file.good()) break;
			if (line.size() == 0) continue;
			if (line[0] != '@') continue;
			FastQ newread;
			if (line.back() == '\r') line.pop_back();
			while (line.back() == ' ') line.pop_back();
			newread.seq_id = line.substr(1);
			std::getline(file, line);
			if (line.back() == '\r') line.pop_back();
			newread.sequence = line;
			std::getline(file, line);
			std::getline(file, line);
			if (line.back() == '\r') line.pop_back();
			while (line.back() == ' ') line.pop_back();
			if (includeQuality)
			{
				newread.quality = line;
				assert(newread.quality.size() == newread.sequence.size());
			}
			f(newread);
		} while (file.good());
	}
	template <typename F>
	static void streamFastqFastaFromStream(std::istream& file, bool includeQuality, F f)
	{
		std::string line;
		std::getline(file, line);
		do
		{
			if (line.size() == 0)
			{
				std::getline(file, line);
				continue;
			}
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
				line.clear();
				std::getline(file, line);
				if (line.size() == 0) continue;
				if (line[0] == '>') break;
				if (line.back() == '\r') line.pop_back();
				while (line.back() == ' ') line.pop_back();
				newread.sequence += line;
			} while (file.good());
			if (includeQuality)
			{
				for (size_t i = 0; i < newread.sequence.size(); i++)
				{
					newread.quality += '!';
				}
			}
			f(newread);
		} while (file.good());
	}
	template <typename F>
	static void streamFastqFastqFromFile(std::string filename, bool includeQuality, F f)
	{
		std::ifstream file {filename};
		streamFastqFastqFromStream(file, includeQuality, f);
	}
	template <typename F>
	static void streamFastqFastaFromFile(std::string filename, bool includeQuality, F f)
	{
		std::ifstream file {filename};
		streamFastqFastaFromStream(file, includeQuality, f);
	}
	template <typename F>
	static void streamFastqFastqFromGzippedFile(std::string filename, bool includeQuality, F f)
	{
		zstr::ifstream file { filename };
		streamFastqFastqFromStream(file, includeQuality, f);
	}
	template <typename F>
	static void streamFastqFastaFromGzippedFile(std::string filename, bool includeQuality, F f)
	{
		zstr::ifstream file { filename };
		streamFastqFastaFromStream(file, includeQuality, f);
	}
	template <typename F>
	static void streamBam(std::string filename, bool includeQuality, F f)
	{
		// https://github.com/pezmaster31/bamtools/wiki/Using-the-API
		BamTools::BamReader reader;
		if (!reader.Open(filename))
		{
			throw std::runtime_error { "Cannot open bam" };
		}
		if (!reader.IsOpen())
		{
			throw std::runtime_error { "Cannot open bam" };
		}
		BamTools::BamAlignment aln;
		while (reader.GetNextAlignment(aln))
		{
			// only primary
			if (!aln.IsPrimaryAlignment()) continue;
			if (aln.IsSupplementary()) continue;
			FastQ newread;
			newread.seq_id = aln.Name;
			newread.sequence = aln.QueryBases;
			if (includeQuality)
			{
				newread.quality = aln.Qualities;
			}
			f(newread);
		}
	}
	template <typename F>
	static void streamFastqFromFile(std::string filename, bool includeQuality, F f)
	{
		bool gzipped = false;
		std::string originalFilename = filename;
		if (filename.size() > 3 && filename.substr(filename.size()-3) == ".gz")
		{
			gzipped = true;
			filename = filename.substr(0, filename.size()-3);
		}
		bool fastq = false;
		bool fasta = false;
		if (filename.size() > 4 && (filename.substr(filename.size()-4) == ".bam" || filename.substr(filename.size()-4) == ".sam"))
		{
			streamBam(originalFilename, includeQuality, f);
			return;
		}
		if (filename.size() > 6 && filename.substr(filename.size()-6) == ".fastq") fastq = true;
		if (filename.size() > 3 && filename.substr(filename.size()-3) == ".fq") fastq = true;
		if (filename.size() > 6 && filename.substr(filename.size()-6) == ".fasta") fasta = true;
		if (filename.size() > 3 && filename.substr(filename.size()-3) == ".fa") fasta = true;
		if (fasta)
		{
			if (gzipped)
			{
				streamFastqFastaFromGzippedFile(originalFilename, includeQuality, f);
				return;
			}
			else
			{
				streamFastqFastaFromFile(originalFilename, includeQuality, f);
				return;
			}
		}
		if (fastq)
		{
			if (gzipped)
			{
				streamFastqFastqFromGzippedFile(originalFilename, includeQuality, f);
				return;
			}
			else
			{
				streamFastqFastqFromFile(originalFilename, includeQuality, f);
				return;
			}
		}
		std::cerr << "Input sequence file cannot be read: " << originalFilename << std::endl;
		std::abort();
	}
	FastQ reverseComplement() const;
	std::string seq_id;
	std::string sequence;
	std::string quality;
};

std::vector<FastQ> loadFastqFromFile(std::string filename, bool includeQuality = true);

#endif
