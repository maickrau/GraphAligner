#ifndef FastqLoader_H
#define FastqLoader_H

#include <exception>
#include <string>
#include <vector>
#include <zstr.hpp> //https://github.com/mateidavid/zstr
extern "C"
{
#include <htslib/sam.h>
}

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
		// modified from https://github.com/samtools/htslib/blob/develop/samples/read_bam.c
		bam1_t *bamdata = nullptr;
		bamdata = bam_init1();
		if (bamdata == nullptr)
		{
			throw std::runtime_error { "Cannot initialize bam data" };
		}
		samFile *infile = NULL;
		infile = sam_open(filename.data(), "r");
		if (infile == nullptr)
		{
			bam_destroy1(bamdata);
			throw std::runtime_error { "Cannot open bam file" };
		}
		sam_hdr_t *in_samhdr = nullptr;
		in_samhdr = sam_hdr_read(infile);
		if (in_samhdr == nullptr)
		{
			sam_close(infile);
			bam_destroy1(bamdata);
			throw std::runtime_error { "Cannot read bam header" };
		}
		while (true)
		{
			FastQ newread;
			int ret_r = sam_read1(infile, in_samhdr, bamdata);
			if (ret_r < 0) break;
			// only primary
			if (bamdata->core.flag & BAM_FSECONDARY) continue;
			if (bamdata->core.flag & BAM_FSUPPLEMENTARY) continue;
			newread.seq_id = bam_get_qname(bamdata);
			uint8_t* data = bam_get_seq(bamdata);
			newread.sequence.resize(bamdata->core.l_qseq);
			for (size_t i = 0; i < newread.sequence.size(); i++)
			{
				newread.sequence[i] = seq_nt16_str[bam_seqi(data, i)];
			}
			if (includeQuality)
			{
				newread.quality.resize(newread.sequence.size());
				for (size_t i = 0; i < newread.sequence.size(); i++)
				{
					newread.quality[i] = bam_get_qual(bamdata)[i]+33;
				}
			}
			f(newread);
		}
		sam_hdr_destroy(in_samhdr);
		bam_destroy1(bamdata);
		sam_close(infile);
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
