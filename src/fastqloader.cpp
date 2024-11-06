#include <algorithm>
#include <fstream>
#include "fastqloader.h"
#include "CommonUtils.h"

std::vector<FastQ> loadFastqFromFile(std::string filename, bool includeQuality)
{
	std::vector<FastQ> result;
	FastQ::streamFastqFromFile(filename, includeQuality, [&result](FastQ& fq) {
		result.emplace_back(std::move(fq));
	});
	return result;
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

std::string FastQ::lowercase(std::string str)
{
	for (size_t i = 0; i < str.size(); i++)
	{
		if (str[i] >= 'A' && str[i] <= 'Z')
		{
			str[i] += 'a'-'A';
		}
	}
	return str;
}
