#ifndef FastqLoader_H
#define FastqLoader_H

#include <string>
#include <vector>

class FastQ {
public:
	FastQ reverseComplement() const;
	std::string seq_id;
	std::string sequence;
	std::string quality;
};

std::vector<FastQ> loadFastqFromFile(std::string filename, bool includeQuality = true);

#endif
