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
	static std::string reverseComplement(std::string str);
};

std::vector<FastQ> loadFastqFromFile(std::string filename);

#endif
