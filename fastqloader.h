#include <string>
#include <vector>

class FastQ {
public:
	std::string sequence;
	std::string quality;
};

std::vector<FastQ> loadFastqFromFile(std::string filename);
