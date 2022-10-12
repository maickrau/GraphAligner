#ifndef DNAString_h
#define DNAString_h

#include <string>
#include <vector>

class DNAString
{
public:
	DNAString();
	DNAString(const std::string& str);
	DNAString& operator=(const std::string& str);
	size_t size() const;
	DNAString reverseComplement() const;
	std::string substr(size_t start, size_t length) const;
	std::string toString() const;
	void rewindIterators(size_t size) const;
private:
	size_t addString(const std::string& str);
	void buildFromString(const std::string& str);
	std::vector<uint64_t> storage;
	size_t realSize;
	size_t lastCharOffset;
	mutable size_t lastSubstringEnd;
	mutable size_t lastIndex;
	mutable size_t lastOffset;
};

#endif
