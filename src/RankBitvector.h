#ifndef RankBitvector_h
#define RankBitvector_h

#include <cstdint>
#include <vector>
#include <iostream>

// "Optimized Succinct Data Structures for Massive Data", Simon Gog, Matthias Petri 2014.
// Section 3: A CACHE FRIENDLY RANK IMPLEMENTATION FOR UNCOMPRESSED BITVECTORS
class RankBitvector
{
public:
	RankBitvector();
	RankBitvector(size_t size);
	void buildRanks();
	void resize(size_t newSize);
	size_t size() const;
	bool get(size_t index) const;
	void set(size_t index, bool value);
	size_t rankOne(size_t index) const;
	size_t rankZero(size_t index) const;
	void clear();
	void save(std::ostream& stream) const;
	void load(std::istream& stream);
	bool operator==(const RankBitvector& other) const;
	bool operator!=(const RankBitvector& other) const;
	bool operator[](size_t index) const;
private:
	std::vector<uint64_t> values;
	bool ranksBuilt;
	uint64_t realSize;
};

#endif
