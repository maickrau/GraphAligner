#ifndef ByteStuff_h
#define ByteStuff_h

#include <cstdint>
#include <tuple>

namespace ByteStuff
{
	void precalculateByteStuff();
	std::tuple<uint8_t, uint8_t, int8_t> VPVNChange(size_t scorediff, size_t sign, size_t low, size_t high);
	extern std::vector<std::tuple<uint8_t, uint8_t, int8_t>> precalcedVPVNChanges;
}

#endif
