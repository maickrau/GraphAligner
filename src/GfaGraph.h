#ifndef GfaGraph_h
#define GfaGraph_h

#include <istream>
#include <ostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "DNAString.h"

class NodePos
{
public:
	NodePos();
	NodePos(int id, bool end);
	int id;
	bool end;
	NodePos Reverse() const;
	bool operator==(const NodePos& other) const;
	bool operator!=(const NodePos& other) const;
};

namespace std 
{
	template <> 
	struct hash<NodePos>
	{
		size_t operator()(const NodePos& x) const
		{
			return hash<int>()(x.id) ^ hash<bool>()(x.end);
		}
	};
	template <> 
	struct hash<std::pair<NodePos, NodePos>>
	{
		size_t operator()(const std::pair<NodePos, NodePos>& x) const
		{
			// simple hashing with hash<NodePos>()(x.first) ^ hash<NodePos>()(x.second) collides each edge formed like (x -> x+1)
			// instead: 
			// https://stackoverflow.com/questions/682438/hash-function-providing-unique-uint-from-an-integer-coordinate-pair
			// https://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
			// and arbitrarily ignore directionality
			size_t pairing = .5 * (x.first.id + x.second.id) * (x.first.id + x.second.id + 1) + x.second.id;
			return hash<size_t>()(pairing);
		}
	};
}

class GfaGraph
{
public:
	GfaGraph();
	static GfaGraph LoadFromFile(std::string filename);
	static GfaGraph LoadFromStream(std::istream& stream);
	std::string OriginalNodeName(int nodeId) const;
	size_t totalBp() const;
	std::vector<DNAString> nodes;
	std::vector<std::tuple<NodePos, NodePos, size_t>> edges;
	std::vector<std::string> originalNodeName;
private:
};

#endif