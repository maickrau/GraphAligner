#ifndef GfaGraph_h
#define GfaGraph_h

#include <istream>
#include <ostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>

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
	static GfaGraph LoadFromFile(std::string filename, bool allowVaryingOverlaps=false, bool warnAboutMissingNodes=false);
	static GfaGraph LoadFromStream(std::istream& stream, bool allowVaryingOverlaps=false, bool warnAboutMissingNodes=false);
	void SaveToFile(std::string filename) const;
	void SaveToStream(std::ostream& stream) const;
	void AddSubgraph(const GfaGraph& subgraph);
	GfaGraph GetSubgraph(const std::unordered_set<int>& ids) const;
	GfaGraph GetSubgraph(const std::unordered_set<int>& nodes, const std::unordered_set<std::pair<NodePos, NodePos>>& edges) const;
	std::string OriginalNodeName(int nodeId) const;
	void confirmDoublesidedEdges();
	std::unordered_map<int, std::string> nodes;
	std::unordered_map<NodePos, std::vector<NodePos>> edges;
	std::unordered_map<std::pair<NodePos, NodePos>, size_t> varyingOverlaps;
	size_t edgeOverlap;
	std::unordered_map<int, std::string> tags;
	std::unordered_map<int, std::string> originalNodeName;
private:
	void numberBackToIntegers();
};

#endif