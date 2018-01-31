#ifndef GfaGraph_h
#define GfaGraph_h

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
			return hash<NodePos>()(x.first) ^ hash<NodePos>()(x.second);
		}
	};
}

class GfaGraph
{
public:
	static GfaGraph LoadFromFile(std::string filename);
	void SaveToFile(std::string filename) const;
	void AddSubgraph(const GfaGraph& subgraph);
	GfaGraph GetSubgraph(const std::unordered_set<int>& ids) const;
	GfaGraph GetSubgraph(const std::unordered_set<int>& nodes, const std::unordered_set<std::pair<NodePos, NodePos>>& edges) const;
	std::unordered_map<int, std::string> nodes;
	std::unordered_map<NodePos, std::vector<NodePos>> edges;
	int edgeOverlap;
private:
	GfaGraph();
};

#endif