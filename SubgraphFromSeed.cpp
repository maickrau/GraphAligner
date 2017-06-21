#include <functional>
#include <algorithm>
#include "vg.pb.h"
#include "ThreadReadAssertion.h"


// Program to find Dijkstra's shortest path using
// priority_queue in STL
#include <bits/stdc++.h>
using namespace std;
# define INF 0x3f3f3f3f


// iPair ==>  Integer Pair
typedef pair<int, int> iPair;

// This class represents a directed graph using
// adjacency list representation
class Graph
{
    int V;    // No. of vertices

    // In a weighted graph, we need to store vertex
    // and weight pair for every edge
    std::vector<list< pair<int, int>>> adj;

public:
    Graph(int V);  // Constructor

    // function to add an edge to graph
    void addEdge(int u, int v, int w);

    // prints shortest path from s
    std::vector<int> shortestPath(std::vector<int> s, int length);
};


vg::Graph ExtractSubgraph(const vg::Graph& graph, vg::Alignment seed, int length)
{
	Graph dijkstraGraph(graph.node_size());
	Graph reverseDijkstraGraph(graph.node_size());
	std::map<int, int> ids;
	for (int i = 0; i < graph.node_size(); i++)
	{
		ids[graph.node(i).id()] = i;
	}
	for (int i = 0; i < graph.edge_size(); i++)
	{
		dijkstraGraph.addEdge(ids[graph.edge(i).from()], ids[graph.edge(i).to()], graph.node(ids[graph.edge(i).to()]).sequence().size());
		reverseDijkstraGraph.addEdge(ids[graph.edge(i).to()], ids[graph.edge(i).from()], graph.node(ids[graph.edge(i).from()]).sequence().size());
	}
	std::vector<int> startNodes;
	for (int i = 0; i < seed.path().mapping_size(); i++)
	{
		startNodes.push_back(ids[seed.path().mapping(i).position().node_id()]);
	}
	auto forward = dijkstraGraph.shortestPath(startNodes, length);
	auto backward = reverseDijkstraGraph.shortestPath(startNodes, length);
	vg::Graph result;
	std::set<int> merged;
	merged.insert(forward.begin(), forward.end());
	merged.insert(backward.begin(), backward.end());
	for (int i = 0; i < graph.node_size(); i++)
	{
		if (merged.count(i) > 0)
		{
			auto newNode = result.add_node();
			newNode->set_sequence(graph.node(i).sequence());
			newNode->set_name(graph.node(i).name());
			newNode->set_id(graph.node(i).id());
		}
	}
	for (int i = 0; i < graph.edge_size(); i++)
	{
		if (merged.count(ids[graph.edge(i).from()]) > 0 && merged.count(ids[graph.edge(i).to()]) > 0)
		{
			auto newEdge = result.add_edge();
			newEdge->set_from(graph.edge(i).from());
			newEdge->set_to(graph.edge(i).to());
			newEdge->set_from_start(graph.edge(i).from_start());
			newEdge->set_to_end(graph.edge(i).to_end());
			newEdge->set_overlap(graph.edge(i).overlap());
		}
	}
	return result;
}


// Allocates memory for adjacency list
Graph::Graph(int V)
{
	this->V = V;
	adj.resize(V);
}

void Graph::addEdge(int u, int v, int w)
{
	adj[u].push_back(make_pair(v, w));
	adj[v].push_back(make_pair(u, w));
}

 //http://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-using-priority_queue-stl/
// Prints shortest paths from src to all other vertices
std::vector<int> Graph::shortestPath(std::vector<int> srcs, int length)
{
    // Create a priority queue to store vertices that
    // are being preprocessed. This is weird syntax in C++.
    // Refer below link for details of this syntax
    // http://geeksquiz.com/implement-min-heap-using-stl/
	priority_queue< iPair, vector <iPair> , greater<iPair> > pq;

    // Create a vector for distances and initialize all
    // distances as infinite (INF)
	vector<int> dist(V, INF);

    // Insert source itself in priority queue and initialize
    // its distance as 0.
	for (auto src : srcs)
	{
		pq.push(make_pair(0, src));
		dist[src] = 0;
	}

    /* Looping till priority queue becomes empty (or all
      distances are not finalized) */
	while (!pq.empty())
	{
        // The first vertex in pair is the minimum distance
        // vertex, extract it from priority queue.
        // vertex label is stored in second of pair (it
        // has to be done this way to keep the vertices
        // sorted distance (distance must be first item
        // in pair)
		int u = pq.top().second;
		pq.pop();

		if (dist[u] > length) break;

        // 'i' is used to get all adjacent vertices of a vertex
		list< pair<int, int> >::iterator i;
		for (i = adj[u].begin(); i != adj[u].end(); ++i)
		{
            // Get vertex label and weight of current adjacent
            // of u.
			int v = (*i).first;
			int weight = (*i).second;

            //  If there is shorted path to v through u.
			if (dist[v] > dist[u] + weight)
			{
                // Updating distance of v
				dist[v] = dist[u] + weight;
				pq.push(make_pair(dist[v], v));
			}
		}
	}

	std::vector<int> result;
	for (size_t i = 0; i < dist.size(); i++)
	{
		if (dist[i] <= length) result.push_back(i);
	}

	return result;
}
