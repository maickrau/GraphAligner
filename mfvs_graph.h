/**
 * Copyright 2012 Alexandre Blondin Masse
 *
 * This file is part of the MFVS project.
 *
 * MFVS is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or any later version.

 * MFVS is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <stack>

namespace mfvs {

using namespace std;

/**************************
 * For representing edges *
 **************************/
struct Edge {
    int source;
    int target;
    Edge() {source = -1; target = -1;}
};

/**********************************
 * Comparison operators for edges *
 **********************************/
bool operator==(const Edge& e1, const Edge& e2);
bool operator<(const Edge& e1, const Edge& e2);
bool operator>(const Edge& e1, const Edge& e2);
bool operator<=(const Edge& e1, const Edge& e2);
bool operator>=(const Edge& e1, const Edge& e2);

/**
 * Class Graph.
 *
 * Handles basic operations on directed graphs and in particular
 * operations involved in finding minimum feedback vertex sets
 * (see {@link http://en.wikipedia.org/wiki/Feedback_vertex_set}).
 *
 * @author   Alexandre Blondin Masse
 * @version  1.0
 */
class Graph {

    public:
        /*************
         * Constants *
         *************/
        static string TYPE_ARCS_LIST;
        static string TYPE_ADJACENCY_LIST;


        /*******************************
         * Constructors and destructor *
         *******************************/

        /**
         * Returns an instance of a directed graph having at most
         * the given number of vertices
         *
         * @param maxSize  the maximum number of vertices of
         *                 this graph
         */
        Graph(int maxSize);

        /**
         * Returns an instance of a directed graph by loading
         * its content from the given file
         *
         * @param filename  the name of the file from which
         *                  the graph will be loaded
         * @param type      the type of file, must be one of
         *                  the values "arcs list" and
         *                  "adjacency list"
         */
        Graph(const string& filename, const string& type);

        /**
         * Deletes this graph
         */
        ~Graph();

        /*******************
         * Basic accessors *
         *******************/

        /**
         * Returns the number of vertices of this graph.
         *
         * @return  the number of vertices of the graph
         */
        int numVertices() const;

        /**
         * Returns the number of edges of this graph.
         *
         * @return  the number of edges of the graph
         */
        int numEdges() const;

        /**
         * Tells whether the given vertex is in this graph
         *
         * @param vertex  the vertex to be checked
         * @return        true if the given vertex belongs to this
         *                graph
         */
        bool hasVertex(int vertex) const;

        /**
         * Tells whether the given edge is in this graph
         *
         * @param source  the source vertex of the edge
         * @param target  the target vertex of the edge
         * @return        true if the edge (source, target)
         *                belongs to this graph
         */
        bool hasEdge(int source, int target) const;

        /**
         * Returns all the vertices of this graph
         *
         * @return  a vector containing all the vertices
         */
        set<int> const& vertices() const;

        /**
         * Returns all the edges of this graph
         *
         * @return  a vector of pointers to vectors
         *          representing the adjacency list
         */
        vector<Edge> edges() const;

        /**
         * Returns the incoming neighbors of the given vertex v,
         * i.e. the vertices u such that (u,v) is an edge of
         * the graph
         *
         * @param vertex  the vertex for which the incoming
         *                neighbors are computed
         * @return        a vector containing the incoming
         *                neighbors of the vertex
         */
        vector<int> const& incomingNeighbors(int vertex) const;

        /**
         * Returns the outgoing neighbors of the given vertex v,
         * i.e. the vertices u such that (v,u) is an edge of
         * the graph
         *
         * @param vertex  the vertex for which the outgoing
         *                neighbors are computed
         * @return        a vector containing the outgoing
         *                neighbors of the vertex
         */
        vector<int> const& outgoingNeighbors(int vertex) const;

        /**
         * Displays a short description of this graph to
         * the screen. More precisely, displays the number
         * of vertices and edges of the graph.
         *
         * @param edges  if true, displays a list of the edges
         *               as well
         */
        void print(bool edges) const;

        /**
         * Same as calling print(false), i.e. displays this
         * graph without the list of its edges.
         */
        void print() const;

        /*******************
         * Basic modifiers *
         *******************/
        /**
         * Adds the given vertex to the graph if it is allowed.
         * If the vertex is already in the graph, it does nothing.
         *
         * @param vertex  the vertex to be added to the graph.
         *                Must be between 0 and maxSize.
         */
        void addVertex(int vertex);

        /**
         * Adds the given edge to the graph if it is allowed.
         * If the edge is already in the graph, it does nothing.
         *
         * @param source  the source vertex of the edge to be
         *                added to the graph. Must be a vertex
         *                of the graph.
         * @param target  the target vertex of the edge to be
         *                added to the graph. Must be a vertex
         *                of the graph.
         */
        void addEdge(int source, int target);

        /**
         * Deletes the given vertex from the graph, if it exists.
         * If the vertex is not in the graph, it does nothing.
         *
         * @param vertex  the vertex to be deleted.
         */
        void deleteVertex(int vertex);

        /**
         * Deletes the given vertices from the graph, if they
         * exist.
         *
         * @param vertices  the vertices to be deleted.
         */
        void deleteVertices(const vector<int>& vertices);

        /**
         * Deletes the edge (source, target) from the graph, if
         * it exists. If the edge is not in the graph, it does
         * nothing.
         *
         * @param source  the source vertex of the edge to be
         *                removed from the graph
         * @param target  the target vertex of the edge to be
         *                removed from the graph
         */
        void deleteEdge(int source, int target);

        /**
         * Deletes the given edges from the graph, if they
         * exist.
         *
         * @param edges  a vector containing the edges to
         *               be removed from the graph
         */
        void deleteEdges(const vector<Edge>& edges);

        /**
         * Merges the given vertex in this graph.
         *
         * Merging a vertex consists in deleting it while
         * preserving the edges between the neighbors. More
         * precisely, for every incoming neighbor u and every
         * outgoing neighbor w of the merged vertex v, we
         * add the edge (u, w) in the resulting graph.
         *
         * @param vertex  the vertex to be merged
         */
        void mergeVertex(int vertex);

        /**
         * Merges the given vertices in this graph.
         *
         * Merging a set of vertices consists in transforming
         * them into a unique vertex while preserving the
         * edge connections. More precisely, the resulting
         * vertex v will have the incoming neighbor u if
         * u is an incoming neighbor of any of the merged
         * vertices and every will have the outgoing neighbor
         * w if w is an outgoing neighbor of any merged vertices.
         * The resulting merged vertex is the first occuring in the
         * given vector.
         *
         * @param vertex  the set of vertices to be merged
         */
        void mergeVertices(vector<int> vertices);

        /**
         * Merges the two given vertices in this graph.
         * See mergeVertex(const vector<int>) for more
         * details.
         *
         * @param vertex1  the first vertex to be merged
         * @param vertex2  the second vertex to be merged
         */
        void mergeVertices(int vertex1, int vertex2);

        /*************
         * Subgraphs *
         *************/

        /**
         * Returns the subgraph induced by the given
         * vertices.
         *
         * Given a directed graph G = (V,E) and a subset
         * U of V, the subgraph of G induced by U is the
         * graph whose set of vertices is U and whose set
         * of edges is E intersected with U x U, i.e. edges
         * belonging to G whose end points are both in U.
         *
         * @param vs  the vertices whose induced subgraph
         *            is queried
         * @return    a graph corresponding to the subgraph
         *            induced by the given vertices
         */
        Graph subgraph(vector<int> vs) const;

        /*************************************
         * Degrees, sources, sinks and loops *
         *************************************/

        /**
         * Returns the in-degree of the given vertex, i.e.
         * its number of incoming neighbors.
         *
         * @param vertex  the vertex whose in-degree is computed
         * @return        a non negative integer representing the
         *                in-degree of the vertex, -1 if the vertex
         *                does not belong to the graph
         */
        int inDegree(int vertex) const;

        /**
         * Returns the out-degree of the given vertex, i.e.
         * its number of outgoing neighbors.
         *
         * @param vertex  the vertex whose out-degree is computed
         * @return        a non negative integer representing the
         *                out-degree of the vertex, -1 if the vertex
         *                does not belong to the graph
         */
        int outDegree(int vertex) const;

        /**
         * Returns the total degree of the given vertex, i.e.
         * the sum of its in-degree and out-degree.
         *
         * @param vertex  the vertex whose degree is computed
         * @return        the total degree of the given vertex
         */
        int degree(int vertex) const;

        /**
         * Returns a vertex having minimum total degree.
         *
         * @return  a vertex whose total degree is minimum
         */
        int minDegreeVertex() const;

        /**
         * Returns a vertex having maximum total degree.
         *
         * @return  a vertex whose total degree is maximum
         */
        int maxDegreeVertex() const;

        /**
         * Tells whether the given vertex is a source of this
         * graph, i.e. it has null in-degree.
         *
         * @param vertex  the vertex to be checked
         * @return        true if the given vertex is a source
         *                of this graph
         */
        bool isSource(int vertex) const;

        /**
         * Tells whether the given vertex is a sink of this
         * graph, i.e. it has null out-degree.
         *
         * @param vertex  the vertex to be checked
         * @return        true if the given vertex is a sink
         *                of this graph
         */
        bool isSink(int vertex) const;

        /**
         * Tells whether the given vertex has a self-loop in
         * this graph.
         *
         * @param vertex  the vertex to be checked
         * @return        true if the given vertex has a
         *                self-loop
         */
        bool hasLoop(int vertex) const;

        /**
         * Returns the sources of this graph, i.e. the
         * vertices having null in-degree.
         *
         * @return  a vector containing the sources of this
         *          graph
         */
        vector<int> sources() const;

        /**
         * Returns the sinks of this graph, i.e. the
         * vertices having null out-degree.
         *
         * @return  a vector containing the sinks of this
         *          graph
         */
        vector<int> sinks() const;

        /**
         * Returns the vertices of this graph having a
         * self-loop.
         *
         * @return  a vector containing the vertices having
         *          a self-loop in this graph
         */
        vector<int> loops() const;


        /*********************************
         * Strongly connected components *
         *********************************/

        /**
         * Returns the strongly connected components of
         * this graph as a vector of vectors. Based of
         * Tarjan's algorithm.
         *
         * @return  a vector of vectors containing each
         *          strongly connected component
         */
        vector< vector<int> > stronglyConnectedComponents();

        /**
         * Returns a vector containing, for each vertex, the
         * number of the strongly connected component it
         * belongs to.
         *
         * @return  a vector giving the strongly connected
         *          component of each vertex
         */
        vector<int> vertexToStronglyConnectedComponentNumber();



        /***********************************
         * Extracting the grounding kernel *
         ***********************************/
        
        /**
         * Returns the grounding kernel of this graph, i.e.
         * the subgraph obtained by removing recursively
         * every sink.
         *
         * @return  a subgraph corresponding to the grounding
         *          kernel of this graph
         */
        Graph groundingKernel() const;

        /*************************
         * Contraction operators *
         *************************/

        /****************
         * Operator IN0 *
         ****************/

        /**
         * Removes all sources from this graph. This is
         * the IN0 operator mentioned in Lin and Jou.
         *
         * @return  a vector containing all removed vertices
         */
        vector<int> in0();

        /*****************
         * Operator OUT0 *
         *****************/

        /**
         * Removes all sinks from this graph. This is
         * the OUT0 operator mentioned in Lin and Jou.
         *
         * @return  a vector containing all removed vertices
         */
        vector<int> out0();

        /*****************
         * Operator LOOP *
         *****************/

        /**
         * Removes all self-loop vertices from this graph.
         * This is the LOOP operator mentioned in Lin and Jou.
         *
         * @return  a vector containing all removed vertices
         */
        vector<int> loop();

        /****************
         * Operator IN1 *
         ****************/

        /**
         * Merges all non self-looped vertices having in-degree
         * equal to one with their unique predecessor. This is
         * the IN1 operator mentioned in Lin and Jou.
         *
         * @return  a vector containing all merged vertices
         */
        vector<int> in1();

        /*****************
         * Operator OUT1 *
         *****************/

        /**
         * Merges all non self-looped vertices having out-degree
         * equal to one with their unique successor. This is
         * the OUT1 operator mentioned in Lin and Jou.
         *
         * @return  a vector containing all merged vertices
         */
        vector<int> out1();

        /****************
         * Operator PIE *
         ****************/

        /**
         * Returns all acyclic edges of this graph.
         *
         * An edge (u,v) is called acyclic if the vertices
         * u and v do not belong to the same strongly
         * connected component.
         *
         * @return  a vector containing the acyclic edges
         */
        vector<Edge> acyclicEdges();

        /**
         * Returns the pi-edges of this graph.
         *
         * An edge (u,v) is called a pi-edge if (v,u) is also
         * an edge of the graph.
         *
         * @return  a vector containing the pi-edges
         */
        vector<Edge> piEdges();

        /**
         * Returns all pseudo-acyclic edges of this graph.
         *
         * An edge (u,v) is called pseudo-acyclic if it is
         * an acyclic edge of this graph after having removed
         * the pi-edges.
         *
         * @return  a vector containing the pseudo
         *          acyclic edges
         */
        vector<Edge> pseudoAcyclicEdges();

        /**
         * Removes all acyclic and pseudo-acyclic edges (see
         * acyclicEdges and pseudoAcyclicEdges for more details.
         *
         * @return  a vector containing the removed
         *          acyclic and pseudo-acyclic edges
         */
        vector<Edge> pie();

        /*****************
         * Operator CORE *
         *****************/

        /**
         * Checks whether the given vertex is a pi-vertex
         *
         * A vertex is called pi-vertex if all its incident
         * edges are pi-edges.
         *
         * @param vertex  the vertex to be checked
         * @return        true if the given vertex is a
         *                pi-vertex
         */
        bool isPiVertex(int vertex);

        /**
         * Returns all pi-vertices.
         *
         * A vertex is called pi-vertex if all its incident
         * edges are pi-edges.
         *
         * @return  a vector containing the pi-vertices of this graph
         */
        vector<int> piVertices();

        /**
         * Checks whether the given vertices form a clique
         * of this graph
         *
         * A set {v1,v2,...vk} of vertices of a directed graph is
         * called a clique if vi is a neighbor of vj for any i,j
         * in {1,2,...,k} and vi does not have a self-loop for all
         * i in{1,2,...,k}.
         *
         * @param vs  the vertices to be checked
         * @return    true if the given vertices form a clique
         *            of this graph, false otherwise
         */
        bool isClique(const vector<int>& vs) const;

        /**
         * Removes a clique of this graph containing at least
         * one core vertex v and return all neighbors of v, if
         * such a clique exists.
         *
         * A vertex is called a core vertex if all its neighbors
         * together with itself form a clique of the graph.
         *
         * @return  the neighbors of the found core vertex 
         */
        vector<int> core();

        /*****************
         * Operator DOME *
         *****************/

        /**
         * Returns the pi-predecessors of the given vertex.
         *
         * @param vertex  the vertex whose pi-predecessors
         *                are queried
         * @return        the set of pi-predecessors of the
         *                given vertex;
         */
        vector<int> piPredecessors(int vertex) const;

        /**
         * Returns the pi-successors of the given vertex.
         *
         * @param vertex  the vertex whose pi-successors
         *                are queried
         * @return        the set of pi-successors of the
         *                given vertex;
         */
        vector<int> piSuccessors(int vertex) const;

        /**
         * Returns the non-pi-predecessors of the given vertex.
         *
         * @param vertex  the vertex whose non-pi-predecessors
         *                are queried
         * @return        the set of non-pi-predecessors of the
         *                given vertex;
         */
        vector<int> nonPiPredecessors(int vertex) const;

        /**
         * Returns the non-pi-successors of the given vertex.
         *
         * @param vertex  the vertex whose non-pi-successors
         *                are queried
         * @return        the set of non-pi-successors of the
         *                given vertex;
         */
        vector<int> nonPiSuccessors(int vertex) const;

        /**
         * Checks whether the given edge is dominated
         *
         * @param source  the source vertex of the edge
         * @param target  the target vertex of the edge
         * @return        true if the given edge is dominated
         */
        bool isDominated(int source, int target) const;

        /**
         * Returns the dominated edges of this graph
         *
         * @return  a vector containing the dominated edges
         *          of this graph
         */
        vector<Edge> dominatedEdges() const;

        /**
         * Removes all dominated edges of this graph
         *
         * @return  a vector containing all removed edges
         */
        vector<Edge> dome();

        /*************
         * Reduction *
         *************/

        /**
         * Reduces the current graph using the given
         * operator
         *
         * @param applyIn0   true if the IN0  operator must be applied
         * @param applyOut0  true if the OUT0 operator must be applied
         * @param applyLoop  true if the LOOP operator must be applied
         * @param applyIn1   true if the IN1  operator must be applied
         * @param applyOut1  true if the OUT1 operator must be applied
         * @param applyPie   true if the PIE  operator must be applied
         * @param applyCore  true if the CORE operator must be applied
         * @param applyDome  true if the DOME operator must be applied
         * @param verbose    true if the trace must be output
         * @return           a vector containing the partial solution
         */
        vector<int> reduce(bool applyIn0,  bool applyOut0, bool applyLoop,
                           bool applyIn1,  bool applyOut1, bool applyPie,
                           bool applyCore, bool applyDome, bool verbose);

        /**
         * Reduces the current graph until it is not reducible
         * anymore
         *
         * @param verbose  if true, then the trace is displayed
         *                 to the screen
         * @return         a vector containing the partial solution
         */
        vector<int> reduce(bool verbose=false);

        /***********************************
         * Cycles and feedback vertex sets *
         ***********************************/

        /**
         * Returns a shortest cycle of this graph.
         *
         * @return             a vector containing, in the right 
         *                     order, the vertices forming the
         *                     shortest cycle
         */
        vector<int> shortestCycle() const;

        /**
         * Checks whether the given graph is acyclic
         *
         * A graph is called acyclic if it contains no
         * cycle
         * 
         * @return  true if the graph is indeed acyclic
         */
        bool isAcyclic();

        /**
         * Tells whether the given vector contains a feedback
         * vertex set of this graph
         *
         * A set of vertices U is called feedback vertex set
         * of a graph G if the graph obtained by deleting U
         * from G is acyclic
         *
         * @param fvs  the vertices to be checked
         * @return     true if the given vertices form a
         *             feedback vertex set of this graph
         */
        bool isFeedbackVertexSet(const vector<int>& fvs) const;

        /**************************
         * Upper bound algorithms *
         **************************/

        /**
         * Computes an upper bound value for the minimum feedback
         * vertex problem set using the eight contraction operators.
         *
         * The integer u is an upper bound value for this problem
         * if every minimum feedback vertex set of the graph is
         * at most of size u.
         *
         * @param verbose  if true, then the trace is displayed
         *                 to the screen
         * @return         an upper bound for the minimum feedback
         *                 vertex set problem
         */
        int upperBoundValue(bool verbose=false) const;

        /**
         * Computes a feedback vertex set U of this graph using
         * the eight contraction operators. Therefore, any minimum
         * feedback vertex set of this graph has size at most
         * the size of U.
         *
         * @param verbose  if true, then the trace is displayed
         *                 to the screen
         * @return         a vector containing a feedback vertex
         *                 set of this graph
         */
        vector<int> upperBoundSolution(bool verbose=false) const;

        /**************************
         * Lower bound algorithms *
         **************************/

        /**
         * Computes a lower bound value for the minimum feedback
         * vertex set problem using the eight contraction operators.
         *
         * The integer l is a lower bound value for this problem
         * if every minimum feedback vertex set of the graph is
         * at least of size l.
         *
         * @param verbose  if true, then the trace is displayed
         *                 to the screen
         * @return         a lower bound for the minimum feedback
         *                 vertex set problem
         */
        int lowerBoundValue(bool verbose=false) const;

        /********************************
         * Minimum feedback vertex sets *
         ********************************/

        /**
         * Finds a minimum feedback vertex set of this graph
         * using a branch-and-bound algorithm based on the
         * eight contraction operators of Lin and Jou.
         *
         * @return  a vector containing a minimum feedback
         *          vertex set of this graph
         */
        vector<int> minimumFeedbackVertexSet(bool verbose=false) const;

    private:
	    vector< vector<int> > _outgoingNeighbors; // The list of outgoing neighbors
	    vector< vector<int> > _incomingNeighbors; // The list of incoming neighbors
        set<int> _vertices;                       // The set of vertices of the graph
        vector<bool> _hasVertex;                  // A map indicating which vertices belong to the graph
	    int _maxVertices;                         // The maximum number of allowed vertices
	    int _numVertices;                         // The current number of vertices
	    int _numEdges;                            // The current number of edges

        // Private functions for constructors
        void GraphFromAdjacencyList(const string& filename);
        void GraphFromArcsList(const string& filename);

        // Private function for constructors
        void _init(int maxSize);

        // Private function for cycle search
        vector<int> _shortestCycle(int vertex, int maxLength);
                                                   // Computes the shortest cycle starting from
                                                   // the given vertex, if its length does not
                                                   // exceed the given max length

        // Private functions and attributes for Tarjan's algorithm
        int _tarjanCurrentIndex;                   // The current index in the search
        vector<int> _tarjanIndex;                  // The search index of each vertex
        vector<int> _tarjanAncestor;               // The lowest ancestor of each vertex
        stack<int> _tarjanStack;                   // The stack used for Tarjan's algorithm
        vector<bool> _tarjanInStack;               // Tells if the vertex is in the stack
        int _currentComponent;                     // The current strongly connected component
        vector<int> _sccsByNum;                    // A vector associating each vertex with
                                                   // its corresponding strongly connected
                                                   // component
        vector< vector<int> > _sccs;               // The strongly connected components
        void _tarjan(int vertex, bool byNumber);   // The recursive function for computing the
                                                   // strongly connected components

        // Private functions for branch-and-bound algorithm finding
        // a minimum feedback vertex set
        vector<int> _mfvs(vector<int> solution, vector<int> bestSolution, int lowerBound,
                          int level, bool reducible, bool verbose);
                                                   // Computes a minimum feedback vertex
                                                   // set based on a branch-and-bound
                                                   // algorithm
};

}

#endif
