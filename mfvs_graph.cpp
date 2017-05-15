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

#include <stdio.h>
#include <string.h>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>
#include <queue>
#include <set>
#include <stack>

#include "mfvs_graph.h"
#include "mfvs_utils.h"

using namespace std;

namespace mfvs {

/****************************
 * Constants initialization *
 ****************************/
string Graph::TYPE_ARCS_LIST = string("arcs list");
string Graph::TYPE_ADJACENCY_LIST = string("adjacency list");

/**********************************
 * Comparison operators for edges *
 **********************************/
bool operator==(const Edge& e1, const Edge& e2) {
    return e1.source == e2.source && e1.target == e2.target;
}

bool operator<(const Edge& e1, const Edge& e2) {
    if (e1.source < e2.source) return true;
    else if (e1.source == e2.source && e1.target < e2.target) return true;
    else return false;
}

bool operator>(const Edge& e1, const Edge& e2) {
    return e2 < e1;
}

bool operator<=(const Edge& e1, const Edge& e2) {
    return e1 < e1 || e1 == e2;
}

bool operator>=(const Edge& e1, const Edge& e2) {
    return e1 > e1 || e1 == e2;
}



/*******************************
 * Constructors and destructor *
 *******************************/

Graph::Graph(int maxSize) {
    if (maxSize > 0) {
        _init(maxSize);
    }
}

Graph::Graph(const string& filename, const string& type) {
    if (type == Graph::TYPE_ARCS_LIST) {
        GraphFromArcsList(filename);  
    } else if (type == Graph::TYPE_ADJACENCY_LIST)
        GraphFromAdjacencyList(filename);  
}

void Graph::GraphFromAdjacencyList(const string& filename) {
    string line;
    ifstream inputFile;
    inputFile.open(filename.c_str());
    if (inputFile.is_open()) {
        getline(inputFile, line);
        istringstream iss(line);
        int maxSize;
        iss >> maxSize;
        _init(maxSize);
        for (int i = 0; i < maxSize; ++i) addVertex(i);
        int currentVertex = 0;
        while (!inputFile.eof()) {
            vector<string> tokens;
            getline(inputFile, line);
            tokenize(line, tokens);
            for (vector<string>::iterator it = tokens.begin(); it != tokens.end(); ++it) {
                istringstream iss2(*it);
                int neighbor;
                iss2 >> neighbor;
                addEdge(currentVertex, neighbor);
            }
            ++currentVertex;
        }
        inputFile.close();
    }
}

void Graph::GraphFromArcsList(const string& filename) {
    //TODO
}

Graph::~Graph() {
}

void Graph::_init(int maxSize) {
    _outgoingNeighbors.resize(maxSize);
    _incomingNeighbors.resize(maxSize);
    _hasVertex.resize(maxSize);
    for (int i = 0; i < maxSize; ++i) _hasVertex[i] = false;
    _maxVertices = maxSize;
    _numVertices = 0;
    _numEdges    = 0;
}


/*******************
 * Basic accessors *
 *******************/
int Graph::numVertices() const {
    return _numVertices;
}

int Graph::numEdges() const {
    return _numEdges;
}

bool Graph::hasVertex(int vertex) const {
    if (vertex < 0 || vertex >= _maxVertices) return false;
    else return _hasVertex[vertex];
}

bool Graph::hasEdge(int source, int target) const {
    if (!hasVertex(source) || !hasVertex(target)) {
        return false;
    } else {
        return find(outgoingNeighbors(source).begin(), outgoingNeighbors(source).end(), target)
               != outgoingNeighbors(source).end();
    }
}

set<int> const& Graph::vertices() const {
    return _vertices;
}

vector<Edge> Graph::edges() const {
    vector<Edge> allEdges;
    for (int i = 0; i < _maxVertices; ++i) {
        for (vector<int>::const_iterator jt = _outgoingNeighbors[i].begin();
             jt != _outgoingNeighbors[i].end(); ++jt) {
            Edge e;
            e.source = i;
            e.target = *jt;
            allEdges.push_back(e);
        }
    }
    return allEdges;
}

vector<int> const& Graph::incomingNeighbors(int vertex) const {
    if (hasVertex(vertex)) return _incomingNeighbors[vertex];
    else throw "The vertex does not exist.";
}

vector<int> const& Graph::outgoingNeighbors(int vertex) const {
    if (hasVertex(vertex)) return _outgoingNeighbors[vertex];
    else throw "The vertex does not exist.";
}

void Graph::print(bool edges) const {
    cout << "Graph of " << numVertices() << " vertices and " << numEdges() << " edges";
    if (edges) {
        cout << ": ";
        for (set<int>::iterator it = _vertices.begin(); it != _vertices.end(); ++it) {
            const vector<int> neighbors = outgoingNeighbors(*it);
            for (int j = 0; j < int(neighbors.size()); ++j) {
                cout << "(" << *it << "," << neighbors[j] << ") ";
            }
        }
    }
}

void Graph::print() const {
    print(false);
}



/*******************
 * Basic modifiers *
 *******************/
void Graph::addVertex(int vertex) {
    if (vertex >= 0 && vertex < _maxVertices && !hasVertex(vertex)) {
        _vertices.insert(vertex);
        _hasVertex[vertex] = true;
        ++_numVertices;
    }
}

void Graph::addEdge(int source, int target) {
    if (hasVertex(source) && hasVertex(target) && !hasEdge(source, target)) {
        _outgoingNeighbors[source].push_back(target);
        _incomingNeighbors[target].push_back(source);
        ++_numEdges;
    }
}

void Graph::deleteVertex(int vertex) {
    if (hasVertex(vertex)) {
        vector<int> outN = outgoingNeighbors(vertex);
        vector<int> inN  = incomingNeighbors(vertex);
        for (int i = 0; i < int(outN.size()); ++i) deleteEdge(vertex, outN[i]);
        for (int i = 0; i < int(inN.size());  ++i) deleteEdge(inN[i], vertex);
        _vertices.erase(vertex);
        _hasVertex[vertex] = false;
        --_numVertices;
    }
}

void Graph::deleteVertices(const vector<int>& vertices) {
    for (vector<int>::const_iterator it = vertices.begin(); it != vertices.end(); ++it)
        deleteVertex(*it);
}

void Graph::deleteEdge(int source, int target) {
    if (hasEdge(source, target)) {
        _outgoingNeighbors[source].erase(remove(_outgoingNeighbors[source].begin(),
                                         _outgoingNeighbors[source].end(), target));
        _incomingNeighbors[target].erase(remove(_incomingNeighbors[target].begin(),
                                         _incomingNeighbors[target].end(), source));
        --_numEdges;
    }
}

void Graph::deleteEdges(const vector<Edge>& edges) {
    for (vector<Edge>::const_iterator it = edges.begin(); it != edges.end(); ++it) {
        deleteEdge(it->source, it->target);
    }
}

void Graph::mergeVertex(int vertex) {
    if (hasVertex(vertex)) {
        vector<int> inN = incomingNeighbors(vertex);
        vector<int> outN = outgoingNeighbors(vertex);
        for (vector<int>::iterator it = inN.begin(); it != inN.end(); ++it)
            for (vector<int>::iterator jt = outN.begin(); jt != outN.end(); ++jt)
                addEdge(*it, *jt);
        deleteVertex(vertex);
    } else throw "The vertex does not exist.";
}

void Graph::mergeVertices(vector<int> vs) {
    if (vs.size() > 0) {
        int v = vs[0];
        for (vector<int>::iterator it = vs.begin(); it != vs.end(); ++it) {
            if (*it != v) {
                vector<int> inN = incomingNeighbors(*it);
                vector<int> outN = outgoingNeighbors(*it);
                vector<int>::iterator kt;
                for (vector<int>::iterator jt = inN.begin(); jt != inN.end(); ++jt) {
                    kt = find(vs.begin(), vs.end(), *jt);
                    if (*kt == v || kt == vs.end()) addEdge(*jt, v);
                }
                for (vector<int>::iterator jt = outN.begin(); jt != outN.end(); ++jt) {
                    kt = find(vs.begin(), vs.end(), *jt);
                    if (*kt == v || kt == vs.end()) addEdge(v, *jt);
                }
                deleteVertex(*it);
            }
        }
    }
}

void Graph::mergeVertices(int vertex1, int vertex2) {
    vector<int> vs = vector<int>();
    vs.push_back(vertex1);
    vs.push_back(vertex2);
    mergeVertices(vs);
}



/*************
 * Subgraphs *
 *************/

Graph Graph::subgraph(vector<int> vs) const {
    Graph h(_maxVertices);
    for (vector<int>::iterator it = vs.begin(); it != vs.end(); ++it)
        h.addVertex(*it);
    for (vector<int>::iterator it = vs.begin(); it != vs.end(); ++it) {
        if (hasVertex(*it)) {
            vector<int> outN = outgoingNeighbors(*it);
            for (vector<int>::iterator jt = outN.begin(); jt != outN.end(); ++jt)
                if (h.hasVertex(*jt)) h.addEdge(*it, *jt);
        }
    }
    return h;
}



/******************************
 * Degrees, sources and sinks *
 ******************************/

int Graph::inDegree(int vertex) const {
    if (hasVertex(vertex)) return int(incomingNeighbors(vertex).size());
    else return -1;
}

int Graph::outDegree(int vertex) const {
    if (hasVertex(vertex)) return int(outgoingNeighbors(vertex).size());
    else return -1;
}

int Graph::degree(int vertex) const {
    if (hasVertex(vertex)) return inDegree(vertex) + outDegree(vertex);
    else return -1;
}

int Graph::minDegreeVertex() const {
    int v = -1;
    int minDegree = 2 * numVertices() + 1;
    for (int i = 0; i < _maxVertices; ++i) {
        if (hasVertex(i) && degree(i) < minDegree) {
            v = i;
            minDegree = degree(i);
        }
    }
    return v;
}

int Graph::maxDegreeVertex() const {
    int v = -1;
    int maxDegree = -1;
    for (int i = 0; i < _maxVertices; ++i) {
        if (hasVertex(i) && degree(i) > maxDegree) {
            v = i;
            maxDegree = degree(i);
        }
    }
    return v;
}

bool Graph::isSource(int vertex) const {
    if (hasVertex(vertex)) return inDegree(vertex) == 0;
    else return false;
}

bool Graph::isSink(int vertex)  const{
    if (hasVertex(vertex)) return outDegree(vertex) == 0;
    else return false;
}

bool Graph::hasLoop(int vertex) const {
    if (!hasVertex(vertex)) return false;
    else {
        vector<int> outN = outgoingNeighbors(vertex);
        return find(outN.begin(), outN.end(), vertex) != outN.end();
    }
}

vector<int> Graph::sources() const {
    vector<int> vertices = vector<int>();
    for (set<int>::iterator it = _vertices.begin(); it != _vertices.end(); ++it) {
        if (isSource(*it)) vertices.push_back(*it);
    }
    return vertices;
}

vector<int> Graph::sinks() const {
    vector<int> vertices = vector<int>();
    for (set<int>::iterator it = _vertices.begin(); it != _vertices.end(); ++it) {
        if (isSink(*it)) vertices.push_back(*it);
    }
    return vertices;
}

vector<int> Graph::loops() const {
    vector<int> vertices = vector<int>();
    for (set<int>::iterator it = _vertices.begin(); it != _vertices.end(); ++it) {
        if (hasLoop(*it)) vertices.push_back(*it);
    }
    return vertices;
}


/*********************************
 * Strongly connected components *
 *********************************/
vector< vector<int> > Graph::stronglyConnectedComponents() {
    _tarjanCurrentIndex = 0;
    _tarjanIndex = vector<int>(_maxVertices);
    _tarjanAncestor = vector<int>(_maxVertices);
    _tarjanStack = stack<int>();
    _tarjanInStack = vector<bool>(_maxVertices);
    _sccs = vector< vector<int> >();
    for (int i = 0; i < _maxVertices; ++i) {
        _tarjanIndex[i] = -1;
        _tarjanInStack[i] = false;
    }
    for (set<int>::iterator it = vertices().begin(); it != vertices().end(); ++it) {
        if (_tarjanIndex[*it] == -1) _tarjan(*it, false);
    }
    return _sccs;
}

vector<int> Graph::vertexToStronglyConnectedComponentNumber() {
    _tarjanCurrentIndex = 0;
    _tarjanIndex = vector<int>(_maxVertices);
    _tarjanAncestor = vector<int>(_maxVertices);
    _tarjanStack = stack<int>();
    _tarjanInStack = vector<bool>(_maxVertices);
    _currentComponent = 0;
    _sccsByNum = vector<int>(_maxVertices);
    for (int i = 0; i < _maxVertices; ++i) {
        _tarjanIndex[i] = -1;
        _tarjanInStack[i] = false;
    }
    for (set<int>::iterator it = vertices().begin(); it != vertices().end(); ++it) {
        if (_tarjanIndex[*it] == -1) _tarjan(*it, true);
    }
    return _sccsByNum;
}

void Graph::_tarjan(int vertex, bool byNumber) {
    _tarjanIndex[vertex] = _tarjanCurrentIndex;
    _tarjanAncestor[vertex] = _tarjanCurrentIndex;
    ++_tarjanCurrentIndex;
    _tarjanStack.push(vertex);
    _tarjanInStack[vertex] = true;
    vector<int> outN = outgoingNeighbors(vertex);
    for (vector<int>::iterator it = outN.begin(); it != outN.end(); ++it) {
        if (_tarjanIndex[*it] == -1) {
            _tarjan(*it, byNumber);
            _tarjanAncestor[vertex] = min(_tarjanAncestor[vertex], _tarjanAncestor[*it]);
        } else if (_tarjanInStack[*it]) {
            _tarjanAncestor[vertex] = min(_tarjanAncestor[vertex], _tarjanIndex[*it]);
        }
    }
    if (_tarjanAncestor[vertex] == _tarjanIndex[vertex]) {
        if (byNumber) {
            int u;
            do {
                u = _tarjanStack.top();
                _sccsByNum[u] = _currentComponent;
                _tarjanStack.pop();
                _tarjanInStack[u] = false;
            } while (u != vertex);
            ++_currentComponent;
        } else {
            vector<int> scc = vector<int>();
            int u;
            do {
                u = _tarjanStack.top();
                scc.push_back(u);
                _tarjanStack.pop();
                _tarjanInStack[u] = false;
            } while (u != vertex);
            _sccs.push_back(scc);
        }
    }
}



/***********************************
 * Extracting the grounding kernel *
 ***********************************/

Graph Graph::groundingKernel() const {
    Graph h(*this);
    vector<int> out0Vertices;
    do {
        out0Vertices = h.out0();
        h.deleteVertices(out0Vertices); 
    } while (int(out0Vertices.size()) != 0);
    return h;
}



/*************************
 * Contraction operators *
 *************************/

/****************
 * Operator IN0 *
 ****************/
vector<int> Graph::in0() {
    vector<int> in0Vertices = sources();
    deleteVertices(in0Vertices);
    return in0Vertices;
}

/*****************
 * Operator OUT0 *
 *****************/
vector<int> Graph::out0() {
    vector<int> out0Vertices = sinks();
    deleteVertices(out0Vertices);
    return out0Vertices;
}

/*****************
 * Operator LOOP *
 *****************/
vector<int> Graph::loop() {
    vector<int> loopVertices = loops();
    deleteVertices(loopVertices);
    return loopVertices;
}

/****************
 * Operator IN1 *
 ****************/

vector<int> Graph::in1() {
    if (int(loops().size()) == 0) {
        vector<int> in1Vertices = vector<int>();
        vector<int> vs(int(vertices().size()));
        copy(vertices().begin(), vertices().end(), vs.begin());
        for (vector<int>::iterator it = vs.begin(); it != vs.end(); ++it) {
            if (hasVertex(*it) && inDegree(*it) == 1 && !hasLoop(*it)) {
                in1Vertices.push_back(*it);
                mergeVertex(*it);
            }
        }
        return in1Vertices;
    } else
        return vector<int>();
}

/*****************
 * Operator OUT1 *
 *****************/

vector<int> Graph::out1() {
    if (int(loops().size()) == 0) {
        vector<int> out1Vertices = vector<int>();
        vector<int> vs(int(vertices().size()));
        copy(vertices().begin(), vertices().end(), vs.begin());
        for (vector<int>::iterator it = vs.begin(); it != vs.end(); ++it) {
            if (hasVertex(*it) && outDegree(*it) == 1 && !hasLoop(*it)) {
                out1Vertices.push_back(*it);
                mergeVertex(*it);
            }
        }
        return out1Vertices;
    } else
        return vector<int>();
}

/****************
 * Operator PIE *
 ****************/

vector<Edge> Graph::acyclicEdges() {
    vector<Edge> ae;
    vector<Edge> allEdges = edges();
    vector<int> vertexToSCC = vertexToStronglyConnectedComponentNumber();
    for (vector<Edge>::iterator it = allEdges.begin(); it != allEdges.end(); ++it) {
        if (vertexToSCC[it->source] != vertexToSCC[it->target]) {
            ae.push_back(*it);
        }
    }
    return ae;
}

vector<Edge> Graph::piEdges() {
    vector<Edge> es;
    vector<Edge> allEdges = edges();
    for (vector<Edge>::iterator it = allEdges.begin(); it != allEdges.end(); ++it) {
        if (hasEdge(it->target, it->source)) {
            es.push_back(*it);
        }
    }
    return es;
}

vector<Edge> Graph::pseudoAcyclicEdges() {
    Graph h(*this);
    vector<Edge> pies = piEdges();
    h.deleteEdges(acyclicEdges());
    h.deleteEdges(pies);    
    return h.acyclicEdges();
}

vector<Edge> Graph::pie() {
    if (int(loops().size()) == 0) {
        vector<Edge> aes = acyclicEdges();
        vector<Edge> paes = pseudoAcyclicEdges();
        deleteEdges(aes);
        deleteEdges(paes);
        vector<Edge> es(aes);
        for (vector<Edge>::iterator it = paes.begin(); it != paes.end(); ++it) {
            es.push_back(*it);
        }
        return es;
    } else
        return vector<Edge>();
}

/*****************
 * Operator CORE *
 *****************/

bool Graph::isPiVertex(int vertex) {
    vector<int> outN = outgoingNeighbors(vertex);
    bool valid = true;
    vector<int>::iterator it = outN.begin();
    while (valid && it != outN.end()) {
        valid = hasEdge(*it, vertex);
        ++it;
    }
    return valid;
}

vector<int> Graph::piVertices() {
    vector<int> vs;
    for (set<int>::iterator it = vertices().begin(); it != vertices().end(); ++it) {
        if (isPiVertex(*it)) vs.push_back(*it);
    }
    return vs;
}

bool Graph::isClique(const vector<int>& vs) const {
    bool clique = true;
    vector<int>::const_iterator it = vs.begin();
    do {
        vector<int>::const_iterator jt = vs.begin();
        do {
            clique = (*it == *jt && !hasEdge(*it, *jt)) || (*it != *jt && hasEdge(*it, *jt));
            ++jt;
        } while (clique && jt != vs.end());
        ++it;
    } while (clique && it != vs.end());
    return clique;
}

vector<int> Graph::core() {
    if (int(loops().size()) == 0) {
        vector<int> piVs = piVertices();
        for (vector<int>::iterator it = piVs.begin(); it != piVs.end(); ++it) {
            vector<int> potentialClique = outgoingNeighbors(*it);
            potentialClique.push_back(*it);
            if (isClique(potentialClique)) {
                deleteVertices(potentialClique);
                potentialClique.pop_back();
                return potentialClique;
            }
        }
    }
    return vector<int>();
}

/*****************
 * Operator DOME *
 *****************/

vector<int> Graph::piPredecessors(int vertex) const {
    vector<int> piPreds;
    vector<int> preds = incomingNeighbors(vertex);
    for (vector<int>::iterator it = preds.begin(); it != preds.end(); ++it)
        if (hasEdge(vertex, *it)) piPreds.push_back(*it);
    return piPreds;
}

vector<int> Graph::piSuccessors(int vertex) const {
    vector<int> piSucc;
    vector<int> succ = outgoingNeighbors(vertex);
    for (vector<int>::iterator it = succ.begin(); it != succ.end(); ++it)
        if (hasEdge(*it, vertex)) piSucc.push_back(*it);
    return piSucc;
}

vector<int> Graph::nonPiPredecessors(int vertex) const {
    vector<int> nonPiPreds;
    vector<int> preds = incomingNeighbors(vertex);
    for (vector<int>::iterator it = preds.begin(); it != preds.end(); ++it)
        if (!hasEdge(vertex, *it)) nonPiPreds.push_back(*it);
    return nonPiPreds;
}

vector<int> Graph::nonPiSuccessors(int vertex) const {
    vector<int> nonPiSucc;
    vector<int> succ = outgoingNeighbors(vertex);
    for (vector<int>::iterator it = succ.begin(); it != succ.end(); ++it)
        if (!hasEdge(*it, vertex)) nonPiSucc.push_back(*it);
    return nonPiSucc;
}

bool Graph::isDominated(int source, int target) const {
    if (!hasEdge(source, target)) return false;
    else if (hasEdge(source, target) && hasEdge(target, source)) return false;
    else {
        vector<int> pTarget = incomingNeighbors(target);
        vector<int> pSource = nonPiPredecessors(source);
        sort(pTarget.begin(), pTarget.end());
        sort(pSource.begin(), pSource.end());
        if (includes(pTarget.begin(), pTarget.end(), pSource.begin(), pSource.end()))
            return true;
        vector<int> sSource = outgoingNeighbors(source);
        vector<int> sTarget = nonPiSuccessors(target);
        sort(sSource.begin(), sSource.end());
        sort(sTarget.begin(), sTarget.end());
        if (includes(sSource.begin(), sSource.end(), sTarget.begin(), sTarget.end()))
            return true;
        else
            return false;
    }
}

vector<Edge> Graph::dominatedEdges() const {
    vector<Edge> allEdges = edges();
    vector<Edge> des;
    for (vector<Edge>::iterator it = allEdges.begin(); it != allEdges.end(); ++it)
        if (isDominated(it->source, it->target)) des.push_back(*it);
    return des;
}

vector<Edge> Graph::dome() {
    if (int(loops().size()) == 0) {
        vector<Edge> des = dominatedEdges();
        deleteEdges(des);
        return des;
    } else
        return vector<Edge>();
}

/*************
 * Reduction *
 *************/
vector<int> Graph::reduce(bool applyIn0,  bool applyOut0, bool applyLoop,
                          bool applyIn1,  bool applyOut1, bool applyPie,
                          bool applyCore, bool applyDome, bool verbose) {
    vector<int> solution;
    if (applyIn0) {
        vector<int> in0Vertices = in0();
        if (verbose) cout << "IN0  : " << int(in0Vertices.size())
                     << " vertex(vertices) has(have) been removed." << endl;
    }
    if (applyOut0) {
        vector<int> out0Vertices = out0();
        if (verbose) cout << "OUT0 : " << int(out0Vertices.size())
                     << " vertex(vertices) has(have) been removed." << endl;
    }
    if (applyLoop) {
        vector<int> loopVertices = loop();
        for (vector<int>::iterator it = loopVertices.begin(); it != loopVertices.end(); ++it)
            solution.push_back(*it);
        if (verbose) cout << "LOOP : " << int(loopVertices.size())
                     << " vertex(vertices) has(have) been removed." << endl;
    }
    if (applyIn1) {
        vector<int> in1Vertices = in1();
        if (verbose) cout << "IN1  : " << int(in1Vertices.size())
                     << " vertex(vertices) has(have) been merged." << endl;
    }
    if (applyOut1) {
        vector<int> out1Vertices = out1();
        if (verbose) cout << "OUT1 : " << int(out1Vertices.size())
                     << " vertex(vertices) has(have) been merged." << endl;
    }
    if (applyPie) {
        vector<Edge> pieEdges = pie();
        if (verbose) cout << "PIE  : " << int(pieEdges.size())
                     << " edge(s) has(have) been removed." << endl;
    }
    if (applyCore) {
        vector<int> coreVertices = core();
        for (vector<int>::iterator it = coreVertices.begin(); it != coreVertices.end(); ++it)
            solution.push_back(*it);
        if (verbose) cout << "CORE : " << int(coreVertices.size())
                     << " vertex(vertices) has(have) been removed." << endl;
    }
    if (applyDome) {
        vector<Edge> domeEdges = dome();
        if (verbose) cout << "DOME : " << int(domeEdges.size())
                     << " edge(s) has(have) been removed." << endl;
    }
    return solution;
}

vector<int> Graph::reduce(bool verbose) {
    int n, m;
    vector<int> solution;
    do {
        n = numVertices();
        m = numEdges();
        vector<int> partialSolution = reduce(true, true, true, true, true,
                                             true, true, true, verbose);
        for (vector<int>::iterator it = partialSolution.begin();
             it != partialSolution.end(); ++it)
            solution.push_back(*it);
    } while (n != numVertices() || m != numEdges());
    return solution;
}

/***********************************
 * Cycles and feedback vertex sets *
 ***********************************/

vector<int> Graph::shortestCycle() const {
    Graph h(*this);
    if (h.isAcyclic()) {
        throw "No shortest cycle ! This graph is acyclic";
    } else {
        int n;
        do {
            n = h.numVertices();
            h.in0();
            h.out0();
        } while (n != h.numVertices());
        vector<int> shortestCycle;
        vector<int> currentCycle;
        int maxLength = h.numVertices();
        for (set<int>::iterator it = vertices().begin(); it != vertices().end(); ++it) {
            currentCycle = h._shortestCycle(*it, maxLength);
            if (int(currentCycle.size() - 1)  < maxLength) {
                shortestCycle = currentCycle;
                maxLength = int(currentCycle.size()) - 1;
            }
        }
        return shortestCycle;
    }
}

/**
 * Returns the shortest cycle starting from the given vertex
 * if its length does not exceed the given maximum length.
 *
 * @param vertex     the vertex from which the cycle must start
 * @param maxLength  the maximum allowed length for the cycle
 * @return           a vector containing in the right order the
 *                   vertices that form the shortest cycle
 */
vector<int> Graph::_shortestCycle(int vertex, int maxLength) {
    queue< list<int> > paths;
    list<int> path;
    path.push_back(vertex);
    paths.push(path);
    bool cycleFound = false;
    while (!cycleFound && !paths.empty()) {
        path = paths.front();
        paths.pop();
        cycleFound = int(path.size()) > 1 && path.front() == path.back();
        if (!cycleFound && int(path.size()) - 1 < maxLength) {
            vector<int> outN = outgoingNeighbors(path.back());
            for (vector<int>::const_iterator it = outN.begin();
                 it != outN.end(); ++it) {
                path.push_back(*it);
                paths.push(path);
                path.pop_back();
            }
        }
    }
    vector<int> cycle(int(path.size()));
    copy(path.begin(), path.end(), cycle.begin());
    return cycle;
}

bool Graph::isAcyclic() {
    return int(loops().size()) == 0 &&
           int(stronglyConnectedComponents().size()) == numVertices();
}

bool Graph::isFeedbackVertexSet(const vector<int>& fvs) const {
    Graph h(*this);
    h.deleteVertices(fvs);
    return h.isAcyclic();
}

/**************************
 * Upper bound algorithms *
 **************************/

int Graph::upperBoundValue(bool verbose) const {
    return int(upperBoundSolution(verbose).size());
}

vector<int> Graph::upperBoundSolution(bool verbose) const {
    vector<int> solution;
    Graph h(*this);
    while (h.numVertices() > 0) {
        if (verbose) {
            cout << endl;
            h.print(false);
            cout << endl;
        }
        vector<int> partialSolution = h.reduce(verbose);
        for (vector<int>::iterator it = partialSolution.begin();
             it != partialSolution.end(); ++it) {
            solution.push_back(*it);
        }
        if (h.numVertices() > 0) {
            int v = h.minDegreeVertex();
            if (v == -1) {
                h.print(true);
                cout << endl;
            }
            h.mergeVertex(v);
        }
    }
    return solution;
}

/**************************
 * Lower bound algorithms *
 **************************/
int Graph::lowerBoundValue(bool verbose) const {
    int lb = 0;
    Graph h(*this);
    while (h.numVertices() > 0) {
        if (verbose) {
            cout << endl;
            h.print(false);
            cout << endl;
        }
        lb += int(h.reduce(verbose).size());
        if (h.numVertices() > 0) {
            vector<int> shortestCycle = h.shortestCycle();
            h.deleteVertices(shortestCycle);
            lb += 1;
        }
    }
    return lb;
}

/********************************
 * Minimum feedback vertex sets *
 ********************************/
vector<int> Graph::minimumFeedbackVertexSet(bool verbose) const {
    Graph h(*this);
    vector<int> vs(h.numVertices());
    copy(h.vertices().begin(), h.vertices().end(), vs.begin());
    return h._mfvs(vector<int>(), vs, 0, 0, true, verbose);
}

vector<int> Graph::_mfvs(vector<int> solution, vector<int> bestSolution, int lowerBound,
                         int level, bool reducible, bool verbose) {
    Graph h(*this);
    vector<int> partialSolution;
    if (reducible) {
        // Reducing the graph
        partialSolution = h.reduce(false);
        for (vector<int>::iterator it = partialSolution.begin(); it != partialSolution.end(); ++it)
            solution.push_back(*it);
        // Partitioning according to strongly connected components
        vector< vector<int> > sccs = h.stronglyConnectedComponents();
        if (int(sccs.size()) > 1) {
            for (vector< vector<int> >::iterator it = sccs.begin(); it != sccs.end(); ++it) {
                Graph scc = h.subgraph(*it);            
                vector<int> sccSolution = scc._mfvs(vector<int>(), *it, 0, 0, false, verbose);
                solution.insert(solution.end(), sccSolution.begin(), sccSolution.end());
            }
            return solution;
        }
    }
    if (verbose) {
        cout << "Level : " << level << endl;
        h.print(false);
        cout << endl;
        cout << "u = " << h.upperBoundValue() << endl;
        cout << "l = " << h.lowerBoundValue() << endl;
    }
    // Lower bound and initial solution
    int lb = h.lowerBoundValue(false);
    if (level == 0) {
        lowerBound = lb;
        bestSolution = h.upperBoundSolution(false);
        for (vector<int>::iterator it = solution.begin(); it != solution.end(); ++it)
            bestSolution.push_back(*it);
    }
    // Bounding the search
    if (h.numVertices() == 0) return solution;
    else if (int(solution.size()) + lb > int(bestSolution.size())) return bestSolution;
    // Branching
    int v = h.maxDegreeVertex();
    Graph leftGraph(h);
    leftGraph.deleteVertex(v);
    solution.push_back(v);
    vector<int> leftSolution = leftGraph._mfvs(solution, bestSolution, lowerBound,
                                               level + 1, true, verbose);
    if (int(leftSolution.size()) < int(bestSolution.size())) bestSolution = leftSolution;
    // Rebounding the search
    if (lb == int(bestSolution.size())) return bestSolution;
    Graph rightGraph(h);
    rightGraph.mergeVertex(v);
    solution.pop_back();
    vector<int> rightSolution = rightGraph._mfvs(solution, bestSolution, lowerBound,
                                                 level + 1, true, verbose);
    if (int(rightSolution.size()) < int(bestSolution.size())) bestSolution = rightSolution;
    return bestSolution;
}

}