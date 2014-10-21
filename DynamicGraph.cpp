/* -------------------------------------------------
      _       _     ___                            
 __ _| |_ _ _(_)___/ __| __ __ _ _ _  _ _  ___ _ _ 
/ _` |  _| '_| / -_)__ \/ _/ _` | ' \| ' \/ -_) '_|
\__, |\__|_| |_\___|___/\__\__,_|_||_|_||_\___|_|  
|___/                                          
    
gtrieScanner: quick discovery of network motifs
Released under Artistic License 2.0
(see README and LICENSE)

Pedro Ribeiro - CRACS & INESC-TEC, DCC/FCUP

----------------------------------------------------
Graphs Implementation with Adj. Matrix and Adj. List

Last Update: 11/02/2012
---------------------------------------------------- */

#include "DynamicGraph.h"
#include "GraphUtils.h"
#include <stdio.h>
#include <algorithm>

DynamicGraph::DynamicGraph(RepType _r) {
  _rtype = _r;
  _init();
}

DynamicGraph::~DynamicGraph() {
  _delete();
}

// ------------------------------
// Graph Creation
// ------------------------------

void DynamicGraph::_init() {
  _num_nodes = _num_edges = 0;

  _adjM             = NULL;
  _adjOut           = NULL;
  _adjIn            = NULL;
  _neighbours       = NULL;  
  _in               = NULL;
  _out              = NULL;
  _num_neighbours   = NULL;
  _array_neighbours = NULL;
}

void DynamicGraph::_delete() {
  int i;

  if (_adjM!=NULL) {
    for (i=0; i<_num_nodes; i++)
      if (_adjM[i]!=NULL) delete[] _adjM[i];
    delete[] _adjM;
  }
  if (_adjIn!=NULL) delete[] _adjIn;
  if (_adjOut!=NULL) delete[] _adjOut;
  if (_neighbours!=NULL) delete[] _neighbours;

  if (_in!=NULL) delete[] _in;
  if (_out!=NULL) delete[] _out;
  if (_out!=NULL) delete[] _num_neighbours;

  if (_array_neighbours!=NULL) {
    for (i=0; i<_num_nodes; i++)
      if (_array_neighbours[i]!=NULL) delete[] _array_neighbours[i];
    delete[] _array_neighbours;
  }
}

void DynamicGraph::zero() {
  int i,j;
  _num_edges = 0;

  for (i=0; i<_num_nodes;i++) {
    _in[i] = 0;
    _out[i] = 0;
    _num_neighbours[i] = 0;
    _adjIn[i].clear();
    _adjOut[i].clear();
    _neighbours[i].clear();
    if (_rtype == MATRIX)
      for (j=0; j<_num_nodes;j++)
        _adjM[i][j]=false;
  }
}

void DynamicGraph::createGraph(int n, GraphType t) {
  int i;

  _num_nodes = n;
  _type = t;

  _delete();

  if (_rtype == MATRIX)
  {
    _adjM = new bool*[n];  
    for (i=0; i<n; i++) _adjM[i] = new bool[n];
  }
  _adjIn      = new vector<int>[n];
  _adjOut     = new vector<int>[n];
  _neighbours = new vector<int>[n];

  _in             = new int[n]; 
  _out            = new int[n];
  _num_neighbours = new int[n];

  zero();
}

void DynamicGraph::prepareGraph() {
  if (_rtype == MATRIX)
    return;
}

void DynamicGraph::addEdge(int a, int b) {
  _adjOut[a].push_back(b);
  _out[a]++;

  _adjIn[b].push_back(a);
  _in[b]++;

  _num_edges++;

  if (!hasEdge(b, a)) {
    _neighbours[a].push_back(b);
    _num_neighbours[a]++;
    _neighbours[b].push_back(a);
    _num_neighbours[b]++;
  }

  if (_rtype == MATRIX)
  {
    if (_adjM[a][b]) return;
    _adjM[a][b] = true;
  }
}

void DynamicGraph::rmEdge(int a, int b) {
  fprintf(stderr, "rmEdge is deprecated for the time being\n");
}

bool DynamicGraph::hasEdge(int a, int b) {
  if (_rtype == MATRIX)
    return _adjM[a][b];
  return false;
}

void DynamicGraph::sortNeighbours() {
  int i;
  for (i=0; i<_num_nodes; i++)
    sort(_neighbours[i].begin(), _neighbours[i].begin()+_neighbours[i].size());
}

void DynamicGraph::sortNeighboursArray() {
  int i;
  for (i=0; i<_num_nodes; i++)
    qsort(_array_neighbours[i], _num_neighbours[i], sizeof(int), GraphUtils::int_compare);
}

void DynamicGraph::makeArrayNeighbours() {
  int i,j;
  vector<int>:: iterator ii;
  _array_neighbours = new int*[_num_nodes];  
  for (i=0; i<_num_nodes; i++) {
    _array_neighbours[i] = new int[_neighbours[i].size()];
    for (ii=_neighbours[i].begin(), j=0; ii!=_neighbours[i].end(); ++ii, j++)
      _array_neighbours[i][j] = *ii;
    _neighbours[i].clear();
  }
}

void DynamicGraph::makeVectorNeighbours() {
  int i,j;
  vector<int>:: iterator ii;

  for (i=0; i<_num_nodes; i++)
    for (j=0; j<_num_neighbours[i]; j++)
      _neighbours[i].push_back(_array_neighbours[i][j]);

  if (_array_neighbours!=NULL) {
    for (i=0; i<_num_nodes; i++)
      if (_array_neighbours[i]!=NULL) delete[] _array_neighbours[i];
    delete[] _array_neighbours;
  }
}
