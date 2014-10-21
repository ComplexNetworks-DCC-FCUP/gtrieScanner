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
Partially Abstract Base Graph Class

Last Update: 11/02/2012
---------------------------------------------------- */

#include "Graph.h"

// -----------------------------------------------

Graph::Graph() {
  _init();
}

Graph::Graph(int n) {
  _init();
  _createGraph(n);
}

Graph::~Graph() {
  _delete();
}


// Must come back and make it faster :)
Graph::Graph(bool a, Graph &g) {
  int i,j;

  _init();
  _createGraph(g.numNodes());

  for (i=0; i<g.numNodes(); i++)
    for (j=0; j<g.numNodes(); j++)
    addEdge(i, j, g.getEdge(i,j));
}


// -----------------------------------------------
// Graph creation
// -----------------------------------------------

void Graph::_init() {
  _adjL = NULL;
  _adjLI = NULL;
  _neighbours = NULL;
  _in_ds = NULL;
  _out_ds = NULL;
  _io_ds = NULL;
  _in = NULL;
  _out = NULL;
  _num_nodes = _num_edges = 0;  
}

void Graph::_createGraph(int num_nodes) {

  _num_nodes = num_nodes;
  _num_edges = 0;

  _delete();

  _adjL = new list<iPair>[num_nodes];
  _adjLI = new list<iPair>[num_nodes];
  _neighbours = new list<int>[num_nodes];
  _in_ds = new int[num_nodes];
  _out_ds = new int[num_nodes];
  _io_ds = new iPair[num_nodes];

  _in = new int[num_nodes];
  _out = new int[num_nodes];

  zero();
}

void Graph::_delete() {
  if (_adjL!=NULL) delete[] _adjL;
  if (_adjLI!=NULL) delete[] _adjLI;
  if (_neighbours=NULL) delete[] _neighbours;
  if (_in_ds!=NULL) delete[] _in_ds;
  if (_out_ds!=NULL) delete[] _out_ds;
  if (_io_ds!=NULL) delete[] _io_ds;
}

// -----------------------------------------------
// Graph manipulation
// -----------------------------------------------

void Graph::zero() {
  int i,j;

  _num_edges=0;
  for (i=0; i<_num_nodes;i++) {
    _in[i] = 0;
    _out[i] = 0;
    _adjL[i].clear();
    _adjLI[i].clear();
    _neighbours[i].clear();
  }

}

// -----------------------------------------------

int *Graph::makeOutDegreeSequence() {
  int i;
  for (i=0; i<numNodes(); i++)
    _out_ds[i]=nodeEdges(i);
  sort(_out_ds, _out_ds+numNodes());
  return _out_ds;
}

int *Graph::makeInDegreeSequence() {
  int i;
  for (i=0; i<numNodes(); i++)
    _in_ds[i]=nodeInEdges(i);
  sort(_in_ds, _in_ds+numNodes());
  return _in_ds;
}

iPair *Graph::makeIODegreeSequence() {
  int i;
  for (i=0; i<numNodes(); i++) {
    _io_ds[i].first = nodeInEdges(i);
    _io_ds[i].second = nodeEdges(i);
  }
  //sort(_io_ds, _io_ds+numNodes());
  return _io_ds;
}


// -----------------------------------------------

int Graph::readFile(char *file_name, bool undir) {
  int i;
  FILE *f = fopen(file_name, "r");

  if (!f) {
    fprintf(stderr,"Could not open \"%s\"!\n",file_name);
    return 0;
  }
  
  int size = 0, a, b, c, max=0;
  std::vector<int> va, vb, vc;

  while (fscanf(f,"%d %d %d", &a, &b, &c)==3) {
      va.push_back(a);
      vb.push_back(b);
      vc.push_back(c);
      if (a>max) max=a;
      if (b>max) max=b;
      size++;
  }

  fclose(f);

  _createGraph(max);

  printf("!! %d %d\n", max, size);

  if (undir) {
    for (i=0; i<size; i++) {
      if (!getEdge(va[i]-1,vb[i]-1))
	addEdge(va[i]-1,vb[i]-1,vc[i]);
      if (!getEdge(vb[i]-1,va[i]-1))
	addEdge(vb[i]-1,va[i]-1,vc[i]);
    }
  } else {
    for (i=0; i<size; i++)
      if (!getEdge(va[i]-1,vb[i]-1))
	addEdge(va[i]-1,vb[i]-1,vc[i]);
  }
  

  return 1;
}

// -----------------------------------------------
// Collect edges
// -----------------------------------------------

void Graph::vecEdges(vEdges &v) {
  int i, count = 0;
  list<iPair>::iterator ii;
  
  v.resize(numEdges());

  for (i=0; i<numNodes(); i++) 
    for(ii=_adjL[i].begin(); ii!=_adjL[i].end(); ii++) {
      v[count].first = i;
      v[count++].second = (*ii).first;
    }
  
}

list<iPair> Graph::adjL(int v) {
  return _adjL[v];
}
