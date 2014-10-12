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

#ifndef _GRAPHMATRIX_
#define _GRAPHMATRIX_

#include "Graph.h"

class GraphMatrix : public Graph {
 private:
  GraphType _type;

  int _num_nodes;
  int _num_edges;

  int *_in;
  int *_out;
  int *_num_neighbours;

  bool **_adjM;
  int  **_array_neighbours;
  vector<int> *_adjOut;
  vector<int> *_adjIn;
  vector<int> *_neighbours;

  void _init();
  void _delete();
  void _removeVector(vector<int> &v, int b);

 public:
  GraphMatrix();
  ~GraphMatrix();

  bool **adjacencyMatrix() {return _adjM;}

  void createGraph(int n, GraphType t);

  GraphType type() {return _type;}

  void zero();

  int numNodes() {return _num_nodes;}
  int numEdges() {return _num_edges;}

  void addEdge(int a, int b); // add edge from a to b
  void rmEdge(int a, int b);  // remove edge from a to b

  bool hasEdge(int a, int b) {return _adjM[a][b];} 
  bool isConnected(int a, int b)  {return _adjM[a][b] || _adjM[b][a];}

  int nodeOutEdges(int a) {return _out[a];}
  int nodeInEdges(int a)  {return _in[a];}
  int numNeighbours(int a) {return _num_neighbours[a];}
  void sortNeighbours();
  void sortNeighboursArray();
  void makeArrayNeighbours();
  void makeVectorNeighbours();

  vector<int> *neighbours(int a) {return &_neighbours[a];}
  int **matrixNeighbours()       {return _array_neighbours;}
  int *arrayNeighbours(int a)    {return _array_neighbours[a];}
  int *arrayNumNeighbours()      {return _num_neighbours;}
  vector<int> *outEdges(int a)   {return &_adjOut[a];}
  vector<int> *inEdges(int a)    {return &_adjIn[a];}
};

#endif


