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
Pedro Paredes - CRACS & INESC-TEC, DCC/FCUP

----------------------------------------------------
Graphs Implementation with 

Last Update: 27/09/2014
---------------------------------------------------- */

#ifndef _DYNGRAPH_
#define _DYNGRAPH_

#include "Graph.h"

typedef enum{MATRIX, BSLIST} RepType;

class DynamicGraph : public Graph {
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

 public:
  DynamicGraph();
  DynamicGraph(RepType _r);
  ~DynamicGraph();

  RepType _rtype;

  void createGraph(int n, GraphType t);

  GraphType type() {return _type;}

  void zero();
  void prepareGraph();

  int numNodes() {return _num_nodes;}
  int numEdges() {return _num_edges;}

  void addEdge(int a, int b); // add edge from a to b

  bool hasEdge(int a, int b);

  int nodeOutEdges(int a) {return _out[a];}
  int nodeInEdges(int a)  {return _in[a];}
  int numNeighbours(int a) {return _num_neighbours[a];}
  void sortNeighbours();
  void sortNeighboursArray();
  void makeArrayNeighbours();
  void makeVectorNeighbours();

  vector<int> *neighbours(int a) {return &_neighbours[a];}
  int *arrayNeighbours(int a)    {return _array_neighbours[a];}
  int *arrayNumNeighbours()      {return _num_neighbours;}
  vector<int> *outEdges(int a)   {return &_adjOut[a];}
  vector<int> *inEdges(int a)    {return &_adjIn[a];}
};

#endif
