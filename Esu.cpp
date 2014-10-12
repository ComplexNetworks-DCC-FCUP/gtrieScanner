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
Esu implementation

Last Update: 11/02/2012
---------------------------------------------------- */

#include "Esu.h"
#include "Isomorphism.h"
#include "Random.h"

// Class static variables
int     Esu::_motif_size = 0;
int     Esu::_graph_size = 0;
int     Esu::_next = 0;
int    *Esu::_current = NULL;
int    *Esu::_ext = NULL;
Graph  *Esu::_g = NULL;
double *Esu::_prob;
GraphTree *Esu::_sg;


/*! Recursively extend a partial subgraph
    \param n the current position in the constructed subgraph 
    \param size the subgraph size
    \param next number of elements in the list "ext"
    \param ext of the nodes that can be used to extend the subgraph */
void Esu::_go(int n, int size, int next, int *ext) {
  _current[size++] = n;

  if (size==_motif_size) {
    char s[_motif_size*_motif_size+1];
    Isomorphism::canonicalStrNauty(_g, _current, s);
    _sg->incrementString(s);

    if (Global::show_occ) {
      fprintf(Global::occ_file, "%s:", s);
      for (int i=0; i<size; i++)
	fprintf(Global::occ_file, " %d", _current[i]+1);
      fputc('\n', Global::occ_file);
    }

  } else {
    int i,j;
    int next2 = next;
    int ext2[_graph_size];

    for (i=0;i<next;i++) ext2[i] = ext[i];

    int *v  = _g->arrayNeighbours(_current[size-1]);
    int num = _g->numNeighbours(_current[size-1]);


    for (i=0;i<num;i++) {
      if (v[i]<=_current[0]) continue;
      for (j=0;j+1<size;j++)
	if (_g->isConnected(v[i],_current[j])) break;
      if (j+1==size) ext2[next2++]=v[i];      
    }

    while (next2>0) {      
      next2--;
      _go(ext2[next2], size, next2, ext2);
    }
  }
  
}

/*! Recursively extend a partial subgraph (sampling version)
    \param n the current position in the constructed subgraph 
    \param size the subgraph size
    \param next number of elements in the list "ext"
    \param ext of the nodes that can be used to extend the subgraph */
void Esu::_goSample(int n, int size, int next, int *ext) {

  _current[size++] = n;

  if (size==_motif_size) {
    char s[_motif_size*_motif_size+1];
    Isomorphism::canonicalStrNauty(_g, _current, s);
    _sg->incrementString(s);
  } else {
    int i,j;
    int next2 = next;
    int ext2[_graph_size];

    for (i=0;i<next;i++) ext2[i] = ext[i];

    int *v  = _g->arrayNeighbours(_current[size-1]);
    int num = _g->numNeighbours(_current[size-1]);
    for (i=0;i<num;i++) {
      if (v[i]<=_current[0]) continue;
      for (j=0;j+1<size;j++)
	if (_g->isConnected(v[i],_current[j])) break;
      if (j+1==size) ext2[next2++]=v[i];      
    }
    while (next2>0) {
      next2--;
      if (Random::getDouble()<=_prob[size])
	_goSample(ext2[next2], size, next2, ext2);
    }
  }
  
}

/*! Make a complete k-census of a Graph
    \param g the graph to be explored
    \param k the size of the subgraphs
    \param sg The GraphTree where the results should be stored */
void Esu::countSubgraphs(Graph *g, int k, GraphTree *sg) {
  int i, v[1];

  _motif_size = k;
  _graph_size = g->numNodes();
  _current = new int[k];
  _ext = new int[_graph_size];
  _next = 0;
  _g = g;
  _sg = sg;

  sg->zeroFrequency();

   for (i=0; i<_graph_size; i++)
    _go(i, 0, 0, v);
  
  delete[] _current;
  delete[] _ext;
}

/*! Make a complete k-census of a Graph (sampling version)
    \param g the graph to be explored
    \param k the size of the subgraphs
    \param sg The GraphTree where the results should be stored
    \param p array of probabilities to use in sampling*/
void Esu::countSubgraphsSample(Graph *g, int k, GraphTree *sg, double *p) {
  int i, v[1];

  _motif_size = k;
  _graph_size = g->numNodes();
  _current = new int[k];
  _ext = new int[_graph_size];
  _next = 0;
  _g = g;
  _sg = sg;
  _prob = p;

  sg->zeroFrequency();

  for (i=0; i<_graph_size; i++)
    if (Random::getDouble()<=_prob[0])
      _goSample(i, 0, 0, v);
  
  delete[] _current;
  delete[] _ext;
}
