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

#ifndef _ESU_
#define _ESU_

#include "Common.h"
#include "GraphTree.h"

/*! This class implements the ESU subgraph enumeration algorithm */
class Esu {
 private:
  static int *_current;
  static int * _ext;
  static int _next;
  static int _graph_size;
  static int _motif_size;
  static Graph * _g;
  static GraphTree *_sg;
  static double *_prob;

  static void _go(int n, int size, int next, int *ext);
  static void _goSample(int n, int size, int next, int *ext);

 public:
  static void countSubgraphs(Graph *g, int k, GraphTree *sg);
  static void countSubgraphsSample(Graph *g, int k, GraphTree *sg, double *p);
  
};

#endif
