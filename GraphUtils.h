/* -------------------------------------------------
      _       _     ___                            
 __ _| |_ _ _(_)___/ __| __ __ _ _ _  _ _  ___ _ _ 
/ _` |  _| '_| / -_)__ \/ _/ _` | ' \| ' \/ -_) '_|
\__, |\__|_| |_\___|___/\__\__,_|_||_|_||_\___|_|  
|___/                                          
    
gtrieScanner: quick discovery of network motifs
Released under Artistic License 2.0
(see README and LICENSE)

----------------------------------------------------
Graph Utilities

Last Update: 11/02/2012
---------------------------------------------------- */

#ifndef _GRAPHUTILS_
#define _GRAPHUTILS_

#include "Graph.h"

class GraphUtils {
 private:
  static int *_degree;
  static int **_ds;
  static int *_neighbours;
  
 public:

  // Compare two integers
  static int int_compare(const void *a, const void *b);

  // Read file 's', with direction 'dir' to graph 'g'
  static void readFileTxt(Graph *g, const char *s, bool dir, bool weight);

  // Convert adjacency matrix to graph of 'size' nodes
  static void strToGraph(Graph *g, const char *s, int size, bool dir);

  static void orderGraph(Graph *old, Graph *g);
  static int _compare_int(const void *a, const void *b);
  static int _compare_int_descending(const void *a, const void *b);
  static int _compare_degree(const void *a, const void *b);
  
};

#endif
