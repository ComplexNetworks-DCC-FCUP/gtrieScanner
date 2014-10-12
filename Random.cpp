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
Randomization methods

Last Update: 11/02/2012
---------------------------------------------------- */

#include "Random.h"

// Initialize pseudo-random generator with seed 's'
void Random::seed(int s) {
  srandom(s);
}

// Pseudo-Random number between 'a' and 'b' (inclusive)
int Random::getInteger(int a, int b) {
  double aux = random() / (double)RAND_MAX;  
  return a + aux*(b-a+1);
}

 // Pseudo-Random number between 0 and 1
double Random::getDouble() {
  return random() / (double)RAND_MAX;  
}

  // Randomize 'g' network with 'num' exchanges per edge and 'tries' attempts per edge
void Random::markovChainPerEdge(Graph *g, int num, int tries) {
  int i, j, k, n, edges, nodes = g->numNodes();
  int a, b, c, d, aux;
  vector<int> *v, *u;
  vector<int>:: iterator ii;

  a=b=c=d=0;
  for (n=0; n<num; n++)
    for (i=0; i<nodes; i++) {
      v = g->outEdges(i);
      edges = v->size();
      a = i;
      for (j=0; j<edges; j++) {
	b = (*v)[j];
	for (k=0; k<tries; k++) {
	  c = getInteger(0, nodes-1);
	  if (a==c || b==c || g->hasEdge(c,b)) continue;
	  aux = g->nodeOutEdges(c); 
	  if (aux==0) continue;
	  u = g->outEdges(c);
	  d = getInteger(1, aux);
	  d = (*u)[d-1];
	  if (a==d || b==d || g->hasEdge(a,d)) continue;
	  break;
	}
	if (k<tries) { // Found an edge to swap!
	  fflush(stdout);
	  g->rmEdge(a,b);  g->rmEdge(c,d);
	  g->addEdge(a,d); g->addEdge(c,b);
	  if (g->type() == UNDIRECTED) {
	    g->rmEdge(b,a);  g->rmEdge(d,c);
	    g->addEdge(d,a); g->addEdge(b,c);
	  }
	}
      }
    }
}


