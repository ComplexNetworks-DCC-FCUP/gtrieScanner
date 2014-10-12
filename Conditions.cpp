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
Symmetry Breaking Conditions Implementation

Last Update: 11/02/2012
---------------------------------------------------- */

#include "Conditions.h"
#include <algorithm>
#include <stdlib.h>

int Conditions::_subgraph_size;

// -------------------------------------------------

int pos;
GMap *fastf;
Graph *fastg;
VVsmallNode *fastvv;
bool *fasts;
bool **adj;

int int_compare(const void *a, const void *b) {
  return (*((int *)b)) - (*((int *)a));
}

GMap::GMap(int nf, int nr) {
  int i;

  numf = nf;
  numr = nr;
  f = new int[nf];
  r = new int[nr];

  for (i=0; i<nf; i++)
    f[i] = INVALID;

  for (i=0; i<nr; i++)
    r[i] = INVALID;
}

GMap::~GMap() {
  delete[] f;
  delete[] r;
}

// -------------------------------------------------

void Conditions::findAutomorphisms(Graph *G, VVsmallNode *vv) {
  int i,j,k;
  int g;

  _subgraph_size = G->numNodes();
  GMap *f = new GMap(_subgraph_size, _subgraph_size);

  int sequence[_subgraph_size][_subgraph_size];

  for (i=0; i<_subgraph_size; i++)
    for (j=0; j<_subgraph_size; j++)
      if (G->isConnected(i,j))
	sequence[i][j] = G->numNeighbours(j);
      else
	sequence[i][j] = 0;

  for (i=0; i<_subgraph_size; i++)
    qsort(sequence[i], _subgraph_size, sizeof(int), int_compare);

  bool support[_subgraph_size*_subgraph_size];
  for (i=0; i<_subgraph_size; i++)
    for (j=0; j<_subgraph_size; j++) {
      for (k=0; k<_subgraph_size; k++)
	if (sequence[i][k] != sequence[j][k]) break;      
      if (k<_subgraph_size) support[i*_subgraph_size+j]=false;
      else support[i*_subgraph_size+j]=true;
    }

  vv->clear();
  fastf  = f;
  fastg  = G;
  fastvv = vv;
  fasts  = support;
  adj = G->adjacencyMatrix();
  for (g=0; g<_subgraph_size; g++)
    if (support[g*_subgraph_size]) {
      f->add(0,g);
      pos = 1;
      isomorphicExtensions();
      f->remove(0);
    }

  delete f;
}

void Conditions::isomorphicExtensions() {
  int i, j, cand[_subgraph_size], ncand;
  int *v, num;

  if (pos == _subgraph_size) {

    smallNode *v = new smallNode[_subgraph_size];
    for (i=0; i<fastf->numf; i++) v[i]=fastf->f[i];
    fastvv->push_back(v);

  } else {    

    list<iPair>::iterator ii;
    int n, m;
    int flag;

    int count[_subgraph_size];

    ncand=0;
    for (i=0; i<_subgraph_size; i++) count[i]=0;

    for (i=0; i<_subgraph_size; i++)     // For all nodes of H already mapped
      if (fastf->f[i]!=INVALID) {        // find their not mapped neighbours
	v = fastg->arrayNeighbours(i);    
	num = fastg->numNeighbours(i);
	for (j=0; j<num; j++)
	  if (fastf->f[v[j]]==INVALID) {
	    if (count[v[j]]==0) 
	      cand[ncand++]=v[j];	    
	    count[v[j]]++;
	  }	
      }

   // Find most constrained neighbour 'm' (with more mapped neighbours)
    m = 0;
    for (i=1; i<ncand; i++)    
      if (count[i]>count[m])  // Later: add more restraining conditions??
	m = i;
    m = cand[m];

    ncand=0;
    bool already[_subgraph_size];
    for (i=0; i<_subgraph_size; i++) already[i]=false;

    for (i=0; i<_subgraph_size; i++)  // For all nodes of G already mapped
      if (fastf->f[i]!=INVALID) {         // find their not mapped neighbours 
	v = fastg->arrayNeighbours(fastf->f[i]);    
	num = fastg->numNeighbours(fastf->f[i]);    
	for (j=0; j<num; j++)
	  if (!already[v[j]] && fastf->r[v[j]]==INVALID && fasts[m*_subgraph_size+v[j]])  {
	    cand[ncand++]=v[j];
	    already[v[j]]=true;
	  }
      }
        
    for (i=0; i<ncand; i++) {
      n = cand[i];

      flag = false;

      for (j=0; j<_subgraph_size; j++)
	if (fastf->f[j]!=INVALID) {
	  if      (adj[m][j] != adj[n][fastf->f[j]])      {flag=true; break;}
	  else if (adj[j][m] != adj[fastf->f[j]][n])      {flag=true; break;}
	}

      if (!flag) {	
	fastf->add(m, n);
	pos++;
	isomorphicExtensions();
	pos--;
	fastf->remove(m);
      }
    }
  }
}

void Conditions::symmetryConditions(Graph *g, list<iPair> *cond) {
  int i, j, k, vvsize, size = g->numNodes();
  iPair p;
  VVsmallNode vv;

  cond->clear();
  findAutomorphisms(g, &vv);

  vvsize = vv.size();
  bool broken[vvsize];
  for (i=0; i<vvsize; i++) broken[i]=false;
  
  for (i=0; i<size; i++) {
    vvsize = vv.size();
    for (j=0; j<vvsize; j++)
      if (!broken[j] && vv[j][i] != i) 
	break;
      
      // There are still nodes not fixed
      if (j<vvsize)
	for (k=i+1; k<size; k++) 
	  for (j=0; j<int(vv.size()); j++)
	    if (!broken[j] && vv[j][i] == k) {
	      p.first  = i;
	      p.second = k;
	      //printf("%d < %d\n", p.first, p.second);
	      cond->push_back(p);
	      break;
	    }
      
      // Reduce set of automorphisms to set that fix 'i'
      for (j=0; j<vvsize; j++)	
	if (vv[j][i] != i)
	  broken[j]=true;
  }

  for (j=0; j<vvsize; j++)
    delete [] vv[j];
  
}

// -------------------------------------------------
