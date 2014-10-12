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

#ifndef _CONDITIONS_
#define _CONDITIONS_

#include "Common.h"
#include "Graph.h"
#include "GTrie.h"

#include <iostream>

class GMap {
 public:
  int numf, numr;
  int *f;
  int *r;

  GMap(int ng, int nr); 
  ~GMap();
  
  void add(int n, int v);
  void remove(int n);  
};

class LessConditions {
 public:
  int size;

  vector<iPair> *v;

  LessConditions(int s);
  ~LessConditions();

  void add(int root, int a, int b);
  void clear();
  bool pass(int root, GMap &f);
};

class Conditions {
 private:
    static int _subgraph_size;

 public:
    static void findAutomorphisms(Graph *g, VVsmallNode *vv);
    static void isomorphicExtensions();
    static void symmetryConditions(Graph *g, list<iPair> *cond);
};

inline
void GMap::add(int n, int v) {
  f[n]=v;
  r[v]=n;
}

inline
void GMap::remove(int n) {
  r[f[n]]=INVALID;
  f[n]=INVALID;
}

inline
bool LessConditions::pass(int root, GMap &f) {
  int a,b;
  vector<iPair>::iterator ii;

  for(ii=v[root].begin(); ii!=v[root].end(); ii++) {
    a = f.f[ii->first];
    b = f.f[ii->second];    
    if (a!=INVALID && b!=INVALID && a>b) {      
      //putchar('*'); fflush(stdout);
      return false;
    }
  }

  return true;
}

inline 
void LessConditions::add(int root, int a, int b) {
  iPair p;
  p.first = a;
  p.second = b;
  v[root].push_back(p);
}


#endif

