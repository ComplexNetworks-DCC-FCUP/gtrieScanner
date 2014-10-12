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
Isomorphism Utilities

Last Update: 11/02/2012
---------------------------------------------------- */

#ifndef _ISOMORPHISM_
#define _ISOMORPHISM_

#include "Common.h"
#include "Graph.h"

#define MAXN MAX_MOTIF_SIZE
#define WORKSPACE_SIZE MAXN*160

#include "nauty/nauty.h"

class Isomorphism {
 private:  
  static bool dir;
  static setword workspace[WORKSPACE_SIZE];
  static int n,m,lab[MAXN], ptn[MAXN], orbits[MAXN];
  static set *gv;
  static graph g[MAXN*MAXM];

  static void _goCan(int x, int pos, const char *in, 
		     char *perm, char *used,
		     char *current, char *best, int size);

  static void _goCan2(int pos, const char *in, int *perm, bool *used, char *best, int size, char *current);
    
 public:
  static void initNauty(int size, bool dir);
  static void finishNauty();

  static void canonicalStrNauty(Graph *g, int *v, char *s);

  static void canonicalNauty(const char *in, char *out, int size);
  static void canonicalBigger(const char *in, char *out, int size);
  static void canonicalBigger2(const char *in, char *out, int size);
  static void canonicalBasedNauty(const char *in, char *out, int size); // GT Canon
  static void canonicalBasedNautyReverse(const char *in, char *out, int size);
  static void canonicalRandom(const char *in, char *out, int size);
};

#endif


