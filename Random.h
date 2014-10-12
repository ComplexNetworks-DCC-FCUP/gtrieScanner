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

#ifndef _RANDOM_
#define _RANDOM_

#include "Common.h"
#include "Graph.h"

class Random {

 public:
  static void seed(int s);             // Initialize pseudo-random generator with seed 's'
  static int getInteger(int a, int b); // Pseudo-Random number between 'a' and 'b' (inclusive)
  static double getDouble();           // Pseudo-Random number between 0 and 1
    
  // Randomize 'g' network with 'num' exchanges per edge and 'tries' attempts per edge
  static void markovChainPerEdge(Graph *g, int num, int tries);
};

#endif
