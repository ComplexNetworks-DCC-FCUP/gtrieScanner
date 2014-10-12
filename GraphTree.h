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
Graph (0-1) Tree Implementation

Last Update: 11/02/2012
---------------------------------------------------- */

#ifndef _GRAPHTREE_
#define _GRAPHTREE_

#include "Common.h"
#include "GTrie.h"

class GraphTreeNode   {
 public:
  GraphTreeNode();
  ~GraphTreeNode();

  int frequency;
  GraphTreeNode *zero, *one;

  void zeroFrequency();
  void incrementString(int pos, char *s);
  void setString(int pos, char *s, int f);
  void addString(int pos, char *s, int f);
  void showFrequency(int pos, char *s);

  void populateGTrie(GTrie *gt, int size, int pos, char *s);
  void populateMap(mapStringInt *m, int size, int pos, char *s);
  void populateGTrieNauty(GTrie *gt, int size, int pos, char *s, bool dir);

  bool equal(GraphTreeNode *gt,   int pos, char *s);
  bool equal(GTrie *gt, int size, int pos, char *s);

  int countGraphs();
  double countOccurrences();
};

class GraphTree   {
 private:

 public:
  GraphTreeNode *root;

  GraphTree();
  ~GraphTree();

  void zeroFrequency();
  void incrementString(char *s);
  void setString(char *s, int f);
  void addString(char *s, int f);
  void showFrequency(int maxsize);

  bool equal(GraphTree *gt, int maxsize);  
  bool equal(GTrie *gt,     int maxsize);  
  void populateGTrie(GTrie *gt, int maxsize);
  void populateMap(mapStringInt *m, int maxsize);
  void populateGTrieNauty(GTrie *gt, int maxsize, bool dir);

  int countGraphs();
  double countOccurrences();
};

#endif
