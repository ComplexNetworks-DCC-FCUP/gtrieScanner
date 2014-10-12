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

#include "GraphTree.h"
#include "Isomorphism.h"

GraphTreeNode::GraphTreeNode() {
  frequency = 0;
  zero = one = NULL;
}

GraphTreeNode::~GraphTreeNode() {
  if (zero != NULL) delete zero;
  if ( one != NULL) delete  one;
}

void GraphTreeNode::zeroFrequency() {
  frequency = 0;
  if (zero != NULL) zero->zeroFrequency();
  if (one  != NULL)  one->zeroFrequency();
}

void GraphTreeNode::incrementString(int pos, char *s) {
  if (s[pos]==0) frequency++;
  else {
    if (s[pos]=='1') {
      if (one==NULL) one = new GraphTreeNode();
      one->incrementString(pos+1, s);
    } else if (s[pos]=='0') {
      if (zero==NULL) zero = new GraphTreeNode();
      zero->incrementString(pos+1, s);
    }
  }
}

void GraphTreeNode::setString(int pos, char *s, int f) {
  if (s[pos]==0) frequency=f;
  else {
    if (s[pos]=='1') {
      if (one==NULL) one = new GraphTreeNode();
      one->setString(pos+1, s, f);
    } else if (s[pos]=='0') {
      if (zero==NULL) zero = new GraphTreeNode();
      zero->setString(pos+1, s, f);
    }
  }
}

void GraphTreeNode::addString(int pos, char *s, int f) {
  if (s[pos]==0) frequency+=f;
  else {
    if (s[pos]=='1') {
      if (one==NULL) one = new GraphTreeNode();
      one->addString(pos+1, s, f);
    } else if (s[pos]=='0') {
      if (zero==NULL) zero = new GraphTreeNode();
      zero->addString(pos+1, s, f);
    }
  }
}

void GraphTreeNode::showFrequency(int pos, char *s) {
  if (zero == NULL && one == NULL) {
    s[pos]=0;
    printf("%s: %d\n", s, frequency);
  } else {
    if (zero != NULL) {
      s[pos]='0';
      zero->showFrequency(pos+1, s);
    }
    if (one != NULL) {
      s[pos]='1';
      one->showFrequency(pos+1, s);
    }
  }
}


bool GraphTreeNode::equal(GraphTreeNode *gt, int pos, char *s) {
  if (zero == NULL && gt->zero!=NULL) return false;
  if (zero != NULL && gt->zero==NULL) return false;
  if (one  == NULL && gt->one !=NULL) return false;
  if (one  != NULL && gt->one ==NULL) return false;

  if (zero == NULL && one == NULL) {
    s[pos]=0;
    if (frequency != gt->frequency) 
      printf("NOT EQUAL TREE %s: %d != %d\n", s, frequency, gt->frequency);
    return (frequency == gt->frequency);
  } else {
    bool fzero, fone;
    fzero = fone = true;
    if (zero != NULL) {
      s[pos]='0';
      fzero = zero->equal(gt->zero, pos+1, s);
    }
    if (one != NULL) {
      s[pos]='1';
      fone = one->equal(gt->one, pos+1, s);
    }
    return fzero && fone;
  }
}


bool GraphTreeNode::equal(GTrie *gt, int size, int pos, char *s) {
  if (zero == NULL && one == NULL) {
    s[pos]=0;
    int aux = gt->frequencyGraphString(size, s);
    if (frequency != aux)
      printf("NOT EQUAL GTRIE %s: %d != %d\n", s, frequency, aux);
    return (frequency == aux);
  } else {
    bool fzero, fone;
    fzero = fone = true;
    if (zero != NULL) {
      s[pos]='0';
      fzero = zero->equal(gt, size, pos+1, s);
    }
    if (one != NULL) {
      s[pos]='1';
      fone = one->equal(gt, size, pos+1, s);
    }
    return fzero && fone;
  }
}

void GraphTreeNode::populateGTrie(GTrie *gt, int size, int pos, char *s) {
  if (zero == NULL && one == NULL) {
    s[pos]=0;
    gt->insertGraphString(size, s);
  } else {
    if (zero != NULL) {
      s[pos]='0';
      zero->populateGTrie(gt, size, pos+1, s);
    }
    if (one != NULL) {
      s[pos]='1';
      one->populateGTrie(gt, size, pos+1, s);
    }
  }
}

void GraphTreeNode::populateMap(mapStringInt *m, int size, int pos, char *s) {
  if (zero == NULL && one == NULL) {
    s[pos]=0;
    char s2[size*size+1];
    Isomorphism::canonicalBasedNauty(s, s2, size);
    (*m)[s2] = frequency;
  } else {
    if (zero != NULL) {
      s[pos]='0';
      zero->populateMap(m, size, pos+1, s);
    }
    if (one != NULL) {
      s[pos]='1';
      one->populateMap(m, size, pos+1, s);
    }
  }
}


void GraphTreeNode::populateGTrieNauty(GTrie *gt, int size, int pos, char *s, bool dir) {
  if (frequency>0) {
    s[pos]=0;
    gt->insertGraphNautyString(size, s, dir, 1);
  }
   
  if (zero != NULL) {
    s[pos]='0';
    zero->populateGTrieNauty(gt, size, pos+1, s, dir);
  }
  if (one != NULL) {
    s[pos]='1';
    one->populateGTrieNauty(gt, size, pos+1, s, dir);
  }
}

double GraphTreeNode::countOccurrences() {
  double aux = frequency;
  if (zero != NULL) aux += zero->countOccurrences();
  if (one != NULL)  aux +=  one->countOccurrences();
  return aux;
}

int GraphTreeNode::countGraphs() {
  int aux = 0;
  if (frequency>0) aux++;
  if (zero != NULL) aux += zero->countGraphs();
  if (one != NULL)  aux +=  one->countGraphs();
  return aux;
}

// -----------------------------------

GraphTree::GraphTree() {
  root = new GraphTreeNode();  
}

GraphTree::~GraphTree() {
  if (root) delete root;
}

void GraphTree::zeroFrequency() {
  root->zeroFrequency();
}

void GraphTree::incrementString(char *s) {
  root->incrementString(0, s);
}

void GraphTree::setString(char *s, int f) {
  root->setString(0, s, f);
}

void GraphTree::addString(char *s, int f) {
  root->addString(0, s, f);
}

void GraphTree::showFrequency(int maxsize) {
  char s[maxsize*maxsize+1];  
  root->showFrequency(0, s);
}

bool GraphTree::equal(GraphTree *gt, int maxsize) {
  char s[maxsize*maxsize+1];  
  return root->equal(gt->root, 0, s);
}

bool GraphTree::equal(GTrie *gt, int maxsize) {
  char s[maxsize*maxsize+1];  
  return root->equal(gt, maxsize, 0, s);
}

void GraphTree::populateGTrie(GTrie *gt, int maxsize) {
  char s[maxsize*maxsize+1];  
  root->populateGTrie(gt, maxsize, 0, s);
}

void GraphTree::populateMap(mapStringInt *m, int maxsize) {
  char s[maxsize*maxsize+1];  
  root->populateMap(m, maxsize, 0, s);
}

void GraphTree::populateGTrieNauty(GTrie *gt, int maxsize, bool dir) {
  char s[maxsize*maxsize+1];  
  root->populateGTrieNauty(gt, maxsize, 0, s, dir);
}

double GraphTree::countOccurrences() {
  return root->countOccurrences();
}

int GraphTree::countGraphs() {
  return root->countGraphs();
}
