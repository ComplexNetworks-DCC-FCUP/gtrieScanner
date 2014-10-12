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
G-Trie Implementation and associated methods

Last Update: 11/02/2012
---------------------------------------------------- */

#include "GTrie.h"
#include "GraphMatrix.h"
#include "GraphTree.h"
#include "Isomorphism.h"
#include "GraphUtils.h"
#include "Conditions.h"
#include "Random.h"
#include "Error.h"
#include <iostream>
#include <string.h>

#include <assert.h>

int *GTrieNode::mymap;
bool *GTrieNode::used;
bool **GTrieNode::adjM;
int **GTrieNode::fastnei;
int *GTrieNode::numnei;
int GTrieNode::glk;
int GTrieNode::numNodes;
bool GTrieNode::isdir;
double *GTrieNode::prob;

list< list<iPair> >::const_iterator jj, jjend;
list<iPair>::const_iterator kk, kkend;

GTrieNode::GTrieNode(int d) {

  depth     = d;

  is_graph  = false;
  frequency = 0;
  nconn=0;

  if (d!=0) {
    in   = new bool[d];
    out  = new bool[d];
    conn = new int[d];
  } else {
    in = out = NULL;
    conn = NULL;
  }

  total_in = total_out = total_edges = 0;

  child.clear();
  cond.clear();
  this_node_cond.clear();

  cond_ok = false;
  cond_this_ok = false;
}

GTrieNode::~GTrieNode() {
  if (depth!=0) {
    if (in  != NULL) delete[] in;
    if (out != NULL) delete[] out;
    if (conn != NULL) delete[] conn;
  }
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    delete *ii;
}

void GTrieNode::showAsText(FILE *f) {
  int i;
  for (i=0; i<depth; i++) fprintf(f,"   ");
  fprintf(f,"[");
  for (i=0; i<depth; i++) 
    fprintf(f,"%s", out[i]?"1":"0");
  fprintf(f,"][");
  for (i=0; i<depth; i++) 
    fprintf(f,"%s", in[i]?"1":"0");
  fprintf(f,"] ");

  if (cond_ok) fprintf(f,"{}");
  else {    
    list< list<iPair> >:: iterator jj;
    list<iPair>:: iterator kk;
    for (jj=cond.begin(); jj!=cond.end(); jj++) {
      fputc('|', f);
      for (kk=(*jj).begin(); kk!=(*jj).end(); kk++) {
	if (kk!=(*jj).begin()) fputc(' ', f);
	fprintf(f,"%d<%d", kk->first, kk->second);
      }
      fputc('|', f);
    }
  }

  fputc('+', f);

  if (cond_this_ok) fprintf(f,"{}");
  else {
    list< list<int> >:: iterator jjj;
    list<int>:: iterator kkk;
    for (jjj=this_node_cond.begin(); jjj!=this_node_cond.end(); jjj++) {
      fputc('|', f);
      if ((*jjj).size()==0) fprintf(f,"{}");
      else {
	for (kkk=(*jjj).begin(); kkk!=(*jjj).end(); kkk++) {
	  if (kkk!=(*jjj).begin()) fputc(',', f);
	  fprintf(f,"%d", *kkk);	
	}
	fprintf(f,"<%d", depth-1);
      }
      fputc('|', f);
    }
  }
      
  fputc(' ', f);
        
  fprintf(f,"%s\n", is_graph?"isGraph":"");

  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    (*ii)->showAsText(f);
}

int GTrieNode::frequencyGraph(Graph *g) {
  if (g->numNodes()==depth)
    return frequency;
  list<GTrieNode *>::iterator ii;
  int i;
  for(ii=child.begin(); ii!=child.end(); ii++) {
    for (i=0; i<=depth; i++)
      if (g->hasEdge(depth, i) != (*ii)->out[i] ||
	  g->hasEdge(i, depth) != (*ii)->in[i])
	break;
    if (i>depth) break;
  }
  if (ii!=child.end())
    return (*ii)->frequencyGraph(g);
  else
    return -1;
}

void GTrieNode::insert(Graph *g) {
  int i;


  if (g->numNodes() == depth) {
    is_graph = true;
  } else {
    list<GTrieNode *>::iterator ii;
    for(ii=child.begin(); ii!=child.end(); ii++) {
      for (i=0; i<=depth; i++)
	if (g->hasEdge(depth, i) != (*ii)->out[i] ||
	    g->hasEdge(i, depth) != (*ii)->in[i])
	  break;
      if (i>depth) break;
    }
    if (ii==child.end()) {
      GTrieNode *gt = new GTrieNode(depth+1);
      for (i=0; i<=depth; i++) {
      	gt->out[i] = g->hasEdge(depth, i);	
      	gt->in[i]  = g->hasEdge(i, depth);
	if (g->isConnected(depth, i)) 
	  gt->conn[(gt->nconn)++] = i;

	if (gt->in[i])  { (gt->total_in)++;  (gt->total_edges)++;}
	if (gt->out[i]) { (gt->total_out)++; (gt->total_edges)++;}
      }
      gt->insert(g);
      child.push_back(gt);
    } else
      (*ii)->insert(g);
  }
}


void GTrieNode::insertCond(Graph *g, list<iPair> *cond) {
  int i;
  GTrieNode *gt;

  if (g->numNodes() == depth) {
    is_graph = true;
  } else {
    list<GTrieNode *>::iterator ii;
    for(ii=child.begin(); ii!=child.end(); ii++) {
      for (i=0; i<=depth; i++)
	if (g->hasEdge(depth, i) != (*ii)->out[i] ||
	    g->hasEdge(i, depth) != (*ii)->in[i])
	  break;
      if (i>depth) break;
    }
    if (ii==child.end()) {
      gt = new GTrieNode(depth+1);
      for (i=0; i<=depth; i++) {
      	gt->out[i] = g->hasEdge(depth, i);	
      	gt->in[i]  = g->hasEdge(i, depth);
	if (g->isConnected(depth, i)) 
	  gt->conn[(gt->nconn)++] = i;
	if (gt->in[i])  { (gt->total_in)++;  (gt->total_edges)++;}
	if (gt->out[i]) { (gt->total_out)++; (gt->total_edges)++;}
      }
      child.push_back(gt);
    } else {
      gt = (*ii);
    }

    gt->insertConditionsFiltered(cond);
    gt->insertCond(g, cond);

  }
}

// Assume a and b are in ascending order
bool GTrieNode::_is_intset_included(list<int> a, list<int> b) {
  list<int>::iterator aa, bb;

  aa = a.begin();
  bb = b.begin();

  while (aa!=a.end() && bb!=b.end()) {
    if (*aa == *bb) {
      aa++;
      bb++;
    } else {
      bb++;
    }      
  }

  if (aa == a.end()) return true;
  else return false;
}


// Assume a and b are in ascending order
bool GTrieNode::_is_pairset_included(list<iPair> a, list<iPair> b) {
  list<iPair>::iterator aa, bb;

  aa = a.begin();
  bb = b.begin();

  while (aa!=a.end() && bb!=b.end()) {
    if ((aa->first == bb->first) && (aa->second == bb->second)) {
      aa++;
      bb++;
    } else {
      bb++;
    }      
  }

  if (aa == a.end()) return true;
  else return false;
}

// Insert Symmetry conditions pertaining to this node
void GTrieNode::insertConditionsFiltered(list<iPair> *conditions) {

  // Already has "empty set" of conditions
  if (cond_ok && cond_this_ok) return;

  list<iPair> aux;
  list<int>   aux_this_node;
  list<iPair>::iterator jj, kk;
  iPair p;


  /*for (jj=conditions->begin(); jj!=conditions->end(); jj++)
    printf("%d < %d  |  ", jj->first, jj->second);
    putchar('\n');*/
  

  // A bit slow but works (search a<b and b<c: remove a<c)
  for (jj=conditions->begin(); jj!=conditions->end(); jj++)
    for (kk=conditions->begin(); kk!=conditions->end(); kk++)
      if (jj->second == kk->first) {
	p.first  = jj->first;
	p.second = kk->second;
	conditions->remove(p); // does not invalidade jj and kk
      }

  /*for (jj=conditions->begin(); jj!=conditions->end(); jj++)
    printf("%d < %d  |  ", jj->first, jj->second);
    putchar('\n');*/

  for (jj=conditions->begin(); jj!=conditions->end(); jj++) {
    if (max(jj->first, jj->second)<=depth-1) // HERE!
      aux.push_back(*jj);
    if (jj->second == depth-1)
      aux_this_node.push_back(jj->first);
 } 

  // Deal with ancestor conditions
  if (!cond_ok) {

    if (aux.size()==0) {
      cond_ok = true;
      cond.clear();
    } else {

      list< list<iPair> >::iterator mm;
      bool is_contained = false;
      // A bit slow, but saves a lot of conditions
      for (mm=cond.begin(); mm!=cond.end();) {
	if (_is_pairset_included(*mm, aux)) {
	  is_contained = true;
	  break;
	}
	else if (_is_pairset_included(aux, *mm)) {
	  //printf("Erasing old condition (ancestors) !!\n");
	  mm = cond.erase(mm);
	  if (mm == cond.end()) break;
	}
	else mm++;
      }

      if (!is_contained) cond.push_back(aux);
      //else printf("avoid insertion (ancestors)!\n");
    }
  }

  return;

 // Deal with this node conditions {
  if (!cond_this_ok) {
    
    if (aux_this_node.size()==0) {
      cond_this_ok = true;
      this_node_cond.clear();
      } else {
    
      list< list<int> >::iterator mm;

      bool is_contained = false;
      // A bit slow, but saves a lot of conditions
      for (mm=this_node_cond.begin(); mm!=this_node_cond.end(); mm++) {
        if (_is_intset_included(*mm, aux_this_node)) {
            is_contained = true;
            break;
        }
        if (_is_intset_included(aux_this_node, *mm)) {
          //printf("Erasing old condition (this_node) !!\n");
          mm = this_node_cond.erase(mm);
       }
      }

      if (!is_contained) this_node_cond.push_back(aux_this_node);
      //else printf("avoid insertion (this_node)!\n");      
    }
  }

   
}

int GTrieNode::countNodes() {
  int aux=1;
  
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    aux += (*ii)->countNodes();

  return aux;
}

void GTrieNode::zeroFrequency() {
  frequency = 0;
  
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    (*ii)->zeroFrequency();
}

void GTrieNode::showFrequency() {
  if (is_graph) printf("%d \n", frequency);
  
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    (*ii)->showFrequency();
}

int GTrieNode::maxDepth() {
  int aux = 0;
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    aux=max(aux, 1+(*ii)->maxDepth());
  return aux;
}



void GTrieNode::goCondUndir() {  
  int i, j, ci, mylim, glaux;
  int ncand;
  int *p;

  mylim = INT_MAX;
  if (!cond_ok) {
    i = 1;
    list< list<iPair> >::const_iterator jj, jjend;
    list<iPair>::const_iterator kk, kkend;
    for (jj=cond.begin(), jjend=cond.end(); jj!=jjend; ++jj) {
      glaux = -1;
      for (kk=(*jj).begin(), kkend=(*jj).end(); kk!=kkend; ++kk)
	if (kk->second<glk && mymap[kk->first] > mymap[kk->second])
	  break;
	else if (kk->second==glk && mymap[kk->first]>glaux)
	  glaux = mymap[kk->first];
      if (kk==kkend) {
	i = 0;
	if (glaux < mylim) mylim=glaux;
      }
    }
    if (i) return;
  }
  if (mylim == INT_MAX) mylim = 0;
    
  ncand=0;
  j=ci=INT_MAX;
  for (i=0; i<nconn; i++) {
    glaux = numnei[mymap[conn[i]]];
    if (glaux<j) {
      ci = mymap[conn[i]];
      j = glaux;
    }
  }

  glaux = j;
  ncand = ci;
  for (p=&fastnei[ncand][j-1], ci= glaux-1; ci>=0; ci--, p--) {    
    i = *p;
    if (i<mylim) break;
    if (used[i]) continue;
    mymap[glk] = i;
    
    bool *b = &adjM[i][0];
    for (j=0; j<glk; j++)
      if (out[j] != *(b+mymap[j]))
	break;
    if (j<glk) continue;
    
    if (is_graph) {
      frequency++;
      if (Global::show_occ) {
	for (int k = 0; k<=glk; k++)
	  for (int l = 0; l<=glk; l++)
	    fputc(adjM[mymap[k]][mymap[l]]?'1':'0', Global::occ_file);
	fputc(':', Global::occ_file);
	for (int k = 0; k<=glk; k++)
	  fprintf(Global::occ_file, " %d", mymap[k]+1);
	fputc('\n', Global::occ_file);
      }
    }

    used[i]=true;
    glk++;
    list<GTrieNode *>::const_iterator ii, iiend;
    for(ii=child.begin(), iiend = child.end(); ii!=iiend; ++ii)
      (*ii)->goCondUndir();
    glk--;
    used[i]=false;
  }
}


void GTrieNode::goCondDir() {  
  int i, j, ci, mylim, glaux;
  int ncand;
  int *p;

  mylim = INT_MAX;
  if (!cond_ok) {
    i = 1;
    list< list<iPair> >::const_iterator jj, jjend;
    list<iPair>::const_iterator kk, kkend;
    for (jj=cond.begin(), jjend=cond.end(); jj!=jjend; ++jj) {
      glaux = -1;
      for (kk=(*jj).begin(), kkend=(*jj).end(); kk!=kkend; ++kk)
	if (kk->second<glk && mymap[kk->first] > mymap[kk->second])
	  break;
	else if (kk->second==glk && mymap[kk->first]>glaux)
	  glaux = mymap[kk->first];
      if (kk==kkend) {
	i = 0;
	if (glaux < mylim) mylim=glaux;
      }
    }
    if (i) return;
  }
  if (mylim == INT_MAX) mylim = 0;
    
  ncand=0;
  j=ci=INT_MAX;
  for (i=0; i<nconn; i++) {
    glaux = numnei[mymap[conn[i]]];
    if (glaux<j) {
      ci = mymap[conn[i]];
      j = glaux;
    }
  }

  glaux = j;
  ncand = ci;
  for (p=&fastnei[ncand][j-1], ci= glaux-1; ci>=0; ci--, p--) {    
    i = *p;
    if (i<mylim) break;
    if (used[i]) continue;
    mymap[glk] = i;

    for (j=0; j<glk; j++)
      if (in[j] != adjM[mymap[j]][i])
	break;
    if (j<glk) continue;
    bool *b = &adjM[i][0];
    for (j=0; j<glk; j++)
      if (out[j] != *(b+mymap[j]))
	break;
    if (j<glk) continue;
    
    if (is_graph) {
      frequency++;
      if (Global::show_occ) {
	for (int k = 0; k<=glk; k++)
	  for (int l = 0; l<=glk; l++)
	    fputc(adjM[mymap[k]][mymap[l]]?'1':'0', Global::occ_file);
	fputc(':', Global::occ_file);
	for (int k = 0; k<=glk; k++)
	  fprintf(Global::occ_file, " %d", mymap[k]+1);
	fputc('\n', Global::occ_file);
      }
    }

    used[i]=true;
    glk++;
    list<GTrieNode *>::const_iterator ii, iiend;
    for(ii=child.begin(), iiend = child.end(); ii!=iiend; ++ii)
      (*ii)->goCondDir();
    glk--;
    used[i]=false;
  }
}



int GTrieNode::countGraphsApp() {
  int aux=0;
  if (is_graph && frequency>0) aux++;
  
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    aux += (*ii)->countGraphsApp();

  return aux;
}

int GTrieNode::countGraphs() {
  int aux=0;
  if (is_graph) aux=1;
  
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    aux += (*ii)->countGraphs();

  return aux;
}

int GTrieNode::countGraphPaths() {
  int aux=0;
  if (is_graph) aux=depth;
  
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    aux += (*ii)->countGraphPaths();

  return aux;
}

double GTrieNode::countOccurrences() {
  double aux=0;
  if (is_graph && frequency>0) aux+=frequency;
  
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    aux += (*ii)->countOccurrences();

  return aux;
}

void GTrieNode::writeToFile(FILE *f) {
  int i, bits, aux, nbytes;

  int nchild       = child.size();
  nbytes=0;
  while (nchild>0) {
    nbytes++;
    nchild/=BASE_FORMAT;
  }

  int first_byte   = nbytes<<1;

  // first bit: isgraph ?
  if (is_graph) BIT_SET(first_byte,0);
  else          BIT_CLEAR(first_byte,0);
  fputc(BASE_FIRST+first_byte,f);

  nchild       = child.size();
  while (nchild>0) {
    fputc(BASE_FIRST+(nchild%BASE_FORMAT),f);
    nchild/=BASE_FORMAT;
  }

  // Connections: outgoing
  for (i=0, bits=0, aux=0; i<depth; i++, bits++) {
    if (bits==BASE_BITS) {
      fputc(BASE_FIRST+aux, f);
      aux=0;
      bits=0;
    }
    if (out[i]) BIT_SET(aux, bits);
  }
  fputc(BASE_FIRST+aux, f);


  // Connections: ingoing
  for (i=0, bits=0, aux=0; i<depth; i++, bits++) {
    if (bits==BASE_BITS) {
      fputc(BASE_FIRST+aux, f);
      aux=0;
      bits=0;
    }
    if (in[i]) BIT_SET(aux, bits);
  }
  fputc(BASE_FIRST+aux, f);

  // Previous Conditions
  if (cond_ok) aux = 0;
  else         aux = cond.size();
  fputc(BASE_FIRST+aux, f);

  if (aux>0) {
    list< list<iPair> >::const_iterator jj, jjend;
    list<iPair>::const_iterator kk, kkend;
    for (jj=cond.begin(), jjend=cond.end(); jj!=jjend; ++jj) {      
      for (kk=(*jj).begin(), kkend=(*jj).end(); kk!=kkend; ++kk) {
	fputc(BASE_FIRST+1+(kk->first),f);
	fputc(BASE_FIRST+1+(kk->second),f);
      }
      fputc(BASE_FIRST, f);
    }
  }
  
  fputc('\n', f);

  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    (*ii)->writeToFile(f);

}

void GTrieNode::readFromFile(FILE *f) {
  int nchilds, i, j, pos, bits, ncond;
  iPair p;
  char aux, buf[MAX_BUF];
  GTrieNode *c;

  if (fgets(buf, MAX_BUF, f)) {
    aux = buf[0]-BASE_FIRST;
    if (BIT_VALUE(aux, 0)) is_graph=true;
    else                   is_graph=false;
    int nbytes = aux>>1;

    pos = 1;

    for (i=0, j=1, nchilds=0; i<nbytes; i++, j*=BASE_FORMAT) {
      aux = buf[pos++]-BASE_FIRST;
      nchilds += int(aux)*j;
    }

    // Connections: outgoing
    aux = buf[pos++]-BASE_FIRST;
    for (i=0, bits=0; i<depth; i++, bits++) {
      if (bits==BASE_BITS) {
	aux = buf[pos++]-BASE_FIRST;
	bits=0;
      }
      if (BIT_VALUE(aux, bits)) out[i] = true;
      else                      out[i] = false;
    }

    // Connections: ingoing
    aux = buf[pos++]-BASE_FIRST;
    for (i=0, bits=0; i<depth; i++, bits++) {
      if (bits==BASE_BITS) {
	aux = buf[pos++]-BASE_FIRST;
	bits=0;
      }
      if (BIT_VALUE(aux, bits)) in[i] = true;
      else                      in[i] = false;
    }

    // Previous Conditions
    aux = buf[pos++]-BASE_FIRST;
    if (aux==0) cond_ok = true;
    else        cond_ok = false;

    if (aux>0) {
      ncond = aux;
      for (i=0; i<ncond; i++) {
	list<iPair> newcond;
	while(1) {
	  aux = buf[pos++]-BASE_FIRST-1;
	  if (aux<0) break;
	  p.first = aux;
	  aux = buf[pos++]-BASE_FIRST-1;
	  p.second = aux;
	  newcond.push_back(p);
	}
	cond.push_back(newcond);
      }
    }
    
    if (buf[pos]!='\n') {
      fprintf(stderr, "ERROR: [%s] !%d!%c!\n", buf, pos, buf[pos]);
      fprintf(stderr,"[%d](%s) |", depth, is_graph?"X":" ");
      for (i=0, bits=0; i<depth; i++, bits++)
	fprintf(stderr,"%s", out[i]?"1":"0");
      fprintf(stderr,"|\n");
    }

    // Conn and nconn variables (was missing)
    for (i=0; i<depth; i++)
      if (out[i] || in[i]) 
      conn[(nconn)++] = i;
    
    
    for (i=0; i<nchilds; i++) {
      c = new GTrieNode(depth+1);
      c->readFromFile(f);
      child.push_back(c);
    }

  }

}

void GTrieNode::populateGraphTree(GraphTree *tree, char *s, int size) {
  int i, pos=depth-1;
  
  for (i=0;i<depth;i++) {
    s[pos*size+i]=out[i]?'1':'0';
    s[i*size+pos]=in[i]?'1':'0';
  }

  if (is_graph)
    tree->setString(s, frequency);

  list<GTrieNode *>::const_iterator ii, iiend;
  for(ii=child.begin(), iiend = child.end(); ii!=iiend; ++ii)
    (*ii)->populateGraphTree(tree, s, size);
}


void GTrieNode::populateMap(mapStringInt *m, char *s, int size) {
  int i, pos=depth-1;

  for (i=0;i<depth;i++) {
    s[pos*size+i]=out[i]?'1':'0';
    s[i*size+pos]=in[i]?'1':'0';
  }

  if (is_graph && frequency>0) (*m)[s]=frequency;

  list<GTrieNode *>::const_iterator ii, iiend;
  for(ii=child.begin(), iiend = child.end(); ii!=iiend; ++ii)
    (*ii)->populateMap(m, s, size);
}

// -------------------------------------

GTrie::GTrie() {
  _root = new GTrieNode(0);
  _root->cond_ok = _root->cond_this_ok = true;
}

GTrie::~GTrie() {
  delete _root;
}

void GTrie::insertGraphCond(Graph *g, list<iPair> *cond) {
  _root->insertCond(g, cond);
}

void GTrie::insertGraph(Graph *g) {
  _root->insert(g);
}

void GTrie::showAsText(FILE *f) {
  _root->showAsText(f);
}

int GTrie::frequencyGraphString(int size, const char *s) {
  char larger[size*size+1];
  Isomorphism::canonicalBasedNauty(s, larger, size);

  Graph *g = new GraphMatrix();
  GraphUtils::strToGraph(g, larger, size, DIRECTED); // May change later
  int aux = _root->frequencyGraph(g);
  delete g;
  
  return aux;
}


void GTrie::insertGraphString(int size, const char *s) {
  char larger[size*size+1];
  Isomorphism::canonicalBigger(s, larger, size);


  Graph *g = new GraphMatrix();
  GraphUtils::strToGraph(g, larger, size, DIRECTED); // May change later
  list<iPair> *cond = new list<iPair>;
  Conditions::symmetryConditions(g, cond);

  insertGraphCond(g, cond);

  delete cond;
  delete g;
}

// Populate g-trie with subgraphs of 'size' read from file 's'
void GTrie::readSubgraphs(int size, bool dir, char *s) {
  char buf[MAX_BUF];

  FILE *f = fopen(s, "r");
  if (!f) Error::msg(NULL);
  while (fscanf(f, "%s", buf)==1) {
    insertGraphNautyString(size, buf, dir, 1);
  }
  fclose(f);
  cleanConditions();  
}

void GTrie::insertGraphNautyString(int size, const char *s, bool dir, int label) {
  char larger[size*size+1];

  if (label==2)
    Isomorphism::canonicalBigger(s, larger, size);
  else if (label==3)
    Isomorphism::canonicalRandom(s, larger, size);
  else if (label==4)
    Isomorphism::canonicalBasedNautyReverse(s, larger, size);
  else
    Isomorphism::canonicalBasedNauty(s, larger, size);

  Graph *g = new GraphMatrix();

  if (dir)  GraphUtils::strToGraph(g, larger, size, DIRECTED);
  else      GraphUtils::strToGraph(g, larger, size, UNDIRECTED);

  g->makeArrayNeighbours();

  list<iPair> *cond = new list<iPair>;
  Conditions::symmetryConditions(g, cond);  

  insertGraphCond(g, cond);

  delete cond;
  delete g;
}

double GTrie::compressionRate() {
  int nodes  = _root->countNodes()-1;
  int graphs = _root->countGraphPaths();
  return 1 - (double)(nodes) / (double)(graphs);
}

int GTrie::maxDepth() {
  return _root->maxDepth();
}

void GTrie::showFrequency() {
  return _root->showFrequency();
}

void GTrie::census(Graph *g) {
  int i;
  int subgraph_size = maxDepth();
  int num_nodes = g->numNodes();

  _root->zeroFrequency();
  
  GTrieNode::mymap = new int[subgraph_size];
  GTrieNode::used  = new bool[num_nodes];
  GTrieNode::numNodes = num_nodes;
  GTrieNode::fastnei  = g->matrixNeighbours();
  GTrieNode::adjM     = g->adjacencyMatrix();
  GTrieNode::numnei   = g->arrayNumNeighbours(); 
  
  if (g->type() == DIRECTED) GTrieNode::isdir = true;
  else                       GTrieNode::isdir = false;
  for (i=0; i<num_nodes; i++)
    GTrieNode::used[i]=false;

  GTrieNode *c = *(_root->child.begin());
  list<GTrieNode *>::iterator ii;
  GTrieNode::glk=1;

  for (i = 0; i<num_nodes; i++) {
    GTrieNode::mymap[0] = i;
    GTrieNode::used[i]=true;
    if (g->type() == DIRECTED)
      for(ii=c->child.begin(); ii!=c->child.end(); ii++)
	(*ii)->goCondDir();
    else
      for(ii=c->child.begin(); ii!=c->child.end(); ii++)
	(*ii)->goCondUndir();
    GTrieNode::used[i]=false;
  }
        
  delete [] GTrieNode::mymap;
  delete [] GTrieNode::used;
}

double GTrie::countOccurrences() {
  return _root->countOccurrences();
}

int GTrie::countGraphsApp() {
  return _root->countGraphsApp();
}

int GTrie::countGraphs() {
  return _root->countGraphs();
}

void GTrie::writeToFile(char *s) {
  FILE *f;
  
  f=fopen(s,"w");
  if (f!=NULL) {
    fputs("GTRIEFORMAT VERSION 1\n", f);
    _root->writeToFile(f);
    fclose(f);
  } else
    Error::msg("Unable to open g-trie output file \"%s\")", s);    
}

void GTrie::readFromFile(char *s) {
  FILE *f;
  char *dummy, buf[MAX_BUF];
  
  f=fopen(s,"r");
  if (!f) Error::msg(NULL);
  dummy=fgets(buf, MAX_BUF, f);
  _root->readFromFile(f);
  fclose(f);
}

void GTrie::populateGraphTree(GraphTree *tree, int size) {
  char s[size*size+1];
  s[size*size]=0;
  _root->populateGraphTree(tree, s, size);
}

void GTrie::populateMap(mapStringInt *m, int size) {
  char s[size*size+1];
  s[size*size]=0;
  _root->populateMap(m, s, size);
}


void GTrie::censusSample(Graph *g, double *p) {
  int i;
  int subgraph_size = maxDepth();
  int num_nodes = g->numNodes();

  _root->zeroFrequency();
  
  GTrieNode::mymap = new int[subgraph_size];
  GTrieNode::used  = new bool[num_nodes];
  GTrieNode::numNodes = num_nodes;
  GTrieNode::fastnei  = g->matrixNeighbours();
  GTrieNode::adjM     = g->adjacencyMatrix();
  GTrieNode::numnei   = g->arrayNumNeighbours(); 
  GTrieNode::prob     = p;

  if (g->type() == DIRECTED) GTrieNode::isdir = true;
  else                       GTrieNode::isdir = false;
  for (i=0; i<num_nodes; i++)
    GTrieNode::used[i]=false;


  list<GTrieNode *>::iterator ii, iiend;
  GTrieNode::glk=0;
  for(ii=_root->child.begin(), iiend = _root->child.end(); ii!=iiend; ii++)
    if (Random::getDouble()<=GTrieNode::prob[0]) {
      (*ii)->goCondSample();
    }
}


void GTrieNode::goCondSample() {
  int i, j, ci, mylim, glaux;
  int ncand;

  if (!cond_ok) {
    list< list<iPair> >::const_iterator jj, jjend;
    list<iPair>::const_iterator kk, kkend;
    for (jj=cond.begin(), jjend=cond.end(); jj!=jjend; ++jj) {
      for (kk=(*jj).begin(), kkend=(*jj).end(); kk!=kkend; ++kk)
	if ( mymap[kk->first]>mymap[kk->second])
	  break;
      if (kk==kkend) break;
    }
    if (jj==jjend) return;
  }

  ncand=0;
  j=ci=INT_MAX;
  if (nconn == 0) ncand = numNodes;
  else {
    for (i=0; i<nconn; i++) {
      glaux = numnei[mymap[conn[i]]];
      if (glaux<j) {
	ci = mymap[conn[i]];
	j = glaux;
      }
    }
    ncand = j;
  }
  int cand[ncand];

  ncand=0;
  if (nconn == 0) { // We are at a node with no connections to ancestors
    for (i=numNodes-1; i>=0; i--)
      if (!used[i]) 
	cand[ncand++]=i;
  } else {
    for (i=0; i<j; i++) {
      glaux = fastnei[ci][i];
      if (!used[glaux])
	cand[ncand++]=glaux;
    }
  }

  if (cond_this_ok) mylim = 0;
  else {
    list< list<int> >::const_iterator jjj;
    list<int>::const_iterator kkk;
    mylim = INT_MAX;
    for (jjj=this_node_cond.begin(); jjj!=this_node_cond.end(); ++jjj) {
      glaux = -1;
      for (kkk=(*jjj).begin(); kkk!=(*jjj).end(); ++kkk)
	if (mymap[*kkk]>glaux) glaux = mymap[*kkk];
      if (glaux<mylim) mylim=glaux;
    }
  }

  for (ci=ncand-1; ci>=0; ci--) {
    i = cand[ci];
    if (i<mylim) break;
    mymap[glk] = i;

    if (isdir) {  
      for (j=0; j<glk; j++)
	if (in[j]  != adjM[mymap[j]][i] ||
	    out[j] != adjM[i][mymap[j]])
	  break;
    } else {
      for (j=0; j<glk; j++)
	if (in[j]  != adjM[mymap[j]][i])
	  break;
    }
    if (j<glk) continue;

    if (is_graph) {      
      if (Global::show_occ) {
	for (int k = 0; k<=glk; k++)
	  for (int l = 0; l<=glk; l++)
	    fputc(adjM[mymap[k]][mymap[l]]?'1':'0', Global::occ_file);
	fputc(':', Global::occ_file);
	for (int k = 0; k<=glk; k++)
	  fprintf(Global::occ_file, " %d", mymap[k]+1);
	fputc('\n', Global::occ_file);
      }
      frequency++;
    }

    used[i]=true;
    glk++;
    list<GTrieNode *>::const_iterator ii, iiend;
    for(ii=child.begin(), iiend = child.end(); ii!=iiend; ++ii)
      if (Random::getDouble() <= prob[glk]) {
	(*ii)->goCondSample();
      }
    glk--;
    used[i]=false;
  }
}

void GTrieNode::clean(int a, int b) {
  list< list<iPair> >::iterator ii;
  list<iPair> ::iterator jj;

  if (cond.size()>0) {

    for (ii=cond.begin(); ii!=cond.end(); ii++) {
      for (jj=ii->begin(); jj!=ii->end();)
	  if (jj->first==a && jj->second==b)
	    jj = ii->erase(jj);
	  else
	    jj++;
      if (ii->size()==0) {
	ii->clear();
	cond_ok = true;
      }
    }    
  }

  list<GTrieNode *>::iterator cc;
  for(cc=child.begin(); cc!=child.end(); cc++)
    (*cc)->clean(a, b);	

}

void GTrieNode::cleanConditions() {
  int a, b;
  list< list<iPair> >::iterator ii;
  list<iPair> ::iterator jj, kk;

  if (cond.size()>0) {
    for (jj=cond.begin()->begin(); jj!=cond.begin()->end();jj++) {
      a=jj->first;
      b=jj->second;
      for (ii=cond.begin(); ii!=cond.end(); ii++) {
	for (kk=ii->begin(); kk!=ii->end(); kk++)
	  if (kk->first==a && kk->second==b) break;
	if (kk==ii->end()) break;
      }
      if (ii==cond.end())  {
	list<GTrieNode *>::iterator ii;
	for(ii=child.begin(); ii!=child.end(); ii++)
	  (*ii)->clean(a, b);	
      }
    }
  }

  list<GTrieNode *>::iterator cc;
  for(cc=child.begin(); cc!=child.end(); cc++)
    (*cc)->cleanConditions();
}


void GTrie::cleanConditions() {
  _root->cleanConditions();
}

// -------------------------------------

