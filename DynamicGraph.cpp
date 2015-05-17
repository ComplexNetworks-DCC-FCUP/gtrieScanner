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
Graphs Implementation with Adj. Matrix and Adj. List

Last Update: 11/02/2012
---------------------------------------------------- */

#include "DynamicGraph.h"
#include "GraphUtils.h"
#include <stdio.h>
#include <algorithm>

#define TRIE_MOD 15
#define TRIE_ORD 4
#define HYBRID_NUMBER 22

#ifdef PRINT_CALLS
FILE* fn;
#endif

struct DynamicGraph::l_list {
  int value;
  l_list* next;
};

struct DynamicGraph::a_trie {
  bool end;
  a_trie* childs[TRIE_MOD + 1];
};

DynamicGraph::a_trie* DynamicGraph::new_trie() {
  a_trie* tmp = new a_trie();
  int i;
  
  tmp->end = false;
  for (i = 0; i <= TRIE_MOD; i++) {
    tmp->childs[i] = NULL;
  }

  return tmp;
}

void DynamicGraph::delete_trie(a_trie* cur) {
  if (cur == NULL)
    return;

  int i;
  for (i = 0; i <= TRIE_MOD; i++) {
    delete_trie(cur->childs[i]);
    if (cur->childs[i] != NULL)
      delete cur->childs[i];
  }
}

DynamicGraph::DynamicGraph(RepType _r) {
  #ifdef PRINT_CALLS
  fn = fopen("hasEdgeCalls.txt", "w");
  #endif

  _cstatus = false;
  _rtype = _r;
  _init();
}

DynamicGraph::DynamicGraph(RepType _r, bool _cs) {
  #ifdef PRINT_CALLS
  fn = fopen("hasEdgeCalls.txt", "w");
  #endif

  _cstatus = _cs;
  _rtype = _r;
  _init();
}

DynamicGraph::~DynamicGraph() {
  _delete();
}

// ------------------------------
// Graph Creation
// ------------------------------

void DynamicGraph::_init() {
  _num_nodes = _num_edges = 0;

  ready = false;
  _adjM             = NULL;
  _adjOut           = NULL;
  _adjIn            = NULL;
  _neighbours       = NULL;
  cache             = NULL;
  trie              = NULL;
  _hashM            = NULL;
  _in               = NULL;
  _out              = NULL;
  _maxL             = NULL;
  _minL             = NULL;
  _num_neighbours   = NULL;
  _array_neighbours = NULL;
}

void DynamicGraph::_delete() {
  int i;

  if (_adjM != NULL) {
    for (i = 0; i < _num_nodes; i++)
      if (_adjM[i] != NULL)
        delete[] _adjM[i];
    delete[] _adjM;
  }
  
  if (_adjIn!=NULL) delete[] _adjIn;
  if (_adjOut!=NULL) delete[] _adjOut;
  if (_neighbours!=NULL) delete[] _neighbours;

  if (_hashM != NULL) {
    int j;
    for (i = 0; i < _num_nodes; i++) {
      for (j = 0; j < _sqrt_nodes; j++) {
        l_list* cur = _hashM[i][j], *prev;
        while (cur != NULL) {
          prev = cur;
          cur = cur->next;
          delete prev;
        }
      }

      delete[] _hashM[i];
    }
    delete[] _hashM;
  }
    
  if (cache != NULL) {
    for (i = 0; i < _num_nodes; i++)
      delete[] cache[i];
    delete[] cache;
  }

  if (_in != NULL) delete[] _in;
  if (_out != NULL) delete[] _out;
  if (_out != NULL) delete[] _num_neighbours;
  if (_maxL != NULL) delete[] _maxL;
  if (_minL != NULL) delete[] _minL;
  if (trie != NULL) {
    for (i = 0; i < _num_nodes; i++) {
      delete_trie(trie[i]);
      delete trie[i];
    }
    delete[] trie;
  }

  if (_array_neighbours != NULL) {
    for (i = 0; i < _num_nodes; i++)
      if (_array_neighbours[i] != NULL)
        delete[] _array_neighbours[i];
    delete[] _array_neighbours;
  }

  ready = false;
}

void DynamicGraph::zero() {
  int i,j;
  _num_edges = 0;
  ready = false;

  for (i = 0; i < _num_nodes; i++) {
    _in[i] = 0;
    _out[i] = 0;
    _num_neighbours[i] = 0;

    if (_cstatus)
      for (j = 0; j <= _log_nodes; j++)
        cache[i][j] = -1;
    _adjIn[i].clear();
    _adjOut[i].clear();
    _neighbours[i].clear();
    if (_rtype == MATRIX)
      for (j=0; j<_num_nodes;j++)
        _adjM[i][j]=false;
  }
}

void DynamicGraph::createGraph(int n, GraphType t) {
  _delete();
  
  int i;
  _num_nodes = n;
  _type = t;
  
  int tmp_log = (int)log2((double) n);
  _log_nodes = 1;
  while (_log_nodes < tmp_log)
    _log_nodes *= 2;
  _log_nodes--;

  tmp_log = (int)sqrt(n) + 1;
  _sqrt_nodes = 1;
  while (_sqrt_nodes < tmp_log)
    _sqrt_nodes *= 2;
  _sqrt_nodes--;

  if (_rtype == MATRIX) {
    _adjM = new bool*[n];  
    for (i = 0; i < n; i++)
      _adjM[i] = new bool[n];
  }
  
  _adjIn      = new vector<int>[n];
  _adjOut     = new vector<int>[n];
  _neighbours = new vector<int>[n];

  if (_cstatus) {
    cache           = new int*[n];
    for (i = 0; i < n; i++)
      cache[i] = new int[_log_nodes + 1];
  }

  _in             = new int[n];
  _out            = new int[n];
  _num_neighbours = new int[n];

  zero();
}

void DynamicGraph::prepareGraph() {
  ready = true;

  _maxL = new int[_num_nodes];
  _minL = new int[_num_nodes];
  int i, j;
  for (i = 0; i < _num_nodes; i++) {
    sort(_adjOut[i].begin(), _adjOut[i].end());
    sort(_adjIn[i].begin(), _adjIn[i].end());

    _maxL[i] = _minL[i] = -1;

    if (_out[i])
    {
      _maxL[i] = _adjOut[i][_out[i] - 1];
      _minL[i] = _adjOut[i][0];
    }
  }
  
  if (_rtype == MATRIX || _rtype == LINEAR)
    return;
  if (_rtype == BSLIST || _rtype == INTER)
    return;
  if (_rtype == HASH || _rtype == HYBRID) {
    _hashM = new l_list**[_num_nodes];

    int i, j;
    for (i = 0; i < _num_nodes; i++) {
      _hashM[i] = new l_list*[_sqrt_nodes + 1];

      for (j = 0; j <= _sqrt_nodes; j++)
        _hashM[i][j] = NULL;
    }

    for (i = 0; i < _num_nodes; i++) {
      for (j = 0; j < _out[i]; j++) {
        l_list* n_node = new l_list();
        n_node->value = _adjOut[i][j];
        n_node->next = _hashM[i][_adjOut[i][j] & _sqrt_nodes];
        _hashM[i][_adjOut[i][j] & _sqrt_nodes] = n_node;
      }
    }
  }
  if (_rtype == TRIE || _rtype == HYBRID) {
    int i, j;
    trie = new a_trie*[_num_nodes];
    
    for (i = 0; i < _num_nodes; i++) {
      trie[i] = new_trie();

      for (j = 0; j < _out[i]; j++) {
        int num = _adjOut[i][j];
        a_trie* cur = trie[i];
        
        while (num) {
          if (cur->childs[num & TRIE_MOD] == NULL)
            cur->childs[num & TRIE_MOD] = new_trie();
          cur = cur->childs[num & TRIE_MOD];
          num >>= TRIE_ORD;
        }
        cur->end = true;
      }
    }
  }
}

void DynamicGraph::addEdge(int a, int b) {
  _adjOut[a].push_back(b);
  _out[a]++;

  _adjIn[b].push_back(a);
  _in[b]++;

  _num_edges++;

  if (_rtype == MATRIX)
    _adjM[a][b] = true;

  if (!hasEdge(b, a)) {
    _neighbours[b].push_back(a);
    _num_neighbours[b]++;

    _neighbours[a].push_back(b);
    _num_neighbours[a]++;
  }
}

void DynamicGraph::rmEdge(int a, int b) {
  fprintf(stderr, "rmEdge is deprecated for the time being\n");
}

bool DynamicGraph::hasEdge(int a, int b) {
  #ifdef PRINT_CALLS
  fprintf(fn, "%d %d\n", a, b);
  #endif

  if (!ready) {
    int i;
    for (i = 0; i < _out[a]; i++)
      if (_adjOut[a][i] == b)
        return true;
    return false;
  }

  if (_cstatus && cache[a][b & _log_nodes] == b)
    return true;

//  if (b < _minL[a] || b > _maxL[a])
//      return false;

  if (_rtype == MATRIX)
    return _adjM[a][b];
  else if (_rtype == LINEAR) {
    int i;
    for (i = 0; i < _out[a]; i++)
      if (_adjOut[a][i] == b)
        return true;

    return false;
  }
  else if (_rtype == BSLIST) {
    int lo = 0, hi = _out[a] - 1, med, vl;
    if (hi < 0)
      return false;
    
    while (lo <= hi) {
      med = (lo + hi) >> 1;
      vl = _adjOut[a][med];
      
      if (vl < b)
        lo = med + 1;
      else if (vl > b)
        hi = med - 1;
      else if (vl == b)
        return true;
    }

    return false;
  }
  else if (_rtype == HASH || (_rtype == HYBRID && _out[a] <= HYBRID_NUMBER)) {
    l_list* cur = _hashM[a][b & _sqrt_nodes];
    while (cur != NULL) {
      if (cur->value == b) {
        if (_cstatus)
          cache[a][b & _log_nodes] = b;
        return true;
      }

      cur = cur->next;
    }
    return false;
  }
  else if (_rtype == TRIE || _rtype == HYBRID) {
    a_trie* cur = trie[a];

    while (cur != NULL && b) {
      cur = cur->childs[b & TRIE_MOD];
      b >>= TRIE_ORD;
    }

    if (_cstatus && cur != NULL && cur->end)
      cache[a][b & _log_nodes] = b;

    return cur != NULL && cur->end;
  }
  else if (_rtype == INTER) {
    int lo = 0, hi = _out[a] - 1, med, vl;
    if (hi < 0 || b < _adjOut[a][0] || b > _adjOut[a][hi])
      return false;
    
    while (_adjOut[a][lo] <= b && _adjOut[a][hi] >= b) {
      med = (lo == hi) ? lo : lo + ((b - _adjOut[a][lo]) * (hi - lo)) / (_adjOut[a][hi] - _adjOut[a][lo]);
      vl = _adjOut[a][med];
      
      if (vl < b)
        lo = med + 1;
      else if (vl > b)
        hi = med - 1;
      else if (vl == b)
        return true;
    }

    return false;
  }

  return false;
}

void DynamicGraph::sortNeighbours() {
  int i;
  for (i = 0; i < _num_nodes; i++)
    sort(_neighbours[i].begin(), _neighbours[i].begin()+_neighbours[i].size());
}

void DynamicGraph::sortNeighboursArray() {
  int i;
  for (i = 0; i < _num_nodes; i++)
    qsort(_array_neighbours[i], _num_neighbours[i], sizeof(int), GraphUtils::int_compare);
}

void DynamicGraph::makeArrayNeighbours() {
  int i,j;
  vector<int>:: iterator ii;
  _array_neighbours = new int*[_num_nodes];  
  for (i = 0; i < _num_nodes; i++) {
    _array_neighbours[i] = new int[_neighbours[i].size()];
    for (ii = _neighbours[i].begin(), j = 0; ii != _neighbours[i].end(); ii++, j++)
      _array_neighbours[i][j] = *ii;
    _neighbours[i].clear();
  }
}

void DynamicGraph::makeVectorNeighbours() {
  int i,j;
  vector<int>:: iterator ii;

  for (i = 0; i < _num_nodes; i++)
    for (j = 0; j < _num_neighbours[i]; j++)
      _neighbours[i].push_back(_array_neighbours[i][j]);

  if (_array_neighbours!=NULL) {
    for (i = 0; i < _num_nodes; i++)
      if (_array_neighbours[i] != NULL)
        delete[] _array_neighbours[i];
    delete[] _array_neighbours;
  }
}
