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

struct DynamicGraph::l_list {
  int value;
  l_list* next;
};

DynamicGraph::DynamicGraph() {
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
  hybrid_ch         = NULL;
  _hashM            = NULL;
  _in               = NULL;
  _out              = NULL;
  _hash_out         = NULL;
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
      if (_hash_out != NULL) {
        for (j = 0; j < _hash_out[i]; j++) {
          l_list* cur = _hashM[i][j], *prev;
          while (cur != NULL) {
            prev = cur;
            cur = cur->next;
            delete prev;
          }
        }
      }

      delete[] _hashM[i];
    }
    delete[] _hashM;
  }
    
  if (_in != NULL) delete[] _in;
  if (_out != NULL) delete[] _out;
  if (_hash_out != NULL) delete[] _hash_out;
  if (_num_neighbours != NULL) delete[] _num_neighbours;

  if (hybrid_ch != NULL) delete[] hybrid_ch;

  if (_array_neighbours != NULL) {
    for (i = 0; i < _num_nodes; i++)
      if (_array_neighbours[i] != NULL)
        delete[] _array_neighbours[i];
    delete[] _array_neighbours;
  }

  ready = false;
}

void DynamicGraph::_deleteAux() {
  int i;

  if (_adjM != NULL) {
    for (i = 0; i < _num_nodes; i++)
      if (_adjM[i] != NULL)
        delete[] _adjM[i];
    delete[] _adjM;
  }

  if (_hashM != NULL) {
    int j;
    for (i = 0; i < _num_nodes; i++) {
      if (_hash_out != NULL) {
        for (j = 0; j < _hash_out[i]; j++) {
          l_list* cur = _hashM[i][j], *prev;
          while (cur != NULL) {
            prev = cur;
            cur = cur->next;
            delete prev;
          }
        }
      }
   
      delete[] _hashM[i];
    }
    delete[] _hashM;
  }

  if (hybrid_ch != NULL) delete[] hybrid_ch;

  ready = false;
}

void DynamicGraph::zero() {
  int i,j;
  _num_edges = 0;
  ready = false;

  for (i = 0; i < _num_nodes; i++) {
    _in[i] = 0;
    _out[i] = 0;
    _hash_out[i] = 0;
    _num_neighbours[i] = 0;

    _adjIn[i].clear();
    _adjOut[i].clear();
    _neighbours[i].clear();
  }
}

void DynamicGraph::createGraph(int n, GraphType t) {
  _delete();
  
  int i;
  _num_nodes = n;
  _type = t;
  
  int tmp_log = 4 * (int)log2((double) n);
  _log_nodes = 1;
  while (_log_nodes < tmp_log)
    _log_nodes *= 2;
  _log_nodes--;
  
  _adjIn      = new vector<int>[n];
  _adjOut     = new vector<int>[n];
  _neighbours = new vector<int>[n];

  _in             = new int[n];
  _out            = new int[n];
  _num_neighbours = new int[n];
  _hash_out       = new int[n];

  zero();
}

void DynamicGraph::prepareGraph() {
  _deleteAux();
  ready = true;

  hybrid_ch = new int[_num_nodes];

  int i, j;
  long long int total = 0;
  for (i = 0; i < _num_nodes; i++) {
    sort(_adjOut[i].begin(), _adjOut[i].end());
    sort(_adjIn[i].begin(), _adjIn[i].end());

    total += _out[i];
  }

  int tmp_log = 2 * total / _num_nodes;
  _sqrt_nodes = 1;
  while (_sqrt_nodes < tmp_log)
    _sqrt_nodes *= 2;
  _sqrt_nodes--;

  _adjM = new bool*[_num_nodes];

  vector<pair<int, int> > vs;
  for (i = 0; i < _num_nodes; i++) {
    hybrid_ch[i] = 1;
    vs.push_back(pair<int, int>(_out[i], i));
  }

  sort(vs.begin(), vs.end());
  reverse(vs.begin(), vs.end());

  for (i = 0; i < min(5 * _sqrt_nodes / 2, _num_nodes); i++) {
    int cur_node = vs[i].second;
    hybrid_ch[cur_node] = 0;
    _adjM[cur_node] = new bool[_num_nodes];
    for (j = 0; j < _num_nodes; j++)
      _adjM[cur_node][j] = false;

    for (j = 0; j < _out[cur_node]; j++)
      _adjM[cur_node][_adjOut[cur_node][j]] = true;
  }

  _hashM = new l_list**[_num_nodes];

  for (i = 0; i < _num_nodes; i++) {
    if (!hybrid_ch[i])
      continue;

    int tmp_log = 5 * _out[i] / 2;
    _hash_out[i] = 1;
    while (_hash_out[i] < tmp_log)
      _hash_out[i] *= 2;
    _hash_out[i]--;

    _hashM[i] = new l_list*[_hash_out[i] + 1];

    for (j = 0; j <= _hash_out[i]; j++)
      _hashM[i][j] = NULL;
  }

  for (i = 0; i < _num_nodes; i++) {
    if (!hybrid_ch[i])
      continue;

    for (j = 0; j < _out[i]; j++) {
      l_list* n_node = new l_list();
      n_node->value = _adjOut[i][j];
      n_node->next = _hashM[i][_adjOut[i][j] & _hash_out[i]];
      _hashM[i][_adjOut[i][j] & _hash_out[i]] = n_node;
    }
  }
}

void DynamicGraph::addEdge(int a, int b) {
  _adjOut[a].push_back(b);
  _out[a]++;

  _adjIn[b].push_back(a);
  _in[b]++;

  _num_edges++;

  if (!hasEdge(b, a)) {
    _neighbours[b].push_back(a);
    _num_neighbours[b]++;

    _neighbours[a].push_back(b);
    _num_neighbours[a]++;
  }
}

void DynamicGraph::rmEdge(int a, int b) {
  ready = false;

  if (!hasEdge(a, b)) return;

  _removeVector(_adjOut[a], b);
  _out[a]--;

  _removeVector(_adjIn[b], a);
  _in[b]--;

  _num_edges--;

  if (!hasEdge(b, a)) {
    _removeVector(_neighbours[a], b);
    _num_neighbours[a]--;
    _removeVector(_neighbours[b], a);
    _num_neighbours[b]--;
  }
}

bool DynamicGraph::hasEdge(int a, int b) {
  if (!ready) {
    int i;
    for (i = 0; i < _out[a]; i++)
      if (_adjOut[a][i] == b)
        return true;
    return false;
  }

  if (!hybrid_ch[a])
    return _adjM[a][b];
  else if (_out[a] < 3) {
    int i;
    for (i = 0; i < _out[a]; i++)
      if (_adjOut[a][i] == b)
        return true;

    return false;
  }
  else {
    l_list* cur = _hashM[a][b & _hash_out[a]];
    while (cur != NULL) {
      if (cur->value == b)
        return true;

      cur = cur->next;
    }
    return false;
  }

  return false;
}

void DynamicGraph::_removeVector(vector<int> &v, int b) {
  int i, s = v.size();
  for (i=0; i<s; i++)
    if (v[i] == b) break;
  if (i<s) v.erase(v.begin()+i);
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
