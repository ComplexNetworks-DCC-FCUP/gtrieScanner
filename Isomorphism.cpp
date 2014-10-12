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

#include "Isomorphism.h"
#include "GraphUtils.h"
#include "GraphMatrix.h"
#include "Conditions.h"

// Static variables
setword Isomorphism::workspace[WORKSPACE_SIZE];
int Isomorphism::n,Isomorphism::m;
set *Isomorphism::gv;
graph Isomorphism::g[MAXN*MAXM];
int Isomorphism::lab[MAXN];
int Isomorphism::ptn[MAXN];
int Isomorphism::orbits[MAXN];
bool Isomorphism::dir;

DEFAULTOPTIONS(options);
statsblk(stats);
graph mm[MAXN*MAXM];


void Isomorphism::initNauty(int size, bool directed) {
  n = size;  
  m = (n + WORDSIZE - 1) / WORDSIZE;
  nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
  
  dir = directed;

  options.getcanon = TRUE;
  options.writeautoms = FALSE;
  if (dir) options.digraph = TRUE;
  else     options.digraph = FALSE;
}

void Isomorphism::finishNauty() {
  nauty_freedyn();
  naugraph_freedyn();
  nautil_freedyn();
}


void Isomorphism::canonicalStrNauty(Graph *myg, int *v, char *s) {
  int i, j, aux;

  static bool **adjM = myg->adjacencyMatrix();

  for (i=0; i<n; i++) {
    gv = GRAPHROW(g,i,m);
    EMPTYSET(gv,m);
    for (j=0; j<n; j++)
      if (adjM[v[i]][v[j]]) ADDELEMENT(gv,j);
  }

  nauty(g,lab,ptn,NULL,orbits,&options,&stats,
	workspace,WORKSPACE_SIZE,m,n,mm);

  aux=0;
  for (i=0; i<n; i++) {
    gv = GRAPHROW(mm,i,m); 
    for (j=0; j<n; j++)
         s[aux++] = ISELEMENT(gv,j)?'1':'0';
  }
  s[aux]=0;
}

void Isomorphism::_goCan(int x, int pos, const char *in, 
	   char *perm, char *used,
	   char *current, char *best, int size) {
  int i, j;

  if (strcmp(current, best)<0) return;

  if (pos == size) {
    strcpy(best, current);
  } else {

    int nedges[size], max=-1;
    for (i=0; i<size; i++)
      if (!used[i]) {
	nedges[i]=0;
	for (j=0; j<size; j++)
	  if (!used[j] && in[i*size+j]=='1')
	    nedges[i]++;
	if (nedges[i]>max) max=nedges[i];
      }
    

    for (i=0; i<size; i++)
      if (!used[i] /*&& nedges[i]==max*/) {
	used[i]=1;
	perm[pos]=i;
	current[pos*size+pos] = in[i*size+i];
	for (j=0; j<pos; j++) {
	  current[pos*size+j] = in[i*size+perm[j]]; 
	  current[j*size+pos] = in[perm[j]*size+i]; 
	}
	_goCan(i, pos+1, in, perm, used, current, best, size);
	current[pos*size+pos] = '2';
	for (j=0; j<pos; j++) {
	  current[pos*size+j] = '2';
	  current[j*size+pos] = '2';
	}
	used[i]=0;
      }
  }
}

void Isomorphism::canonicalNauty(const char *in, char *best, int size) {
  int i, j, aux;

  for (i=0; i<n; i++) {
    gv = GRAPHROW(g,i,m);
    EMPTYSET(gv,m);
    for (j=0; j<n; j++)
      if (in[i*size+j]=='1') ADDELEMENT(gv,j);
  }

  nauty(g,lab,ptn,NULL,orbits,&options,&stats,
	workspace,WORKSPACE_SIZE,m,n,mm);

  aux=0;
  for (i=0; i<n; i++) {
    gv = GRAPHROW(mm,i,m); 
    for (j=0; j<n; j++)
         best[aux++] = ISELEMENT(gv,j)?'1':'0';
  }
  best[aux]=0;
}


void Isomorphism::canonicalBasedNauty(const char *in, char *best, int size) {

  int i, j, k, min_i, ss, mymap[size];
  int degree[size], last_degree[size], total_degree[size];
  char used[size], articulation[size];
  int nbfs, bfs[size];

  for (i=0; i<size; i++) {
    degree[i] = 0;
    used[i]=0;
    for (j=0; j<size; j++) {
      if (in[i*size+j]=='1') degree[i]++;
      if (in[i+size*j]=='1') degree[i]++;
    }
    total_degree[i] = last_degree[i] = degree[i];
  }

  // Find the node with the smaller number of connections
  for (ss = size-1; ss>=0; ss--) {

    if (ss>3) {
      // Find articulation points (cannot be removed)
      for (i=0; i<size; i++) {
	if (used[i]) continue;
	used[i]=2;
	j=0; while (used[j]) j++;
	nbfs=1; bfs[0]=j; used[j]=2; j=0;
	while (j<nbfs) {
	  for (k=0; k<size; k++)
	    if (!used[k] && (in[k*size+bfs[j]]=='1' || in[k+size*bfs[j]]=='1')) {
	      used[k]=2;
	      bfs[nbfs++]=k;
	    }
	  j++;
	}
	if (nbfs!=ss) articulation[i]=1;
	else articulation[i]=0;
	for (j=0; j<size; j++)
	  if (used[j]==2) used[j]=0;
      }
    } else for (i=0; i<size; i++) articulation[i]=0;

    min_i = -1;
    for (i=0; i<size; i++) {
      if (used[i] || articulation[i]) continue;
      if (min_i<0 || degree[i] < degree[min_i]) min_i=i;
      else if (degree[i] == degree[min_i]) {
	if (last_degree[i] < last_degree[min_i]) min_i=i;
	else if (last_degree[i] == last_degree[min_i]) {
	  if (total_degree[i] < total_degree[min_i]) min_i = i;
	  else if (total_degree[i] == total_degree[min_i]) {
	    for (j=ss+1; j<size; j++)
	      if (in[i*size+mymap[j]] != in[min_i*size+mymap[j]]) {
		if (in[i*size+mymap[j]]=='0') min_i=i;
		break;
	      }
	    if (j==size) {
	      for (j=ss+1; j<size; j++)
		if (in[i+size*mymap[j]] != in[min_i+size*mymap[j]]) {
		  if (in[i+size*mymap[j]]=='0') min_i=i;
		  break;
		}	      
	    }
	  }
	}
      }
    }
    
    if (ss!=0 && degree[min_i]==0) {
      fprintf(stderr, "[%s]\n", in);
      fprintf(stderr, "!!! OH NO %d [canonicalBasedNauty]\n", ss); fflush(stdout);
      exit(1);
    }
    
    for (i=0; i<size; i++) {
      last_degree[i] = degree[i];
      if (in[i*size+min_i]=='1') degree[i]--;
      if (in[i+size*min_i]=='1') degree[i]--;
    }
    mymap[ss] = min_i;
    used[min_i]=1;
  }

  best[size*size]=0;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      best[i*size+j] = in[mymap[i]*size+mymap[j]];  
}


void Isomorphism::canonicalBigger(const char *in, char *best, int size) {
  int i,j;
  int nedges[size], max;
  char perm[size], current[size*size+1], used[size];

  for (i=0, max=0; i<size; i++) {
    nedges[i]=0;
    for (j=0; j<size; j++)
      nedges[i]+=in[i*size+j]-'0';
    if (nedges[i] > nedges[max]) max = i;
  }


  for (i=0; i<size*size; i++)
    best[i] = in[size*size-1-i];    
  best[size*size]=0;

  for (i=0; i<size; i++) used[i]=0;
  for (i=0; i<size*size; i++) current[i]='2';

  current[size*size]='\0';
  for (i=0; i<size; i++)
    if (nedges[i]==nedges[max]) {
      used[i]=1;
      perm[0]=i;
      current[0]=in[i*size+i];
      _goCan(i, 1, in, perm, used, current, best, size);
      current[0]=i;
      used[i]=0;
    }
}

void Isomorphism::_goCan2(int pos, const char *in, int *perm, bool *used, char *best, int size, char *current) {
  int i, j;

  if (strcmp(current,best)<0) return;

  if (pos==size) {
    strcpy(best, current);
    return;
  }
  
  // At the beggining everyone is a possible candidate
  bool candidate[size];
  for (i=0; i<size; i++)
    if (used[i]) candidate[i]=false;
    else         candidate[i]=true;


  for (i=0; i<pos; i++) {
    for (j=0; j<size ;j++)
      if (candidate[j] && in[perm[i]*size+j]=='1') break;
    if (j<size) // at least one neighbour
      for (j=0; j<size ;j++)
	if (candidate[j] && in[perm[i]*size+j]=='0') candidate[j]=false;
  }

  // Among those that are still candidates, select the
  // biggest neighborhood on those candidates
  int degree[size];
  int max=-1;
  for (i=0; i<size; i++) {
    degree[i]=0;
    if (candidate[i]) {
      for (j=0; j<size; j++)
	if (i!=j && candidate[j] && in[j*size+i]=='1')
	  degree[i]++;
      if (degree[i]>max) max=degree[i];
    }
  }

  for (i=0; i<size; i++)
    if (candidate[i] && degree[i]==max) {
      perm[pos]=i;
      used[i]=true;
      current[pos*size+pos] = in[i*size+i];
      for (j=0; j<pos; j++) {
	current[pos*size+j] = in[i*size+perm[j]]; 
	current[j*size+pos] = in[perm[j]*size+i]; 
      }

      _goCan2(pos+1, in, perm, used, best, size, current);

      current[pos*size+pos] = '2';
      for (j=0; j<pos; j++) {
	current[pos*size+j] = '2';
	current[j*size+pos] = '2';
      }
      used[i]=false;
    }

}

void Isomorphism::canonicalBigger2(const char *in, char *best, int size) {
  int i;

  int perm[size];
  bool used[size];

  for (i=0; i<size; i++) used[i]=false;

  char current[size*size+1];
  for (i=0; i<size*size; i++) current[i]='2';
  current[size*size]=0;

  strcpy(best, in);

  _goCan2(0, in, perm, used, best, size, current);
}


void Isomorphism::canonicalRandom(const char *in, char *best, int size) {
  int n, i, j, aux;
  bool used[size];
  int v[size];
  int ncand, cand[size];
  bool adj[size][size];

  for (i=0; i<size; i++) {
    used[i]      = false;
    for (j=0; j<size; j++)
      if (in[i*size+j]=='1' || in[j*size+i]=='1') adj[i][j]=true;
      else adj[i][j]=false;
  }


  for (n=0; n<size; n++) {
    
    ncand = 0;
    for (i=0; i<size; i++) {
      if (used[i]) continue;
      for (j=0;j<n;j++)
	if (adj[i][v[j]]) break;
      if (j<n) cand[ncand++]=i;
    }

    if (ncand==0)
      for (i=0; i<size; i++)
	if (!used[i]) cand[ncand++]=i;

    aux = random()%ncand;

    v[n]=cand[aux];
    used[cand[aux]]=true;
  }
  

  best[size*size]=0;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      best[i*size+j] = in[v[i]*size+v[j]];
}

void Isomorphism::canonicalBasedNautyReverse(const char *in, char *best, int size) {

  int n, i, j, max_i, conn, max_conn;
  bool used[size], flag;
  int v[size];
  int degree_levels = size;
  int degree[degree_levels][size];
  bool adj[size][size];

  for (i=0; i<size; i++) {
    used[i]      = false;
    degree[0][i] = 0;
    for (j=0; j<size; j++)
      if (in[i*size+j]=='1'|| in[j*size+i]=='1') {
	adj[i][j]=true;
	degree[0][i]++;
      }
      else adj[i][j]=false;
  }
 
  for (j=1; j<degree_levels; j++)
    for (n=0; n<size; n++) {
      degree[j][n]=0;
      for (i=0; i<size; i++)
	if (adj[n][i]) degree[j][n]+=degree[j-1][i];
    }

  for (n=0; n<size; n++) {
    max_i = max_conn = INT_MAX;
    for (i=0;i<size;i++) {
      if (used[i]) continue;
      conn = 0;

      for (j=0;j<n;j++)
	if (adj[i][v[j]]) conn++;
      
      flag = false;

      if (max_conn==INT_MAX) flag=true;
      else if (conn>0 && conn < max_conn) flag = true;
      else if (max_conn==0 && conn>0) flag = true;
      else if (conn == max_conn) {
	for (j=0; j<degree_levels; j++)
	    if (degree[j][i] != degree[j][max_i]) {
	      if (degree[j][i] < degree[j][max_i]) flag=true;
	      break;
	    }
	if (j==degree_levels) {
	  for (j=0; j<n; j++)	  
	    if (adj[i][v[j]] != adj[max_i][v[j]]) {
	      if (!adj[i][v[j]]) flag = true;
	      break;
	    }
	}
      }

      if (flag) {
	max_i      = i;
	max_conn   = conn;
      }
    }
    used[max_i]=true;
    v[n] = max_i;
  }

  best[size*size]=0;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      best[i*size+j] = in[v[i]*size+v[j]];
}

