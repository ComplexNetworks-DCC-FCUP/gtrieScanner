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
Command Line Client functions

Last Update: 11/02/2012
---------------------------------------------------- */

#ifndef _CMDLINE_
#define _CMDLINE

#include "GraphMatrix.h"
#include "GraphUtils.h"
#include "Error.h"
#include "Common.h"
#include "Esu.h"
#include "GTrie.h"
#include "Isomorphism.h"
#include "Timer.h"
#include "Random.h"

class CmdLine {
 private:
  static char graph_file[MAX_BUF];
  static char gtrie_file[MAX_BUF];
  static char subgraphs_file[MAX_BUF];
  static char output_file[MAX_BUF];
  static char occ_file[MAX_BUF];

  static bool dir;
  static bool occurrences;
  static bool create;

  static int motif_size;
  static int random_number;
  static int random_seed;
  static int random_exchanges;
  static int random_tries;

  static double time_original;
  static double time_random;

  static MethodType method;
  static FormatType format;
  static OutputType output;
  static Graph *g;

  static FILE *f_output;
  static FILE *f_occ;

  static GTrie *gt_original;
  static GTrie *gt;

  static GraphTree sg_original;

  static time_t t_start;
  
  static void about();
  static void defaults();
  static void parse_cmdargs(int argc, char **argv);
  static void run_esu(Graph *g, GraphTree *sg);
  static void run_gtrie(Graph *g, GraphTree *sg);
  static void run_subgraphs(Graph *g, GraphTree *sg);

  static MethodType str_to_method(char *s);
  static FormatType str_to_format(char *s);
  static OutputType str_to_output(char *s);

  static int compare_results(const void *a, const void *b);

  static void prepare_graph();
  static void prepare_files();
  static void compute_original();
  static void compute_results();
  static void show_results(ResultType *res, int nres);

  static void create_gtrie();

 public:
  static void init(int argc, char **argv);
  static void finish();
  static void decide_action();  
};

#endif
