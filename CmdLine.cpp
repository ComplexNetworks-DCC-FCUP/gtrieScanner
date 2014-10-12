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

Last Update: 23/02/2012
---------------------------------------------------- */

#include "CmdLine.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <cmath>

// ----------------------------------------------

// Variable declarations

char CmdLine::graph_file[MAX_BUF];
char CmdLine::gtrie_file[MAX_BUF];
char CmdLine::subgraphs_file[MAX_BUF];
char CmdLine::output_file[MAX_BUF];
char CmdLine::occ_file[MAX_BUF];

bool CmdLine::dir;
bool CmdLine::occurrences;
bool CmdLine::create;

int  CmdLine::motif_size;
int  CmdLine::random_number;
int  CmdLine::random_seed;
int  CmdLine::random_exchanges;
int  CmdLine::random_tries;

double CmdLine::time_original;
double CmdLine::time_random;

MethodType CmdLine::method;
FormatType CmdLine::format;
OutputType CmdLine::output;
Graph *CmdLine::g;

FILE *CmdLine::f_output;
FILE *CmdLine::f_occ;

GTrie *CmdLine::gt;
GTrie *CmdLine::gt_original;

GraphTree CmdLine::sg_original;

time_t CmdLine::t_start;

// ----------------------------------------------

// Create a g-trie from a list of subgraphs
void CmdLine::create_gtrie() {

  // Check motif size
  if (motif_size<MIN_MOTIF_SIZE || motif_size>MAX_MOTIF_SIZE)
    Error::msg("Invalid subgraph size (%d <= SIZE <= %d)", MIN_MOTIF_SIZE, MAX_MOTIF_SIZE);

  // Create the g-trie
  printf("Creating G-Trie\n");
  gt_original = new GTrie();
  Timer::start(0);
  gt_original->readSubgraphs(motif_size, dir, subgraphs_file);
  Timer::stop(0);
  printf("Creation time: %.2f\n", Timer::elapsed(0));
  printf("Nr %d-subgraphs in g-trie: %d\n", motif_size, gt_original->countGraphs()); 
  printf("Compression rate = %.2f%%\n\n", gt_original->compressionRate()*100);

  // Write g-trie to file
  gt_original->writeToFile(output_file);
  printf("G-Trie written to file \"%s\"\n", output_file);

  delete gt_original;
}

// ----------------------------------------------

// Run ESU algorithm on graph 'g' and store results on GraphTree 'gt'
void CmdLine::run_esu(Graph *g, GraphTree *sg) {
  Esu::countSubgraphs(g, motif_size, sg);
}

// Run SUBGRAPHS method on graph 'g' and store results on GraphTree 'gt'
void CmdLine::run_subgraphs(Graph *g, GraphTree *sg) {
  printf("Creating G-Trie\n");
  gt_original = new GTrie();
  Timer::start(0);
  gt_original->readSubgraphs(motif_size, dir, subgraphs_file);
  Timer::stop(0);
  printf("Creation time: %.2f\n", Timer::elapsed(0));
  printf("Compression rate = %.2f%%\n\n", gt_original->compressionRate()*100);

  gt_original->census(g);
  gt_original->populateGraphTree(sg, motif_size);
}

// Run GTRIES method on graph 'g' and store results on GraphTree 'gt'
void CmdLine::run_gtrie(Graph *g, GraphTree *sg) {
  printf("Reading G-Trie\n");
  gt_original = new GTrie();
  Timer::start(0);
  gt_original->readFromFile(gtrie_file);
  Timer::stop(0);
  printf("Creation time: %.2f\n", Timer::elapsed(0));
  printf("Compression rate = %.2f%%\n\n", gt_original->compressionRate()*100);

  gt_original->census(g);
  gt_original->populateGraphTree(sg, motif_size);
}

// ----------------------------------------------

// Compare two different motif results (for sorting)
int CmdLine::compare_results(const void *a, const void *b) {
  ResultType *r1 = (ResultType *)a;
  ResultType *r2 = (ResultType *)b;

  if (r1->z_score > r2->z_score) return -1;
  if (r1->z_score < r2->z_score) return +1;
  if (r1->f_original > r2->f_original) return -1;
  if (r1->f_original < r2->f_original) return +1;
  return 0;
}

// ----------------------------------------------

// Prepare files
void CmdLine::prepare_files() {

  // Check output format
  if (output == NOOUTPUT)
    Error::msg("No valid output format specified");    

  if (!strcmp(output_file,  DEFAULT_RESULTS)) {
    if (output == TEXT) strcat(output_file, ".txt");
    else if (output == HTML) strcat(output_file, ".html");
  }

  f_output = fopen(output_file, "w");
  if (f_output==NULL)
    Error::msg("Unable to open output file \"%s\"", output_file);    

  if (occurrences) {
    f_occ = fopen(occ_file, "w");
    if (f_occ==NULL)
      Error::msg("Unable to open occurrences file \"%s\"", occ_file);
  }
}

// Prepare the original graph for computation
void CmdLine::prepare_graph() {

  // Check motif size
  if (motif_size<MIN_MOTIF_SIZE || motif_size>MAX_MOTIF_SIZE)
    Error::msg("Invalid subgraph size (%d <= SIZE <= %d)", MIN_MOTIF_SIZE, MAX_MOTIF_SIZE);

  // Check if graph file name was given
  if (!strcmp(graph_file, INVALID_FILE))
    Error::msg("No graph file specified");

  // Check graph format
  if (format == NOFORMAT)
    Error::msg("No valid graph format specified");

  // Read the graph file
  g = new GraphMatrix();

  // Use simple or simple_weight text format
  if (format == SIMPLE)
    GraphUtils::readFileTxt(g, graph_file, dir, false);
  else if  (format == SIMPLE_WEIGHT)
    GraphUtils::readFileTxt(g, graph_file, dir, true);
  else printf("??");

  // sort and create array of neighbours
  g->sortNeighbours();
  g->makeArrayNeighbours();

  // Print chosen parameters
  printf("Subgraph Size: %d\n", motif_size);
  printf("Graph File: %s\n", graph_file);
  printf("%s, %d Nodes, %d Edges\n", dir?"Directed":"Undirected",g->numNodes(), dir?g->numEdges():g->numEdges()/2);
}

// Count subgraphs on original network
void CmdLine::compute_original() {

  // Print method description
  if (method == ESU)
    printf("Method:     ESU on original network\n");
  else if (method == GTRIE)
    printf("Method:     GTRIE with file containing complete g-trie\n");
  else if (method == SUBGRAPHS)
    printf("Method:     GTRIE with subgraphs read from file\n");

  // Compute frequency
  Global::show_occ = occurrences;
  Global::occ_file = f_occ;
  puts("\nCounting subgraph frequency on 'ORIGINAL NETWORK'");
  Timer::start(0);
  if (method == ESU) run_esu(g, &sg_original);
  else if (method == SUBGRAPHS) run_subgraphs(g, &sg_original);
  else if (method == GTRIE) run_gtrie(g, &sg_original);
  Timer::stop(0);  
  printf("%d subgraphs, ",   sg_original.countGraphs());
  printf("%.0f occurrences\n", sg_original.countOccurrences());
  printf("Time elapsed: %.6fs\n\n", Timer::elapsed(0));
  Global::show_occ = false;

  time_original = Timer::elapsed(0);
}

// Compute random networks and output results
void CmdLine::compute_results() {
  int i, j;
  mapStringInt:: const_iterator ii, iiend;

  // Create map and init results
  mapStringInt m_original;  
  sg_original.populateMap(&m_original, motif_size);
  ResultType res[m_original.size()];
  for (ii=m_original.begin(), iiend=m_original.end(), i=0; ii!=iiend; ii++, i++) {
    res[i].s = strdup((ii->first).c_str());    
    res[i].f_original = ii->second;
    res[i].z_score = res[i].avg_random = res[i].dev_random = 0;
  }

  // Do we have random networks to compute?
  if (random_number > 0) {
    gt = new GTrie;
    sg_original.populateGTrieNauty(gt, motif_size, dir);
    
    mapStringInt m_count[random_number];

    // Generate all random networks
    printf("Computing random networks: ");
    time_random = 0;
    for (i=0; i<random_number; i++) {      

      // Create new random network from previous one
      g->makeVectorNeighbours();
      Random::markovChainPerEdge(g, random_exchanges, random_tries);
      g->sortNeighbours();
      g->makeArrayNeighbours();

      // Compute census
      Timer::start(0);
      gt->census(g);
      Timer::stop(0);
      time_random += Timer::elapsed(0);
      gt->populateMap(&m_count[i], motif_size);
      fputc('.', stdout);      
    }
    fputc('\n', stdout);
    time_random /= (double)random_number;
    printf("Avg time per random network: %.6fs\n\n", time_random);

    // Compute significance
    for (ii=m_original.begin(), iiend=m_original.end(), i=0; ii!=iiend; ii++, i++) {
      // Average frequency
      double avg = 0;
      for (j=0; j<random_number; j++)
	avg += m_count[j][ii->first];
      avg /= random_number;

      // Standard deviation
      double dev=0;
      for (j=0; j<random_number; j++)
	dev += double(m_count[j][ii->first]-avg)*double(m_count[j][ii->first]-avg)/double(random_number-1);
      dev = sqrt(dev);

      // zscore
      double zscore = (ii->second - avg)/dev;

      res[i].avg_random = avg;
      res[i].dev_random = dev;
      res[i].z_score    = zscore;
    }
  }

  // Sort results
  qsort(res, m_original.size(), sizeof(ResultType), compare_results);

  // print results
  show_results(res, m_original.size());
  printf("Results written to file \"%s\"\n", output_file);
  if (occurrences)
    printf("Occurences on original network written to file \"%s\"\n", occ_file);

  for (i=0; i<(int)m_original.size(); i++)
    free(res[i].s);
}

// ----------------------------------------------

// Print results
void CmdLine::show_results(ResultType *res, int nres) {
  int i, j, k;
  bool html = (output==HTML)?true:false;

  if (html) fprintf(f_output, "<html><head>\n"
		    "<title>%s results</title>\n"
		    "<style type=\"text/css\">"
		    "p,body,td,th,li  {font-size: 13px; font-family: arial, geneva, \"sans serif\";color: black; font-weight: normal;}\n"
		    "body, tr {background-color: #fbfbfb;}\n"
		    "hr {color: black;}\n"
		    "table {background-color: black;}\n"
		    "th {color:white;}\n"
		    "td {text-align: right;}\n"
		    ".odd  {background-color: #d1d1d1;}\n"
		    ".even {background-color: #f2f2f2;}\n"
		    ".hd   {background-color: #444444; color: white;}\n"
		    ".pre  {white-space: pre;}\n"
		    "</style>"
		    "</head>\n<body>",
		    PROGRAM_NAME);

  fprintf(f_output, "%s%s Results%s\n", html?"<h1>":"", PROGRAM_NAME, html?"</h1>":"");
  
  if (html) fprintf(f_output, "<hr>\n");
  else      fprintf(f_output, "%s\n", SEPARATOR);

  fprintf(f_output, "%sGeneral Information%s\n", html?"<h2>":"", html?"</h2>":"\n");

  time_t t_end = time(0);
  struct tm *tm_start = localtime(&t_start);
  fprintf(f_output, "%sStart of Computation:%s %02dh%02dm%02ds %02d/%02d/%02d\n",
	  html?"<li><b>":"",html?"</b>":"",
	  tm_start->tm_hour, tm_start->tm_min, tm_start->tm_sec,
	  tm_start->tm_mday, tm_start->tm_mon+1, 1900+tm_start->tm_year);
  struct tm *tm_end   = localtime(&t_end);
  fprintf(f_output, "%sEnd of Computation:%s %02dh%02dm%02ds %02d/%02d/%02d\n",
	  html?"<li><b>":"",html?"</b>":"",
	  tm_end->tm_hour, tm_end->tm_min, tm_end->tm_sec,
	  tm_end->tm_mday, tm_end->tm_mon+1, 1900+tm_end->tm_year);

  if (html) fprintf(f_output, "<br>&nbsp;\n");
  else      fprintf(f_output, "\n");

  fprintf(f_output, "%sSubgraph Size:%s %d\n", html?"<li><b>":"", html?"</b>":"", motif_size);
  fprintf(f_output, "%sGraph File:%s \"%s\"\n", html?"<li><b>":"", html?"</b>":"",graph_file);
  fprintf(f_output, "%sDirected:%s %s\n", html?"<li><b>":"", html?"</b>":"",dir?"YES":"NO");
  fprintf(f_output, "%sNr Nodes:%s %d\n", html?"<li><b>":"", html?"</b>":"",g->numNodes());
  fprintf(f_output, "%sNr Edges:%s %d\n", html?"<li><b>":"", html?"</b>":"",dir?g->numEdges():g->numEdges()/2);
  
  if (html) fprintf(f_output, "<br>&nbsp;\n");
  else      fprintf(f_output, "\n");

  fprintf(f_output, "%sMethod:%s ", html?"<li><b>":"", html?"</b>":"");
  if (method == ESU)
    fprintf(f_output, "ESU on original network\n");
  else if (method == GTRIE)
    fprintf(f_output, "GTRIE with file containing complete g-trie\n");
  else if (method == SUBGRAPHS)
    fprintf(f_output, "GTRIE with subgraphs read from file\n");

  fprintf(f_output, "%sDifferent Types of Subgraphs Found [Original Network]:%s %d\n", html?"<li><b>":"", html?"</b>":"", sg_original.countGraphs());
  fprintf(f_output, "%sSubgraph Occurrences Found [Original Network]:%s %.0f\n", html?"<li><b>":"", html?"</b>":"", sg_original.countOccurrences());
  fprintf(f_output, "%sTime for computing census on original network%s: %.6fs\n", html?"<li><b>":"", html?"</b>":"", time_original);
  fprintf(f_output, "%sAverage time for census on random network%s: %.6fs\n", html?"<li><b>":"", html?"</b>":"", time_random);

  if (html) fprintf(f_output, "<br>&nbsp;\n");
  else      fprintf(f_output, "\n");

  fprintf(f_output, "%sNumber of random networks:%s %d\n", html?"<li><b>":"", html?"</b>":"", random_number);
  fprintf(f_output, "%sRandom seed:%s %d\n", html?"<li><b>":"", html?"</b>":"", random_seed);
  fprintf(f_output, "%sExchanges per edge:%s %d\n", html?"<li><b>":"", html?"</b>":"", random_exchanges);
  fprintf(f_output, "%sNumber of tries per exchange:%s %d\n", html?"<li><b>":"", html?"</b>":"", random_tries);

  if (html) fprintf(f_output, "<hr>\n");
  else      fprintf(f_output, "%s\n", SEPARATOR);

  fprintf(f_output, "%sMotif Analysis Results%s\n", html?"<h2>":"", html?"</h2>":"");

  char adj[motif_size*motif_size+motif_size];

  if (html) fprintf(f_output, "<table cellpadding=\"3\" cellspacing=\"2\">\n<tr class=\"hd\"><th colspan=\"2\">Subgraph</th><th>Org. Frequency</th><th>Z-score</th><th>Rnd. Frequency</th></tr>\n");
  else fprintf(f_output, "\nGraph%*s   Org_Freq |  Z-score |    Rnd_Avg +/-    Rnd_Dev\n\n", motif_size>5?motif_size-5:0, "");
  for (i=0; i<nres; i++) {

    for (j=0, k=0; res[i].s[j]; j++, k++) {
      if (j>0 && j%motif_size==0) adj[k++] = '\n';
      adj[k] = res[i].s[j];
    }
    adj[k]=0;

    if (html) 
      fprintf(f_output, "<tr class=\"%s\"><td><img src=\"http://www.dcc.fc.up.pt/gtries/graph.php?%swidth=75&height=75&adj=%s\"></td><td class=\"pre\">%s</td><td>%d</td><td>%.2f</td><td>%.2f +/- %.2f</td></tr>\n",
	      (i%2)?"odd":"even", dir?"dir&":"",
	      res[i].s, adj, res[i].f_original,
	      res[i].z_score, res[i].avg_random, res[i].dev_random);
    else
      fprintf(f_output, "%s%*s %10d | %8.2f | %10.2f +/- %10.2f\n\n",
	      adj, motif_size<5?5-motif_size:0, "", res[i].f_original,
	      res[i].z_score, res[i].avg_random, res[i].dev_random);
  }
  if (html) fprintf(f_output, "</table>\n");


}

// ----------------------------------------------

// Execute program with the chosen parameters
void CmdLine::decide_action() {

  t_start = time(0);

  if (create) {
    create_gtrie();
  } else {
      // Check method
      if (method == NOMETHOD)
	Error::msg("No valid method specified");    
      
      prepare_graph();
      prepare_files();
      compute_original();
      compute_results();
  }

}

// ----------------------------------------

// Initialize everything
void CmdLine::init(int argc, char **argv) {
  about();
  defaults();
  parse_cmdargs(argc, argv);
  Isomorphism::initNauty(motif_size, dir);
}

// Finish everything
void CmdLine::finish() { 
  if (!create) {
    fclose(f_output);
    if (occurrences) fclose(f_occ);
  }
  Isomorphism::finishNauty();
}

// About mesage
void CmdLine::about() {
  puts(SEPARATOR);
  printf("%s (version %s)\n", PROGRAM_NAME, VERSION );
  puts(SEPARATOR);
}

// ----------------------------------------

// Load default values
void CmdLine::defaults() {

  strcpy(graph_file, INVALID_FILE);
  dir        = false;
  motif_size = -1;

  method = ESU;

  strcpy(output_file, DEFAULT_RESULTS);
  strcpy(occ_file, DEFAULT_OCC);

  random_number    = 0;
  random_seed      = time(0);
  random_exchanges = 3;
  random_tries     = 10;

  create = false;
  format = SIMPLE_WEIGHT;
  output = TEXT;
  occurrences = false;
}

// ----------------------------------------------

// Convert string to method type
MethodType CmdLine::str_to_method(char *s) {
  if      (!strcmp(s, "esu"))       return ESU;
  else if (!strcmp(s, "gtrie"))     return GTRIE;
  else if (!strcmp(s, "subgraphs")) return SUBGRAPHS;
  else return NOMETHOD;
}

// Convert string to graph format type
FormatType CmdLine::str_to_format(char *s) {
  if      (!strcmp(s, "simple"))        return SIMPLE;
  else if (!strcmp(s, "simple_weight")) return SIMPLE_WEIGHT;
  else return NOFORMAT;
}

// Convert string to graph format type
OutputType CmdLine::str_to_output(char *s) {
  if      (!strcmp(s, "txt"))  return TEXT;
  else if (!strcmp(s, "html")) return HTML;
  else return NOOUTPUT;
}

// Parse all command line arguments
void CmdLine::parse_cmdargs(int argc, char **argv) {
    
  for (int i=1; i<argc; i++) {

    // Create g-trie ?
    if (!strcmp("-c",argv[i]) || !strcmp("--create",argv[i])) {
      create=true;
      strcpy(subgraphs_file, argv[++i]);
    }
    
    // Graph file
    else if (!strcmp("-g",argv[i]) || !strcmp("--graph",argv[i])) {
      strcpy(graph_file, argv[++i]);
    }

    // Output file
    else if (!strcmp("-o",argv[i]) || !strcmp("--output",argv[i])) {
      strcpy(output_file, argv[++i]);
    }

    // Output file
    else if (!strcmp("-oc",argv[i]) || !strcmp("--occurrences",argv[i])) {
      occurrences = true;
      strcpy(occ_file, argv[++i]);
    }

    // Output format
    else if (!strcmp("-t",argv[i]) || !strcmp("--type",argv[i])) {
      output = str_to_output(argv[++i]);
    }

    // Size of motifs to consider
    else if (!strcmp("-s",argv[i]) || !strcmp("--size",argv[i])) {
      motif_size = atoi(argv[++i]);
    }

    // Method for set of subgraphs
    else if (!strcmp("-m",argv[i]) || !strcmp("--method",argv[i])) {
      method = str_to_method(argv[++i]);
      if (method == SUBGRAPHS) strcpy(subgraphs_file, argv[++i]);
      else if (method == GTRIE) strcpy(gtrie_file, argv[++i]);
    }

    // Graph format
    else if (!strcmp("-f",argv[i]) || !strcmp("--format",argv[i])) {
      format = str_to_format(argv[++i]);
    }
    
    // Directed Graph
    else if (!strcmp("-d",argv[i]) || !strcmp("--directed",argv[i])) {
      dir=true;
    }

    // Undirected Graph
    else if (!strcmp("-u",argv[i]) || !strcmp("--undirected",argv[i])) {
      dir=false;
    }

    // Number of random networks
    else if (!strcmp("-r",argv[i]) || !strcmp("--random",argv[i])) {
      random_number = atoi(argv[++i]);
    }

    // Random seed
    else if (!strcmp("-rs",argv[i]) || !strcmp("--rseed",argv[i])) {
      random_seed = atoi(argv[++i]);
    }

    // Number of exchanges per edge
    else if (!strcmp("-re",argv[i]) || !strcmp("--rexchanges",argv[i])) {
      random_exchanges = atoi(argv[++i]);
    }

    // Number of tries per node
    else if (!strcmp("-rt",argv[i]) || !strcmp("--rtries",argv[i])) {
      random_tries = atoi(argv[++i]);
    }

  }

  // If no random seed given, initialize with time
  // (not an optimal choice, but present here for portability)
  if (random_seed<0) Random::seed(time(NULL));
  else               Random::seed(random_seed);
}
