#!/usr/bin/python

####
#
# Reads a list of edges file, removes
# duplicates and self loops and writes
# back to file
#
####

import sys
import networkx as nx

if len(sys.argv) != 4 and len(sys.argv) != 3:
  print "Usage: cleaner.py <in-file> <out-file> <type>\n(type: 0 - undirected, 1 - directed)"
  print "\tor"
  print "Usage: cleaner.py <in/out-file> <type>\n(type: 0 - undirected, 1 - directed)"
  exit()

infile = sys.argv[1]
outfile = sys.argv[2 if len(sys.argv) == 4 else 1]
directed = sys.argv[3 if len(sys.argv) == 4 else 2] == '1'

if not(directed):
  G = nx.read_edgelist(infile, create_using=nx.Graph())
else:
  G = nx.read_edgelist(infile, create_using=nx.DiGraph())

fl = open(outfile, "w")

ls = []
for fr, nei in G.adjacency_iter():
  for to, _ in nei.items():
    if (directed and int(fr) != int(to)) or (not(directed) and int(fr) < int(to)):
      ls.append((int(fr), int(to)))

ls.sort()
for i in ls:
  fl.write("%d %d\n" % (i[0], i[1]))

fl.close()
