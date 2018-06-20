#!/usr/bin/python

import sys

fname = sys.argv[1] # input the outtree file from phylip that has the tree info for all isolates.
# get one isolate from each leafset: nearby branched isolates in the tree

fh = open(fname, 'r')
fout = open(fname+".independent.txt", 'w')

alllines = ""
for line in fh:
    line = line.rstrip('\n')
    alllines = alllines + line

leftbrackets = alllines.split('(')
for lb in leftbrackets:
    if ')' in lb:
        rightbracidx = lb.find(')') # find index for the first right bracket
        if rightbracidx > -1: # if right bracket doesn't exist, the index will return -1
            leafset = lb[0:rightbracidx] 
            leaves = leafset.split(',')
            strains = []
            branchlen = []
            for l in leaves:
                content = l.split(':')
                strains.append(content[0])
                branchlen.append(float(content[1]))
            maxlen = max(branchlen) # find isolate farthest from the root of the tree
            maxidx = branchlen.index(maxlen) # find index of the strain with max branchlen
            fout.write(strains[maxidx]+"\n")
fh.close()
fout.close()
