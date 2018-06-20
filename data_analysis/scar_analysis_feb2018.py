#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 12})
from matplotlib.ticker import NullFormatter
import random
import numpy as np
from scipy.stats.stats import pearsonr
import scipy.stats as ss
from sklearn.neighbors import KernelDensity
import math
import scipy.spatial.distance
from scar_analysis_functions_feb2018 import *

print "scipy:", scipy.version.full_version
print "matplotlib:", matplotlib.__version__
print "numpy:", np.version.full_version

genes, gposdict = getGenePos() # get start and end positions of the H37Rv genes
# get essential genes
essgenes = getEssGenes()
protgenes, compgenes = getH37Rv_proteins("H37Rv_genpept.gp") # GenPept file of reference genome H37Rv
print "There are ", len(genes), " genes of which ", len(protgenes), " are protein coding."
ess_protcoding = set(essgenes).intersection(set(protgenes))
print "There are ", len(essgenes), " ess genes of which ", len(ess_protcoding), " are protein coding"

# read scar genes and find how many have low complexity
scarfh = open("scargenes.txt", 'r')
allscargenes = [] # list of scargenes
for line in scarfh:
    line = line.rstrip("\n")
    allscargenes.append(line)
scarfh.close()
scargenes = list(set(allscargenes))
print "Genes that have scars:", len(scargenes)

# plot gene quadrant for scar indels
getScarIndels(gposdict, scargenes, compgenes) # get scar indels
#exit()

snpfile = "all_qscore200_SNPpos.txt"
indelfile = "all_indels_feb2018.txt"
FIGPREFIX = "all_indels_"
#indelfile = "frameshift_indels_only_feb2018.txt"
#FIGPREFIX = "frameshift_indels_"

indelpos = []
fh = open(indelfile, 'r')
for line in fh:
    contents = line.split()
    indelpos.append(int(contents[0]))
fh.close()
numindelpos = len(indelpos)

# get indel positions that are far away from other indels, so that the decline in LC/SE can be investigated.
dist = 100 # if two mutations/indels are within this distance, they are clustered together
#snppos, snpclusters, snporphans, snpmutdist = getClusters(snpfile, dist) # inter-snp distance of 250 leaves no orphans
indelpos, indelclusters, indelorphans, indelmutdist = getClusters(indelfile, dist)
#print "There are ", len(snppos), " snps in ", len(snpclusters), " snp clusters and ", len(snporphans), " orphans."
#print "There are ", len(indelpos), " indels in ", len(indelclusters), " indel clusters and ", len(indelorphans), " orphans."
#orphandist = [] # distances between indel positions that are not near other indels
#for i in range(len(indelorphans)-1):
#    orphandist.append(indelorphans[i+1]-indelorphans[i])
#print "distance between orphan indels: min, max, and avg:", min(orphandist), max(orphandist), float(sum(orphandist))/len(orphandist)

# read the LC and H scores computed for the H37Rv genome
LCscores = readScoreFile("flankingregion_H37Rv_LCscore_w10.out")
SEscores = readScoreFile("flankingregion_H37Rv_SEscore_w10.out")
meanLC = np.mean(np.array(LCscores))
stdLC = np.std(np.array(LCscores))
meanSE = np.mean(np.array(SEscores))
stdSE = np.std(np.array(SEscores))

LCthres = 0.55 #meanLC - 2 * stdLC # 0.55 is the lowest LC value in the -7 to +14 window around an orphan indel position
SEthres = 0.931 #meanSE - 2 * stdSE # 0.931 is the lowest H value in the -7 to +14 window around an orphan indel position

print "LC mean, sd, and thres:", meanLC, stdLC, LCthres 
print "SE mean, sd, and thres:", meanSE, stdSE, SEthres 

# get a sample of non-indel pos that are 100 bases away from indel pos on both sides
nonindel_sample = getNonIndelPos(len(LCscores), indelpos, 10000, 100) # len(LCscores) is genome length
print "Num of pos in nonindel sample:", len(nonindel_sample)
# get complexity scores in window -windowdist/2 to +windowdist/2 for the nonindel positions
windowdist = 40
LC_nonindel_itrap = plotComplexNbg(nonindel_sample, LCscores, windowdist/2, "Average $LC$", "nonindel_LC")
H_nonindel_itrap = plotComplexNbg(nonindel_sample, SEscores, windowdist/2, "Average $H$", "nonindel_H")

LCindeltrap = plotComplexNbg(indelorphans, LCscores, windowdist/2, "Average $LC$", "LC") # plot LC scores of indels for neighboring pos as far as windowdist/2
SEindeltrap = plotComplexNbg(indelorphans, SEscores, windowdist/2, "Average $H$", "H") # plot SE scores of indels for neighboring pos as far as windowdist/2

print "LC indel trap:", LCindeltrap[windowdist/2-7:windowdist/2+16] # the window -7 to +15 has the indeltrap signature
print "SE indel trap:", SEindeltrap[windowdist/2-7:windowdist/2+16] # the window -7 to +15 has the indeltrap signature
print "LC nonindel itrap:", LC_nonindel_itrap[windowdist/2-7:windowdist/2+16] # the window -7 to +15 has the _nonindeltrap signature
print "H nonindel itrap:", H_nonindel_itrap[windowdist/2-7:windowdist/2+16] # the window -7 to +15 has the _nonindeltrap signature

# randomly sample nonidl pos and compare to indel positions
#compareIndelNonindel(LCscores, SEscores, d, indelpos, FIGPREFIX)

geneindels, geneindelpos = getGenewiseIndels( indelfile, gposdict, genes ) # get num of indels in each gene
geneswithindels = [g for g in geneindels if geneindels[g] > 0]
print "Num genes with indels:", len(geneswithindels)
peppegenes = [] # names of PE PPE genes
num_allindels = 0 # total number of indels
num_peppeindels = 0 # total number of indels in PE PPE genes
for g in geneindels:
    num_allindels += geneindels[g]
    if g.startswith("PE") or g.startswith("PPE"):
        peppegenes.append(g)
        num_peppeindels += geneindels[g]

print "Total num of genic indels:", num_allindels, " Num indels in PE PPE:", num_peppeindels
peppe_protcoding = set(peppegenes).intersection(set(protgenes))
non_peppe_non_ess_genes = []
for pg in protgenes:
    if pg not in essgenes and pg not in peppegenes:
        non_peppe_non_ess_genes.append(pg)
    else:
        continue
print "There are ", len(genes), " genes and ", len(gposdict), " gpos."
print "There are ", len(peppegenes), " PE PPE genes of which ", len(peppe_protcoding), " are protein coding"
print "There are ", len(non_peppe_non_ess_genes), " genes that are neither essential nor PE PPE." 
#plotIndelDist_fig1a(indelfile, genes, geneindels, essgenes, peppegenes, FIGPREFIX)
#exit()
plotIndelDist(indelfile, geneindels, essgenes, peppegenes, FIGPREFIX)
exit()
plotIndelQuad(geneindelpos, gposdict, compgenes, ess_protcoding, peppe_protcoding, scargenes, FIGPREFIX)
exit()

print "Num genes with indels:", len(geneswithindels)
indel_ess_genes = set(geneswithindels).intersection(set(essgenes))
print "Intersection of ess genes and indel genes:", len(indel_ess_genes), "\n", indel_ess_genes 

#exit()
allgenes_avgLC, allgenes_avgSE, allgenes_stdLC, allgenes_stdSE = getGeneComplexity(gposdict, LCscores, SEscores)
allgenes_LCfrac = getLowComplexFrac(gposdict, LCscores, LCthres)
allgenes_SEfrac = getLowComplexFrac(gposdict, SEscores, SEthres)

# plot complexity features (avg complexity score and low complex frac) vs num indels in the gene
#plotGenewiseComplexity(geneindels, allgenes_avgLC, FIGPREFIX+"numindels_vs_avgLC_pergene_feb2018.pdf", "Mean LC of Gene", scargenes, essgenes)
plotGenewiseComplexity(geneindels, gposdict, allgenes_avgLC, FIGPREFIX+"numindels_vs_avgLC_pergene_ALLgenes2_feb2018.pdf", "Mean $LC$ of Gene", scargenes, essgenes)
plotGenewiseComplexity(geneindels, gposdict, allgenes_avgSE, FIGPREFIX+"numindels_vs_avgSE_pergene_ALLgenes2_feb2018.pdf", "Mean $H$ of Gene", scargenes, essgenes)
plotGenewiseComplexity(geneindels, gposdict, allgenes_LCfrac, FIGPREFIX+"numindels_vs_fracLC_pergene_ALLgenes2_feb2018.pdf", "Fraction of Gene With Below Threshold $LC$", scargenes, essgenes)
plotGenewiseComplexity(geneindels, gposdict, allgenes_SEfrac, FIGPREFIX+"numindels_vs_fracSE_pergene_ALLgenes2_feb2018.pdf", "Fraction of Gene With Below Threshold $H$", scargenes, essgenes)

# plot SE/LC scores for the scar genes
plotComplexity(["Rv3902c"], gposdict, LCscores, SEscores, meanLC, meanSE, "essgene") # only do 60/461 ess genes
plotComplexity(["espI"], gposdict, LCscores, SEscores, meanLC, meanSE, "scargene") # only do 60/461 ess genes
#plotComplexity(uniq_scargenes, gposdict, LCscores, SEscores, meanLC, meanSE, "scargene")
