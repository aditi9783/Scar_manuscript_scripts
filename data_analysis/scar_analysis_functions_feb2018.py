#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 10})
from matplotlib.ticker import NullFormatter
import random
import numpy as np
from scipy.stats.stats import pearsonr
import scipy.stats as ss
#from sklearn.neighbors import KernelDensity
import math
import scipy.spatial.distance

# start of getH37Rv_proteins #################
def getH37Rv_proteins(fname): # read genpept file and extract names of protein coding and complement genes
    fh = open(fname, 'r')
    protnames = [] # names of protein coding genes only (no tRNA genes)
    compgenes = [] # complement genes
    flag = 1
    gname = ''
    locustag = ''
    for line in fh:
        if line.startswith("LOCUS"):
            flag = 1
        if flag == 1:
            if "/gene=" in line:
                content = line.split('="')
                gname = content[1].rstrip('"\n')
            elif "/locus_tag" in line:
                content = line.split('="')
                locustag = content[1].rstrip('"\n')
            elif "/coded_by" in line:
                gname_or_ltag = '' # get gname in format used at other places in this code
                if gname == '': # only locus tag present
                    gname_or_ltag = locustag
                else:
                    gname_or_ltag = gname
                if "complement" in line: # complement gene
                    compgenes.append(gname_or_ltag)
                protnames.append(gname_or_ltag)
                flag = 0
                gname = ''
                locustag = ''
    return protnames, compgenes
# end of getH37Rv_proteins #################

# start of getClusters ####################
def getClusters( fname, d ): 
    clusters = [] # 2d list, where each elt is a list of positions that have mutations (snp/indel) within the distance d
    fh = open(fname, 'r')
    allmutpos = [] # all mutation positions
    for line in fh:
        line = line.rstrip("\n")
        contents = line.split()
        pos = int(contents[0]) # get the position number
        allmutpos.append(pos)
    fh.close()

    prevpos = 0
    thiscluster = []
    orphans = [] # mutations that are not clustered with other muts
    intermutdist = [] # distance between two consecutive mutations
    sorted_ipos = sorted(allmutpos)
    for pos in sorted_ipos:
        mutdist = pos-prevpos
        intermutdist.append(mutdist)
        if mutdist < d: # the consecutive mut positions are within the threshold distance
            thiscluster.append(pos)
        else: # mutations are far apart
            if len(thiscluster) > 1: # more than one mutation was in a cluster
                clusters.append(thiscluster[:])
            elif len(thiscluster) == 1:
                orphans.extend(thiscluster)
            thiscluster = [pos] # re-initialize this cluster
        prevpos = pos
    return allmutpos, clusters, orphans, intermutdist
# end of getClusters ####################

# start of readScoreFile ###############
def readScoreFile(scorefile):
    sfh = open(scorefile, 'r')
    score = [] # index is the genome position-1 (thus, index 0 is pos 1), and value is the score for that position
    for line in sfh:
        if line.startswith("Mean"):
            break
        line = line.rstrip("\n")
        content = line.split()
        score.append(float(content[1]))
    sfh.close()
    return score
# end of readScoreFile #################

# start of getScores ###################
def getScores(pos, score, w, prefix): # return complexity scores for pos in mutpos
    posscore = []
    for p in pos:
        for i in range(p-1-w, p+w): # for positions in the window p-w to p+w)
            if i >= 0 and i < len(score): # i is within the genome
                posscore.append(score[i])
    avg = sum(posscore)/len(posscore)
    plt.hist(posscore, bins=50)
    plt.savefig(prefix+"_scoredist_feb2018.pdf")
    plt.close()
    return avg
# end of getScores ###################

# start of mutComplexScore ##############
def mutComplexScore( mutpos, score, w, prefix ): # read in the mutation positions and complexity score (SE or LC) and report scores for mutpos
    # get complexity score for the mutation positions
    avgmutscore = getScores(mutpos, score, w, prefix+"mutpos")
    nummut = len(mutpos)
    N = 5 # generate N random samples of same size as nummut from the genome and gather their avg complexity scores
    avgrandscore = []
    for i in range(N):
        randpos = random.sample(range(1, len(score)+1), nummut) 
        avgrandscore.append(getScores(randpos, score, w, prefix+"randpos"+str(i)))
    print "Average mut pos score:", avgmutscore
    print "Avg rand pos scores:", avgrandscore
    
# end of mutComplexScore ##############

# start of checkRegion ##########################
def checkRegion( pos, posrange ): # check if this position is in the posrange 
    flag = 0
    for tup in posrange:
        if pos >= tup[0] and pos <= tup[1]: # pos is within the repeat region 
            flag = 1 
            break
    return flag
# end of checkRegion ##########################

# start of coOccur ####################
def coOccur(poslist, posclusters): # see if the pos in list overlap the position range in posclusters
    posranges = [] # get the position ranges in the posclusters
    p_in_clusters = [] # positions that lie withing the regions defined by the clusters
    for c in posclusters:
        posranges.append([c[0], c[-1]])
    genomecovered = 0 # num of bases in the genome that fall within the clusters
    for pr in posranges:
         genomecovered += pr[1]-pr[0]
    #    print pr, pr[1]-pr[0]
    for p in poslist:
        flag = checkRegion(p, posranges)
        if flag == 1: # pos lies in the posranges)
            p_in_clusters.append(p)
    print "Genome covered by clusters:", genomecovered
    numpos = len(poslist)
    num_in_cluster = len(p_in_clusters)
    print "Total pos:", len(poslist), " Num in clusters:", len(p_in_clusters), float(num_in_cluster)/numpos
    
# end of coOccur ####################

# start of getGenePos ####################
def getGenePos(): # get gene start and end positions
    genes = [] # list of tuples of gene start, gene end, gene name such that start < end (can't identify complement genes)
    gposdict = {} # key: gene name, value: tuple of start and end pos for the gene.
    genefile = "genes_loci.txt" # get gene start and end positions for the reference genome H37Rv
    genefh = open(genefile, 'r') # H37Rv genes
    genicregion = [0 for i in range(0,4411529)] # positions that lie in a genic region will be assigned 1
    for line in genefh: # genes are already sorted by position in this genefile
        line = line.rstrip("\n")
        content = line.split()
        start = int(content[1])
        end = int(content[2])
        gname = content[3]
        gname_elts = gname.split("_")
        ltag = gname_elts.pop()
        genename = "_".join(gname_elts)

        if start < end:
            for gp in range(start, end+1):
                genicregion[gp] = 1
            if genename == "": # only locus tag for this gene
                gposdict[ltag] = [start, end]
                genes.append([ start, end, ltag ])
            else:
                gposdict[genename] = [start, end]
                genes.append([ start, end, genename ])
        else:
            for gp in range(end, start+1):
                genicregion[gp] = 1
            if genename == "": # only locus tag for this gene
                gposdict[ltag] = [end, start]
                genes.append([ end, start, ltag ])
            else:
                gposdict[genename] = [end, start]
                genes.append([ end, start, genename ])
    print "Total length of genic region:", sum(genicregion)
    genefh.close()
    return genes, gposdict
# end of getGenePos #####################

# start of checkGene ##########################
def checkGene( pos, genes ): # check if this position is in genic region
    glist = []
    for tup in genes:
        if pos >= tup[0] and pos <= tup[1]: # pos is within gene
            glist.append(tup[2])
    return glist
# end of checkGene ##########################

# start of getGenesLowComplexPos ##############
def getGenesLowComplexPos( genes, scores, thres ):
    genes_low_complex = [] # genes with positions that have low complexity
    for i in range(len(scores)):
        if scores[i] < thres: # score below thres. Low complexity pos.
            pos = i + 1 # since index starts from 0, genome pos is index + 1
            glist = checkGene( pos, genes )
            genes_low_complex.extend(glist)
    return set(genes_low_complex) # rendering a set removes duplicates
# end of getGenesLowComplexPos ##############

# start of writeToFile ###################
def writeToFile(ds, fname): # write list/set data structure to file
    fout = open(fname, 'w')
    for val in ds:
        fout.write(val+"\n")
    fout.close()
# end of writeToFile ###################

# start of getCommonGenes ##############
def getCommonGenes( set1, set2 ): # set 1 has gene name or locus tag, set 2 has genename_locus tag. 
    common = [] # list of common genes
    for g2 in set2:
        content = g2.split("_") # PE_PGRS_Rv gives three elements when split like that
        ltag = content.pop() # remove last element
        gname = "_".join(content) # join remaining elements.
        if gname in set1 or ltag in set1:
            common.append(g2)
    return common
# end of getCommonGenes ##############

# start of plotComplexity ############
def plotComplexity(glist, gpos, LCscore, SEscore, meanLC, meanSE, prefix):
    for gname in glist:
        start, end = gpos[gname]
        LCgene = [] 
        SEgene = []
        for i in range(start-1, end): # index in LC and SE scores is from 0, but genome pos starts from 1. Thus adjust for that.
            LCgene.append(LCscore[i])
            SEgene.append(SEscore[i])
        # plot the LC and SE scores for the gene, mainly to see if there are dips in the complexity
        avgLCgene = np.mean(np.array(LCgene))
        avgSEgene = np.mean(np.array(SEgene))
        f, axarr = plt.subplots(2, sharex=True)
        axarr[0].plot(range(0,len(LCgene)), LCgene, label="LC")
        axarr[1].plot(range(0,len(SEgene)), SEgene, label="SE")
        axarr[0].set_ylabel("LC")
        axarr[1].set_ylabel("H")
        axarr[0].set_title(gname)
        axarr[0].axhline(y=avgLCgene, linestyle='-', color='red')
        axarr[1].axhline(y=avgSEgene, linestyle='-', color='red')
        axarr[0].axhline(y=meanLC, linestyle='-', color='k')
        axarr[1].axhline(y=meanSE, linestyle='-', color='k')
        f.savefig(prefix+"_LC_SE_scores_"+gname+"_feb2018.pdf")
        plt.close()
# end of plotComplexity #############

# start of getEssGenes ##############
def getEssGenes(): # read DeJesus mBio 2017 data on gene function calls, and return the essential genes
    essgenes = []
    fname = "H37Rv_geneannotation_DeJesusMBIO2017.txt"
    fh = open(fname, 'r')
    for line in fh:
        if line.endswith("ES\n"): # essential gene
            content = line.split()
            gname = ""
            if content[1] == "-": # only locus tag present
                gname = content[0]
            else:
                gname = content[1]
            essgenes.append(gname)
    return essgenes
# end of getEssGenes ##############

# start of getGeneComplexity ######
def getGeneComplexity( gpos, LCscore, SEscore): # compute avg gene complexity for each gene
    LCavg = {} # key: gene name, val: avg LC score for that gene
    SEavg = {} # key: gene name, val: avg SE score for that gene
    LCstd = {} # key: gene name, val: stdev of LC score for that gene
    SEstd = {} # key: gene name, val: stdev of SE score for that gene
    for gname in gpos:
        start, end = gpos[gname]
        LCgene = []
        SEgene = []
        for i in range(start-1, end): # index in LC and SE scores is from 0, but genome pos starts from 1. Thus adjust for that.
            LCgene.append(LCscore[i])
            SEgene.append(SEscore[i])
        LCavg[gname] = np.mean(np.array(LCgene))
        SEavg[gname] = np.mean(np.array(SEgene))
        LCstd[gname] = np.std(np.array(LCgene))
        SEstd[gname] = np.std(np.array(SEgene))
    return LCavg, SEavg, LCstd, SEstd
# end of getGeneComplexity ######

# start of getLowComplexSpan #####
def getLowComplexFrac(gpos, scores, meanscore): # for each gene, determine the fraction of genome with score < meanscore
    genefrac = {} # key: gene name, val: frac of gene that have score less than the genome average
    for gname in gpos:
        start, end = gpos[gname]
        glen = end-start+1
        belowavgN = 0 # num positions with below average complexity
        for i in range(start-1, end): # index in LC and SE scores is from 0, but genome pos starts from 1. Thus adjust for that.
            if scores[i] < meanscore:
                belowavgN += 1 
        lowcomplexfrac = float(belowavgN)/glen
        genefrac[gname] = lowcomplexfrac
    return genefrac
# end of getLowComplexSpan #####

# start of getGenewiseIndels ####
def getGenewiseIndels(indelfile, gpos, genes):
    geneindel = {} # key: gene name, value: num indels in that gene
    geneindelpos = {} # key: gene name, value: list of genome positions that have indels in that gene
    for g in gpos: # initialize geneindel
        geneindel[g] = 0
        geneindelpos[g] = []
    num_genicindels = 0
    num_intergenicindels = 0

    fh = open(indelfile, 'r')
    numallindels = 0
    numindels_multiplegenes = 0 # sometimes same indel occurs in multiple genes due to overlapping gene boundried. Count how many times this happens
    peppe_indelpos = [] # list of all positions in pe ppe genes that have indels
    for line in fh:
        line = line.rstrip("\n")
        contents = line.split()
        ipos = int(contents[0]) 
        numallindels += 1
        glist = checkGene(ipos, genes) # check which genes have this indel position
        if len(glist) > 0: # there is at least one gene that has this position
            num_genicindels += 1
            numindels_multiplegenes += len(glist)-1
            for gn in glist:
                geneindel[gn] += 1 # add num of indels to this gene
                geneindelpos[gn].append(ipos)
                if gn.startswith("PE") or gn.startswith("PPE"):
                    peppe_indelpos.append(ipos)    
        else: # intergenic indels
            num_intergenicindels += 1 
    print "Total number of unique indels:", numallindels
    print "Num unique genic indels: ", num_genicindels, "\tNum unique intergenic indels: ", num_intergenicindels
    print "Num times an indel overcounted in absolute count (not uniq counts) due to occurence in multiple genes:", numindels_multiplegenes
    uniq_peppe_indelpos = set(peppe_indelpos)
    print "Num unique indels in pe ppe genes:", len(uniq_peppe_indelpos)
    return geneindel, geneindelpos
# end of getGenewiseIndels ####

# start of plotScatter ########
def plotScatter(xarr, yarr, ckeyarr, alphaarr, xlabeltext, ylabeltext, figname):
    f, ax = plt.subplots()
    ax.margins(0.05)
    for i in range(len(xarr)):
        #ax.scatter(scar_x, scar_y, c='cyan', marker="o", s=slist, lw=0, label="Scar genes") # 
        ax.scatter(xarr[i], yarr[i], c=ckeyarr[i], marker="o", lw=0, alpha=alphaarr[i]) # 
    ax.set_xlabel(xlabeltext)
    ax.set_ylabel(ylabeltext)
    if xlabeltext != "Gene Length":
        ax.set_ylim([0.0,1.0])
        ax.set_xlim([-2, 90])
    f.savefig(figname)
    plt.close()
# end of plotScatter ########

# start of plotGenewiseComplexity ###
def plotGenewiseComplexity(gindels, gposdict, gfeature, figname, ylabeltext, scargenes, essgenes):
    x = [] # num indels
    y = [] # feature
    mkey = []
    peppe_x = []
    peppe_y = []
    scar_x = []
    scar_y = []
    ess_x = []
    ess_y = []
    allbutppe_x = []
    allbutppe_y = []
    noness_nonppe_nonscar_x = []
    noness_nonppe_nonscar_y = []
    peppelen_x = []
    peppelen_y = []
    esslen_x = []
    esslen_y = []
    nonppe_noness_len_x = []
    nonppe_noness_len_y = []
    for g in gposdict:
        if g.startswith("PE") or g.startswith("PPE"):
            glen = abs(gposdict[g][1] - gposdict[g][0])
            #x.append(glen)
            #x.append(gindels[g])
            #y.append(gindels[g])
            #y.append(gfeature[g])
            peppe_x.append(gindels[g])
            peppe_y.append(gfeature[g])
            peppelen_y.append(gindels[g])
            peppelen_x.append(glen)
        elif g in scargenes:
            glen = abs(gposdict[g][1] - gposdict[g][0])
            #x.append(glen)
            x.append(gindels[g])
            #y.append(gindels[g])
            y.append(gfeature[g])
            scar_x.append(gindels[g])
            scar_y.append(gfeature[g])
            #scar_y.append(gindels[g])
            #scar_x.append(glen)
            #allbutppe_x.append(gindels[g])
            #allbutppe_y.append(gfeature[g])
            nonppe_noness_len_y.append(gindels[g])
            nonppe_noness_len_x.append(glen)
        elif g in essgenes:
            glen = abs(gposdict[g][1] - gposdict[g][0])
            #x.append(glen)
            #x.append(gindels[g])
            #y.append(gindels[g])
            #y.append(gfeature[g])
            ess_x.append(gindels[g])
            ess_y.append(gfeature[g])
            esslen_y.append(gindels[g])
            esslen_x.append(glen)
            #allbutppe_x.append(gindels[g])
            #allbutppe_y.append(gfeature[g])
            #if "LC" in ylabeltext and gfeature[g] > 0.7:
            #    print g, gfeature[g]
        else:
            glen = abs(gposdict[g][1] - gposdict[g][0])
            x.append(gindels[g])
            y.append(gfeature[g])
            # only select 20% of non ess, non PE PPE, non scar genes for less cluttered plot
            if random.uniform(0,1) <= 0.2:
                noness_nonppe_nonscar_x.append(gindels[g])
                noness_nonppe_nonscar_y.append(gfeature[g])
                nonppe_noness_len_y.append(gindels[g])
                nonppe_noness_len_x.append(glen)
            allbutppe_x.append(gindels[g])
            allbutppe_y.append(gfeature[g])

    print "\nAll genes except essential and PE PPE:", len(x), " Pearson corr and pval between ", ylabeltext, " and num indels/gene:", pearsonr(x,y)
    print "PE PPE genes:", len(peppe_x), " Pearson corr and pval between ", ylabeltext, " and num indels/gene:", pearsonr(peppe_x,peppe_y)
    #print "Scar genes:", len(scar_x), " Pearson corr and pval between ", ylabeltext, " and num indels/gene:", pearsonr(scar_x,scar_y)
    print "Essential genes:", len(ess_x), " Pearson corr and pval between ", ylabeltext, " and num indels/gene:", pearsonr(ess_x,ess_y)
    #print "All genes except for PE PPE genes:", len(allbutppe_x), " Pearson corr and pval between ", ylabeltext, " and num indels/gene:", pearsonr(allbutppe_x,allbutppe_y)
    print "PE PPE genes:", len(peppelen_x), " Pearson corr and pval between genelen and num indels/gene:", pearsonr(peppelen_x,peppelen_y)
    print "Essential genes:", len(esslen_x), " Pearson corr and pval between genelen and num indels/gene:", pearsonr(esslen_x,esslen_y)
    print "Non PPE Non Essential genes:", len(nonppe_noness_len_x), " Pearson corr and pval between genelen and num indels/gene:", pearsonr(nonppe_noness_len_x,nonppe_noness_len_y)
    # Welch t-test to compare population means
    print " Number of noness nonPPE and non scar genes randomly selected: ", len(noness_nonppe_nonscar_y)
    meantest = ss.ttest_ind(scar_y,allbutppe_y,axis=None,equal_var=False)
    print "Welch t-test (unequal variances) for comparing the diff in -", ylabeltext, "- between ", len(allbutppe_y), " scar and (non ess non PPE non scar), t statistic and p value:", meantest
    
    xarr = [ess_x, peppe_x, x]
    yarr = [ess_y, peppe_y, y]
    ckeyarr = ['b', 'r', 'g']
    alphaarr = [0.8, 0.8, 0.1]
    xlabeltext = "Number of Unique Indels in Gene"
    plotScatter(xarr, yarr, ckeyarr, alphaarr, xlabeltext, ylabeltext, figname)
    xarr2 = [scar_x, allbutppe_x]
    yarr2 = [scar_y, allbutppe_y]
    ckeyarr2 = ['k', 'g']
    alphaarr2 = [1.0, 0.1]
    plotScatter(xarr2, yarr2, ckeyarr2, alphaarr2, xlabeltext, ylabeltext, "scargenes_"+figname)
    xarr3 = [ess_x, peppe_x, noness_nonppe_nonscar_x[0:500]]
    yarr3 = [ess_y, peppe_y, noness_nonppe_nonscar_y[0:500]]
    ckeyarr3 = ['b', 'r', 'g']
    alphaarr3 = [0.8, 0.8, 0.4]
    plotScatter(xarr3, yarr3, ckeyarr3, alphaarr3, xlabeltext, ylabeltext, "randomsample_500_nonppe_noness_nonscar_genes_"+figname)
    xarr4 = [esslen_x, peppelen_x, nonppe_noness_len_x[0:500]]
    yarr4 = [esslen_y, peppelen_y, nonppe_noness_len_y[0:500]]
    ckeyarr4 = ['b', 'r', 'g']
    alphaarr4 = [0.8, 0.8, 0.4]
    plotScatter(xarr4, yarr4, ckeyarr4, alphaarr4, "Gene Length", "Number of Unique Indels in Gene", "allindels_genelen_vs_numuniqindels_peppe_ess_500rest_ALLgenes2_feb2018.pdf")
    return
# end of plotGenewiseComplexity ###

# start of plotIndelDist_fig1a #########
def plotIndelDist_fig1a(indelfile, genes, geneindels, essgenes, peppegenes, FIGPREFIX): # all ess and peppe, not just protein coding ones
    # FIG1 : distribution of insertions and deletions based on indel size
    geneinslen = {}
    genedellen = {}
    for gtype in ["ess", "peppe", "rest", "intergenic"]:
        geneinslen[gtype] = [0 for i in range(0,6)] # number of insertions of size given by the index
        genedellen[gtype] = [0 for i in range(0,6)] # number of deletions of size given by the index
    insertions = [0 for i in range(0,6)] # number of insertions of size given by the index
    deletions = [0 for i in range(0,6)] # number of deletions of size given by the index
    largeins = 0 # number of insertions of size > 5
    largedel = 0 # number of deletions of size > 5
    fh = open(indelfile, 'r')
    for line in fh:
        line = line.rstrip("\n")
        contents = line.split()
        ipos = int(contents[0])
        glist = checkGene(ipos, genes) # check which genes have this indel position
        indelsize = int(contents[2])

        if contents[1] == "+": # insertion
            if indelsize < 6:
                if len(glist) > 0: # genic indel
                    for gn in glist:
                        if gn in essgenes:
                            geneinslen["ess"][indelsize] += 1
                        elif gn in peppegenes:
                            geneinslen["peppe"][indelsize] += 1
                        else:
                            geneinslen["rest"][indelsize] += 1
                else: # intergenic indel
                    geneinslen["intergenic"][indelsize] += 1
            else:
                largeins += 1
        elif contents[1] == "-": # deletion
            if indelsize < 6:
                if len(glist) > 0: # genic indel
                    for gn in glist:
                        if gn in essgenes:
                            genedellen["ess"][indelsize] += 1
                        elif gn in peppegenes:
                            genedellen["peppe"][indelsize] += 1
                        else:
                            genedellen["rest"][indelsize] += 1
                else: # intergenic indel
                    genedellen["intergenic"][indelsize] += 1
            else:
                largedel += 1
    fh.close()
    xlabels = ["> -5", "-5", "-4", "-3", "-2", "-1", "0", "+1", "+2", "+3", "+4", "+5", "> +5"]
    yvals = {}
    print "rest:"
    print genedellen["rest"], geneinslen["rest"]
    for gtype in ["ess", "peppe", "rest", "intergenic"]:
        yvals[gtype] = [0] # value for >5 len deletion
        for i in range(1,len(genedellen[gtype])):
            yvals[gtype].append(genedellen[gtype][len(genedellen[gtype])-i])
        yvals[gtype].append(0) # append the value for 0-size deletions
        for i in range(1,len(geneinslen[gtype])):
            yvals[gtype].append(geneinslen[gtype][i])
        yvals[gtype].append(0) # append the value for >5 size insertions
    yvallen = len(yvals["ess"])
    yvals_extremes = [largedel]
    for i in range(0, yvallen-2):
        yvals_extremes.append(0)
    yvals_extremes.append(largeins)
    print xlabels
    cscheme = {}
    for gtype in ["ess", "peppe", "rest", "intergenic"]:
        print gtype, yvals[gtype]
    print "extreme:", yvals_extremes
    f, ax = plt.subplots()
    plt.bar(range(len(xlabels)), yvals["intergenic"], align='center', color='gray', alpha=0.8, label="Intergenic") # 0s at start and end for the extreme values
    plt.bar(range(len(xlabels)), yvals["rest"], bottom=np.array(yvals["intergenic"]), align='center', color='green', alpha=0.8, label="NENP")
    plt.bar(range(len(xlabels)), yvals["peppe"], bottom=np.array(yvals["rest"])+np.array(yvals["intergenic"]), align='center', color='red', alpha=0.8, label="PE PPE")
    plt.bar(range(len(xlabels)), yvals["ess"], bottom=np.array(yvals["peppe"])+np.array(yvals["rest"])+np.array(yvals["intergenic"]), align='center', color='blue', alpha=0.8, label="Essential")
    plt.bar(range(len(xlabels)), yvals_extremes,  align='center', color='black', alpha=0.8, label="Indels >5 nt")
    plt.xticks(range(len(xlabels)), xlabels) 
    plt.xlabel("Indel Length")
    #plt.xlabel("(deletions)\t\t\tIndel Length\t\t\t(insertions)")
    plt.ylabel("Frequency")
    plt.legend()
    #f.text(0.5, 0.04, '(deletions)\t\t(insertions)', ha='center')
    plt.tight_layout()
    f.savefig(FIGPREFIX+"_indellendist_feb2018_2.pdf")
    plt.close()

# end of plotIndelDist_fig1a ###

# start of plotIndelDist #########
def plotIndelDist(indelfile, geneindels, essgenes, peppegenes, FIGPREFIX): # all ess and peppe, not just protein coding ones
    # FIG1 : distribution of insertions and deletions based on indel size
    insertions = [0 for i in range(0,6)] # number of insertions of size given by the index
    largeins = 0 # number of insertions of size > 5
    deletions = [0 for i in range(0,6)] # number of deletions of size given by the index
    largedel = 0 # number of deletions of size > 5
    fh = open(indelfile, 'r')
    for line in fh:
        line = line.rstrip("\n")
        contents = line.split()
        indelsize = int(contents[2])
        if contents[1] == "+": # insertion
            if indelsize < 6:
                insertions[indelsize] += 1
            else:
                largeins += 1
        elif contents[1] == "-": # deletion
            if indelsize < 6:
                deletions[indelsize] += 1
            else:
                largedel += 1
    fh.close()
    xlabels = ["> -5", "-5", "-4", "-3", "-2", "-1", "0", "+1", "+2", "+3", "+4", "+5", "> +5"]
    yvals = [largedel]
    for i in range(1,len(deletions)):
        yvals.append(deletions[len(deletions)-i])
    yvals.append(0) # append the value for 0-size deletions
    for i in range(1, len(insertions)):
        yvals.append(insertions[i])
    yvals.append(largeins)
    print "deletions:", deletions
    print "ins:", insertions
    print xlabels
    print yvals
    f, ax = plt.subplots()
    plt.bar(range(len(xlabels)), yvals, align='center', alpha=0.5)
    plt.xticks(range(len(xlabels)), xlabels) 
    plt.xlabel("Indel Length")
    plt.ylabel("Frequency")
    #f.text(0.5, 0.04, '(deletions)\t\t(insertions)', ha='center')
    plt.tight_layout()
    f.savefig(FIGPREFIX+"_indellendist_feb2018.pdf")
    plt.close()

    # FIG 2: Frequency of indels per gene in essential, PE PPE, and the rest of the genes
    maxindelnum = 10 # max num indels/gene that I want to plot
    essgeneindels = [0 for i in range(0,maxindelnum+1)]
    peppegeneindels = [0 for i in range(0,maxindelnum+1)]
    restgeneindels = [0 for i in range(0,maxindelnum+1)]
    ess_manyidl = 0 # num of ess genes that has more than maxindelnum indels
    peppe_manyidl = 0 # num of peppe genes that has more than maxindelnum indels
    rest_manyidl = 0 # num of rest genes that has more than maxindelnum indels
    peppe_and_ess = set(essgenes).intersection(set(peppegenes))
    print "Num of genes in geneindels: ", len(geneindels)
    print "Number of indels in pe ppe genes that are also essential genes:"
    print "Num of PE PPE genes ", len(peppegenes)
    #for g in peppe_and_ess:
    #    print g, geneindels[g]
    for g in geneindels:
        if geneindels[g] <= maxindelnum:
            if g in peppegenes:
                peppegeneindels[geneindels[g]] += 1
            elif g in essgenes:
                essgeneindels[geneindels[g]] += 1
            else:
                restgeneindels[geneindels[g]] += 1
        else:
            if g in peppegenes:
                peppe_manyidl += 1
            elif g in essgenes:
                ess_manyidl += 1
                print "Essential genes with ", geneindels[g], " indels:", g 
            else:
                rest_manyidl += 1
    print "Num_idl Ess PE_PPE Rest -- Num genes"
    for i in range(maxindelnum+1):
        print i, essgeneindels[i], peppegeneindels[i], restgeneindels[i]
    print "Num of peppe, ess, and rest of the genes that have indels more than ", maxindelnum, ":", peppe_manyidl, ess_manyidl, rest_manyidl
    f, ax = plt.subplots()
    barw = 0.2 # bar width
    plt.bar([i for i in range(0,maxindelnum+1)],[float(v)/459 for v in essgeneindels], width=barw, color='b', align='center', label="Essential") 
    plt.bar([i-barw for i in range(0,maxindelnum+1)],[float(v)/165 for v in peppegeneindels], width=barw, color='r', align='center', label="PE PPE")
    plt.bar([i+barw for i in range(0,maxindelnum+1)],[float(v)/(4110-459-165) for v in restgeneindels], width=barw, color='g', align='center', label="NENP")
    plt.xticks(range(0,maxindelnum+1))
    plt.xlabel("Num Indels")
    plt.ylabel("Fraction of Genes")
    plt.legend(fontsize=10)
    f.savefig(FIGPREFIX+"_numindels_vs_genetype_frac_feb2018.pdf")
    plt.close()
# end of plotIndelDist ###########

# start of plotIndelQuad ##########
def plotIndelQuad(geneindelpos, gposdict, compgenes, essgenes, peppegenes, scargenes, FIGPREFIX): # ess and peppe here are protein coding only
    nquad = 10.0 # num of quadrants
    quadb = 1.0/nquad # boundries for each quad
    quads = [quadb*i for i in range(0, int(nquad)+1)]
    numquads = len(quads)  
    essidlquads = [0 for i in range(numquads)] # initialize num indels in gene quad for ess genes
    peppeidlquads = [0 for i in range(numquads)] # initialize num indels in gene quad for peppe genes
    restidlquads = [0 for i in range(numquads)] # initialize num indels in gene quad for rest genes
    checkgene = "fadD32"
    print "for ", checkgene, ": gstart, gend, and indel pos:", gposdict[checkgene], "\n", geneindelpos[checkgene]
    for g in geneindelpos:
        if len(geneindelpos[g]) > 0: # at least one indel in this gene
            gstart, gend = gposdict[g]
            glen = float(abs(gstart-gend))
            for ipos in geneindelpos[g]:
                posfrac = float(ipos-gstart)/glen # subtract gstart to get indel position in the gene index not genome index
                posquad = 0
                if g == checkgene:
                    print ipos, ipos-gstart, posfrac
                for i in range(0,numquads):
                    #print g, scar, posfrac, i, quads[i]
                    if g == checkgene:
                        print "\t", quads[i]
                    if posfrac > quads[i] and posfrac <= quads[i+1]: # the quad that is just below the position fraction
                        if g in compgenes: # complement gene
                            posquad = numquads - i - 2 # -2 because quads are numbered 0-5 and not 1-6
                        else:
                            posquad = i
                        if g == checkgene:
                            print "\t\tposfrac just above quad", i
                            print "posquad: ", posquad
                        break
                if g in essgenes:
                    essidlquads[posquad] += 1
                    if posquad < 2:
                        print "Ess gene: ", g, ipos, ipos-gstart, posquad, gposdict[g], glen
                elif g in peppegenes:
                    peppeidlquads[posquad] += 1
                #elif g in scargenes:
                else:
                    restidlquads[posquad] += 1
    print "ess genes indel quads:", essidlquads, sum(essidlquads)
    print "peppe genes indel quads:", peppeidlquads, sum(peppeidlquads)
    print "rest genes indel quads:", restidlquads, sum(restidlquads)
    #print "scar genes indel quads:", restidlquads, sum(restidlquads)
    quad_xticks = []
    prev = "0.0"
    for q in quads:
        quad_xticks.append(prev+"-"+str(q))
        prev = str(q)
    f, ax = plt.subplots()
    ax.margins(0.05)
    barw = 0.2 # bar width
    #plt.bar([i for i in range(numquads)],[float(v)/sum(essidlquads) for v in essidlquads], width=barw, color='g', align='edge', label="Essential") 
    #plt.bar([i-barw for i in range(numquads)],[float(v)/sum(peppeidlquads) for v in peppeidlquads], width=barw, color='r', align='edge', label="PE PPE")
    #plt.bar([i+barw for i in range(numquads)],[float(v)/sum(restidlquads) for v in restidlquads], width=barw, color='b', align='edge', label="Rest")
    plt.plot([i for i in range(numquads-1)],[float(v)/sum(essidlquads) for v in essidlquads[:-1]], 'ob-', label="Essential") 
    plt.plot([i for i in range(numquads-1)],[float(v)/sum(peppeidlquads) for v in peppeidlquads[:-1]], 'or-', label="PE PPE")
    plt.plot([i for i in range(numquads-1)],[float(v)/sum(restidlquads) for v in restidlquads[:-1]], 'og-', label="NENP")
    #plt.plot([i for i in range(numquads-1)],[float(v)/sum(restidlquads) for v in restidlquads[:-1]], 'oc-', label="Scar")
    plt.xticks([i for i in range(numquads-1)],quad_xticks[1:], fontsize=10, rotation=30)
    plt.xlabel("Gene Quadrant")
    plt.ylabel("Fraction of Indels in a Quadrant")
    plt.legend(ncol=3, loc=2, fontsize=10)
    plt.tight_layout()
    f.savefig(FIGPREFIX+"_ess_PEPPE_remaining_genes_indels_vs_genequadrant_feb2018.pdf")
    plt.close()
            
# end of plotIndelDist #########

# start of indelComplexityScore ####
def indelComplexityScore(indelpos, nonidlpos, scores, ytext, iternum, filetype):
    idlscores = [scores[i-1] for i in indelpos if i < len(scores)] # complexity scores for indel positions
    nonidlscores = [scores[i-1] for i in nonidlpos] # complexity scores for the non indel positions
    #print "score for first 10 sites :", nonidlscores[0:10]
    print "Mean score for indel pos:", sum(idlscores)/len(idlscores)
    print "Mean score for nonindel pos:", sum(nonidlscores)/len(nonidlscores)
    print "Num indel pos:", len(idlscores), " Num non indel pos:", len(nonidlscores)
    meantest = ss.ttest_ind(idlscores,nonidlscores,axis=None,equal_var=False)
    print "Welch t-test (unequal variances) for comparing the diff in means, t statistic and p value:", meantest
    #plt.close() # close any opened plot windows
    plt.hist(idlscores, bins=50, color='b', normed=True, alpha=0.5, label="Indel")
    plt.hist(nonidlscores, bins=50, color='r', normed=True, alpha=0.5, label="Non-indel")
    plt.xlabel(ytext+" scores")
    plt.legend()
    plt.savefig(filetype+"_"+ytext+"_iter"+str(iternum)+"_indel_nonindl_scores_hist_feb2018.pdf")
    plt.close()
    return
    X1 = np.array(idlscores)[:, np.newaxis]
    X2 = np.array(nonidlscores)[:, np.newaxis]
    X_plot = np.linspace(-2, 2, 500)[:, np.newaxis]

    fig, ax = plt.subplots()
    kde1 = KernelDensity(kernel='gaussian', bandwidth=0.05).fit(X1)
    kde2 = KernelDensity(kernel='gaussian', bandwidth=0.05).fit(X2)
    log_dens1 = kde1.score_samples(X_plot)
    log_dens2 = kde2.score_samples(X_plot)
    ax.fill(X_plot[:, 0], np.exp(log_dens1), fc='blue', alpha=0.3, label="Indel")
    ax.fill(X_plot[:, 0], np.exp(log_dens2), fc='red', alpha=0.3, label="Non-indel")
    #ax.text(-3.5, 0.31, "Gaussian Kernel Density")

    ax.plot(X1[:, 0], np.zeros(X1.shape[0]) -0.03, '+b', markersize=4)
    ax.plot(X2[:, 0], np.zeros(X2.shape[0]) -0.1, '+r', markersize=4)
    if ytext == "H":
        ax.set_xlim(0.5,1.2)
    else:
        ax.set_xlim(-0.2,1.2)
    #ax.set_ylim(-0.02, 0.34)
    ax.set_ylabel('Normalized Density')
    ax.set_xlabel(ytext+" scores")
    plt.legend()
    fig.savefig(filetype+"_"+ytext+"_iter"+str(iternum)+"_indel_nonindel_density_feb2018.pdf")

# end of indelComplexityScore ####

# start of plotComplexNbg ####
def plotComplexNbg(indelorphans, scores, nbgpos, ytext, scoretype):
    nbgscores = [] # list of lists, each row is scores in neighboring region of an indel
    fout = open(scoretype+"_w"+str(nbgpos)+"_indeltrap_scores.csv", 'w')
    windowrange = [str(x) for x in range(-nbgpos, nbgpos+1)]
    #print windowrange
    fout.write(",".join(windowrange)+"\n")
    for ipos in indelorphans:
        wleft = ipos-nbgpos
        wright = ipos+nbgpos
        thisipos = [] # complexity scores in window of this indel pos
        for i in range(wleft-1, wright): # pos 1 has index 0 in scores, thus wleft-1 till wright
            thisipos.append(scores[i])
        #print ipos, scores[ipos-1], thisipos
        nbgscores.append(thisipos[:])
        thisipos_str = [str(x) for x in thisipos]
        fout.write(",".join(thisipos_str)+"\n")
    fout.close()
    # compute avg and sd of score at each pos in the window
    nbgscore_arr = np.array(nbgscores)
    colwisemean = np.mean(nbgscore_arr, axis=0)
    colwisestd = ss.sem(nbgscore_arr, axis=0)
    #colwisestd = np.std(nbgscore_arr, axis=0)
    #for i in range(len(colwisemean)):
    #    if i+25 > len(colwisemean)-1:
    #        break
    #    print i, colwisemean[i], sum(colwisemean[i:i+25])/len(colwisemean[i:i+25])
    
    # indels have low scores from 8 bases upstream to 14 bases downstream of an indel pos. get average score for this window
    avgwindowscore = sum(colwisemean[nbgpos-8:nbgpos+15])/len(colwisemean[nbgpos-8:nbgpos+15])
    
    print "Avg score in the itrap region: ", scoretype, avgwindowscore
    xvals = range(-nbgpos, nbgpos+1) 
    #print xvals
    plt.figure()
    plt.errorbar(xvals, colwisemean, yerr=colwisestd, mfc='k', mec='k')
    if scoretype.endswith("H"):
        plt.ylim([0.930,0.945])
    elif scoretype.endswith("LC"):
        plt.ylim([0.54,0.601])
    plt.ylabel(ytext)
    plt.savefig("Indelpos_neighbors_"+scoretype+"_scores_sem_feb2018.pdf")
    plt.close()
    return colwisemean
# end of plotComplexNbg ####

# start of absDist ###
def absDist(pos, scores, itrap): # find sum of absolute differences between each elt in vectors
    dist = []
    ilen = len(itrap)  
    flank = ilen/2 # region on both sides of a pos that should be looked at
    print "indeltrap len: ", ilen
    for p in pos:
        #print p
        vec2 = scores[p-flank-1:p+flank] # -1 bec p is from index 1 and scores is from index 0
        if len(vec2) == len(itrap): # some positions at the end may not have 20 bases flanking region
            d = 0.0
            for i in range(len(itrap)):
                d += abs(itrap[i] - vec2[i])
                #print itrap[i], vec2[i], abs(itrap[i] - vec2[i])
            dist.append(d)
            #dist.append(scipy.spatial.distance.mahalanobis(itrap, vec2))
    return dist 
# end of absDist ###

# start of matchPattern ###
def matchPattern(pos, scores, itrap):
    itrap_pattern = []
    for i in range(1,len(itrap[1:])):
        if itrap[i] > itrap[i-1]:
            itrap_pattern.append("i") # "i"ncrease in LC/H value from the prev pos
        elif itrap[i] < itrap[i-1]:
            itrap_pattern.append("d") # "d"ecrease in LC/H value from the prev pos
        else:
            itrap_pattern.append("s") # "s"ame LC/H value as in prev pos
    #print itrap
    #print len(itrap_pattern), itrap_pattern
    upstream_flank = 8
    downstream_flank = 17
    counter = 0

    dist = [] # number of positions at which the pattern does not match the vector

    for p in pos:
        counter += 1
        vec = scores[p-upstream_flank-1:p+downstream_flank-1]
        vec_pattern = []
        for i in range(1,len(vec[1:])):
            if vec[i] > vec[i-1]:
                vec_pattern.append("i") # "i"ncrease in LC/H value from the prev pos
            elif vec[i] < vec[i-1]:
                vec_pattern.append("d") # "d"ecrease in LC/H value from the prev pos
            else:
                vec_pattern.append("s") # "s"ame LC/H value as in prev pos
        #print len(vec), vec
        #print itrap_pattern
        #print vec_pattern
        d = 0
        for i in range(len(vec_pattern)):
            if vec_pattern[i] == "s":
                continue
            elif vec_pattern[i] == itrap_pattern[i]:
                continue
            else:
                d += 1
        dist.append(d)
        if counter == 1000:
            return dist
# end of matchPattern ###

# start of predictIndel ###
def predictIndel(indelpos, indelorphans, nonidlpos, scores, indeltrap, scoretype):
    ipos_TP = [] # indel positions that are correctly predicted by indeltrap
    iorphan_TP = [] # indel orphans that are correctly predicted by indeltrap. This is more imp because indeltrap came from this data.
    #matchPattern(indelorphans, scores, indeltrap)
    #return

    iorphan_dist = matchPattern(indelorphans, scores, indeltrap) # distance between indeltrap and windows with indelorphan in center
    ipos_dist = matchPattern(indelpos, scores, indeltrap) # distance between indeltrap and windows with indel pos in center
    nonidl_dist = matchPattern(nonidlpos, scores, indeltrap) # distance between indeltrap and windows with nonidl pos in center
    #iorphan_dist = absDist(indelorphans, scores, indeltrap) # distance between indeltrap and windows with indelorphan in center
    plt.hist(iorphan_dist, bins=20, histtype='stepfilled', normed=True, alpha=0.3, label="iorphan")
    #plt.hist(ipos_dist, bins=20, histtype='stepfilled', normed=True, alpha=0.3, label="ipos")
    plt.hist(nonidl_dist, bins=20, histtype='stepfilled', normed=True, alpha=0.3, label="nonidl")
    #plt.title("iorphan_dist")
    plt.legend()
    plt.savefig(scoretype+"_predictIndel_dist_hist.pdf")
    plt.close()
    return

    #ipos_dist = absDist(indelpos, scores, indeltrap) # distance between indeltrap and windows with indel pos in center
    plt.hist(ipos_dist, bins=50)
    plt.title("ipos_dist")
    plt.savefig("ipos_dist_hist.pdf")
    plt.close()

    #nonidl_dist = absDist(nonidlpos, scores, indeltrap) # distance between indeltrap and windows with nonidl pos in center
    plt.hist(nonidl_dist, bins=100)
    plt.title("nonidl_dist")
    plt.savefig("nonidl_dist_hist.pdf")
    plt.close()
# end of predictIndel ###

# start of getScarIndels ###############
def getScarIndels(gposdict, scargenes, compgenes): # get scar indels
    nquad = 10.0 # num of quadrants
    quadb = 1.0/nquad # boundries for each quad
    quads = [quadb*i for i in range(0, int(nquad)+1)]
    numquads = len(quads)  
    scargenequad = [0 for i in range(numquads-1)] # initialize num indels in gene quad for scar genes
    scaridldist = []
    scarindels = [] # each row represents gene quadrants for scars with 2 indels

    fh = open("all_confirmed_scars.txt", 'r')
    gene_indelpos = {} # gene name is key, value is list of gene positions (not genome pos) that have scar indels 

    header = fh.readline()
    for line in fh:
        content = line.split('\t')
        gname = content[0]
        fullname = ""
        for g in scargenes:
            if gname in g:
                fullname = g
        # get numbers (gene positions) for each scar
        pos = [int(s) for s in content[1].split() if s.isdigit()]
        gstart, gend = gposdict[fullname]
        glen = float(abs(gstart-gend))
        compflag = 0
        if fullname in compgenes:
            compflag = 1
        if len(pos) == 2: # only consider scars that have two indels
            scaridldist.append(abs(pos[0]-pos[1]))
            #print fullname, glen, pos, compflag
            thisscar = []
            for ipos in pos:
                posfrac = float(ipos)/glen # scar indels are noted by indel pos in gene, not genome 
                posquad = 0
                for i in range(0,numquads):
                    #print fullname, scar, posfrac, i, quads[i]
                    if posfrac > quads[i] and posfrac <= quads[i+1]: # the quad that is just below the position fraction
                        posquad = i # in all_confirmed_scars.txt, ipos is wrt translation start site in complement genes
                        thisscar.append(posquad)
                        break
                scargenequad[posquad] += 1
                if fullname == "Rv2020c":
                    print fullname, glen, ipos, posquad
            scarindels.append(sorted(thisscar))
    fh.close()
    print "scar gene quads: ", scargenequad
    print "distances between scar indels: ", scaridldist
    sorted_scarindels = sorted(scarindels) #,key=lambda x: (x[0],x[1])) # scar indel quadrants sorted by the quad of 1st indel in gene

    nindels = len(sorted_scarindels)
   
    quad_xticks = []
    prev = "0.0"
    for q in quads:
        quad_xticks.append(prev+"-"+str(q))
        prev = str(q)
    f, ax = plt.subplots()
    ax.margins(0.05)
    xs = np.arange(numquads)
    for si in sorted_scarindels:
        #print si, nindels
        fullarr = [None for i in range(numquads)]
        for posquad in si:
            fullarr[posquad] = nindels
        nindels -= 1
        nparr = np.array(fullarr).astype(np.double)
        mask = np.isfinite(nparr)

        plt.plot(xs[mask], nparr[mask], marker='o', color='cyan')
    plt.xticks([i for i in range(numquads-1)],quad_xticks[1:], fontsize=10)
    plt.xlabel("Gene Quadrant")
    plt.ylabel("Indels")
    f.savefig("scar_wise_indels_in_gene_quad_feb2018.pdf")
    plt.close()
    
    bins = range(1,371,5)
    hist_distance, bins_distance = np.histogram(scaridldist, bins)
    print len(hist_distance), len(bins_distance)
    f, ax = plt.subplots()
    ax.margins(0.04, 0)
    plt.bar(bins[:-1], hist_distance, width = 1.0, align='center')
    #plt.bar(bins[:-1], [math.log(x, 5) for x in hist_distance], width = 1.5, align='center')
    plt.ylim([0,25])
    plt.xlabel("Distance between scar indels (nt)")
    plt.ylabel("Frequency")
    ax.xaxis.set_ticks_position('none')
    f.savefig("distance_bet_scarindels_feb2018.pdf")
    plt.close()

    f, ax = plt.subplots()
    ax.margins(0.05)
    barw = 0.4 # bar width
    plt.plot([i for i in range(len(scargenequad))],[float(v)/sum(scargenequad) for v in scargenequad], 'oc-') 
    plt.xticks([i for i in range(numquads-1)],quad_xticks[1:], fontsize=10)
    plt.xlabel("Gene Quadrant")
    plt.ylim([0.0,0.5])
    plt.ylabel("Fraction of Indels in a Quadrant")
    f.savefig("scarindels_genequadrant_feb2018.pdf")
    plt.close()

# end of getScarIndels ###############

# start of getNonIndelPos ############
def getNonIndelPos(genomeL, indelpos, numpos, d): # return numpos non-indel pos that are at least d nt away from indel pos on both sides
    nonidlpos = []
    allgenomepos = [i+1 for i in range(genomeL)]
    randpos = random.sample(allgenomepos, numpos) 
    for r in randpos:
        flag = 0
        for ipos in indelpos:
            if r > ipos-d and r < ipos+d: # the randomly selected genome position is within dist bases of an indel pos, don't use it 
                flag = 1
                break
        if flag == 0:
            nonidlpos.append(r)
    return nonidlpos
# end of getNonIndelPos ############

# start of compareIndelNonindel #######
def compareIndelNonindel(LCscores, SEscores, d, indelpos, FIGPREFIX): # randomly sample nonindel pos and compare to indel pos
    genomeL = len(LCscores)-1 # all genome pos for which scores are calculated. Last line has mean and SD. thus -1.
    allgenomepos = [i+1 for i in range(genomeL)]
    d = 20
    numindelpos = len(indelpos)
    ##nonidlpos = [x for x in allgenomepos if x not in indelpos]
    for i in range(1,11): # randomly sample genome positions that are non-indel, and compute p value for diff of mean
        print "\n===\niteration ",i
        randpos = random.sample(allgenomepos, numindelpos+2700)
        randnonidlpos = []
        for r in randpos:
            flag = 0
            for ipos in indelpos:
                if r > ipos-d and r < ipos+d: # the randomly selected genome position is within dist bases of an indel pos, don't use it 
                    flag = 1
                    break
            if flag == 0:
                randnonidlpos.append(r)
        nonidlpos = randnonidlpos[:]
        print "Num indel and nonindel pos: ", numindelpos, len(randnonidlpos)
        print "T test for diff in mean for LC indel score and non indel scores:"
        indelComplexityScore(indelpos, randnonidlpos, LCscores, "LC", i, FIGPREFIX)
        print "T test for diff in mean for SE indel score and non indel scores:"
        indelComplexityScore(indelpos, randnonidlpos, SEscores, "H", i, FIGPREFIX)

# end of compareIndelNonindel #######
