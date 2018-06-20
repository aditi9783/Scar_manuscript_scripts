#!/usr/bin/env python

# Reads a file that contains a list of evolutionarily independent non-scarred isolates, selects a random sample of 200 isolates,
# accesses data dir of each of these 200 isoaltes to determine SNPs in that isolate, does the same for scarred gene isolates
# computes pairwise SNP differences for all isolates (scarred isolates for each gene + 200 random isolates),
# writes the pairwise SNP differences as a distance matrix that serves as input to the PHYLIP program for generating trees.

import os, subprocess
import random

# start of findMutations #################
def findMutations( seqlist ):
    mutmatrix = []
    strainids = []
    #fout = open("distmat_strain_ids.out", 'w') # write the strain ids as key for the distance matrix
    counter = 0
    for sq in seqlist:
        seqdir, sd = sq
        strainids.append(sd)
    #    fout.write(str(counter)+"\t"+sd+"\n")
        mutlist = []
        fh = open(seqdir+sd+".vcf", 'r')
        for line in fh:
            if line.startswith("#"):
                    continue
            else:
                content = line.split("\t")
                if float(content[5]) > 200.0: # minimum quality score for calling a mutation
                    if len(content[3]) == 1 and len(content[4]) == 1: # both the ref base and mutated base are of len 1 => no indels, only SNPs 
                        mut = content[1]+content[4]
                        mutlist.append( mut )
        fh.close()
        mutmatrix.append( mutlist )
    #fout.close()
    return mutmatrix, strainids

# end of findMutations ###################

# start of calcDistances #################
def calcDistances( mutmatrix, strainids ):
    distmat = []
    for i in range(len(mutmatrix)):
        distmat.append([0 for x in range(len(mutmatrix))]) # initialize the distances
    
    for i in range(len(mutmatrix)):
        for j in range(i+1, len(mutmatrix)):
            intersection = list( set(mutmatrix[i]) & set(mutmatrix[j]) )
            union = list( set(mutmatrix[i]) | set(mutmatrix[j]) )
            num_shared_mut = len(intersection) # number of mutations common in both strains
            num_all_mut = len(union) # number of unique mutations present in both strains
            num_mut_diff = num_all_mut - num_shared_mut
            distmat[i][j] = num_mut_diff
            #distmat[i][j] = num_shared_mut
            distmat[j][i] = distmat[i][j] # because the matrix is symmetric
    return distmat
    
# end of calcDistances ###################

# start of printPhylipDistMat ############
def printPhylipDistMat( strainids, distmat, genename, outpath ):
    fout = open(outpath+genename+"_and_200independent_distmat_phylip.txt", 'w')
    numstrains = len(strainids)
    fout.write("\t"+str(numstrains)+"\n")

    for i in range(numstrains):
        # Zhang acc numbers are variable len. ENA acc numbers are either 9 or 10 in length. Phylip requires that strain names be of length 10
        padreq = 10 - len(strainids[i])
        for j in range(padreq):
            strainids[i] = strainids[i]+" " # add space so that strain id string is of len 10
        fout.write(strainids[i]+"\t"+"\t".join([str(v) for v in distmat[i]])+"\n")

    fout.close()
# end of printPhylipDistMat ############

# start of getGenesMultipleIsolates ####
def getGenesMultipleIsolates(fname): # get genes that have scars in multiple isolates
    genescar_multisolates = {} # key: name of gene that has scar in multiple isolates, value: empty list in which the scar isolates id will saved
    fh = open(fname, 'r')
    for line in fh:
        contents = line.split()
        numisolates = int(contents[0])
        if numisolates > 1:
            if ".map" in contents[1]: # ignore the last line of file that shows total num isolates with scars
                str1 = contents[1].rstrip(".map\n")
                list1 = str1.split("_")
                gname = list1[-1]
                genescar_multisolates[gname] = [] # initialize empty list for the gene that has scars in multiple isoaltes
    print "Num of genes with scars in multiple isoaltes: ", len(genescar_multisolates)
    return genescar_multisolates
# end of getGenesMultipleIsolates ####

if __name__ == '__main__':
    isolatefiles = ["zhang_isolates_scars_confirmed_noPGRS.out", "glynn_isolates_scars_confirmed_noPGRS.out", "walker_isolates_scars_confirmed_noPGRS.out"]
    genescar_multisolates = getGenesMultipleIsolates("numisolates_scarredgenes.out")
    scarisolates = []
    isolatesseen = []
    for f in isolatefiles:
        fh = open(f, 'r')
        for line in fh:
            content = line.split()
            seqdir = content[0]
            sd = content[1]
            gene = content[2]
            #if gene == "espI": # only select espI isolates for now
            if gene in genescar_multisolates: # this gene has scars in multiple isolates
                if sd not in isolatesseen:
                    #scarisolates.append([seqdir, sd]) # this line is for espI gene
                    genescar_multisolates[gene].append([seqdir, sd])
                    #isolatesseen.append(sd) # some isolates have scars in multiple genes, and they were being added twice. Add each isolate once.
        fh.close()

    #print "Number of isolates with scar in espI:", len(scarisolates) # for espI gene
    # to the espI scar isolates, add 500 independent isolates
    indepisolates = []
    # read file that contains the full paths of data dir for evolutionarily independent isolates
    ifh = open("all_highcov_isolates_numSNPdiff_distmat_phylip.outtree.independent.txt.fullpath.nonscarred", 'r')
    for line in ifh:
        line = line.rstrip("\n")
        seqdir, sd = line.split()
        indepisolates.append([seqdir, sd])

    #allisolates = scarisolates + randsample_indep
    #print "Total num of scar + indep isolates:", len(allisolates)
    #mutmatrix, strainids = findMutations( allisolates )
    #distmat = calcDistances( mutmatrix, strainids )
    #printPhylipDistMat( strainids, distmat ) # print in the distance matrix format read by Phylip

    n_random = 9 # number of times random samples are to be generated

    for i in range(n_random):
        print "\nIter", i
        outpath = "/convergent_evol/rep"+str(i)+"/" # location where the phylip distmat matrices are to be written for generating NJ trees
        randsample_indep = random.sample(indepisolates, 200) 
        for g in genescar_multisolates:
            allisolates = genescar_multisolates[g] + randsample_indep
            print "Number of isolates with scar in ", g, ":", len(genescar_multisolates[g]) 
            print "Total num of scar + indep isolates:", len(allisolates)
            mutmatrix, strainids = findMutations( allisolates )
            distmat = calcDistances( mutmatrix, strainids )
            printPhylipDistMat( strainids, distmat, g, outpath ) # print in the distance matrix format read by Phylip
