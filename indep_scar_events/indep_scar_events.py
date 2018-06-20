#!/usr/bin/python

# Reads a directory containing Newick trees with extension ".outtree" and output number of scar clades in each tree

import os
import re
import numpy as np

def getScarIsolates(fh):
    scarsids = {} # key: color code for 
    numisolates = 0 
    for line in fh:
        numisolates += 1
        line = line.rstrip("\n")
        contents = line.split() # index 1: color of scar in tree, last index: isolate name.
        if contents[1] not in scarsids:
            scarsids[contents[1]] = [contents[-1]] # initialize
        else:
            scarsids[contents[1]].append(contents[-1]) # this scar (denoted of color in tree) has been seen before in another isolate
    return scarsids, numisolates 

def drawSimpleTree(sids, tree): 
    skeys = {} # key: sid, value: index of list that has the key
    for i in range(len(sids)):
        skeys[sids[i]] = i
    for sk in skeys:
        regex_sid = re.escape(sk)
        regex_numid = re.escape(str(skeys[sk]))
        tree = re.sub(regex_sid, regex_numid, tree)
    print tree

def extendClades(tree, cladehits, notseensids, sids_to_del): # see if sids not in clades are adjacent to one and form a clade
    #print "Call to extendClades:"
    #print "\tClades:", cladehits
    #print "\tnotseensids: ", notseensids
    #print "\tNum clades: ", len(cladehits), " Num notseensids: ", len(notseensids)
    # see how many unseen sids belong to a clade (i.e. when placed adjacent to a clade, they show up in a tree)
    for i in range(len(cladehits)):
        cladecontent = list(cladehits[i])
        #print "i: ", i, cladehits[i]
        for j in range(len(notseensids)):
            #print "\tj:  ", j, notseensids[j]
            extend_left = extend_right = "" # initialize
            if cladecontent[0] == "(": # clade is terminal leaves from left side
                extend_left = notseensids[j]+","+cladehits[i]
            else:
                extend_left = notseensids[j]+",("+cladehits[i]
            if cladecontent[-1] == ")": # clade has terminal leaves from right side
                extend_right = cladehits[i]+","+notseensids[j]
            else:
                extend_right = cladehits[i]+"),"+notseensids[j]
            #print "left and right:", extend_left, "\t", extend_right

            if extend_left in tree:
                cladehits[i] = "("+extend_left+")" # update this clade
                sids_to_del.append(j) # save for deletion the index of the sid that was added to clade from unseensids list
                tree, cladehits, notseensids, sids_to_del = extendClades(tree, cladehits, notseensids, sids_to_del) # recursive call to keep updating clades
                return tree, cladehits, notseensids, sids_to_del
            elif extend_right in tree:
                #print "right in tree"
                cladehits[i] = "("+extend_right+")" # update this clade
                sids_to_del.append(j) # save for deletion the index of the sid that was added to clade from unseensids list
                tree, cladehits, notseensids, sids_to_del = extendClades(tree, cladehits, notseensids, sids_to_del) # recursive call to keep updating clades
                return tree, cladehits, notseensids, sids_to_del
                #continue
    # check if any cladehits are now adjoining and can be merged
    newclades = []
    flag = 0
    if len(cladehits) > 1: # at least 2 clades needed to check adjoinment
        step = 1
        i=0
        while i < len(cladehits)-1:
        #    print "i", i, cladehits[i]
            flag = 0
            for j in range(i+1, len(cladehits)):
        #        print "j", j, cladehits[j]
                adjclade = "("+cladehits[i]+","+cladehits[j]+")"
        #        print "adjacent ", adjclade
                if adjclade in tree:
                    newclades.append(adjclade)
                    flag = 1
            if flag == 0: # this clade was separate from other clades
                newclades.append(cladehits[i])
                i = i+1
            else:
                # both ith and jth clade have been added to newclades. Update "i" to skip both i and j
                i = i+2
        #        print "adjacent in tree!!! i is now", i
                #i = min(j+1, len(cladehits)-2) # -2 because i shoud denote second last index of cladehist
        if flag == 0: # second last elt didn't adjoin the last elt, add the last elt to newclades 
            newclades.append(cladehits[-1])
        if len(newclades) != len(cladehits): # some clades were adjacent           
        #    print "new: ", newclades, "\tcladehits:", cladehits
            cladehits = newclades
            tree, cladehits, notseensids, sids_to_del = extendClades(tree, cladehits, notseensids, sids_to_del)
        else:
            return (tree, cladehits, notseensids, sids_to_del)
    return (tree, cladehits, notseensids, sids_to_del)

def getClades(outtree, stypesids):
    treestr = ""
    scarclades = {} # sections of tree (clades) that have the scarred isolates. Key: sids in this clade, value: clade with isolates
    with open(outtree) as fh: # remove all white spaces at ends of line and concatenate all lines
        treestr = "".join(line.strip() for line in fh)
    fh.close()
    # remove all digits and special characters from the tree string
    tree_nodigits = re.sub("(\d+\.\d+)", "", treestr)
    tree_onlysids = re.sub("[\:\.]", "", tree_nodigits)
    #print tree_onlysids
    #regex = r'(\(*[A-Z\,0-9]*\(*[A-Z\,0-9]*\(*[A-Z\,0-9]*\([A-Z\,0-9]+\)[A-Z\,0-9]*\)*[A-Z\,0-9]*\)*[A-Z\,0-9]*\)*\,*)'
    regex = r'(\([A-Z\,0-9]+\))'
    clades = re.findall(regex, tree_onlysids)
    terminal_leaves_sids = []
    cladehits = []

    for c in clades: # get clades that have at least one scar sid
        for i in range(len(stypesids)):
            if stypesids[i] in c: # if this strain id is in this clade
                cladehits.append(c)
                break

    # get those clades that have only scar sids. 
    allcladesids = []
    incompleteclades = [] # clades in which one more sids are not scar sid
    for ch in cladehits:
        #print "\nClade:\n", ch
        allidsinclade = re.findall(r'([A-Za-z0-9\-]+)', ch)
        #print "clade", ch, "has sids", allidsinclade
        flag = 0
        for cid in allidsinclade:
            if cid not in stypesids: # if this clade sid is not a scar sid
                flag = 1
                break
        if flag == 0: # all sids in clade are scar ids
        #    print len(allidsinclade), "sids in this clade: ", allidsinclade 
            allcladesids.extend(allidsinclade)
            #drawSimpleTree(allidsinclade, ch) # draw a tree with sids denoted as alphabets for easy reading
        else: # there are some sids in this clade that are non-scar. all sids then should go to singleton set (i.e. NOT A SCAR CLADE)
            incompleteclades.append(ch)
    notseensids = [x for x in stypesids if x not in allcladesids] # scar sids not found in scar clades      
    scarclades = [x for x in cladehits if x not in incompleteclades]
    #print "\n", len(notseensids), " sids not found in terminal leaf groups. Their tree segments are:"
    #for s in notseensids:
    #    regex2 = r'.{0,50}'+ re.escape(s) + r'.{0,50}'
    #    treesegment = re.findall(regex2, tree_onlysids)
    #    print s, " in tree segment: ", treesegment
    sids_to_del = [] #unseen sids that will extend existing clades and thus need to be deleted from the singleton set
    tree, finalclades, allsingletons, sids_to_del = extendClades(tree_onlysids, scarclades, notseensids, sids_to_del)
    #print "sids to del:", sids_to_del
    #print "original singletons:", allsingletons
    finalsingletons = [allsingletons[idx] for idx in range(len(allsingletons)) if idx not in sids_to_del] # singletons that extended clades are to be excluded from this list
    indepevents = len(finalclades)+len(finalsingletons)
    print "Total number of evolutionarily independent occurences of this scar: ", len(finalclades)+len(finalsingletons)
    print "clades: ", finalclades
    print "singletons: ", finalsingletons
    return indepevents

############ MAIN ###############
dirpath = "../convergent_evol/" # directory where tree files are stored
reps = ["_main"] + [str(i) for i in range(9)] # replicates are named "main", and then from 0-8: total 10 reps
ievents = []
for r in reps: # looking at 0-8 reps
    reppath = dirpath+"rep"+r+"/"
    allf = os.listdir(reppath)
    n_multisids = 0
    n_singlesidperscar = 0
    n_multisidperscar = 0

    all_evol_events = 0 # indep evol events for all scars in all genes

    for f in allf:
        if f.endswith(".outtree"):
            contents = f.split("_")
            gname = contents[0]
            #if gname != "espB":
            #    continue
            scarf = open(dirpath+"allscarisolates_ornament_"+gname+".map", 'r') # read ornament file to get isolates with scars, and isolates that have the same scar
            scarsids, numisolates = getScarIsolates(scarf)
            if numisolates > 1: # only those genes that had scars in more than one isolate
                n_multisids += 1
                print "\n\n", gname
                for stype in scarsids:
                    fillregex = re.escape("fill:")+ r'(.+)' + re.escape(";stroke") # get the fill color of this scar type
                    scarfill = re.findall(fillregex, stype)
                    print "\tscar type: \"", scarfill[0], "\" is present in ", len(scarsids[stype]), " sids."
                    #print scarsids[stype]
                    if len(scarsids[stype]) > 1: # if multiple isolates have this particular scar
                        n_multisidperscar += 1
                        n_events = getClades(reppath+f, scarsids[stype])
                        all_evol_events += n_events
                    else:
                        n_singlesidperscar += 1

    print "genes with scars in >1 isolate: ", n_multisids
    print "\t amongst those genes, scars found in a single isolate: ", n_singlesidperscar, ", and scars found in multiple isolates: ", n_multisidperscar 
    print "\t amongst scars in multiple isolates, total number of evol indep events of scars: ", all_evol_events

    print "\n\nREP", r, ": Total number of evolutionarily independent scar events:"
    print "REP", r, ": Genes with scar in a single isolate:", 57-26 # 57 total scar genes
    print "REP", r, ": Amongst genes with scars in > 1 isolate (26 genes), num isolates with a unique scar in a single isolate:", n_singlesidperscar
    print "REP", r, ": num evol indep events for scars in these 26 genes that are present in multiple isoaltes:", all_evol_events
    ievent_thisrep = 57-26+n_singlesidperscar+all_evol_events
    print "REP", r, ": Total num of evol indep events for ALL scars in ALL genes:", ievent_thisrep
    ievents.append(ievent_thisrep)

print "Indep events from all reps:", ievents
print "Average:", np.mean(ievents), " Stdev:", np.std(ievents)

