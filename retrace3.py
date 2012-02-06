#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
"""
ReTrace - find branching pathways in metabolic networks
Copyright (C) 2009-2010 Esa Pitkänen
See file COPYING for license notice.
"""

# Default options

DEFAULT_NUMPATHS = [100,1]     # 100 paths at first level, then 1 at second, third, ...
DEFAULT_MAXDEPTH = 3           # Find pathways with max. 3 branches
DEFAULT_TRACED_ATOMS = ["C"]   # Trace only carbons

ALGORITHM_KSP = "SimpleYen"

DIR_CONSTRAINT_FORWARD = ">"
DIR_CONSTRAINT_BACKWARD = "<"
DIR_CONSTRAINT_FORBIDDEN = "-"

import sys, os, datetime, random, re, shutil, getopt, time, math

from tracerutil import *
from htmlexport import *

from metabolism.parser import kegg_ligand_parser as klp
from graph.graph import Graph
from graph.shortestPaths import dijkstra, KSPDijkstra
from ds.comb import comb, cartesian

#import rea      # REA: solver for k (non-simple) shortest paths 
import kspyen   # A simple implementation of Yen's algorithm for k shortest simple paths

class Pathway:

    def __init__(self, src, tgt, rpairs, Z, compositeMap):
        self.src = src
        self.tgt = tgt
        self.rpairs = rpairs
        self.Z = Z
        self.compositeMap = compositeMap

    def __str__(self):
        s = "Pathway %s -> %s (Z = %.1f)\n" % (",".join(self.src), self.tgt, self.Z)
        s += "   Rpairs:         %s\n" % (list(self.rpairs))
        s += "   Composite map:  %s\n" % (self.compositeMap)
        return s

def pairScoreCmp(x, y):
    if x[1]:
        if not y[1]:
            return 1
        else:
            return int(x[1] - y[1])
    elif y[1]:
        return -1
    else:
        return 0

def reverseMap(m):
    rmap = {}
    for atom in m:
        rmap[m[atom]] = atom
    return rmap

def getDirectedAtomMap(ag, pair, dir):
    if dir == "f":
        atomMap = ag.atomMap[pair]
    elif dir == "r":
        atomMap = reverseMap(ag.atomMap[pair])
    else:
        raise RuntimeError, "Unknown dir: %s" % (dir)
    return atomMap        

def getIntermediates(r, mol, pro, atom, ag):
    intermediates = set()

    if mol in ag.edges and pro in ag.edges[mol]:
        molpairs = ag.edges[mol][pro]
    elif pro in ag.edges and mol in ag.edges[pro]:
        molpairs = ag.edges[pro][mol]
    else:
        # This molecule does not map to the product
        return intermediates

    for rp in molpairs:

        if r not in ag.pairToReaction[rp]:
            continue

        dir = pro == ag.pairToSub[rp]
        if dir:
            # rpair in pro -> im direction
            rmap = getDirectedAtomMap(ag, rp, "f")
        else:
            # rpair in im -> pro direction
            rmap = getDirectedAtomMap(ag, rp, "r")

        # Check if an atom of this potential intermediate
        # is mapped to sub
        if atom in rmap:
            #print "*** Necessary intermediate %s-%s and rpair %s" % (mol, rmap[atom], rp)
            # Note that the rpair rp that links the new branch
            # to the existing pathway needs to be in _reverse_ direction
            # that the map found above, hence "not dir"
            intermediates.add((mol, rmap[atom], (rp, not dir)))

    return intermediates

def determineRequiredIntermediates(path, ta, ag, net):
    atom = ta

    intermediates = set()

    # path = [.., pair, tgtatom, tgt]
    pathpos = len(path) - 3    # last rpair in path

    while pathpos != 1:  # while not in path head

        pair, srcatom, tgtatom = path[pathpos].split("-")        # e.g. (A06006_f, 1, 7)
        pair, dir = pair.split("_")                     # e.g. (A06006, f)

        # We go backwards from target
        if dir == "f":
            dir = "r"
        else:
            dir = "f"
        amap = getDirectedAtomMap(ag, pair, dir)

        if atom not in amap:
            # This atom is not mapped to by the current pair
            
            # This is the metabolite on the path
            # where this atom does not originate from.
            sub = path[pathpos - 1].split("-")[0]
            pro = path[pathpos + 1].split("-")[0]

            # Determine which metabolite should be produced to
            # reach this atom: find reactions that 
            # associate with this rpair

            # Find substrates of each reaction from which
            # the missing atom could be transferred from.
            for r in ag.pairToReaction[pair]:
                re = net.reactions[r]

                # Reaction re can be in either direction:
                # determine direction by checking on which side
                # does the pathway substrate sub appear 
                if sub in re.substrates:
                    potentialIntermediates = set(re.substrates)
                else:
                    potentialIntermediates = set(re.products)
                potentialIntermediates.remove(sub)
                if len(potentialIntermediates) == 0:
                    print "Error in atom mappings: reaction %s should have other reactants in addition to %s" % (re, sub)
                    continue

                # Check which of the molecules (exactly one if mappings are correct!)
                # is able to originate the missing atom
                for mol in potentialIntermediates:
                    intermediates.update(getIntermediates(r, mol, pro, atom, ag))

            break  # while pathpos != -1
        else:  # this atom can be traced further back towards sources
            atom = amap[atom]
            pathpos -= 2       # go to previous rpair

    return intermediates

def doFindPathsAndIntermediates(g, tdir, net, sources, unsatisfied, target, ag, pairScores, npaths, tracedSourceAtoms, reactionDirConstraints, tracedAtoms):

    tempfn = "%s/temp_ksp%d.txt" % (tdir, random.randint(10**5, 10**6))

    if npaths == 1:
        (ddist, dpaths) = dijkstra(g, SROOT, target)
        node = target
        path = []
        while node in dpaths:
            path.insert(0, node)
            node = dpaths[node]
        path.insert(0, node)
        paths = [path]
        acyclics, cyclics = 1, 0

    else:
        if ALGORITHM_KSP == "REA":
            # REA was utilized to check amount of cyclic paths in preliminary tests
            cyclicpaths = rea.REAksp(g, npaths, SROOT, target)
        elif ALGORITHM_KSP == "SimpleYen":
            cyclicpaths = kspyen.kspSimpleYen(g, npaths, SROOT, target)
        else:
            print "Unknown algorithm: " % (ALGORITHM_KSP)
            assert(0)

        paths = kspyen.removeCyclicPaths(cyclicpaths)
        acyclics, cyclics = len(paths), len(cyclicpaths) - len(paths)
        print "%d/%d acyclic paths" % (len(paths), len(cyclicpaths))

    results = []

    if paths == None:
        return results, 0, 0

    for path in paths:

        if len(path) < 4:
            # This should not happen but if it happens, lets not abort
            print "Warning: short path (assertion failed):", path
            continue

        # Collect the set of atom nodes which need to be
        # covered by subsequent path queries
        intermediates = set()

        # Traverse path backwards from each unsatisfied atom, 
        # check which atoms cannot be traced back to sources

        for u in unsatisfied:
            tgt, ta = u

            if ag.atomTypes[tgt][ta] not in tracedAtoms:
                continue

            pathend = path[-2].split("-")[0]
            if tgt != pathend:
                intermediates.add((tgt, ta))
                continue

            newIntermediates = determineRequiredIntermediates(path, ta, ag, net)
            intermediates.update(newIntermediates)

        # Extract rpairs from path and return them as result
        res = set()
        
        for i in range(3, len(path) - 1, 2):
            rpair, stuff = path[i].split("_")
            dir = stuff.split("-")[0] == "f"  # True: rpair in forward direction
            res.add((rpair, dir))

        results.append((res, intermediates))

    return results, acyclics, cyclics

# Identifiers for dummy nodes; internal use only

TROOT = "troot"
SROOT = "root"

def hashResult(path):
    s = ""
    items = list(path)
    items.sort()
    for p, dir in items:
        s += "%s_%s*" % (p, dir)
    return hash(s)

def hashPartialResult(path, unsat):
    s = ""
    items = list(path)
    items.sort()
    for p, dir in items:
        s += "%s_%s=" % (p, dir)
    items = list(unsat)
    items.sort()
    for mol, atom in items:
        s += "%s=%d*" % (mol, atom)
    return hash(s)

def pathContainsRpairBothDirections(joinedpath):
    pairs = set()
    for rpair, dir in joinedpath:
        if rpair in pairs:
            return True
        pairs.add(rpair)
    return False

def doFindPathways(g, tdir, net, sources, unsatisfied, target, ag, pairScores, npaths, 
                   tracedSourceAtoms, reactionDirConstraints, tracedAtoms, maxdepth, 
                   greedyfinish, storePartials):
    
    finished = {}
    incomplete = {}
    
    tasks = []
    tasks.append((set(), unsatisfied, 1))

    visited = set()
    
    round = 1

    ndups = 0
    ndepthreached = 0
    npruned = 0

    numAcyclic = 0
    numCyclic = 0

    while len(tasks) > 0:

        # path: partial pathway, [(rpair, dir), (rpair, dir), ...]
        # usnodes: unsatisfied atom nodes
        # depth: current search depth
        path, usnodes, depth = tasks.pop()

        if depth >= len(npaths):
            k = npaths[-1]
        else:
            k = npaths[depth - 1]

        if greedyfinish and depth > 1:
            k = 1 

        status("*** Round %d (depth %d, k=%d): %d tasks, %d complete, %d incomplete, %d dups, %d max depth, %d pruned" % (round, depth, k, len(tasks), len(finished), len(incomplete), ndups, ndepthreached, npruned))
 
        round += 1

        if len(usnodes) == 0:
            h = hashResult(path)
            if h not in finished:
                finished[h] = path
                if h in incomplete:
                    del incomplete[h]
            else:
                ndups += 1
            continue

        if depth > maxdepth:
            ndepthreached += 1
            h = hashResult(path)
            if h not in finished or h not in incomplete:
                incomplete[h] = path
            continue

        # Add unsatisfied->troot edges
        for u in usnodes:
            unode = "%s-%s" % (u[0], u[1])
            g.addEdge(unode, TROOT, 0.0)

        results, acyclics, cyclics = doFindPathsAndIntermediates(g = g, tdir = tdir, net = net, 
                                              sources = sources, unsatisfied = usnodes, 
                                              target = TROOT, ag = ag, 
                                              pairScores = pairScores,
                                              npaths = k, 
                                              tracedSourceAtoms = tracedSourceAtoms, 
                                              reactionDirConstraints = reactionDirConstraints, 
                                              tracedAtoms = tracedAtoms)

        numAcyclic += acyclics
        numCyclic += cyclics

        # Remove extra edges so that we do not find the same
        # path again
        for u in usnodes:
            unode = "%s-%s" % (u[0], u[1])
            if TROOT in g.E[unode]:   # no idea why sometimes this is not true; this is a potential source of nasty bugs: if doFindPathsAndIntermediates modifies g then g might "erode" during search and cause unpredictable behaviour; EDIT: maybe because previously usnodes contained duplicate u[0],u[1]-entries?
                g.removeEdge(unode, TROOT)

        for respath, resunsatisfied in results:
            newpartials = 0

            unsatmols = set()

            # The following code generates the different combinations
            # of RPAIR additions to the pathway
            
            cands = {}
            for vals in resunsatisfied:
                if len(vals) == 3:
                    # unsatisfied node with rpair candidates
                    uMol, uAtom, uRpair = vals
                    unsatmols.add((uMol, uAtom))
                else:
                    # unsatisfied node with no rpair candidates:
                    # this case will remain in the pool of unsatisfied
                    # nodes
                    uMol, uAtom = vals
                    unsatmols.add((uMol, uAtom))
                    continue

                uRpairId, uRpairDir = uRpair
                key = "%s-%s" % (uRpairId, uRpairDir)
                if key not in cands:
                    cands[key] = set()
                cands[key].add("%s-%s" % (uMol, uAtom))
            rpairCov = {}
            for uRpair in cands:
                molatoms = list(cands[uRpair])
                molatoms.sort()
                molatoms = ",".join(molatoms)
                if molatoms not in rpairCov:
                    rpairCov[molatoms] = set()
                rpairCov[molatoms].add(uRpair)
            choices = []
            for node in rpairCov:
                choices.append(list(rpairCov[node]))
            choices = cartesian(*choices)  # unpack args from list

            # Now choices holds all variations how to extend the pathway

            for rps in choices:
                joinedpath = path.union(respath)

                # Add the indicated rpair to pathway to connect the new branch
                # to the rest of the pathway
                for rpairIdDir in rps:
                    rpairId, rpairDir = rpairIdDir.split("-")
                    joinedpath.add((rpairId, rpairDir))

                if pathContainsRpairBothDirections(joinedpath):
                    continue

                h = hashPartialResult(joinedpath, unsatmols)

                if h not in visited:

                    tasks.append((joinedpath, unsatmols, depth + 1))
                    visited.add(h)
                    newpartials += 1

                    if storePartials:
                        h = hashResult(joinedpath)
                        if h not in incomplete or h not in finished:
                            incomplete[h] = joinedpath

        if len(results) == 0:
            # No paths found because connections to unsatisfied metabolites
            # have been removed in pruning.
            h = hashResult(path)
            if h not in finished or h not in incomplete:
                incomplete[h] = path
                npruned += 1

    status("%d completed paths" % (len(finished)))
    status("%d incomplete paths" % (len(incomplete)))

    acyclicFreq = 1.0 * numAcyclic / (numAcyclic + numCyclic)

    return finished.values(), incomplete.values(), round, acyclicFreq

def determineAtomGraphPruning(g, ag, sroot, target, prunedGraphSize):

    for i in range(1, len(ag.atomTypes[target]) + 1):
        g.addEdge("%s-%s" % (target, i), target, 1.0)
        g.addEdge(target, "%s-%s" % (target, i), 10000)
    
    (sd, sp) = dijkstra(g, sroot)
    (td, tp) = dijkstra(g, target)

    for i in range(1, len(ag.atomTypes[target]) + 1):
        g.removeEdge("%s-%s" % (target, i), target)
        g.removeEdge(target, "%s-%s" % (target, i))

    revDists = {}                
    removals = set()
    for u in g.V:
        removals.add(u)
        if u in sd and u in td:
            d = sd[u] + td[u]
            if d in revDists:
                revDists[d].append(u)
            else:
                revDists[d] = [u]
    rdists = revDists.keys()
    rdists.sort()
    count = 0

    ind = 0
    inclusions = set()
    cutval = None
    for d in rdists:
        for u in revDists[d]:
            ind += 1
            if not cutval:
                removals.remove(u)
                if ind < prunedGraphSize:
                    pass
                else:
                    cutval = sd[u] + td[u]
            elif sd[u] + td[u] < cutval:
                removals.remove(u)
            else:
                pass

    print "%d in removals" % (len(removals))

    return removals

def findPathways(tdir, net, sources, tgt, ag, pairScores, npaths, tracedSourceAtoms, reactionDirConstraints, rpairConstraints, tracedAtoms, maxdepth, greedyfinish, prunedGraphSize, storePartials, weighting):

    resultPathways = []

    tracedTgtAtoms = set()
    for i in range(1, len(ag.atomTypes[tgt]) + 1):
        type = ag.atomTypes[tgt][i]
        if type in tracedAtoms:
            tracedTgtAtoms.add(i)

    status("Looking for paths from %s to %s" % (",".join(sources), tgt))

    minTracedSources = sys.maxint
    tracedSrcAtoms = {}
    for src in sources:
        tracedSrcAtoms[src] = set()
        for i in range(1, len(ag.atomTypes[src]) + 1):
            if (src in tracedSourceAtoms and i in tracedSourceAtoms[src]) or (src not in tracedSourceAtoms and ag.atomTypes[src][i] in tracedAtoms):
                tracedSrcAtoms[src].add("%s-%d" % (src, i))

        status("Tracing %d/%d atoms in source %s" % (len(tracedSrcAtoms[src]), ag.atomCount[src], src))
        if len(tracedSrcAtoms[src]) < minTracedSources:
            minTracedSources = len(tracedSrcAtoms[src])

    status("Tracing %d/%d atoms in target %s" % (len(tracedTgtAtoms), ag.atomCount[tgt], tgt))

    g = buildAtomTraceGraph(ag, tracedAtoms, pairScores, weighting, 
                            sources = tracedSrcAtoms, target = None, minCommonAtoms = 1)

    removeConstrainedRpairs(ag, g, rpairConstraints, tracedAtoms)

    status("%d nodes and %d edges in atom graph" % (g.numNodes(), g.numEdges()))

    # Prune the graph if requested

    if prunedGraphSize != None:
        removals = determineAtomGraphPruning(g, ag, SROOT, tgt, prunedGraphSize)
        g2 = buildAtomTraceGraph(ag, tracedAtoms, pairScores, weighting, 
                                 sources = tracedSrcAtoms, target = None, minCommonAtoms = 1, 
                                 skipnodes = removals)
        removeConstrainedRpairs(ag, g2, rpairConstraints, tracedAtoms)

        print "%d nodes and %d edges in the pruned atom graph" % (g2.numNodes(), g2.numEdges())
    else:
        g2 = g

    unsatisfied = set()

    # Initialise the set of unsatisfied nodes with target atoms,
    # set rpair (3rd item) each to None, as we do not need to
    # add an rpair after we have satisfied a target node
    # to connect node to pathway.
    for i in range(1, len(ag.atomTypes[tgt]) + 1):
        unsatisfied.add((tgt, i))

    completePaths, incompletePaths, numRounds, acyclicFreq = doFindPathways(g2, tdir, net, sources, unsatisfied, TROOT, ag,
                                                    pairScores, npaths, tracedSourceAtoms, 
                                                    reactionDirConstraints, tracedAtoms, maxdepth,
                                                    greedyfinish, storePartials)

    return completePaths, incompletePaths, numRounds, acyclicFreq


def parseReactionDirConstraints(reactionDirFn, ag, net):
    reactionDirConstraints = {}
    if reactionDirFn != None and reactionDirFn != "-":
        rdf = open(reactionDirFn)
        for s in rdf:
            if s.startswith("#"):
                continue
            r, dir = s.strip().split()
            if r in reactionDirConstraints:
                print "Error: direction of %s constrained more than once in %s" % (r, reactionDirFn)
                sys.exit(1)
            if dir != DIR_CONSTRAINT_FORWARD and dir != DIR_CONSTRAINT_BACKWARD and dir != DIR_CONSTRAINT_FORBIDDEN:
                print "Error: unknown constraint %s for %s in %s" % (dir, r, reactionDirFn)
            #reactionDirConstraints[r] = dir == ">"
            reactionDirConstraints[r] = dir

    rpairConstraints = {}

    for rpair in ag.pairToReaction:
        rl = ag.pairToReaction[rpair]
        if rl == None:
            continue
        constrained = set()
        sub = ag.pairToSub[rpair]
        pro = ag.pairToPro[rpair]
        # Map reaction-level constraints to rpair level
        for r in rl:
            if r in reactionDirConstraints:
                re = net.reactions[r]
                if reactionDirConstraints[r] == DIR_CONSTRAINT_FORWARD:  
                    # reaction must occur in substrates -> products dir
                    if sub in re.substrates and pro in re.products:
                        constrained.add(DIR_CONSTRAINT_FORWARD)
                    elif sub in re.products and pro in re.substrates:
                        constrained.add(DIR_CONSTRAINT_BACKWARD)
                elif reactionDirConstraints[r] == DIR_CONSTRAINT_BACKWARD:
                    # reaction must occur in products -> substrates dir
                    if sub in re.substrates and pro in re.products:
                        constrained.add(DIR_CONSTRAINT_BACKWARD)
                    elif sub in re.products and pro in re.substrates:
                        constrained.add(DIR_CONSTRAINT_FORWARD)
                elif reactionDirConstraints[r] == DIR_CONSTRAINT_FORBIDDEN:
                    constrained.add(DIR_CONSTRAINT_FORWARD)
                    constrained.add(DIR_CONSTRAINT_BACKWARD)
                else:
                    print "Invalid reactionDirConstraint:", reactionDirConstraints[r]
                    assert(0) # Invalid reactionDirConstraint                    

        if len(constrained) == 1: 
            rpairConstraints[rpair] = constrained.pop()
            assert(rpairConstraints[rpair] == DIR_CONSTRAINT_FORWARD or rpairConstraints[rpair] == DIR_CONSTRAINT_BACKWARD)
        elif len(constrained) == 2:
            rpairConstraints[rpair] = DIR_CONSTRAINT_FORBIDDEN

    return reactionDirConstraints, rpairConstraints

class ReactionScore:
    def __init__(self, id, score, seq1, seq2, evalue, ecs):
        self.id = id
        self.score = score
        self.seq1 = seq1
        self.seq2 = seq2
        self.evalue = evalue
        self.ecs = ecs

def parseScores(scoresfn, net):
    scores = {}
    scoreOrigin = {}

    for r in net.reactions:
        scores[r] = 0

    if scoresfn != None:
        scoresf = open(scoresfn)           # score file: 1st col kegg id, 2nd col score, "?" denotes unknown score
        for s in scoresf:
            if s.startswith("#"):
                continue
            vals = s.strip().split("\t")
            scores[vals[0]] = vals[1]
            scoreOrigin[vals[0]] = ReactionScore(vals[0], vals[1], vals[2], vals[3], vals[4], vals[5])
    return scores, scoreOrigin

def calculateRpairScoresAndConstraints(ag, scores, reactionDirConstraints):
    pairScores = {}
    pairConstraints = {}

    for pair in ag.pairToReaction:
        maxScore = None
        maxReaction = None
        dirConstraints = set()

        pairScores[pair] = []

        if pair not in pairConstraints:
            pairConstraints[pair] = set()

        if ag.pairToReaction[pair]:
            for react in ag.pairToReaction[pair]:

                if react in scores and scores[react] != "?":
                    score = float(scores[react])
                else:
                    score = None
                pairScores[pair].append((react, score))

                if react in reactionDirConstraints:
                    pairConstraints[pair].add("%s%s" % (react, reactionDirConstraints[react]))
                else:
                    pairConstraints[pair].add("%s-" % (react))

        pairScores[pair].sort(pairScoreCmp)

    return pairScores, pairConstraints

def parseSources(sources):
    parsedSources = []
    tracedSourceAtoms = {} # compound -> set of traced index
    for src in sources:
        m = re.search("(\w+)-([\d/]+)", src)
        if m:
            srcmol = m.group(1)
            tracedSourceAtoms[srcmol] = set()
            atomixes = m.group(2).split("/")
            for atom in atomixes:
                tracedSourceAtoms[srcmol].add(int(atom))
        else:
            srcmol = src
        parsedSources.append(srcmol)
    return parsedSources, tracedSourceAtoms

def containsSameReactionInBothDirections(path, ag, net):
    reactions = {}

    for rpair, dir in path:
        sub = ag.pairToSub[rpair]
        for r in ag.pairToReaction[rpair]:
            re = net.reactions[r]
            if sub in re.substrates:
                if r in reactions and reactions[r] == False:
                    return True
                reactions[r] = True
            elif sub in re.products:
                if r in reactions and reactions[r] == True:
                    return True
                reactions[r] = False
            else:
                assert(0)

    return False

def isTargetMappingOk(path, ag, tgt):
    reachednodes = set()

    i = 0
    reactions = set()

    for rpair, dir in path:

        if ag.pairToSub[rpair] == tgt:
            rmap = getDirectedAtomMap(ag, rpair, 'f')
        elif ag.pairToPro[rpair] == tgt:
            rmap = getDirectedAtomMap(ag, rpair, 'r')
        else:
            continue

        if i == 0:
            reactions = set(ag.pairToReaction[rpair])
        else:
            reactions.intersection_update(ag.pairToReaction[rpair])

        i += 1

        for atom in rmap:
            node = "%s-%d" % (tgt, atom)
            if node in reachednodes:
                return False
            reachednodes.add(node)

    if len(reactions) == 0:
        return False

    return True

def processPaths(paths, ag, sources, tgt, tracedAtoms, net, discardInvalidTargetMappings, discardReactionsBothDirections, discardMultipleTargetMappings):
    """Construct Pathway instances from the pathways found by the algorithm and perform some validity checks on pathways if required:
discardInvalidTargetMapping:    discard a pathway where the pathway RPAIRs map to the same target atom more than once (a local property).
discardReactionsBothDirections: discard a pathway where the same reaction is necessarily used in both directions (note that operating in RPAIR representation makes this check a bit involved)
discardMultipleTargetMappings:  discard a pathway where composite mapping refers to the same target atom more than once (a global property).

"""

    processedPaths = []

    numDiscarded = 0

    for path, isComplete in paths:

        if discardInvalidTargetMappings and isTargetMappingOk(path, ag, tgt) == False:
            numDiscarded += 1
            continue

        if discardReactionsBothDirections and containsSameReactionInBothDirections(path, ag, net):
            numDiscarded += 1
            continue

        atoms = []
        queue = []
        visited = set()        
        origin = {}

        for src in sources:
            for sa in range(1, len(ag.atomTypes[src]) + 1):
                if ag.atomTypes[src][sa] in tracedAtoms:
                    node = "%s-%s" % (src, sa)
                    origin[node] = set([node])
                    queue.append(node)

        targetCounts = {}

        discard = False
        while len(queue) > 0:
            subatom = queue.pop(0)

            if subatom in visited:
                continue

            visited.add(subatom)

            sub, atom = subatom.split("-")
            atom = int(atom)

            if sub == tgt:
                if atom not in targetCounts:
                    targetCounts[atom] = 1
                else:
                    targetCounts[atom] += 1
                    discard = True

            # Find out where this atom is mapped to by path rpairs
            for rpair, dir in path:

                rmap = pro = None
                if ag.pairToSub[rpair] == sub:
                    # This rpair maps at least one atom from sub in fwd dir
                    rmap = getDirectedAtomMap(ag, rpair, 'f')
                    pro = ag.pairToPro[rpair]
                elif ag.pairToPro[rpair] == sub:
                    # This rpair maps at least one atom from sub in rev dir
                    rmap = getDirectedAtomMap(ag, rpair, 'r')
                    pro = ag.pairToSub[rpair]

                if rmap != None and atom in rmap:
                    # Rpair maps this atom to some other atom
                    node = "%s-%s" % (pro, rmap[atom])
                    
                    if node not in origin:
                        o = set()
                        o.update(origin[subatom])
                        origin[node] = o
                    else:
                        origin[node].update(origin[subatom])

                    queue.append(node)

        targetAtoms = 0        
        for i in range(1, len(ag.atomTypes[tgt]) + 1):
            if ag.atomTypes[tgt][i] in tracedAtoms:
                targetAtoms += 1

        compositeMap = {}

        targetsReached = 0
        for v in visited:
            mol, atom = v.split("-")
            if mol == tgt:
                targetsReached += 1
                for org in origin[v]:
                    if org not in compositeMap:
                        compositeMap[org] = set([v])
                    else:
                        compositeMap[org].add(v)

        Z = 1.0 * targetsReached / targetAtoms

        # Discard pathways that map target atom more than once?
        if discardMultipleTargetMappings and discard:
            numDiscarded += 1
            continue
        
        p = Pathway(sources, tgt, path, Z, compositeMap)
        processedPaths.append(p)

    return processedPaths, numDiscarded

def calculatePathwayStatistics(paths):

    lens = []
    total = 0.0
    for path in paths:
        l = len(path.rpairs)
        total += l
        lens.append(l)

    if len(lens) == 0:
        return 0, 0

    mean = total / len(lens)
 
    diff = 0.0
    for l in lens:
        diff += (l - mean) * (l - mean)
    diff = diff / len(lens)
    
    return mean, math.sqrt(diff)

def checkSourceTargetExist(sources, tgt, compounds):
    problems = False
    for src in sources:
        smol = src.split("-")[0]
        if smol not in compounds:
            print "Source %s not found in KEGG" % (smol)
            problems = True
    tmol = tgt.split("-")[0]
    if tmol not in compounds:
        print "Target %s not found in KEGG" % (tmol)
        problems = True
    if problems:
        print "Specify sources and target with KEGG identifiers (e.g., C00031)"
        sys.exit(1)

def retrace(tdir, sources, tgt, npaths, maxdepth, scoresfn, rdir, reactionDirFn, 
            tracedAtoms, greedyfinish, prunedGraphSize, storePartials,
            weights, disablePathOutput):

    status("Target dir: %s" % (tdir))
    status("Source(s): %s" % (sources))
    status("Target: %s" % (tgt))
    status("Paths for each missing metabolite: %s" % (",".join(map(str, npaths))))
    status("Scores file: %s" % (scoresfn))
    status("Database: %s" % (rdir))
    status("Use greedy finish: %s" % (greedyfinish))
    status("Store partial results: %s" % (storePartials))
    status("Disable path-specific output: %s" % (disablePathOutput))
    if prunedGraphSize != None:
        status("Prune graph to %s nodes" % (prunedGraphSize))
    else:
        status("Search unpruned atom graph")

    if weights == None or weights.startswith("u"):
        weighting = TRACE_GRAPH_WEIGHTS_UNWEIGHTED
    elif weights.startswith("s"):
        weighting = TRACE_GRAPH_WEIGHTS_RSCORE
    elif weights.startswith("a"):
        weighting = TRACE_GRAPH_WEIGHTS_NATOMS
    else:
        print "Unknown value %s to parameter -w. Try u, s or a." % (weights)
        sys.exit(2)

    if dotAvailable():
        status("Graphviz dot available")
    else:
        status("Warning: Graphviz dot unavailable")

    #if rea.isREAavailable():
    #    status("REA available")
    #else:
    #    status("Error: REA unavailable")
    #    sys.exit(1)

    searchParameters = {"Target dir (-o)" : tdir,
                        "NumPaths (-k)" : ",".join(map(str, npaths)),
                        "MaxDepth (-m)" : maxdepth,
                        "Scores (-c)" : scoresfn,
                        "Database (-d)" : rdir,
                        "GreedyFinish (-g)" : greedyfinish,
                        "PruneToSize (-p)" : prunedGraphSize,
                        "StorePartials (-i)" : storePartials,
                        "Weighting (-w)" : weights,
                        "Disable path-specific output" : disablePathOutput}

    sources = sources.split(",")

    st = datetime.datetime.now()

    f2 = open("%s/compound" % (rdir)) 
    status("Parsing compound...")
    compounds = klp.parse_compound(f2)

    f2 = open("%s/reaction" % (rdir)) 
    status("Parsing reaction...")
    net = klp.parse_network(f2)

    #print net.reactions["R08796"]

    f2 = open("%s/rpair" % (rdir)) 
    status("Parsing rpair...")
    rpairerrorsf = open("%s/rpair-errors.txt" % (tdir), "w")
    ag = klp.parseAtomGraphFromRpair(f2, verbose = False, errfile = rpairerrorsf)

    status("Parsing reaction direction constraints...")
    reactionDirConstraints, rpairConstraints = parseReactionDirConstraints(reactionDirFn, ag, net)

    status("Parsing scores...")
    scores, scoreOrigin = parseScores(scoresfn, net) 

    status("Calculating rpair scores and determining rpair direction constraints...")
    pairScores, ppairConstraints = calculateRpairScoresAndConstraints(ag, scores, reactionDirConstraints)

    sources, tracedSourceAtoms = parseSources(sources)

    checkSourceTargetExist(sources, tgt, compounds)

    startTime = datetime.datetime.now()

    cpuTime1 = time.clock()

    completePaths, incompletePaths, numRounds, acyclicFreq = findPathways(tdir = tdir, net = net, 
                                                  sources = sources, tgt = tgt, 
                                                  ag = ag, pairScores = pairScores, 
                                                  npaths = npaths, tracedSourceAtoms = tracedSourceAtoms, 
                                                  reactionDirConstraints = reactionDirConstraints, 
                                                  rpairConstraints = rpairConstraints,
                                                  tracedAtoms = tracedAtoms, maxdepth = maxdepth, 
                                                  greedyfinish = greedyfinish, 
                                                  prunedGraphSize = prunedGraphSize, 
                                                  storePartials = storePartials,
                                                  weighting = weighting)
    status("Pathfinding done: %d complete and %d incomplete paths found." % (len(completePaths), len(incompletePaths)))

    cpuTime2 = time.clock()  # note that this works only for short runs because of the overflow problem with clock()
    elapsedCpu = cpuTime2 - cpuTime1
    
    compTime = datetime.datetime.now() - startTime

    searchParameters["AcyclicFreq"] = acyclicFreq
    searchParameters["NumRounds"] = numRounds
    searchParameters["Sources"] = ",".join(sources)
    searchParameters["Targets"] = tgt

    paths = []
    for p in completePaths:
        paths.append((p, True))
    for p in incompletePaths:
        paths.append((p, False))
        
    pathways, numDiscarded = processPaths(paths, ag, sources, tgt, tracedAtoms, net,
                                          discardInvalidTargetMappings = False, 
                                          discardReactionsBothDirections = False, 
                                          discardMultipleTargetMappings = False)
    
    status("Discarded %d invalid pathways" % (numDiscarded))
    
    searchParameters["Discarded pathways"] = numDiscarded
    searchParameters["Pathways found"] = len(pathways)
    searchParameters["Computation time, total"] = str(compTime)
    searchParameters["CPU time, total"] = str(elapsedCpu)
    if len(pathways) > 0:
        searchParameters["Computation time, time/path"] = str(compTime / (len(pathways)))
        searchParameters["CPU time, time/path"] = str(elapsedCpu / (len(pathways)))
    else:
        searchParameters["Computation time, time/path"] = "N/A"
        searchParameters["CPU time, time/path"] = "N/A"

    searchParameters["PathwayAvgLen"], searchParameters["PathwayLenStdDev"] = calculatePathwayStatistics(pathways)

    prepareHtmlOutput(tdir)

    writePathwaysToHtml(tdir, sources, [tgt], pathways, compounds, ag, net, scores, searchParameters, scoreOrigin, disablePathOutput)

    status("All done.")

def usage():
    print """Usage: %s -o dir -s sources -t target -acdgkmnpqr
Find branching metabolic pathways connecting source metabolites to target.
See file COPYING for license notice.

Required options:
   -o, --outputdir      Directory where html output is written to.
   -s, --source         Source metabolites
   -t, --target         Target metabolite
   -d, --database       KEGG database directory

Optional options:
   -a, --tracedatoms    Atom types to trace, comma separated (e.g., "C,N") (default: %s)
   -c, --scores         Reaction scores
   -g, --greedy         Use greedy finish (numpaths=1 after first round)
   -i, --incomplete     Store also partial pathways (default: discard partials)
   -k, --numpaths       Number of paths in each step. (default: %s)
                        Give a comma-separated list for level-specific paths, e.g., 50,10,1.
   -m, --maxdepth       Maximum search depth (default: %d)
   -p, --prune          Prune atom graph to specific size (default: no pruning)
   -r, --dirlimits      Reaction direction constraints
   -w, --weights        Atom graph edge weights ((u)niform, (s)cores, (a)toms)
""" % (sys.argv[0], ",".join(DEFAULT_TRACED_ATOMS), ",".join(map(str, DEFAULT_NUMPATHS)), DEFAULT_MAXDEPTH)
    sys.exit(2)

def main():

    # Default parameters
    tdir = sources = tgt = scoresfn = dbdir = reactionDirFn = None 
    npaths = DEFAULT_NUMPATHS
    maxdepth = DEFAULT_MAXDEPTH
    tracedatoms = DEFAULT_TRACED_ATOMS
    greedyfinish = False
    storePartials = False
    prune = None
    weights = None
    disablePathOutput = False

    try:
        optlist, args = getopt.getopt(sys.argv[1:], "o:t:s:k:d:c:m:r:a:gp:iw:", ["outputdir=", "target=", "source=", "numpaths=", "maxdepth=", "scores=", "database=", "dirlimits=", "tracedatoms=", "greedyfinish", "prune=", "incomplete", "weights=", "disable-path-output"])
    except getopt.GetoptError:
        usage()
    for k, v in optlist:
        if k in ["-o", "--outputdir"]:
            tdir = v
        elif k in ["-s", "--source"]:
            sources = v
        elif k in ["-t", "--target"]:
            tgt = v
        elif k in ["-k", "--numpaths"]:
            npaths = map(int, v.split(","))
        elif k in ["-m", "--maxdepth"]:
            maxdepth = int(v)
        elif k in ["-c", "--scores"]:
            scoresfn = v
        elif k in ["-d", "--database"]:
            dbdir = v
        elif k in ["-r", "--dirlimits"]:
            reactionDirFn = v
        elif k in ["-a", "--tracedatoms"]:
            tracedatoms = v.split(",")
        elif k in ["-g", "--greedy"]:
            greedyfinish = True
        elif k in ["-p", "--prune"]:
            prune = int(v)
        elif k in ["-i", "--incomplete"]:
            storePartials = True
        elif k in ["-w", "--weights"]:
            weights = v
        elif k in ["--disable-path-output"]:
            disablePathOutput = True
 
    if not (tdir and sources and tgt and npaths and maxdepth and dbdir):
        usage()

    if not os.path.exists(tdir):
        os.mkdir(tdir)

    if not os.path.exists(dbdir):
        print "Database directory %s does not exist" % (dbdir)
        sys.exit(2)

    if not os.path.exists("%s/reaction" % (dbdir)):
        print "Database directory %s does not contain 'reaction' - not a KEGG database?"
        sys.exit(2)

    if scoresfn != None and not os.path.exists(scoresfn):
        print "Score file %s not found" % (scoresfn)
        sys.exit(2)
    
    retrace(tdir = tdir, sources = sources, tgt = tgt, 
            npaths = npaths, maxdepth = maxdepth, scoresfn = scoresfn, 
            rdir = dbdir, reactionDirFn = reactionDirFn, 
            tracedAtoms = tracedatoms, greedyfinish = greedyfinish, 
            prunedGraphSize = prune, 
            storePartials = storePartials,
            weights = weights,
            disablePathOutput = disablePathOutput)
if __name__ == '__main__':
    main()

