# -*- coding: iso-8859-15 -*-
"""
ReTrace - find branching pathways in metabolic networks
Copyright (C) 2009 Esa Pitkänen
See file COPYING for license notice.

Utility functions for tracepath.py and computePairwiseAtomDistances.py
"""

import sys, datetime

from graph.graph import Graph

def status(s):    
    sys.stdout.write("[%s] %s\n" % (datetime.datetime.now(), s))
    sys.stdout.flush()

TRACE_GRAPH_WEIGHTS_UNWEIGHTED = 0
TRACE_GRAPH_WEIGHTS_RSCORE = 1      # Use reaction scores from file
TRACE_GRAPH_WEIGHTS_NATOMS = 2      # Use 1/n where n is the number of atoms mapped in edge

SCORE_BASIS = 50.0

DUMMY_DIST = 1000.0
NONE_SCORE = 0

def convScore(L):
    """Converts a score to edge weight."""
    maxscore = 0
    for re, score in L:
        if score > maxscore:
            maxscore = score
    return 1.0 / (SCORE_BASIS + maxscore) * 100

def buildAtomTraceGraph(ag, tracedAtoms, pairScores = None, weighting = TRACE_GRAPH_WEIGHTS_UNWEIGHTED,
                        sources = None, target = None, minCommonAtoms = 1, skipnodes = set(),
                        reactionDirConstraints = {}):
    g = Graph()
    #ccc = 0
    for u in ag.edges:
        #ccc += 1
        #print ccc, len(ag.edges)
        types = ag.atomTypes[u]
        for v in ag.edges[u]:
            pairs = ag.edges[u][v]

            # Each molecule pair might have multiple rpairs
            for pair in pairs:

                amtypes = ag.atomMapTypes[pair]
                tracedCount = 0
                for type in amtypes:
                    if type in tracedAtoms:
                        tracedCount += 1
               
                if tracedCount < minCommonAtoms: # Does this rpair carry enough traced atoms?
                    continue 

                # Find out the score of this rpair
                if weighting == TRACE_GRAPH_WEIGHTS_RSCORE:
                    assert(pairScores)
                    if pair in pairScores:
                        #print pairScores[pair]
                        score = convScore(pairScores[pair])
                        #score = 1
                    else:
                        score = convScore(NONE_SCORE)
                elif weighting == TRACE_GRAPH_WEIGHTS_NATOMS:
                    #score = 2**(-tracedCount)     # As described in the 03/2009 manuscript
                    score = 1.0 / tracedCount
                    #score = 1.0
                elif weighting == TRACE_GRAPH_WEIGHTS_UNWEIGHTED:
                    score = 1.0
                else:
                    print "unknown TRACE_GRAPH_SCORE: %d" % (TRACE_GRAPH_SCORE)
                    assert(0)

                # Generate edges for each mapped atom pair of this rpair
                amap = ag.atomMap[pair]

                for src in amap:

                    if types[src] not in tracedAtoms:
                        continue

                    tgt = amap[src]

                    srcnode = "%s-%d" % (u, src)
                    tgtnode = "%s-%d" % (v, tgt)
                    fwdpair = "%s_f-%d-%d" % (pair, src, tgt)
                    revpair = "%s_r-%d-%d" % (pair, src, tgt)

                    if srcnode in skipnodes or tgtnode in skipnodes:
                        continue

                    # forward direction
                    g.addEdge(srcnode, fwdpair, score)
                    g.addEdge(fwdpair, tgtnode, 0.0)

                    # reverse direction
                    g.addEdge(tgtnode, revpair, score)
                    g.addEdge(revpair, srcnode, 0.0)

    # Add source and target root nodes, if necessary
    if sources != None:
        for source in sources: # source = mol
            g.addEdge("root", source, 0.0)
            for i in sources[source]: # i = mol-atom
                g.addEdge(source, i)
                if target != None:
                    g.addEdge(i, target, DUMMY_DIST)
    if target != None:
        for i in range(ag.atomCount[target]):
            g.addEdge("%s-%d" % (target, i + 1), target, 0.0)
            g.addEdge(target, "%s-%d" % (target, i + 1), DUMMY_DIST)
    return g                    

def removeConstrainedRpairs(ag, g, rpc, tracedAtoms):
    for rp in rpc:
        sub = ag.pairToSub[rp]
        pro = ag.pairToPro[rp]
        dir = rpc[rp]
        amap = ag.atomMap[rp]
        types = ag.atomTypes[sub]

        for src in amap:

            if types[src] not in tracedAtoms:
                continue

            tgt = amap[src]


            srcnode = "%s-%d" % (sub, src)
            tgtnode = "%s-%d" % (pro, tgt)

            if dir == "<" or dir == "-":
                fwdpair = "%s_f-%d-%d" % (rp, src, tgt)
                g.removeEdge(srcnode, fwdpair)
                g.removeEdge(fwdpair, tgtnode)
            if dir == ">" or dir == "-":
                revpair = "%s_r-%d-%d" % (rp, src, tgt)
                g.removeEdge(tgtnode, revpair)
                g.removeEdge(revpair, srcnode)
