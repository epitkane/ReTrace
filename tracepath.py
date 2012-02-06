#!/usr/bin/env python
#
# Find atom trace paths from source to target:
#
# 1. Remove from consideration the rpair mappings that are
#    too far away (in the sense of rpair scores) from both
#    source and target.
# 2. Find k shortest individual atom trace paths from source to target
#    (trace paths extracted from rpair mappings)
# 3. Find out how many atoms actually are transferred from source
#    to target via rpair maps used in the trace paths found in step 1
#
# This version can be called from tracer.py.
#
# Changes:
#    * Find trace paths in a graph with all rpairs removed that
#      do not carry at least z atoms/carbons through:
#      Start from z = min(#sourceatoms, #targetatoms), go down to 1 (maybe 2).
#    * Traced atom types can be specified with TRACED_ATOMS
#    * Traced atoms are highlighted in the result.
# TODO Find quickly the maximum number of atoms/carbons that
#      can be transferred from source to target. This is a baseline
#      result for k shortest path finding.
#
# Dependencies:
#    * Program ksp through graph.ksp (implements Yen's k shortest paths)
#    * MolImage JSP at sysdb.cs.helsinki.fi/tomcat/tkt_icom/MolImage/MoleculeImage
#    * KEGG molecule, reaction and rpair pictures at www.cs.helsinki.fi/group/icomic/kegg-pics
#
# See find_atom_path_first.py for the first attempt to do this.
#
# $Log: find_atom_path4.py,v $
# Revision 1.3  2008/06/15 19:41:29  epitkane
# Fixed a problem with reaction scoring leading into incorrect scores for
# most of the reactions. Combined pathway coverage scored now only
# by target atoms.
#
# Revision 1.2  2008/06/13 19:54:47  epitkane
# Fixed an issue with z upper bound -1. Fixed Pathway.__hash__.
#
# Revision 1.1  2008/06/09 10:47:37  epitkane
# First cvs version. Combining of traces into pathways not completely
# implemented: selection of good atom configurations is done, but
# selecting good reaction combinations within the atom configuration
# is not.
#
#

LAPTOP = 0

import sys, os, datetime, random, re

sys.path.append("/group/home/icomic/software/python/packages")

from metabolism.parser import kegg_ligand_parser as klp
from graph.graph import Graph
from graph.ksp import findKShortestPaths, checkKSP, KSP
from graph.shortestPaths import dijkstra
from ds.comb import comb, cartesian

# How many traces are combined at maximum? Values over 6 are bad, really baaaad
MAX_TRACES_TO_COMBINE = 4

# How many best pathways are generated for each network size?
MAX_PATHWAYS_PER_NETWORK_SIZE = 5

# Maximum coverage score difference allowed in generated reaction combinations
COVERAGE_SCORE_MARGIN = 1

# How high coverage score a map has to have to be included
COMBINE_THRESHOLD = 2  # at least one source, one target => at least two atoms

MAX_REACTION_COMBINATION_TO_SHOW = 10

SCORE_BASIS = 500.0
NONE_SCORE = 0

DUMMY_DIST = 1000.0

IGNORE_NONCARBONS = True

TRACED_ATOMS = set(["C"])

#REQN_REACTIONS = 400
#K = 100

#REQN_REACTIONS = 2000
#K = 100

REQN_REACTIONS = 200
#K = 20

MAX_PATH_TO_SHOW = 40

# z: Minimum number of atoms in rpair required
# Algorithm performs a sweep COMMON_ATOM_REQ_LOWER <= z <= COMMON_ATOM_REQ_UPPER
#COMMON_ATOM_REQ_LOWER = 3
#COMMON_ATOM_REQ_UPPER = None   # If None, use min(sourceatoms, targetatoms)
#COMMON_ATOM_REQ_UPPER = 3   # If None, use min(sourceatoms, targetatoms)

#MOLVIEWER="http://sysdb.cs.helsinki.fi/MolImage/MoleculeImage"
#MOLVIEWER="http://localhost:4444/MolImage/MoleculeImage"

MOLVIEWER="http://sysdb.cs.helsinki.fi/tomcat/tkt_icom/MolImage/MoleculeImage"
HTMLBASE = "http://www.cs.helsinki.fi/group/icomic/kegg-pics"

#SCORESFN = "/home/epitkane/tracer/data/treesei_scores_metacyc_to_kegg.txt"
#SCORESFN = "/group/home/icomic/data/reconstruction/scores/treesei_scores_metacyc_to_kegg.txt"

if LAPTOP:
    KEGGDIR = "/home/epitkane/work/data/kegg/2008.05.19"
else:
    KEGGDIR = "/group/home/icomic/data/kegg/ligand/LATEST/"


class AtomConfiguration(dict):
    def __init__(self, src, tgt):
        self.src = src
        self.tgt = tgt

    def getTargetMoleculeAtomPos(self):
        a = set()
        for atom in iter(self):
            a.add("%s-%s" % (self.tgt, atom))
        return a
    
    def getSourceMoleculeAtomPos(self):
        a = set()
        for atom in iter(self):
            a.add("%s-%s" % (self.src, self.get(atom)))
        return a
                  
    def __hash__(self):
        s = "%s=%s=" % (self.src, self.tgt)
        for a in iter(self):
            s += "%s-%s:" % (self.get(a), a)
        return hash(s)

    def __str__(self):
        s = "%s, %s: " % (self.src, self.tgt)
        for a in iter(self):
            s += "%s->%s, " % (self.get(a), a)
        return s.strip(", ")
                  
class Pathway:

    def __init__(self, src, tgt, transferredAtoms, origin, score, rpairs, reactions, substrates, tracedAtoms):
        self.src = src
        self.tgt = tgt
        self.transferredAtoms = transferredAtoms
        self.origin = origin
        self.score = score
        self.rpairs = rpairs
        self.reactions = reactions
        self.substrates = substrates
        self.tracedAtoms = tracedAtoms        

    def __hash__(self):
        "Hash by reactions, reaction order is irrelevant"
        objs = set()
        for p in self.reactions:
            objs.add(p)
        s = "%s=%s=" % (self.src, self.tgt)
        keys = list(objs)
        keys.sort()
        for o in objs:
            s += "%s-" % (o)

        return hash(s)

def cmpTrace(x, y):
    """Sort primarily by the number of transferred atoms, secondarily
    by the number of reactions, thirdly by score."""
    if len(x.transferredAtoms) < len(y.transferredAtoms):
        return 1
    elif len(x.transferredAtoms) > len(y.transferredAtoms):
        return -1
    elif len(x.reactions) < len(y.reactions):
        return -1
    elif len(x.reactions) > len(y.reactions):
        return 1
    else:
        return int(y.score / len(y.reactions) - x.score / len(x.reactions))

def hashReactionSet(rs):
    s = ""
    for r in rs:
        s += "%s===" % (r)
    return hash(s)

def status(s):    
    sys.stdout.write("[%s] %s\n" % (datetime.datetime.now(), s))
    sys.stdout.flush()

def convScore(sc):
    #return 1.0 / (SCORE_BASIS + sc) * 100
    return 1.0

def reverseMap(m):
    rmap = {}
    for atom in m:
        rmap[m[atom]] = atom
    return rmap

def findCutoff(sdists, tdists, reqn):
    totalDists = []
    for u in sdists:
        if u not in tdists:
            continue
        d = sdists[u] + tdists[u] - DUMMY_DIST
        totalDists.append(d)
    totalDists.sort()
    return totalDists[reqn]

def buildAtomTraceGraph(ag, pairScores, source, target, minCommonAtoms, skipnodes = set()):
    g = Graph()
    for u in ag.edges:
        #if u in skipnodes:
        #    continue
        types = ag.atomTypes[u]
        for v in ag.edges[u]:
            #if v in skipnodes:
            #    continue
            pairs = ag.edges[u][v]

            # Each molecule pair might have multiple rpairs
            for pair in pairs:

                amtypes = ag.atomMapTypes[pair]
                tracedCount = 0
                for type in amtypes:
                    if type in TRACED_ATOMS:
                        tracedCount += 1
                
                if tracedCount < minCommonAtoms: # Does this rpair carry enough traced atoms?
                    continue 

                # Find out the score of this rpair
                if pair in pairScores and pairScores[pair] != None:
                    score = convScore(pairScores[pair])
                else:
                    score = convScore(NONE_SCORE)

                # Generate edges for each mapped atom pair of this rpair
                amap = ag.atomMap[pair]

                for src in amap:

                    if types[src] not in TRACED_ATOMS:
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

    #for i in range(ag.atomCount[source]):
    #    g.addEdge(source, "%s-%d" % (source, i + 1))
    for i in range(ag.atomCount[target]):
        g.addEdge("%s-%d" % (target, i + 1), target, 0.0)
        g.addEdge(target, "%s-%d" % (target, i + 1), DUMMY_DIST)
    return g                    

def buildGraph(ag, skipnodes = set()):

    g = Graph()
    for u in ag.edges:
        if u in skipnodes:
            continue
        for v in ag.edges[u]:
            if v in skipnodes:
                continue

            # add an edge for each map between u and v

            pairs = ag.edges[u][v]

            for pair in pairs:
                if pair in ag.carbonCount:
                    score = 1.0 / (ag.carbonCount[pair] + 1)
                else:
                    score = 1.0

                #if pair in pairScores and pairScores[pair] != None:
                #    #score = convScore(pairScores[pair])
                #    score = ag.carbonCount[pair]
                #else:
                #    #score = convScore(NONE_SCORE)
                g.addEdge(u, pair, score)
                g.addEdge(v, pair, score)
                g.addEdge(pair, v, 0.0)
                g.addEdge(pair, u, 0.0)

    return g

def getDirectedAtomMap(ag, pair, dir):
    if dir == "f":
        atomMap = ag.atomMap[pair]
    elif dir == "r":
        atomMap = reverseMap(ag.atomMap[pair])
    else:
        raise RuntimeError, "Unknown dir: %s" % (dir)
    return atomMap        

def writeTracesToFile(tdir, sources, tgt, resultPathways, compounds, scores, ag, maxPath):
    "Write a html result file of individual trace paths"

    of = open("%s/traces-%s-to-%s.html" % (tdir, "-".join(sources), tgt), "w")
    of.write("""
<html><head><title>Atom traces from %s to %s</title></head><body>\n""" % (",".join(sources), tgt))

    tgtname = tgt
    if tgt in compounds:
        tgtname = compounds[tgt].names[0].strip()
    
    srcnames = {}
    for src in sources:
        if src in compounds:
            srcnames[src] = compounds[src].names[0].strip()
        else:
            srcnames[src] = src

    of.write("<p><b>Sources:</b>")
    for src in sources:
        of.write("""%s (<a href=\"http://www.genome.jp/dbget-bin/www_bget?compound+%s\">%s</a>)
<img src=\"%s?mol=%s&show_carbon_indices=true&highlight=%s\">\n
""" % (srcnames[src], src, src, MOLVIEWER, src, "0"))
    of.write("""<p><b>Target:</b>
%s (<a href=\"http://www.genome.jp/dbget-bin/www_bget?compound+%s\">%s</a>)
<img src=\"%s?mol=%s&show_carbon_indices=true&highlight=%s\"></p>\n
""" % (tgtname, tgt, tgt, MOLVIEWER, tgt, "0"))

    # Report a summary of paths

    of.write("<br>\n<table border=\"1\"><tr><th>Path</th><th>#reactions</th><th>Total score</th><th>Avg. score</th><th>Source atoms</th><th>Target atoms</th></tr>\n")
    i = 1
    cmpUse = {}
    allSubs = set()
    for path in resultPathways:
        cmpUse[i] = set()
        for sub in path.substrates:
            subname = sub.split("-")[0]
            allSubs.add(subname)
            cmpUse[i].add(subname)

        totalScore = 0.0
        for react in path.reactions:
            if react in scores and scores[react] != "?":
                totalScore += float(scores[react])

        satoms = []
        tatoms = []
        for atom in path.origin:
            satoms.append(path.origin[atom])
            tatoms.append(atom)

        of.write("<tr><td><a href=\"#path%d\">Path %d</td><td>%d reactions</td><td>%.2f</td><td>%.2f</td><td>%s</td><td>%s</td></tr>\n" % (i, i, len(path.reactions), totalScore, totalScore / len(path.reactions), ",".join(map(str, satoms)), ",".join(map(str, tatoms))))
        i += 1
        #if i == maxPath:
        #    break

    of.write("</table>\n")
    
    # Show metabolite usage table
    #
    
    of.write("<h2>Metabolite usage</h2>\n")
    of.write("<table border=\"1\">\n")
    of.write("<tr><td>&nbsp;</td>")
    keys = list(allSubs)
    keys.sort()
    for k in keys:
        of.write("<td>%s</td>" % (compounds[k].names[0]))
    of.write("\n")

    for i in range(1, len(resultPathways) + 1):
        of.write("<tr><td><a href=\"#path%d\">Path %d</a></td>" % (i, i))
        for j in range(len(keys)):
            if keys[j] in cmpUse[i]:
                mark = "X"
            else:
                mark = "&nbsp;"
            of.write("<td align=\"center\">%s</td>" % (mark))
        of.write("</tr>\n")
    of.write("</table>\n")

    # Show detailed descriptions for each pathway
    # 

    i = 1
    for path in resultPathways:
        #tAtoms, origin, score, rpairs, reactions, substrates, tracedAtoms = path
        of.write("<br><h2 id=\"path%d\">Path %d</h2>\n<p>\n" % (i, i))
        src = path.src

        tgthighlights = []
        srchighlights = []
        for atom in path.origin:
            tgthighlights.append(str(atom))
            srchighlights.append(str(path.origin[atom]))

        of.write("""<p><b>Source:</b>
%s <img src=\"%s?mol=%s&show_carbon_indices=true&highlight=%s\">\n
""" % (srcnames[src], MOLVIEWER, src, ",".join(srchighlights)))
        of.write("""<b>Target:</b>
%s <img src=\"%s?mol=%s&show_carbon_indices=true&highlight=%s\">\n
""" % (tgtname, MOLVIEWER, tgt, ",".join(tgthighlights)))

        for atom in path.origin:
            of.write("<br>%s -> %s\n" % (path.origin[atom], atom))
        of.write("\n<table border=\"1\"><tr><th>Reaction id</th><th>Score</th><th>Substrate</th><th>Reaction</th><th>RPair</th></tr>\n")
        #print list(tAtoms), origin, score, len(reactions), reactions, rpairs
        
        # Show individual reaction steps and substrates for
        # this pathway

        #track = 0
        for ix in range(len(path.reactions)):
            react = path.reactions[ix]
            rpair = path.rpairs[ix]

            #if rpair == "A00021":
            #    print "pair found", i
            #    track = 1

            sub = path.substrates[ix]
            submol = sub.split("-")[0]
            subname = sub
            if submol in compounds:
                subname = compounds[submol].names[0].strip()
            score = "?"
            if react in scores:
                score = scores[react]

            atomHighlights = ",".join(map(str, list(path.tracedAtoms[ix])))

            of.write("""
<tr><td>%s (<a href=\"http://www.genome.jp/dbget-bin/www_bget?compound+%s\">%s</a>)
 -> <a href=\"http://www.genome.jp/dbget-bin/www_bget?rn+%s\">%s</a></td>
<td>%s</td>
<td><img src=\"%s?mol=%s&show_carbon_indices=true&highlight=%s\"></td>
<td><img src=\"%s/reaction_gif_small/%s.gif\"></td>
<td><img src=\"%s/rpair_gif_small/%s.gif\"></td></tr>\n
""" % (subname, submol, sub, react, react, score, 
       MOLVIEWER, submol, atomHighlights,
       HTMLBASE, react,
       HTMLBASE, rpair))

        #if track:
        #    print path.origin
        #    print path.score
        #    print path.transferredAtoms
        #    print path.tracedAtoms

        # Show the product in the end of the table
        of.write("""
<tr><td>%s (<a href=\"http://www.genome.jp/dbget-bin/www_bget?compound+%s\">%s</a>)</td>
<td>&nbsp;</td>
<td><img src=\"%s?mol=%s&show_carbon_indices=true&highlight=%s\"></td>
<td>&nbsp;</td>
<td>&nbsp</td></tr>\n
""" % (tgtname, tgt, tgt, MOLVIEWER, tgt, ",".join(tgthighlights)))

        of.write("</table></p>\n")
        
        i += 1

        #if i == maxPath:
        #    break
    of.write("</body></html>\n")
    of.close()

def writePathwaysToFile(tdir, sources, tgt, combos, compounds, scores, ag, 
                        atomConfigs, resultPathways):
    "Write a html result file of combination pathways"

    # combos: a list of (acHashes, coverageScore, {network size -> [(reactions, totalScore, avgScore)]})

    srcnames = {}
    for src in sources:
        if src in compounds:
            srcnames[src] = compounds[src].names[0].strip()
        else:
            srcnames[src] = src

    tgtname = compounds[tgt].names[0].strip()

    of = open("%s/pathways-%s-to-%s.html" % (tdir, "-".join(sources), tgt), "w")
    of.write("""
<html><head><title>Pathways from %s to %s</title></head><body>\n""" % (",".join(sources), tgt))

    of.write("<p><b>Sources:</b> ")
    for src in sources:
        of.write("""%s (<a href=\"http://www.genome.jp/dbget-bin/www_bget?compound+%s\">%s</a>)
<img src=\"%s?mol=%s&show_carbon_indices=true&highlight=%s\">\n
""" % (srcnames[src], src, src, MOLVIEWER, src, "0"))
    of.write("""<p><b>Target:</b>
%s (<a href=\"http://www.genome.jp/dbget-bin/www_bget?compound+%s\">%s</a>)
<img src=\"%s?mol=%s&show_carbon_indices=true&highlight=%s\"></p>\n
""" % (tgtname, tgt, tgt, MOLVIEWER, tgt, "0"))

    #
    # 1. Write a list of all possible atom map configurations with
    # links to detailed list of different reaction combinations for each configuration
    # - combinedPathway is already sorted by the score of configuration (source and target atom coverage)
    of.write("<table border=\"1\"><tr><th>&nbsp;</th><th>Source atoms</th><th>Target atoms</th><th>Coverage Score</th><th>Network sizes</th></tr>\n")

    i = 1
    for path in combos:    # path: a combination of atom configs (found in resultPathways)
        acHashes, coverageScore, pathways = path  # acHashes: indexes atomConfigs, each atomConfig holds one or more (alternative) trace
        sa = []
        ta = []

        networkSizes = pathways.keys()
        networkSizes.sort(lambda x, y: int(x) - int(y))

        numberOfNetworks = 0
        pruneNetworks = None
        for nr in networkSizes:
            numberOfNetworks += len(pathways[nr])
            if numberOfNetworks > MAX_REACTION_COMBINATION_TO_SHOW:
                pruneNetworks = nr
                break
        if pruneNetworks:
            networkSizes = networkSizes[0:pruneNetworks]

        networkSizesStr = "%d-%d" % (networkSizes[0], networkSizes[-1])
        
        for acHash in acHashes:
            # origin: source, target, target atom positions -> source atom positions
            # traceIxes: list of indecies to resultPathways

            origin, traceIxes = atomConfigs[acHash]

            tgthighlights = []
            srchighlights = []

            for atom in origin:
                sa.append("%s-%d" % (origin.src, origin[atom]))
                ta.append("%s-%d" % (tgt, atom))

        of.write("<tr><td><a href=\"#map%d\">Map %d</a></td><td>%s</td><td>%s</td><td>%.2f</td><td>%s</td></tr>\n" % (i, i, ",".join(sa), ",".join(ta), coverageScore, networkSizesStr))

        i += 1
               
    of.write("</table>\n")
#
# copy'n'paste
#
#

    i = 1
    for path in combos:    # path: a combination of atom configs (found in resultPathways)
        acHashes, coverageScore, pathways = path  # acHashes: indexes atomConfigs, each atomConfig holds one or more (alternative) trace
        combinedReactions = set()
        sa = []
        ta = []

        networkSizes = pathways.keys()
        networkSizes.sort(lambda x, y: int(x) - int(y))

#         numberOfNetworks = 0
#         pruneNetworks = None
#         for nr in networkSizes:
#             numberOfNetworks += len(pathways[nr])
#             if numberOfNetworks > MAX_REACTION_COMBINATION_TO_SHOW:
#                 pruneNetworks = nr
#                 break
#         if pruneNetworks:
#             networkSizes = networkSizes[0:pruneNetworks]

#        networkSizesStr = "%d-%d" % (networkSizes[0], networkSizes[-1])
        
        bestScore = bestScoreNr = None
        for nr in networkSizes:
            if not bestScore or pathways[nr][0][2] > bestScore:
                bestScore = pathways[nr][0][2]
                bestScoreNr = nr

        of.write("""
<h2 id=\"map%d\">Map %d</h2>
<p>
This pathway map contains %d component atom map(s) with total coverage score of <b>%.1f</b>.
<br>
Smallest pathway found for this map contains <b>%d reactions</b> (pathway score is %.2f).
<br>
The best pathway score found for this map is <b>%.2f</b> (pathway has %d reactions).
<br></p>""" % (i, i, len(acHashes), coverageScore, networkSizes[0], pathways[networkSizes[0]][0][2], bestScore, bestScoreNr))
        
        of.write("""
<table border="1">
<tr><th>Source</th><th>Target</th></tr>
""")
        for acHash in acHashes:
            # origin: source, target, target atom positions -> source atom positions
            # traceIxes: list of indecies to resultPathways

            origin, traceIxes = atomConfigs[acHash]

            tgthighlights = []
            srchighlights = []
            for atom in origin:
                tgthighlights.append(str(atom))
                srchighlights.append(str(origin[atom]))

            of.write("""<tr>
<td>%s (<a href=\"http://www.genome.jp/dbget-bin/www_bget?compound+%s\">%s</a>)
<img src=\"%s?mol=%s&show_carbon_indices=true&highlight=%s\">
</td>
<td>%s (<a href=\"http://www.genome.jp/dbget-bin/www_bget?compound+%s\">%s</a>)
<img src=\"%s?mol=%s&show_carbon_indices=true&highlight=%s\">
</td></tr>\n
""" % (srcnames[origin.src], origin.src, origin.src, MOLVIEWER, origin.src, ",".join(srchighlights),
       tgtname, origin.tgt, origin.tgt, MOLVIEWER, origin.tgt, ",".join(tgthighlights)))

            for atom in origin:
                sa.append("%s-%d" % (origin.src, origin[atom]))
                ta.append("%s-%d" % (tgt, atom))

        of.write("</table>")  # close map table

        of.write("""
<p>
The following pathways transfer atoms according to this map.
<br>
</p>
""")

        of.write("<table border=\"1\"><tr><th>#Reactions</th><th>Total score</th><th>Avg. score</th><th>Traces</th><th>Reactions</th></tr>")

        for nr in networkSizes:
            for reactionCombination in pathways[nr]:
                reactions, totalComboScore, avgComboScore, traceCombo = reactionCombination

                traceCombo = list(traceCombo)
                traceCombo.sort(lambda x, y: x - y)
                traceStr = ""
                for trace in traceCombo:
                    traceStr += "<a href=\"traces-%s-to-%s.html#path%d\">%d</a>," % ("-".join(sources), tgt, trace + 1, trace + 1)
                traceStr = traceStr.rstrip(",")

                assert(nr == len(reactions))
                of.write("<tr><td>%d</td><td>%.2f</td><td>%.2f</td><td>%s</td><td>" % (nr, totalComboScore, avgComboScore, traceStr))
                reKeys = list(reactions)
                reKeys.sort()
                sst = ""
                for react in reKeys:
                    sst += "<a href=\"http://www.genome.jp/dbget-bin/www_bget?rn+%s\">%s</a>," % (react, react)
                of.write(sst.rstrip(","))
                of.write("</td></tr>\n")

        of.write("</table>\n</p>\n")

        i += 1
               
#
#
# copy'n'paste END
#
#


    #
    # 2. Write reaction combinations for each atom map configuration
    # - Sort by number of reactions and/or score
    #

    of.write("</body>\n</html>\n")
    of.close()

def findTraces(tdir, sources, tgt, ag, commonAtomReqUpper, commonAtomReqLower, pairScores, npaths, pairScoreReactions, tracedSourceAtoms):

    resultPathways = []

    tracedTgtAtoms = set()
    for i in range(1, len(ag.atomTypes[tgt]) + 1):
        type = ag.atomTypes[tgt][i]
        if type in TRACED_ATOMS:
            tracedTgtAtoms.add(i)

    for src in sources:

        status("Looking for paths from %s to %s" % (src, tgt))

        tracedSrcAtoms = set()
        for i in range(1, len(ag.atomTypes[src]) + 1):
            type = ag.atomTypes[src][i]
            if src in tracedSourceAtoms:
                if i in tracedSourceAtoms[src]:
                    tracedSrcAtoms.add(i)
            elif type in TRACED_ATOMS:
                tracedSrcAtoms.add(i)

        status("Tracing %d/%d atoms in %s, %d/%d atoms in %s" % (len(tracedSrcAtoms), ag.atomCount[src], src, len(tracedTgtAtoms), ag.atomCount[tgt], tgt))

        minNumAtoms = min(len(tracedSrcAtoms), len(tracedTgtAtoms))
        carUp = commonAtomReqUpper
        if carUp > minNumAtoms or carUp == -1:
            carUp = minNumAtoms
        if carUp < commonAtomReqLower:
            carUp = commonAtomReqLower

        # Build a rpair graph to find out rpair maps that are too
        # far away from both source and target nodes

    #    src = "C00031"  # glucose
    #    tgt = "C00251"  # chorismate
    #    tgt = "C00041"  # L-alanine
    #    tgt = "C00135"  # Histidine, 1 carbon gets through
    #    tgt = "C00097"  # cysteine
    #    tgt = "C00049"  # Asparatate
    #    tgt = "C00123"   # Leucine

        #minNumAtoms = 3
        #MINIMUM_NUMBER_OF_COMMON_ATOMS = 1

        #
        # Find trace paths
        #

        status("Looking for a minimum of %d...%d traced atoms" % (commonAtomReqLower, carUp))

        for minCommonAtoms in range(carUp, commonAtomReqLower - 1, -1):

            status("Building atom trace graph for rpairs with z >= %d" % (minCommonAtoms))

            g = buildAtomTraceGraph(ag, pairScores, src, tgt, minCommonAtoms, skipnodes = set())

            (tgtDists, tgtPaths) = dijkstra(g, tgt)

            # seek out atom trace paths for each position of source

            for pos in tracedSrcAtoms:
            #for pos in range(1, ag.atomCount[src] + 1):
            #    if ag.atomTypes[src][pos] not in TRACED_ATOMS:
            #       continue

                st = datetime.datetime.now()

                node = "%s-%d" % (src, pos)
                (srcDists, srcPaths) = dijkstra(g, node)

                if "%s-%d" % (src, pos) not in tgtDists:
                    status("Can't find %s-%d in tgtdists" % (src, pos))
                    continue

                total = 0.0
                for d in srcDists:
                    total += srcDists[d]
                status("%d: %d nodes reachable, %.2f avg dist, %.2f tgt dist" % (pos, len(srcDists), total / len(srcDists), tgtDists["%s-%d" % (src, pos)]))
                #cutoff = total / len(srcDists) * 0.70  # fancy heuristic
                cutoff = findCutoff(srcDists, tgtDists, REQN_REACTIONS)
                status("Cutoff: %d" % (cutoff))
                removals = set()
                for u in g.V:
                    if u not in srcDists or u not in tgtDists or srcDists[u] + tgtDists[u] - DUMMY_DIST > cutoff:
                        removals.add(u)
                g2 = buildAtomTraceGraph(ag, pairScores, src, tgt, minCommonAtoms, removals)
                status("%d nodes in g, %d nodes in g2" % (len(g.V), len(g2.V)))

                # find k shortest paths from source-pos node to target (any position)
                k = npaths
                status("Finding %d shortest paths from %s to %s" % (k, node, tgt))
                tempfn = "%s/temp_ksp%d.txt" % (tdir, random.randint(10**5, 10**6))
                paths = findKShortestPaths(g2, k, node, tgt, tempfn)

                if paths == None:
                    status("KSP failed!")
                    sys.exit(0)

                # Find for each atom map pair path how many atoms actually are
                # transferred on the path
                ix = 1
                track = 0 # DEBUGGGGGGGGGGGGGGGGGGG
                for p in paths:
 #                   if ix == 1:
 #                       print "%d: %s" % (ix, p)
                    pair = p[1].split("-")[0]
                    pair, dir = pair.split("_")
                    amap = getDirectedAtomMap(ag, pair, dir)
                    transferredAtoms = set(amap.keys())  # what atoms are transferred in the first map?
                    origin = AtomConfiguration(src, tgt)   # atom pos in current step -> atom pos in source
                    for atom in transferredAtoms:
                        origin[atom] = atom
#                    if ix == 1:
#                        status("  Initial atoms: %s (%d in total)" % (transferredAtoms, len(transferredAtoms)))
#                        print "Origin: ", origin

                    reactions = []
                    rpairs = []
                    substrates = []
                    score = 0.0
                    tracedAtoms = []
                    # p[1] is the first pair, p[3] next, ...
                    for i in range(1, len(p) - 1, 2):
                        neworigin = AtomConfiguration(src, tgt)

                        pair, srcatom, tgtatom = p[i].split("-")        # e.g. (A06006_f, 1, 7)
                        pair, dir = pair.split("_")                     # e.g. (A06006, f)


                        #reactions.append(ag.pairToReaction[pair])
                        reactions.append(pairScoreReactions[pair])
                        substrates.append(p[i - 1])
                        tracedAtoms.append(transferredAtoms)
                        rpairs.append(pair)
                        if pair in pairScores and pairScores[pair] != None:
                            score += convScore(pairScores[pair])
                        else:
                            score += convScore(NONE_SCORE)

#                        if pair == "A00021":
#                            print "pair found:", p, len(resultPathways)
#                            print "score:", score
#                            track = 1

                        amap = getDirectedAtomMap(ag, pair, dir)
                        mappedAtoms = set()
                        for atom in transferredAtoms:
                            if atom in amap:
                                mappedAtoms.add(amap[atom])             # find out where this atom goes
                                neworigin[amap[atom]] = origin[atom]
                        transferredAtoms = mappedAtoms

                        #if track:
                        #    print "transferred: ", transferredAtoms


                        origin = neworigin
#                        if ix == 1:
#                            print "Origin: ", origin
#                            status("  Step %d: atoms transferred: %s (%d in total)" % (i, transferredAtoms, len(transferredAtoms)))

#                     if ix == 1:
#                         print "TA before removals:", transferredAtoms
#                         print "tgt atom types:", ag.atomTypes[tgt]
#                         print "Traced atoms:", TRACED_ATOMS

                    # remove non-carbons, if necessary
                    removals = set()
                    for i in transferredAtoms:
                        if ag.atomTypes[tgt][i] not in TRACED_ATOMS:
                            removals.add(i)
                            del origin[i]
                    transferredAtoms.difference_update(removals)

#                     if ix == 1:
#                         print "TA after removals:", transferredAtoms
#                         print "Removals:", removals

#                         print ag.atomMap["A00021"]
#                         print ag.atomMapTypes["A00021"]

#                        sys.exit(1)

                    #if track:
                    #    print "Adding:", Pathway(src, tgt, transferredAtoms, origin, score, rpairs, reactions, substrates, tracedAtoms)

                    #status("Pathway carries %d atoms" % (len(transferredAtoms)))
                    resultPathways.append(Pathway(src, tgt, transferredAtoms, origin, score, rpairs, reactions, substrates, tracedAtoms))

                    ix += 1

                status("That took %s time and now we have %d pathways" % (datetime.datetime.now() - st, len(resultPathways)))

    return resultPathways

def combineTraces(resultPathways, sources, tgt, ag):
    "Combine trace paths into trace networks"

    atomConfigs = {}
    for i in range(len(resultPathways)):
        p = resultPathways[i]
        h = hash(p.origin)
        if h not in atomConfigs:
            atomConfigs[h] = (p.origin, [])
        atomConfigs[h][1].append(i)

    # print number of different combinations
    for c in atomConfigs:
        #for i in atomConfigs[c]:
        #    print resultPathways[i].origin, i
        print len(atomConfigs[c][1]), atomConfigs[c][0]
    
    combinedPathways = []  # holds ([path hashes], score) pairs

    # generate all n-combinations of trace paths, n = 1,2,...
    keys = atomConfigs.keys()

    maxTracesToCombine = min(len(atomConfigs), MAX_TRACES_TO_COMBINE)
    
    for i in range(1, maxTracesToCombine):
        status("Finding %d trace combinations..." % (i))
        for c in comb(keys, i):       # c: a combination of i traces
            #print i, c
            # find out how well these traces combined cover the source and target atoms
            coverage = {}
            for pathhash in c:

                origin = atomConfigs[pathhash][0]

                # Score the combined path

                for atom in origin.getTargetMoleculeAtomPos():
                    if atom not in coverage:
                        coverage[atom] = 1
                    else:
                        coverage[atom] += 1
                # Ignore source atoms
                #for atom in origin.getSourceMoleculeAtomPos():
                #    if atom not in coverage:
                #        coverage[atom] = 1
                #    else:
                #        coverage[atom] += 1
                  
            #print "i=%s, c=%s, Coverage: %s" % (i, c, coverage)
            # Score this pathway
            score = 0.0
            # Score positions that have been mapped at least once
            for atom in coverage:
                if coverage[atom] == 1:
                    score += 1
                elif coverage[atom] > 1:
                    score -= 100
                #score += 2 - coverage[atom]   # score = 1, 0, -1, ..., score = 1 when position is mapped once
                #print "Scored mapped pos: ", atom, "score now", score
            # Score positions that have not been mapped
            # Ignore source atoms
#             for mol in sources:
#                 #print "Checking", mol, ", counts:", ag.atomCount[mol], ", types:", ag.atomTypes[mol]
#                 for pos in range(1, ag.atomCount[mol] + 1):
#                     if ag.atomTypes[mol][pos] in TRACED_ATOMS and "%s-%d" % (mol, pos) not in coverage:
#                         score -= 1
#                         #print "Scores not mapped pos:", mol, pos, "score:", score
            for pos in range(1, ag.atomCount[tgt] + 1):
                if ag.atomTypes[tgt][pos] in TRACED_ATOMS and "%s-%d" % (tgt, pos) not in coverage:
                    score -= 1
                    #print "Scores not mapped pos:", tgt, pos, "score:", score
            if score > COMBINE_THRESHOLD:
                combinedPathways.append((c, score))

    return combinedPathways, atomConfigs

def chooseReactionCombinations(combinedPathways, atomConfigs, tracePathways, scores):
    """For the best combined pathways (=atom configurations with a good coverage),
choose good reaction combinations."""
    # combinedPathways: list of ([atom config hashes], coverage score)
    # atomConfigs: atom config hash -> (origin, [trace indecies])

    bestScore = None
    combos = []
    i = 0
    for cp in combinedPathways:                 # combinedPathway: combination of atom configs
        acHashes, coverageScore = cp

        if bestScore:
            assert(coverageScore <= bestScore)
            if coverageScore < bestScore - COVERAGE_SCORE_MARGIN:
                continue   # this atom config does not high enough coverage score so we skip it
        else:
            bestScore = coverageScore

        #status("CombinedPathway: %s acHashes, %s cov. score" % (len(acHashes), coverageScore))

        # for each atom config, pick a trace such that the total score
        # of traces combined is high and/or the total number of reactions is low
        allTraces = []
        for acHash in acHashes:
            origin, traces = atomConfigs[acHash]
            #status(origin)
            allTraces.append(traces)

        #status("allTraces: %s " % (allTraces))

        traceCombos = cartesian(*allTraces)

        #status("traceCombos: %s " % (traceCombos))

        comboScores = {}   # network size -> network
        
        combReactHashes = set()
        for traceCombo in traceCombos:          # traceCombo: one trace from each atom config
            combinedReactions = set()
            for traceix in traceCombo:
                trace = tracePathways[traceix]
                for react in trace.reactions:
                    combinedReactions.add(react)
            # what is this combo's reaction score and number of reactions?
            totalScore = 0.0
            for react in combinedReactions:
                if react in scores and scores[react] != "?":
                    totalScore += float(scores[react])
            nr = len(combinedReactions)
            if nr not in comboScores:
                comboScores[nr] = []
            h = hashReactionSet(combinedReactions)
            if h in combReactHashes:
                continue  # duplicate reaction set
            combReactHashes.add(h)
            comboScores[nr].append((combinedReactions, totalScore, totalScore / len(combinedReactions), traceCombo))
       
        for cs in comboScores:
            #status("  #reactions = %d" % (cs))
            comboScores[cs].sort(lambda x, y: int(y[2] - x[2]))
            comboScores[cs] = comboScores[cs][0:MAX_PATHWAYS_PER_NETWORK_SIZE]
            #for path in comboScores[cs]:
            #    status("    Total score: %s, avg. score: %s, reactions: %s" % (path[1], path[2], path[0]))

        combos.append((acHashes, coverageScore, comboScores))

    return combos

def main():
    if not checkKSP():
        print "Problem with KSP"
        sys.exit(1)

    # params: tdir, source, target, k
    tdir = sys.argv[1]          # directory to store results in
    sources = sys.argv[2]           # source metabolite(s), comma separated
    tgt = sys.argv[3]           # target metabolite
    npaths = int(sys.argv[4])   # number of traces / atom position
    maxPath = int(sys.argv[5])  # max paths to show in results
    commonAtomReqLower = sys.argv[6]  # min z (see header)
    commonAtomReqUpper = sys.argv[7]  # max z (see header)
    scoresfn = sys.argv[8]      # scores file: kegg id -> score

    if not re.match("\d+", commonAtomReqLower):
        commonAtomReqLower = 2
    else:
        commonAtomReqLower = int(commonAtomReqLower)
        if commonAtomReqLower < 1:
            commonAtomReqLower = 1

    if not re.match("\d+", commonAtomReqUpper):
        commonAtomReqUpper = -1  # changed later to match a valid number
    else:
        commonAtomReqUpper = int(commonAtomReqUpper)

    status("Target dir: %s" % (tdir))
    status("Source(s): %s" % (sources))
    status("Target: %s" % (tgt))
    status("Traces / position: %d" % (npaths))
    status("Max paths to show: %d" % (maxPath))
    status("Least number of common atoms in rpair: %d" % (commonAtomReqLower))
    status("Highest number of common atoms in rpair: %d" % (commonAtomReqUpper))
    status("Scores file: %s" % (scoresfn))

    sources = sources.split(",")

    rdir = KEGGDIR
    scoresf = open(scoresfn)           # score file: 1st col kegg id, 2nd col score, "?" denotes unknown score

    st = datetime.datetime.now()

    f = open("%s/rpair" % (rdir))                 # rpair
    
    f2 = open("%s/compound" % (rdir))   # compound
    status("Parsing compound...")
    compounds = klp.parse_compound(f2)
    
    status("Parsing rpair...")
    ag = klp.parseAtomGraphFromRpair(f)

    status("Parsing scores...")
    scores = {}   # reactionId -> score
    for s in scoresf:
        if s.startswith("#"):
            continue
        vals = s.strip().split("\t")
        scores[vals[0]] = vals[1]

    pairScores = {}
    pairScoreReactions = {}
    
    status("Calculating rpair scores...")
    for pair in ag.pairToReaction:
        maxScore = None
        maxReaction = None
        if ag.pairToReaction[pair]:
            for react in ag.pairToReaction[pair]:
                if maxReaction == None:
                    maxReaction = react  # to ensure we have at least one reaction in result list
                if react in scores and scores[react] != "?" and (maxScore == None or maxScore < float(scores[react])):
                    maxScore = float(scores[react])
                    maxReaction = react
        pairScores[pair] = maxScore
        pairScoreReactions[pair] = maxReaction

    #
    # Parse sources
    # 
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
    sources = parsedSources

    #
    # Find traces from each source to the target
    #
    
    status("Finding atom traces from %s to %s..." % (",".join(sources), tgt))
    resultPathways = findTraces(tdir, sources, tgt, ag, commonAtomReqUpper, commonAtomReqLower, pairScores, npaths, pairScoreReactions, tracedSourceAtoms)
    status("Pathfinding done. Result traces found: %d" % (len(resultPathways)))

    # Remove duplicate traces
    uniquePathways = []

    pathHashes = set()
    i = 1
    for path in resultPathways:
        h = hash(path)
        if h in pathHashes:
            continue
        pathHashes.add(h)
        uniquePathways.append(path)

    #
    # Sort the traces first by the number of mapped atoms, second by score,
    # third by the number of reactions. 
    # Note: combineTraces stored indecies to uniquePathways, so it has to be
    # sorted before the call to combineTraces
    uniquePathways.sort(cmpTrace)

    status("Combining %d unique traces into pathways..." % (len(uniquePathways)))
    combinedPathways, atomConfigs = combineTraces(uniquePathways, sources, tgt, ag)
    status("Found %d combined pathways with score > 0" % (len(combinedPathways)))

    #
    # Sort the combined pathways by the score
    #
    combinedPathways.sort(lambda x, y: int(y[1] - x[1]))

    #
    # Pick good reaction combinations for each combined pathway
    # 
    status("Finding good reaction combinations for each trace combination...")
    combos = chooseReactionCombinations(combinedPathways, atomConfigs, uniquePathways, scores)

    status("Writing traces to file...")

    writeTracesToFile(tdir, sources, tgt, uniquePathways, compounds, scores, ag, maxPath)

    status("Writing pathways to file...")

    writePathwaysToFile(tdir, sources, tgt, combos, compounds, scores, ag, 
                        atomConfigs, resultPathways)

    status("All done.")


if __name__ == '__main__':
    main()
