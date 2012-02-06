# -*- coding: iso-8859-15 -*-
"""
ReTrace - find branching pathways in metabolic networks
Copyright (C) 2009 Esa Pitkänen
See file COPYING for license notice.
"""

import sys, os, shutil, subprocess, random

# Set DRAW_MOL_STRUCTS to 1 if you want to have ReTrace draw
# pathway diagrams with molecule structures. You need to provide 
# the KEGG gif directory in MOL_STRUCT_DIR.

DRAW_MOL_STRUCTS = 0
MOL_STRUCT_DIR = "kegg/gif"

# Name of the file written in result directory which contains
# a one-line summary from each query performed.

SUMMARY_FNAME = "summary.txt"

# How many best-scoring reactions to show at maximum per rpair?

REACTIONS_IN_RPAIR = 1

# Reactions with a score under SCORE_THRESHOLD will be reported
# separately

SCORE_THRESHOLD = 50

# Metabolites which are drawn separately in visualization
# Note that these settings have no effect on the analysis itself.

SIDE_METABOLITES_IN_VISUALISATION = set(["C00001",   # water
                                         "C00007",   # O2
                                         "C00011",   # CO2
                                         "C00014",   # NH3
                                         "C00027",   # H2O2
                                         "C00009",   # Orthophosphate
                                         "C00003",   # NAD+
                                         "C00004",   # NADH
                                         "C00005",   # NADPH
                                         "C00006",   # NADP+
                                         "C00080",   # H+ 
                                         "C00013",   # orthophosphate
                                         "C00002",   # ATP
                                         "C00008"])  # ADP

# Individual pathway results will be written to MAPDIR directory

MAPDIR = "paths"

# URL of KEGG DBGET for crosslinking reaction and metabolite data

KEGG_DBGET = "http://www.genome.jp/dbget-bin/www_bget?"

# Graphviz (pathway visualization) options

# Graphviz layout program to use
DOT = "dot"
# Image file type of pathway diagrams
IMAGE_TYPE = "png"

# URL of the web program (ReMatch) to visualize molecules with
# See www.cs.helsinki.fi/group/sysfys for details on ReMatch

MOLVIEWER = "http://sysdb.cs.helsinki.fi/ReMatch/MoleculeImage"
MOLVIEWER_PARAMS = "&mapping=RPAIR"

# Where can we find the javascript code for sortable tables
# and what is the name of the code file

HTML_DIR = "html"
JS_SORTABLE_TABLE = "webtoolkit.sortabletable.js"

# CSS code for sortable tables

HTML_HEAD_SORTABLE_TABLE = """
  <script type="text/javascript" src="%s"></script>
	<style>
		table {
			text-align: left;
			font-size: 12px;
			font-family: verdana;
			background: #c0c0c0;
		}
 
		table thead  {
			cursor: pointer;
		}
 
		table thead tr,
		table tfoot tr {
			background: #c0c0c0;
		}
 
		table tbody tr {
			background: #f0f0f0;
		}
 
		td, th {
			border: 1px solid white;
		}
	</style>
""" % (JS_SORTABLE_TABLE)

HTML_BODY_SORTABLE_TABLE = """
  <script type="text/javascript">
  var t = new SortableTable(document.getElementById('%s'), 100);
  </script>
"""

def dotAvailable():
    """Check whether we can use Graphviz dot tool."""
    ret = subprocess.call("%s -V" % (DOT), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    return ret == 0
              
def getCompositeMapStr(m):
    """Return a string representation of the composite map corresponding to the given atom mapping."""
    compositeMapStr = ""
    pairs = {}
    for src in m:
        srcmol, sa = src.split("-")
        for tgt in m[src]:
            tgtmol, ta = tgt.split("-")
            pair = "%s->%s" % (srcmol, tgtmol)
            if pair not in pairs:
                pairs[pair] = []
            pairs[pair].append("%s->%s" % (sa, ta))
    for pair in pairs:
        pairs[pair].sort()
        compositeMapStr += "%s:[%s], " % (pair, ",".join(pairs[pair]))
    return compositeMapStr.rstrip(", ")

def prepareHtmlOutput(tdir):
    """Copy code for sortable html tables to output directory."""
    mapdir = "%s/%s" % (tdir, MAPDIR)
    if not os.path.exists(mapdir):
        os.mkdir(mapdir)

    # copy javascript for sortable tables to tdir
    shutil.copy("%s/%s" % (HTML_DIR, JS_SORTABLE_TABLE), "%s/%s" % (tdir, JS_SORTABLE_TABLE))
    shutil.copy("%s/%s" % (HTML_DIR, JS_SORTABLE_TABLE), "%s/%s/%s" % (tdir, MAPDIR, JS_SORTABLE_TABLE))

def getMolName(mol, compounds):
    """Return a name for the molecule, if such exists."""
    if mol in compounds:
        name = compounds[mol].names[0].strip()
    else:
        name = mol
    return name

def writePathwaysToHtml(tdir, sources, targets, pathways, compounds, ag, net, scores,
                        searchParameters, scoreOrigin, disablePathOutput):
    """Write html and text files describing ReTrace results."""
    drawPathwayDiagram = dotAvailable()
    
    srcnames = {}
    for src in sources:
        srcnames[src] = getMolName(src, compounds)
    
    tgtnames = {}
    for tgt in targets:
        tgtnames[tgt] = getMolName(tgt, compounds)
            
    of = open("%s/pathways-%s-to-%s.html" % (tdir, "-".join(sources), "-".join(targets)), "w")
    of3 = open("%s/pathways-%s-to-%s.txt" % (tdir, "-".join(sources), "-".join(targets)), "w")

    of3.write("# Pathways from %s (%s) to %s (%s)\n" % (",".join(srcnames.values()), ",".join(srcnames.keys()), ",".join(tgtnames.values()), ",".join(tgtnames.keys())))
    items = searchParameters.keys()
    items.sort()

    if not os.path.exists("%s/%s" % (tdir, SUMMARY_FNAME)):
        summaryf = open("%s/%s" % (tdir, SUMMARY_FNAME), "w")
        summaryf.write("#")
        for p in items:
            summaryf.write("%s " % (p.replace(" ", "")))
        summaryf.write("\n")
        summaryf.close()

    print "Appending results to summary file %s..." % (SUMMARY_FNAME)

    summaryf = open("%s/%s" % (tdir, SUMMARY_FNAME), "a")
    for p in items:
        summaryf.write("%s\t" % (str(searchParameters[p])))
    summaryf.write("\n")
    summaryf.close()

    print "Writing query result files..."

    for p in items:
        of3.write("# %s:\t%s\n" % (p, searchParameters[p]))

    of.write("""
<html>
<head>
  <title>Pathways from %s to %s</title>
  %s
</head>
<body>\n""" % (",".join(sources), ",".join(targets), HTML_HEAD_SORTABLE_TABLE))

    of.write("<p><b>Pathways from %s to %s</b></p>" % ("-".join(sources), "-".join(targets)))

    of.write("<p><b>Sources:</b> ")
    for src in sources:
        of.write("""%s (<a href=\"http://www.genome.jp/dbget-bin/www_bget?compound+%s\">%s</a>)
<img src=\"%s?mol=%s&show_carbon_indices=true&highlight=%s%s\">\n
""" % (srcnames[src], src, src, MOLVIEWER, src, "0", MOLVIEWER_PARAMS))
    of.write("""<p><b>Target:</b>""")
    for tgt in targets:
        of.write("""%s (<a href=\"http://www.genome.jp/dbget-bin/www_bget?compound+%s\">%s</a>)
<img src=\"%s?mol=%s&show_carbon_indices=true&highlight=%s%s\"></p>\n
""" % (tgtnames[tgt], tgt, tgt, MOLVIEWER, tgt, "0", MOLVIEWER_PARAMS))

    of.write("""<table id=\"pathwaytable\" border=\"1\">
<thead>
<tr><th>&nbsp;</th><th>Composite mapping</th><th>Z</th><th>Average score</th><th>Rpairs</th><th>Reactions</th><th>Zero scores</th><th>Scores under threshold</th></tr>
</thead>
<tbody>\n""")
    of3.write("#PathIndex ZoScore AvgReactionScore NumRpairs NumReactions ZeroScores ScoresUnderThr CompositeMap Rpairs Reactions\n")

    srclist = "-".join(sources)
    tgtlist = "-".join(targets)

    i = 1
    for path in pathways:
        fn = "path-%s-to-%s-%d.html" % (srclist, tgtlist, i)

        avgScore = 0
        ms = 0
        ut = 0
        reCount = 0

        #rpairStr = ",".join(path.rpairs)
        reStr = ""
        rpairStr = ""

        for rpair, dir in path.rpairs:
            rpairStr += "%s," % (rpair)
            missing = True
            under = True
            for re in ag.pairToReaction[rpair]:
                reCount += 1
                reStr += "%s," % (re)

                if re not in scores:
                    continue

                sc = scores[re]
                if sc == "0" or sc == "-" or sc == "?":
                    pass
                elif int(sc) < SCORE_THRESHOLD:
                    missing = False
                    avgScore += int(sc)
                else:
                    missing = under = False
                    avgScore += int(sc)

            if missing:
                ms += 1
            if under:
                ut += 1

        if reCount > 0:
            avgScore = 1.0 * avgScore / reCount
        else:
            avgScore = 0
        reStr = reStr.rstrip(",")
        rpairStr = rpairStr.rstrip(",")
        
        compositeMapStr = getCompositeMapStr(path.compositeMap)

        of.write("""<tr>
  <td><a href=\"%s/%s\">Path %d</a></td>
  <td>%s</td>
  <td>%.2f</td>
  <td>%s</td>
  <td>%d</td>
  <td>%d</td>
  <td>%s</td>
  <td>%s</td>
</tr>\n""" % (MAPDIR, fn, i, 
              compositeMapStr,
              path.Z,
              avgScore,
              len(path.rpairs),
              reCount,
              ms, 
              ut))

        of3.write("%d\t%.2f\t%2.f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n" % (i, path.Z, avgScore, 
                                                                    len(path.rpairs), 
                                                                    reCount, ms, ut, 
                                                                    compositeMapStr,
                                                                    rpairStr, reStr))

        i += 1

    of.write("</tbody><tfoot></tfoot></table>\n")
    of.write(HTML_BODY_SORTABLE_TABLE % ("pathwaytable"))
    of.close()

    of3.close()

    if disablePathOutput == False:

        print "Writing results for individual paths..."

        i = 1
        for path in pathways:
            fn = "path-%s-to-%s-%d.html" % (srclist, tgtlist, i)
            writePathToHtml(tdir, i, sources, targets, fn, path, ag, net, 
                            compounds, scores, drawPathwayDiagram, scoreOrigin)
            i += 1



def bfsOrder(reactions, sources, net):
    """Obsolete for now. Used to determine edge directions in pathway diagrams."""
    order = {}
    queue = []
    
    for r in reactions:
        re = net.reactions[r]
        for src in sources:
            if src in re.substrates:
                order[r] = (0, True)  # src -> reaction
                queue.append(r)
            elif src in re.products:
                order[r] = (0, False) # src <- reaction, reverse dir
                queue.append(r)

    next = 1
    while len(queue) > 0:
        r = queue.pop(0)
        re = net.reactions[r]
        #print re

        for sub in re.substrates:
            mol = net.molecules[sub]
            for r2 in mol.consumers:
                if r2 in reactions and r2 not in order:
                    order[r2] = (next, True)
                    next += 1
                    queue.append(r2)
            for r2 in mol.producers:
                if r2 in reactions and r2 not in order:
                    order[r2] = (next, False)
                    next += 1
                    queue.append(r2)

        for pro in re.products:
            mol = net.molecules[pro]
            for r2 in mol.consumers:
                if r2 in reactions and r2 not in order:
                    order[r2] = (next, True)
                    next += 1
                    queue.append(r2)
            for r2 in mol.producers:
                if r2 in reactions and r2 not in order:
                    order[r2] = (next, False)
                    next += 1
                    queue.append(r2)

    #print order
    return order

def getMolDotId(mol):
    if mol in SIDE_METABOLITES_IN_VISUALISATION:
        label = "%s_%d" % (mol, random.randint(10000,99999))
    else:
        label = mol
    return label

def getMolDotStr(mol, molid, compounds, drawMolStructs):
    url = "%scpd:%s" % (KEGG_DBGET, mol)
    if drawMolStructs:
        name = ""
        tooltip = compounds[mol].names[0].strip()
        images = ",tooltip=\"%s\",shape=\"custom\",shapefile=\"%s/%s.gif\"" % (tooltip, MOL_STRUCT_DIR, mol)
    else:
        name = compounds[mol].names[0].strip()
        images = ""
    return "\t%s [label=\"%s\",URL=\"%s\"%s];\n" % (molid, name, url, images)


def getReactionStr(r, compounds):
    rstr = ""
    for sub in r.substrates:
        name = compounds[sub].names[0].strip()
        rstr += "%s + " % (name)
    rstr = "%s <=> " % (rstr.rstrip(" + "))
    for pro in r.products:
        name = compounds[pro].names[0].strip()
        rstr += "%s + " % (name)
    rstr = rstr.rstrip(" + ")
    return rstr

def convertPathToImage(tdir, reactions, imagefn, ag, sources, targets, net, 
                       compounds, scores, mapfn, drawMolStructs):

    #ordering = bfsOrder(reactions, sources, net)
    #for r in ordering:
    #    print r, ordering[r]

    dotfn = "temp%d.dot" % (random.randint(100000,999999))
 
    o = open("%s/%s" % (tdir, dotfn), "w")
    o.write("digraph path {\n" % ())
    for src in sources:
        o.write("\t%s [style=\"filled\", fillcolor=\"green\"];\n" % (src))
    for tgt in targets:
        o.write("\t%s [style=\"filled\", fillcolor=\"yellow\"];\n" % (tgt))

    for re, dir in reactions:

        r = net.reactions[re]
        score = "?"
        if re in scores:
            score = scores[re]
        color = scorecolor(score)

        url = "%srn+%s" % (KEGG_DBGET, re)

        o.write("\t%s [tooltip=\"%s\",shape=\"box\",style=\"filled\",fillcolor=\"%s\",URL=\"%s\"];\n" % (re, getReactionStr(r, compounds), color, url))

        #bfsrank, dir = ordering[re]

        if dir == True:
            lhs, rhs = r.substrates, r.products
        else:
            lhs, rhs = r.products, r.substrates

        for sub in lhs:
            id = getMolDotId(sub)
            molstr = getMolDotStr(sub, id, compounds, drawMolStructs)
            o.write(molstr)
            o.write("\t%s -> %s;\n" % (id, re))

        for pro in rhs:
            id = getMolDotId(pro)
            molstr = getMolDotStr(pro, id, compounds, drawMolStructs)
            o.write(molstr)
            o.write("\t%s -> %s;\n" % (re, id))

    o.write("}\n")
    o.close()

    subprocess.call("%s -Tcmapx -o %s -T%s %s/%s -o %s" % (DOT, mapfn, IMAGE_TYPE, tdir, dotfn, imagefn), shell = True)

    os.remove("%s/%s" % (tdir, dotfn))

def scorecolor(score):
    if score == "0" or score == "?" or score == "-":
        color = "indianred"
    elif int(score) < SCORE_THRESHOLD:
        color = "lightblue"
    else:
        color = "palegreen"
    return color

def cmpreaction(s1, s2):
    score1 = 0.0
    if s1 != "-" and s1 != "?":
        score1 = float(s1)
    score2 = 0.0
    if s2 != "-" and s2 != "?":
        score2 = float(s2)
    return cmp(score1, score2)

def writePathToHtml(tdir, ix, sources, targets, fn, path, ag, net, compounds, scores, drawPathwayDiagram, scoreOrigin):
    of = open("%s/%s/%s" % (tdir, MAPDIR, fn), "w")
    
    imagebase = "path-%s-to-%s-%d" % ("-".join(sources), "-".join(targets), ix)

    # image with simple nodes
    imagefn = "%s/%s/%s.%s" % (tdir, MAPDIR, imagebase, IMAGE_TYPE)

    # image with molecule structure images as nodes
    structimagefn = "%s/%s/%s-struct.%s" % (tdir, MAPDIR, imagebase, IMAGE_TYPE)

    mapfn = "%s/%s/%s.map" % (tdir, MAPDIR, imagebase)
    structmapfn = "%s/%s/%s-struct.map" % (tdir, MAPDIR, imagebase)

    maxreaction = REACTIONS_IN_RPAIR

    reactions = set()
    for rpair, rpairdir in path.rpairs:
        # Determine edge direction for each reaction based on 
        # rpair's direction
        if rpairdir:
            sub = ag.pairToSub[rpair]
        else:
            sub = ag.pairToPro[rpair]

        items = list(ag.pairToReaction[rpair])
        items.sort(lambda x, y: -cmpreaction(scores[x], scores[y]))
        for re in items[0:min(maxreaction, len(items))]:
            r = net.reactions[re]
            if sub in r.substrates:
                ri = (re, True)
            else:
                ri = (re, False)
            reactions.add(ri)

    if drawPathwayDiagram:
        convertPathToImage(tdir, reactions, imagefn, ag, sources, targets, 
                           net, compounds, scores, mapfn, False)
        if DRAW_MOL_STRUCTS:
            convertPathToImage(tdir, reactions, structimagefn, ag, sources, targets, 
                               net, compounds, scores, structmapfn, True)

    # Determine which atoms have been mapped by this pathway in sources and targets
    srchighlight = {}
    tgthighlight = {}
    for srcmol in path.compositeMap:
        src, sa = srcmol.split("-")
        if src not in srchighlight:
            srchighlight[src] = set()
        srchighlight[src].add(sa)
        for tgtmol in path.compositeMap[srcmol]:
            tgt, ta = tgtmol.split("-")
            if tgt not in tgthighlight:
                tgthighlight[tgt] = set()
            tgthighlight[tgt].add(ta)

    of.write("""
<html>
<head>
  <title>Pathway %d from %s to %s</title>
  %s
</head>
<body>\n""" % (ix, ",".join(sources), ",".join(targets), HTML_HEAD_SORTABLE_TABLE))

    of.write("<h2>Path %d from %s to %s</h2>\n" % (ix, ",".join(sources), ",".join(targets)))

    of.write("<p>Back to <a href=\"../pathways-%s-to-%s.html\">all pathways</a></p>\n" 
             % ("-".join(sources), "-".join(targets)))

    of.write("<p><b>Sources:</b> ")
    for src in sources:
        srchlstr = 0
        if src in srchighlight:
            srchlstr = ",".join(list(srchighlight[src]))
        of.write("""%s (<a href=\"http://www.genome.jp/dbget-bin/www_bget?compound+%s\">%s</a>)
<img src=\"%s?mol=%s&show_carbon_indices=true&highlight=%s%s\">\n
""" % (getMolName(src, compounds), src, src, MOLVIEWER, src, 
       srchlstr, MOLVIEWER_PARAMS))
    of.write("""<p><b>Target:</b>""")
    for tgt in targets:
        tgthlstr = "0"
        if tgt in tgthighlight:
            tgthlstr = ",".join(list(tgthighlight[tgt]))

        of.write("""%s (<a href=\"http://www.genome.jp/dbget-bin/www_bget?compound+%s\">%s</a>)
<img src=\"%s?mol=%s&show_carbon_indices=true&highlight=%s%s\"></p>\n
""" % (getMolName(tgt, compounds), tgt, tgt, MOLVIEWER, tgt, 
       tgthlstr, MOLVIEWER_PARAMS))

    of.write("""<p>
Z: %s<br/>\n
Composite map: %s
</p>\n""" % (path.Z, getCompositeMapStr(path.compositeMap)))

    of.write("<p>Showing %d best reactions for each rpair.</p>\n" % (maxreaction))

    of.write("""<table id=\"pathwaytable\" border=\"1\">
<thead>
<tr><th>RPAIR</th><th>Reaction</th><th>Score</th><th>Seq1</th><th>Seq2</th><th>Evalue</th><th>ECs</th><th>Equation</th></tr>
</thead>
<tbody>\n""")
    
    for rpair, rpairdir in path.rpairs:
        items = list(ag.pairToReaction[rpair])
        items.sort(lambda x, y: -cmpreaction(scores[x], scores[y]))
        
        maxi = min(maxreaction, len(items))
        for re in items[0:maxi]:
            r = net.reactions[re]
            score = "0"
            if re in scores:
                score = scores[re]
            if score == "?":
                score = "0"   # otherwise unable to sort html table by score column correctly

            rstr = getReactionStr(r, compounds)


            seq1 = seq2 = evalue = ecstr = ""
            if re in scoreOrigin:
                so = scoreOrigin[re]
                seq1 = so.seq1
                seq2 = "<a href=\"http://www.uniprot.org/uniprot/%s\">%s</a>" % (so.seq2.split("|")[0], so.seq2)
                
                evalue = so.evalue
                ecs = so.ecs.split(",")
                ecstr = ""
                for ec in ecs:
                    ecstr += "<a href=\"http://www.genome.jp/dbget-bin/www_bget?%s\">%s</a>," % (ec, ec)
                ecstr = ecstr.rstrip(",")

            of.write("<tr><td>%s</td>" % (rpair))
            of.write("<td><a href=\"http://www.genome.jp/dbget-bin/www_bget?rn+%s\">%s</a></td><td><font style=\"BACKGROUND-COLOR: %s\">%s</font></td>" % (re, re, scorecolor(score), score))
            of.write("<td>%s</td><td>%s</td><td>%s</td><td>%s</td>\n" % (seq1, seq2, evalue, ecstr))
            of.write("<td>%s</td></tr>\n" % (rstr))

    of.write("</tbody><tfoot></tfoot></table>\n")
    of.write(HTML_BODY_SORTABLE_TABLE % ("pathwaytable"))

    if drawPathwayDiagram:

        # Link to the pathway diagram with molecule structures
        of.write("<p><b>Pathway diagram</b>")
        if DRAW_MOL_STRUCTS:
            of.write("(see <a href=\"%s-struct.%s\">pathway with molecule structures</a> instead)</p>\n" % (imagebase, IMAGE_TYPE))

        # Regular pathway diagram embedded in html
        of.write("<p><img border=0 src=\"%s.%s\" usemap=\"#path\"/></p>\n" % (imagebase, IMAGE_TYPE))
        # Write imagemap to html
        mapf = open(mapfn)
        for s in mapf:
            of.write(s)
        mapf.close()
 
    of.write("</body></html>\n")

    of.close()
