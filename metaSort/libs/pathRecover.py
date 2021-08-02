############################################################################
# Copyright (c) 2011-2014 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import indexGraph
import config
import lostLinks
from log import get_logger
# set log file
logger = get_logger()

def nodesOngraph(allNodes, soap):
    """
    caculate the number of node in the graph
    """
    status = 0 # use this to swith whether or not use pelink
    if not allNodes:
        return set(), set()
    
    for one in allNodes:
        if isinstance(one, basestring):
            nodes = set((int(i) for i in allNodes))
        elif isinstance(one, int):
            nodes = set(allNodes)
        else:
            logger.error('Unknown object type in the allNodes',
                         exit_with_code=1
                         )
        break
    
    pnodes = soap.dPlink.vertices() & nodes 
    dnodes = soap.dDlink.vertices() & nodes
    sharedNodes = pnodes & dnodes
    
    rate = len(pnodes | dnodes)*100.0 / len(nodes)
    rate = '%.1f' %rate
    logger.info('There are ' + str(len(pnodes)) + ' in the plink graph')
    logger.info('There are ' + str(len(dnodes)) + ' in the dlink graph')
    logger.info('There are ' + str(len(sharedNodes)) + ' in both of two graphs')
    logger.info('About ' + rate + '% of the total nodes are on the graph')
    return pnodes, dnodes

def kahnTopoSort(graph):
    """
    Implementation of Kahn's algorithm of sorting graph
    """
    import collections
    indegree = {u:0 for u in graph}
    for u in graph:
        for v in graph[u]:
            indegree[v] += 1
    
    Q = collections.deque() # collect nodes with zero in-degree
    for u in indegree:
        if indegree[u] == 0:
            Q.appendleft(u)
    
    L = [] # list for order of nodes
    while Q:
        u = Q.pop()
        L.append(u)
        for v in graph[u]:
            indegree[v] -= 1
            if indegree[v] == 0:
                Q.appendleft(v)
    
    if len(L) == len(graph):
        return L
    else: # # if there is a cycle,then return an empty list
        return L 

def revCom(string):
    """
    reverse complement the given string
    """
    comp = { 'A' : 'T', 'T' : 'A',
             'C' : 'G', 'G' : 'C',
             'a' : 'T', 't' : 'A',
             'c' : 'G', 'g' : 'C',
             'n' : 'N', 'N' : 'N'
            }
    result = ''
    for i in string[::-1]:
        result += comp[i]
    
    return result

def connectSubgraphs(graph, soap, mappedNodes):
    """
    Find the lost links that has been deleted by SOAPdenovo
    criteria: contig length > 100, there are no links except the lost dlink between two nodes
    the most important there are pelinks to support the lost dlink
    """
    logger.info('Finding new links')
    
    # get the subgraph numbers, noly the links that connect different subgraphs will used
    nodes = set([i for i in graph.vertices() if soap.length[i] >= 100])
    # mark the nodes
    scl = graph.findSubgraph(out='cl')
    group = 0
    for i in mappedNodes - graph.vertices():
        scl[i] = group
        group -= 1
        if soap.length[i] >= 100:
            nodes.add(i)
    
    graphSeqs = {}
    with open(soap.transedContig,'r') as af:
        for line in af:
            info = line.rstrip().split()
            if int(info[0]) in nodes:
                if len(info[-1]) < 200:
                    graphSeqs[int(info[0])] = info[-1]
                    rev = soap.revNode(int(info[0]))
                    graphSeqs[rev] = revCom(info[-1])
                else:
                    graphSeqs[int(info[0])] = info[-1][:100] + info[-1][-100:]
                    revSeq = revCom(info[-1])
                    rev = soap.revNode(int(info[0]))
                    graphSeqs[rev] = revSeq[:100] + revSeq[-100:]
    
    # scan kmers to math the nodes
    links = []
    krange = (i for i in xrange(20,100) if i != soap.kmersize)
    for size in krange:
        kmerleft, kmerright = {}, {}
        for node in nodes:
            if node not in graphSeqs:
                continue
            if len(graphSeqs[node]) > size:
                if graphSeqs[node][:size] not in kmerleft:
                    kmerleft[graphSeqs[node][:size]] = []
                kmerleft[graphSeqs[node][:size]].append(node)
                
                if graphSeqs[node][-size:] not in kmerright:
                    kmerright[graphSeqs[node][-size:]] = []
                kmerright[graphSeqs[node][-size:]].append(node)
        
        # key is kmer string
        for key in kmerright:
            if key in kmerleft:
                for up in kmerright[key]:
                    for down in kmerleft[key]:
                        links.append((up, down, 0-size))
                        nodes.discard(down)
                    nodes.discard(up)
    
    # first filter the newlink by contig length, the error rate of two very short contig overlaps is very high
    links = [link for link in links if soap.length[link[0]]>100 and soap.length[link[1]]>100]
    
    # second the newlink should connect two distinct subgraphs
    newlinks = [link for link in links if scl[link[0]] != scl[link[1]]]
    
    # third: there are no links except the lost dlink between two nodes
    passedLinks = []
    
    # fourth the lost link must have pelink suport
    for link in newlinks:
        if link[1] in graph.plinkNodes(link[0], soap.dPlink):
            passedLinks.append(link)
    
    logger.info('{0:,} new links are created.'.format(len(passedLinks)))
    
    return passedLinks

def algin(left, right, kmersize, idenCut):
    """
    Algin two seqs 
    """
    initLen = kmersize - 1
    for step in xrange(initLen, 5, -1):
        match, mismatch = 0, 0
        subseq1 = left[-step:]
        subseq2 = right[:step]
        for pos, value in enumerate(subseq1):
            if value == subseq2[pos]:
                match += 1
            else:
                mismatch += 1
        if match*1.0/(match+mismatch)>=idenCut:
            return -step
    else:
        return 0

def linkNodes(graph, gNodes, maxStep, node):
    """
    bfs searching on the graph using the node
    """
    pairs = set()
    for path in graph.bfs(node, step = 2*maxStep):
        if set(path[1:]) & gNodes:
            tmp = [i for i in path if i in gNodes]
            for i in xrange(len(tmp)-1):
                pairs.add(tuple([tmp[i], tmp[i+1]]))

    return pairs

def genomeGraph(soap, mappedNodes, candidate, step=3):
    """
    Search the dlink graph and plink graph based on the mappedNodes and candidate nodes
    and filter the false paths to extract the target genome paths
    """
    
    # determine the searching stratagy
    try:
        maxStep = step
    except:
        logger.error('step shold be a int',
                     exit_with_code = 1
                     )
    
    logger.debug('\nSearching max step size: {}'.format(maxStep))
    
    # generate the seed nodes on the graph
    pnodes, dnodes = nodesOngraph(mappedNodes|candidate, soap)
    pnodes = pnodes - dnodes
    
    # search the path 
    dPaths, pathNodes = [], set()
    for node in dnodes:
        for path in soap.dDlink.bfs(node, step = maxStep):
            if len(path) > 1:
                dPaths.append(path)
                pathNodes.update(path)
    
    # build a graph
    tmpg = indexGraph.Graph(soap.index)
    for path in dPaths:
        if len(path)> 1:
            for p in xrange(len(path)-1):
                o, t = path[p], path[p+1]
                ro, rt = soap.revNode(path[p+1]), soap.revNode(path[p])
                tmpg.addEdge(o, t, soap.dDlink[o][t])
                tmpg.addEdge(ro, rt, soap.dDlink[o][t])
    
    ########################### graph trimming ##########################
    logger.info('Trimming graph ...')
    
    # trim the nodes that no contribution to link two known nodes
    gNodes = tmpg.vertices() & (mappedNodes | candidate)
    
    pairs = set()
    for node in gNodes:
        for path in tmpg.bfs(node, step = 2*maxStep):
            if set(path[1:]) & gNodes:
                tmp = [i for i in path if i in gNodes]
                for i in xrange(len(tmp)-1):
                    pairs.add(tuple([tmp[i], tmp[i+1]]))
    
    accept = set()
    for s, t in pairs:
        for path in tmpg.allPaths(s, t, step = 2*maxStep):
            accept.update(path)
    accept.update([soap.revNode(i) for i in accept])
    
    _ = tmpg.deleteNodes(tmpg.vertices() - accept)
    
    # every subgraph should have at least 2 candiates to support the feasibility
    groups = tmpg.findSubgraph('gp')
    filter1 = set()
    for num, value in groups.iteritems():
        if not set(groups[num]) & mappedNodes:
            #if len(set(groups[num]) & candidate) <= 4:
            #    filter1.update(value)
            filter1.update(value)
    _ = tmpg.deleteNodes(filter1)
    
    # link the mappednodes using a increasing stepsize
    addtPath = []
    leftMappedNodes = mappedNodes - tmpg.vertices()
    for k in xrange(maxStep+1, 6):
        for ele in leftMappedNodes:
            for path in soap.dDlink.bfs(ele, step = k):
                inter = set(path) & mappedNodes
                pos = sorted([path.index(i) for i in inter])[-1]
                addtPath.append(path[:pos+1])
    
    for path in addtPath:
        if len(path)> 1:
            for p in xrange(len(path)-1):
                o, t = path[p], path[p+1]
                ro, rt = soap.revNode(path[p+1]), soap.revNode(path[p])
                tmpg.addEdge(o, t, soap.dDlink[o][t])
                tmpg.addEdge(ro, rt, soap.dDlink[o][t])
    
    # create new dlink in case of pelink support number is 1
    newlinks = connectSubgraphs(tmpg, soap, mappedNodes)
    
    # add the robust newlinks to the graph
    for link in newlinks:
        up, down, overlap = link
        tmpg.addEdge(up, down, overlap)
        tmpg.addEdge(soap.revNode(down), soap.revNode(up), overlap)
        
    # searching on the pelink by using pnodes as marker (seeds)
    gNodes = tmpg.vertices() | mappedNodes
    peStep, peOnly = 1, []
    for node in pnodes:
        for path in soap.dPlink.bfs(node, step = peStep):
            if len(set(path)) > 1:
                # in case of self mapped pelinks
                if path[-1] in gNodes:
                    # have support pelinks, check if there are dlinks to link them
                    mean, sd = soap.dPlink[path[0]][path[1]][:2]
                    gapRange = [mean-3*sd, mean+3*sd]
                    paths = tmpg.allPaths(path[0], path[1])
                    # check if the distance betweeen two nodes mathes
                    if not paths:
                        peOnly.append(path)
        
    graphSeqs = soap.fetchSeqs(set([i for k in peOnly for i in k]))
    for o, t in peOnly:
        matches = algin(graphSeqs[o], graphSeqs[t], soap.kmersize, 0.9)
        if matches:
            tmpg.addEdge(o, t, -matches)
            tmpg.addEdge(soap.revNode(t), soap.revNode(o), -matches)
        else: 
            tmpg.addEdge(o, t, soap.dPlink[o][t])
            tmpg.addEdge(soap.revNode(t), soap.revNode(o), soap.dPlink[o][t])
        
    # doing statistic analysis of recovered nodes
    gNodes = tmpg.vertices()
        
    logger.info()
    logger.info('Recovering and searching statistics:')
    logger.info('#'*60)
    logger.info('Number of seeds:                      {0:,}'.format(len(mappedNodes)/2))
    logger.info('Graph size:                           {0:,}'.format(len(gNodes|mappedNodes)/2))
    logger.info('Number of recovered nodes:            {0:,}'.format(len(gNodes-mappedNodes)/2))
    logger.info('  By contig feature:                  {0:,}'.format(len(gNodes & candidate)/2))
    logger.info('  By graph searching:                 {0:,}'.format(len(gNodes - (candidate | mappedNodes))/2))
    logger.info('#'*60)
    
    return tmpg


