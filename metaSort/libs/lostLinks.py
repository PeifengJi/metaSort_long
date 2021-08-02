############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import itertools
import collections
import copy
import config
import indexGraph
from log import get_logger
# set log file
logger = get_logger()

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

def lenBfsNodes(soap, node, lencut, lastNode=0, reverse=0):
    """
    bfs algorithm to traversal the graph under the limited distance
    """
    queue = collections.deque([[node]])
    searched, visitedNodes = set(), set()
    # there is no visitedNodes, no matter the redundant, beacuse if we set vistiedNodes
    # there will be a huge rate of unbalance link in two strand (caused by bubbule)
    visitedNodes = set()
    while queue:
        if len(queue)>1000:
            print 's'
            return set([node]) if not reverse else set([soap.revNode(node)])
        tmp = queue.pop()
        last = tmp[-1]
        if last in visitedNodes:
            searched.update(tmp)
            continue
        visitedNodes.add(last)
        if last in soap.dDlink:
            for second in soap.dDlink[last]:
                longerPath = tmp + [second]
                acclen = sum([soap.length[i]-soap.kmersize for i in longerPath[1:]])
                if acclen >= lencut:
                    searched.update(longerPath) if not lastNode else searched.update(longerPath[-lastNode:])
                    # too much node 
                    if lastNode and len(searched) > 200:
                        print 'passed'
                        return set([node]) if not reverse else set([soap.revNode(node)])
                    elif not lastNode and len(searched) > 500:
                        print 'passed'
                        return set([node]) if not reverse else set([soap.revNode(node)])
                else:
                    queue.appendleft(longerPath)
        else:
            searched.update(tmp) if not lastNode else searched.update(tmp[-lastNode:])
            # too much node 
            if lastNode and len(searched) > 200:
                print 'passed'
                return set([node]) if not reverse else set([soap.revNode(node)])
            elif not lastNode and len(searched) > 500:
                print 'passed'
                return set([node]) if not reverse else set([soap.revNode(node)])
    
    if not reverse:
        return searched
    else:
        return set([soap.revNode(node) for node in searched])

def lenShortest(soap, dlink, start, end, lencut):
    """
    Use bfs algorithm to get the true dirrection and order of two plink nodes
    """
    path, searched, visited = [], [], set()
    queue = collections.deque([[start]])
    poll = set([end, soap.revNode(end)])
    while queue:
        tmp = queue.pop()
        last = tmp[-1]
        if last in visited:
            continue
        visited.add(last)
        if last in dlink:
            for second in dlink[last]:
                if second in poll:
                    return 1
                else:
                    longerPath = tmp + [second]
                    acclen = sum([soap.length[i] for i in longerPath[1:]]) - soap.kmersize*len(longerPath)
                    if acclen <= lencut:
                        queue.appendleft(longerPath)
    return 0

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

def stretchingMatch(soap, dlink, graphSeqs, rkey, lkey, mindis, maxdis):
    """
    Find the lost dlinks by stretching strategy
    """
    rPath = {(rkey,) : 0} # gapsize as value
    lvkey = soap.revNode(lkey)
    lPath = {(lvkey,) : 0}
    extend = {(rkey,) : 1, (lvkey,) : 1}
    steps, maxMatches = [-50, -30, -20, -15, -10], []
    while 1:
        rRout, lRout = {}, {}
        for path in rPath:
            if extend[path]:
                if path[-1] in dlink:
                    for key in dlink[path[-1]]:
                        newid = tuple(list(path) + [key])
                        extendGapLen = rPath[path] + soap.length[key] - soap.kmersize
                        rRout[newid] = extendGapLen
                        extend[newid] = 1 if extendGapLen < maxdis else 0
                else:
                    rRout[path] = rPath[path]
                    extend[path] = 0
            else:
                rRout[path] = rPath[path]
        
        for path in lPath:
            if extend[path]:
                if path[-1] in dlink:
                    for key in dlink[path[-1]]:
                        newid = tuple(list(path) + [key])
                        extendGapLen = lPath[path] + soap.length[key] - soap.kmersize
                        lRout[newid] = extendGapLen
                        extend[newid] = 1 if extendGapLen < maxdis else 0
                else:
                    lRout[path] = lPath[path]
                    extend[path] = 0
            else:
                lRout[path] = lPath[path]
        
        for path in extend.keys():
            if path in rRout or path in lRout:
                continue
            else:
                del extend[path]
        
        if len(extend) > 20000: # if too many paths, exit
            return set()
        
        # compare the rout
        matches = []
        for comb in itertools.product(rRout.iterkeys(), lRout.iterkeys()):
            # gapMin = rRout[comb[0]]+lRout[comb[1]]-soap.kmersize
            gapMax = rRout[comb[0]] + lRout[comb[1]]
            if mindis <= gapMax <= maxdis:
                # filter other obstacle paths to reduce the searching time and memory
                lv = soap.revNode(comb[1][-1])
                if (comb[0][-1], lv) not in matches:
                    score = algin(graphSeqs[comb[0][-1]], graphSeqs[lv], soap.kmersize, 0.9)
                    if score:
                        matches.append((comb[0][-1], lv, score))
        
        if set(extend.values()) == set([0]):
            for comb in itertools.product(rRout.iterkeys(), lRout.iterkeys()):
                lv = soap.revNode(comb[1][-1])
                if (comb[0][-1], lv) not in matches:
                    score = algin(graphSeqs[comb[0][-1]], graphSeqs[lv], soap.kmersize, 0.9)
                    if score:
                        matches.append((comb[0][-1], lv, score))
        
        if matches:
            matches.sort(key=lambda x:x[-1])
            if not maxMatches or maxMatches[-1][-1] > matches[0][-1]:
                maxMatches = [matches[0]]
                for pos in xrange(1, len(matches)):
                    if matches[pos][-1] == maxMatches[-1][-1]:
                        maxMatches.append(matches[pos])
                    else:
                        break
            else:
                return set(maxMatches)
            
            # filter the overlap to filtered noise pathes to reduce searching time and save memory
            if steps:
                conserved = set([ele for ele in matches if ele[-1] <= steps[-1]])
                if conserved:
                    steps.pop()
                    rnodes, lnodes = set(), set()
                    for ele in conserved:
                        rnodes.add(ele[0])
                        lnodes.add(soap.revNode(ele[1]))
                    
                    for ele in (key for key in rRout.keys() if key[-1] in rnodes):
                        del rRout[ele]
                    
                    for ele in (key for key in lRout.keys() if key[-1] in lnodes):
                        del lRout[ele]
        # also converged
        elif not matches and maxMatches:
            return set(maxMatches)
        
        # if no nodes satisfy the situation, go to next cycle
        if set(extend.values()) == set([0]):
            return set() if not maxMatches else set(maxMatches)
        else:
            rPath, lPath = copy.deepcopy(rRout), copy.deepcopy(lRout)

def findLostDlinks(soap):
    """
    Find the lost dlinks that deleted by SOAPdenovo when construct de brujin graph
    Use pelink as a ref to find lost dlinks
    """
    graphSeqs = {}
    with open(soap.transedContig,'r') as af:
        for line in af:
            info = line.rstrip().split()
            rev = soap.revNode(int(info[0]))
            if len(info[-1]) <= 2*soap.kmersize:
                graphSeqs[int(info[0])] = info[-1]
                graphSeqs[rev] = revCom(info[-1])
            else:
                graphSeqs[int(info[0])] = info[-1][:soap.kmersize] + info[-1][-soap.kmersize:]
                revSeq = revCom(info[-1])
                graphSeqs[rev] = revSeq[:soap.kmersize] + '_' + revSeq[-soap.kmersize:]
    
    if not soap.dDlink:
        # load ARC file
        dlink = indexGraph.indexLinks(soap.arcfile, soap.index, soap.kmersize, linkType='d')
    else:
        dlink = copy.deepcopy(soap.dDlink)
    oriNum = dlink.linkNum()
    
    # start to find lost dlink, if there is pelink between two nodes, but no dlink rout, then re-link them
    # first find the links that has overlaps
    visited = set()
    for rkey, lkey in soap.dPlink.edges():
        # first check if there is a dlink rout can link two nodes, if exists pass (do nothing)
        if (rkey, lkey) in visited:
            continue
         # format of plink: 0 miu of gap 1 sigma of gap 3 support 4 score 5 lib id
        ae = soap.dPlink[rkey][lkey]
        if ae[2] < config.SUPPORTNUM:
            continue
        minGap, maxGap = ae[0]-4*ae[1], ae[0]+4*ae[1]
        if minGap < 0:
            visited.add((rkey, lkey))
            visited.add((soap.revNode(lkey), soap.revNode(rkey)))
            if rkey in dlink and lkey in dlink[rkey]:
                continue
            else:
                score = algin(graphSeqs[rkey], graphSeqs[lkey], soap.kmersize, 0.9)
                # update dlink
                if rkey not in [lkey, soap.revNode(lkey)]:
                    dlink.addEdge(rkey, lkey, score)
                    dlink.addEdge(soap.revNode(lkey), soap.revNode(rkey), score)
    
    # handle the gap size larger than zero
    backedLinks = set()
    for rkey, lkey in soap.dPlink.edges():
        # first check if there is a dlink rout can link two nodes, if exists pass (do nothing)
        if (rkey, lkey) in visited:
            continue
        visited.add((rkey, lkey))
        visited.add((soap.revNode(lkey), soap.revNode(rkey)))
        # 0 miu of gap 1 sigma of gap 3 support 4 score 5 lib id
        inf = soap.dPlink[rkey][lkey]
        if inf[2] < config.SUPPORTNUM:
            continue
        minGap, maxGap = inf[0]-4*inf[1], inf[0]+4*inf[1]
        if minGap > 0:
            # to see if there is a d-rout can link the two nodes
            flag = lenShortest(soap, dlink, rkey, lkey, maxGap)
            if not flag:
                # get  all of the node in the maxdis range and compare them
                try:
                    matches = stretchingMatch(soap, dlink, graphSeqs, rkey, lkey, minGap, maxGap)
                except:
                    continue
                if len(matches) == 1:
                    backedLinks.update(matches)
                elif len(matches) > 1:
                    counter = 0
                    for ele in matches:
                        score = algin(graphSeqs[ele[0]], graphSeqs[ele[1]], soap.kmersize, 1)
                        if score:
                            backedLinks.add(ele)
                            counter += 1
                    if not counter:
                        for ele in matches:
                            backedLinks.add(ele)
        
        # update dlinks
        for ele in backedLinks:
            if ele[0] not in [ele[1], soap.revNode(ele[1])]:
                dlink.addEdge(ele[0], ele[1], ele[2])
                dlink.addEdge(soap.revNode(ele[1]), soap.revNode(ele[0]), ele[2])
        
        newNum = dlink.linkNum()
        logger.debug('{0:,} new dlinks are created !'.format(newNum - oriNum))
        
        return  dlink
