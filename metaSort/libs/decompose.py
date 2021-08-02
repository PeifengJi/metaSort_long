############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import collections
import os
import cPickle as pickle
import itertools
import config
import sys
import math
import indexGraph
import copy
from pybloom import BloomFilter
from log import get_logger
import multiprocessing

# set log file
logger = get_logger()

def unwarp_self(worker, node):
    return worker.converge(node)

def child_initialize(graph, contigLen, indegree, outdegree):
    """
    in case of large dataset of graph is passed to the function (which will result only one cpu occupied)
    http://stackoverflow.com/questions/25825995/python-multiprocessing-only-one-process-is-running#
    the huge volume varible should be shared by child threads
    """
    global tmpGraph, length, ind, oud
    tmpGraph = graph
    length = contigLen
    ind = indegree
    oud = outdegree

def converge(node):
    
    #print multiprocessing.current_process()
    paths, visitedNodes, stack = [], set(), [[node]]
    stack = [[node]]
    while stack:
        tmp = stack.pop()
        last = tmp[-1]
        if last in visitedNodes:
            paths.append(tmp)
            continue
        visitedNodes.add(last)
        if last in tmpGraph:
            for sec in tmpGraph[last]:
                longer = tmp + [sec]
                pathLen = sum([length[i] for i in longer[1:-1]])
                if pathLen >= config.CRLEN:
                    paths.append(tmp)
                else:
                    stack.append(longer)
        else:
            paths.append(tmp)
    
    crossPoints = {i:set() for i in tmpGraph[node]}
    for path in paths:
        if len(path) > 1:
            crossPoints[path[1]].update(path[1:])
    
    links = collections.defaultdict(dict)
    combinations = itertools.combinations(crossPoints.keys(), 2)
    for comb in combinations:
        inter = set(crossPoints[comb[0]]).intersection(crossPoints[comb[1]])
        if inter:
            links[comb[0]][comb[1]] = links[comb[1]][comb[0]] = 1
    
    if not links:
        return (node, 1)
    else:
        downNodes = set(crossPoints.keys())
        left = downNodes - set(links.keys())
        if len(left) is not len(downNodes):
            suspects = set()
            for ele in left:
                if (ind[ele], oud[ele]) == (1,0):
                    pass
                else:
                    suspects.add(ele)
                    return (node, 1)
    return (node, 0)

class decomposition():
    """
    Decompose the mixed de bruijn graph
    """
    def __init__(self, soap, threads):
        
        self.soap = soap
        self.cpu = int(threads)
        logger.debug('Para crbranch length is: {0}'.format(config.CRLEN))
        #set output files
        self.ind, self.oud = self.soap.dDlink.degree()
        self.gpsfile = os.path.join(soap.database, 'paritionS')
    
    def covGroup(self, points):
        
        # set paramerters of gmm module
        clusterArgs = []
        for com in self.componentsNu:
            clusterArgs.append([com,
                                self.covariance,
                                self.inits,
                                self.iters,
                                points])
        
        # parallels to run gmm
        pool = multiprocessing.Pool(self.threads)
        results = pool.map(cluster, clusterArgs)
        
        # determine the cluster number by using BIC
        bics = [(r[0],r[1]) for r in results]
        
        # get the optimal cluster number
        min_bic, optimalC, converged = min(results,key=lambda x: x[0])
        if not converged:
            logging.error(("Optimal bic score was reached for non convergent "
                           "cluster number {0}, exiting without clustering "
                           "output").format(optimalC),
                          exit_with_code=1
                          )
        
        gmm = GMM(n_components=optimalC,
                  covariance_type=self.covariance,
                  n_init=self.inits,
                  n_iter=self.iters,
                  random_state=random_seed
                  ).fit(points)
        
        # use the posterior probabilities of each data point to cluster
        clusters = gmm.predict(points)
        
        return clusters
    
    def subgraph(self, graph, initNodes):
        
        # initialize the cluster
        cl = {i:1 for i in initNodes}
        tmp, count = [], 1
        for node in cl:
            if cl[node] is 1:
                tmp.append(node)
                while tmp:
                    last = tmp.pop()
                    cl[last] = count + 1
                    for second in graph[last]:
                        if second in cl and cl[second] is 1:
                            tmp.append(second)
                count += 1
        
        group = collections.defaultdict(set)
        for node in cl:
            group[cl[node]].add(node)
        
        return group
    
    def breakEdeges(self, nodes, sizeLimit):
        
        logger.debug('subgraph length: {0}'.format(len(nodes)/2))
        
        # prepare data
        tmpg = indexGraph.Graph(self.soap.index)
        transedG, gnodes = {}, set()
        tlen, tin, tout = 0, 0, 0
        crpoints = []
        for node in nodes:
            tlen += self.soap.length[node]
            tin += self.ind[node]
            out = self.oud[node]
            tout += out
            if out > 1:
                crpoints.append(node)
            # construct a new graph for decomposition
            if node in self.soap.dDlink:
                soure = node if self.soap.index[node] is 0 else node-1
                gnodes.add(soure)
                if soure not in transedG:
                    transedG[soure] = {}
                for sec in self.soap.dDlink[node]:
                    tmpg.addEdge(node, sec, 1)
                    sink = sec if self.soap.index[sec] is 0 else sec-1
                    transedG[soure][sink] = 1
                    gnodes.add(sink)
        logger.debug('crpoints number: {0}'.format(len(crpoints)))
        
        pool = multiprocessing.Pool(processes = self.cpu, initializer = child_initialize, initargs = (tmpg, self.soap.length, self.ind, self.oud, ))
        result = []
        for node in crpoints:
            result.append(pool.apply_async(converge, (node, )))
        pool.close()
        pool.join()
        
        candidates = set()
        for ele in result:
            res = ele.get()
            if res[1] == 1:
                candidates.add(res[0])
        
        # get the cr-branch nodes
        #candidates = set([node for node in crpoints if self.converge(tmpg, node)])
        logger.debug('Pre canidates: {0}'.format(len([i if self.soap.index[i] is 0 else i-1 for i in candidates])))
        
        # remove tips
        tips = set()
        for key in candidates:
            if self.ind[key] == 1 and self.oud[key] > 1:
                alive, end = set(), set()
                for sec in tmpg[key]:
                    if self.ind[sec] == 1 and self.oud[sec] == 0:
                        end.add(sec)
                    else:
                        alive.add(sec)
                if end and not alive:
                    tips.add(key)
                if end and alive:
                    if len(alive) is 1:
                        tips.add(key)
        
        tips.update([self.soap.revNode(i) for i in tips])
        candidates = candidates - tips
        logger.debug('Tips: {0}'.format(len(tips)))
        
        candidates = set([i if self.soap.index[i] is 0 else i-1 for i in candidates])
        if not candidates:
            return []
        
        # caculate the entropy of each node and sort them
        entropy = collections.defaultdict(int)
        for node in candidates:
            #print tin, tout, tlen
            s = [1-self.ind[node]*1.0/tin, 1-self.oud[node]*1.0/tout, self.soap.length[node]*1.0/tlen]
            entropy[node] = - sum(f*math.log(f, 2) for f in s if f !=0)
        
        order = sorted(entropy, key=lambda x:entropy[x])
        logger.debug('{} candidates'.format(len(order)))
        
        passed = []
        gsize = len(gnodes)
        if gsize > 1000000:
            step = 5000
            logger.debug('Too large the subgraph is, step {0} will be used.'.format(step))
            exclude = BloomFilter(capacity=gsize, error_rate=0.001)
            for pos in xrange(step-1, len(order)-step+1, step):
                for ele in order[:pos]:
                    _ = exclude.add(ele)
                gps = self.subgraph(transedG, (i for i in gnodes if i not in exclude))
                if len(gps) == 1:
                    continue
                count = 0
                for nu in gps:
                    if len(gps[nu]) < sizeLimit:
                        passed.append(gps[nu])
                        for ele in gps[nu]:
                            _ = exclude.add(ele)
                        count += 1
                if count is len(gps) or not gps:
                    break
            else:
                for nu in gps:
                    if len(gps[nu]) >= sizeLimit:
                        passed.append(gps[nu])
        elif gsize > 100000:
            step = 500
            logger.debug('Too large the subgraph is, step {0} will be used.'.format(step))
            exclude = BloomFilter(capacity=gsize, error_rate=0.001)
            for pos in xrange(step-1, len(order)-step+1, step):
                for ele in order[:pos]:
                    _ = exclude.add(ele)
                gps = self.subgraph(transedG, (i for i in gnodes if i not in exclude))
                if len(gps) == 1:
                    continue
                count = 0
                for nu in gps:
                    if len(gps[nu]) < sizeLimit:
                        passed.append(gps[nu])
                        for ele in gps[nu]:
                            _ = exclude.add(ele)
                        count += 1
                if count is len(gps) or not gps:
                    break
            else:
                for nu in gps:
                    if len(gps[nu]) >= sizeLimit:
                        passed.append(gps[nu])
        elif gsize > 10000:
            step = 50
            logger.debug('Too large the subgraph is, step {0} will be used.'.format(step))
            exclude = BloomFilter(capacity=gsize, error_rate=0.001)
            for pos in xrange(step-1, len(order)-step+1, step):
                for ele in order[:pos]:
                    _ = exclude.add(ele)
                gps = self.subgraph(transedG, (i for i in gnodes if i not in exclude))
                if len(gps) == 1:
                    continue
                count = 0
                for nu in gps:
                    if len(gps[nu]) < sizeLimit:
                        passed.append(gps[nu])
                        for ele in gps[nu]:
                            _ = exclude.add(ele)
                        count += 1
                if count is len(gps) or not gps:
                    break
            else:
                for nu in gps:
                    if len(gps[nu]) >= sizeLimit:
                        passed.append(gps[nu])
        else:
            exclude = set()
            for node in order:
                if node in exclude:
                    continue
                exclude.add(node)
                gps = self.subgraph(transedG, gnodes-exclude)
                if len(gps) is 1:
                    continue
                count = 0
                for nu in gps:
                    if len(gps[nu]) < sizeLimit:
                        exclude.update(gps[nu])
                        passed.append(gps[nu])
                        count += 1
                if count is len(gps) or not gps:
                    break
            else:
                for nu in gps:
                    if len(gps[nu]) >= sizeLimit:
                        passed.append(gps[nu])
        
        #aver = sum([len(i) for i in passed])/len(passed)
        #logger.info('It is divided into {0} disconnect components, average number of nodes is: {1}'.format(
        #                    len(passed), aver))
        
        return passed
    
    def partition(self):
        
        # first check if the files exists
        if os.path.exists(self.gpsfile) and os.path.getsize(self.gpsfile):
            logger.info('File: partitionS exists, skip')
        else:
            # find subgraph of the total grpah
            oriGroup = self.soap.dDlink.findSubgraph(out='gp')
            top10ids = sorted(oriGroup, key=lambda x:len(oriGroup[x]), reverse=True)[:10]
            oriNum = [len(oriGroup[i]) for i in top10ids]
            logger.info('\n###################################################################')
            logger.info('Befor disconnection, there are {} groups'.format(len(oriGroup)))
            logger.info('The size of top 10 largest group are: \n{}'.format(oriNum))
            
            logger.info('Disconnecting the graphs ...')
            # decompose the total graph into small subgraphs
            gps, count = {}, 1
            for nu in oriGroup:
                if len(oriGroup[nu]) >= 2*config.SMALLSIZE:
                    passed = self.breakEdeges(oriGroup[nu], config.SMALLSIZE)
                    for ele in passed:
                        gps[count] = ele
                        count += 1
                else:
                    gps[count] = set([i for i in oriGroup[nu] if self.soap.index[i]==0])
                    count += 1
            
            # store the decomposition file
            with open(self.gpsfile, 'wb') as af:
                newg = {}
                for nu in gps:
                    if len(gps[nu])>1:
                        newg[nu] = gps[nu]
                pickle.dump(newg, af)
            
            aft = sorted(newg, key=lambda x:len(newg[x]), reverse=True)[:10]
            aftNum = [len(newg[i]) for i in aft]
            logger.info('\nTotal {} groups are generated'.format(len(newg)))
            logger.info('The size of top 10 largest group are: \n{}'.format(aftNum))
            logger.info('#####################################################################\n')
            logger.info('Cluster data is stored into file: {0}.'.format(self.gpsfile))

class acclibs():
    
    def __init__(self, soap, metaSdir):
        
        self.fn = os.path.join(metaSdir, config.LIBPATHS)
        self.soap = soap
    
    def update(self, links):
        
        counter = 0
        if not os.path.exists(self.fn) or not os.path.getsize(self.fn):
            paths = []
            for path in links.itervalues():
                orders, gaps = [], []
                for pos, align in enumerate(path):
                    node = align[-1] if align[-2] is "1" else self.soap.revNode(align[-1])
                    orders.append(node)
                    if pos > 0:
                        up, down = path[pos-1], path[pos]
                        gaps.append([
                                    down[0]-up[0],
                                    down[0]-up[1],
                                    down[1]-up[0],
                                    down[1]-up[1]
                                     ])
                paths.append([orders, gaps])
        else:
            with open(self.fn, 'rb') as af:
                paths = pickle.load(af)
            
            for path in links.itervalues():
                orders, gaps, refresh = [], [], 2
                for pos, align in enumerate(path):
                    node = align[-1] if align[-2] is "1" else self.soap.revNode(align[-1])
                    orders.append(node)
                    if pos > 0:
                        up, down = path[pos-1], path[pos]
                        gaps.append([
                                    down[0]-up[0],
                                    down[0]-up[1],
                                    down[1]-up[0],
                                    down[1]-up[1]
                                     ])
                for pos, ele in enumerate(paths):
                    if set(ele[0]).issubset(set(orders)):
                        position = pos
                        refresh = 1
                        break
                    elif set(ele[0]).issuperset(set(orders)):
                        refresh = 0
                
                if refresh is 1:
                    paths.pop(position)
                    paths.append([orders, gaps])
                    counter += 1
                elif refresh is 2:
                    paths.append([orders, gaps])
                    counter += 1
            logger.info('Total {0} new paths are added'.format(counter))
           
        # index the node position for the sake of running time
        index = {}
        for pos, path in enumerate(paths):
            for node in path[0]:
                if node not in index:
                    index[node] = set()
                rev = self.soap.revNode(node)
                if rev not in index:
                    index[rev] = set()
                index[node].add(pos)
                index[rev].add(pos)
        
        with open(self.fn, 'wb') as af:
            pickle.dump(paths, af)
            pickle.dump(index, af)
    
    def load(self):
        
        with open(self.fn, 'rb') as af:
            self.paths = pickle.load(af)
            self.index = pickle.load(af)

