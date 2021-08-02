############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import collections
import os
import copy
from log import get_logger

# set log file
logger = get_logger()

class Graph:
    
    def __init__(self, index):
        self.__graphDict = {}
        self.index = index
    
    def __len__(self):
        """
        Get the keys of graph
        """
        return len(self.__graphDict)
    
    def addEdge(self, start, end, addtionInfo):
        if start not in self.__graphDict:
            self.__graphDict[start] = {}
        
        if end not in self.__graphDict[start]:
            self.__graphDict[start][end] = addtionInfo
    
    def __contains__(self, key):
        if not isinstance(key, int):
            logger.error('Please provide a int to represent the node on the graph',
                         exit_with_code=1
                         )
        if key in self.__graphDict:
            return 1
        else:
            return 0
    
    def __getitem__(self, key):
        """
        get the keys of dict of dict, key must be provided
        """
        if not key:
            raise Exception('Please provide a key to graph, like: graph[key]')
            #logger.error('Please provide a key to graph, like: graph[key]',
            #             exit_with_code=1
            #             )
        
        if key in self.__graphDict:
            return self.__graphDict[key]
        else:
            rev = reverseNode(key, self.index)
            if rev in self.__graphDict:
                raise Exception('node ' + str(key) + ' is a append node')
            else:
                raise Exception('node ' + str(key) + ' is not on the graph')
    
    def bfs(self, node, step = 5):
        """
        use BFS algorithm to search the nodes on the graph in a given length range,
        because there may be some subspecies in the path, so we used a step size to get the structure
        """
        queue = collections.deque([[node]])
        searched = []
        # there is no visitedNodes, no matter the redundant, beacuse if we set vistiedNodes
        # there will be a huge rate of unbalance link in two strand (caused by bubbule)
        while queue:
            tmp = queue.pop()
            last = tmp[-1]
            if last in self.__graphDict:
                for second in self.__graphDict[last]:
                    longerPath = tmp + [second]
                    if second in tmp: # a cycle formed
                        searched.append(longerPath)
                    elif len(longerPath) == step + 1:
                        searched.append(longerPath)
                    else:
                        queue.appendleft(longerPath)
            else:
                searched.append(tmp)
        
        return searched   # return a structer link [[1,2,4],[3,5,6],[2,3,4]]
    
    def keys(self):
        """
        Only get the nodes the outdegree larger than one, return a list
        """
        return self.__graphDict.keys()
    
    def iterkeys(self):
        """
        Only get the nodes the outdegree larger than one, return a generator
        """
        return self.__graphDict.iterkeys()
    
    def vertices(self):
        """
        get all the nodes on the graph
        """
        vertex = set()
        for node in self.__graphDict:
            vertex.add(node)
            for second in self.__graphDict[node]:
                vertex.add(second)
        
        return vertex
    
    def findSubgraph(self, out='cl'):
        """
        find the subgraphs of the graph
        NOTE: just suite for dlink graph
        """
        cl = {}
        tmp = []
        count = 1
        undGraph = {}
        # turn the graph into a undirected graph
        for ele in self.__graphDict:
            source = ele-1 if self.index[ele] == 2 else ele
            cl[source] = 1
            if source not in undGraph:
                undGraph[source] = {}
            for node in self.__graphDict[ele]:
                sink = node-1 if self.index[node] == 2 else node
                undGraph[source][sink] = 1
                cl[sink] = 1
        
        for node in cl:
            if cl[node] == 1:
                tmp.append(node)
                while tmp:
                    last = tmp.pop()
                    cl[last] = count + 1
                    for second in undGraph[last]:
                        if cl[second] == 1:
                            tmp.append(second)
                count += 1
        for key in cl.keys(): # you can not use iterkeys(), beacuse we change the size of dict: cl
            if self.index[key] == 0:
                cl[key + 1] = cl[key]
        
        group = {}
        for ele in cl:
            if cl[ele] not in group:
                group[cl[ele]] = []
            group[cl[ele]].append(ele)
        
        if out == 'cl':
            return cl
        elif out == 'gp':
            return group
        elif out == 'b':
            return cl, group
    
    def checkSymetric(self):
        """
        Check if the graph is symmetrical, if not some nodes are not deleted
        """
        #first check if the nodes are ok
        plus, minus = set(), set()
        for node in self.vertices():
            if self.index[node]==0:
                plus.add(node)
            elif self.index[node]==2:
                minus.add(node-1)
        if len(plus) != len(minus):
            print(plus - minus)
            logger.error('Error, the graph is not symmetrical',
                         exit_with_code=1
                         )
        # then links
        # links = set([(s, e) for s in self.__graphDict for e in self.__graphDict[s]])
        links = set()
        for s in self.__graphDict:
            for e in self.__graphDict[s]:
                links.add((s,e))
        # check
        record = set()
        for ele in links:
            if ele not in record:
                revEle = tuple([reverseNode(i, self.index) for i in ele[::-1]])
                if revEle in links:
                    # ok, symetric
                    record.update([ele, revEle])
                else:
                    logger.error('no links between ' + str(revEle[0]) + ' ' + str(revEle[1]),
                                 exit_with_code=1
                                 )
        return 1
    
    def deleteNodes(self, nodes):
        """
        Delete the nodes provided, and refresh the graph. Also output the singltions
        Because when delete some nodes, some very small subgraphs will contain only one node
        """
        if isinstance(nodes, int) or isinstance(nodes, basestring):
            logger.error('The parameter should be a list or set or tuple, please check it',
                         exit_with_code=1
                         )
        
        # reverse complement the node
        pairs = set([k for i in nodes for k in (i, reverseNode(i, self.index))])
        
        preNodes, aftNodes, graph = set(), set(), {}
        # output to a list
        for key in self.__graphDict:
            for second in self.__graphDict[key]:
                edge = [key, second]
                preNodes.update(edge)
                if not pairs.intersection(edge):
                    if key not in graph:
                        graph[key] = {}
                    graph[key][second] = self.__graphDict[key][second]
                    aftNodes.update(edge)
        
        # new indicator
        self.__graphDict = graph
        self.checkSymetric()
        
        # get singltons
        singltons = preNodes - aftNodes - pairs
        singltons.update([reverseNode(i, self.index) for i in singltons])
        
        return singltons
    
    def cutEdge(self, upNode, downNode=None):
        if not downNode:
            del self.__graphDict[upNode]
        else:
            del self.__graphDict[upNode][downNode]
            if len(self.__graphDict[upNode]) == 0:
                del self.__graphDict[upNode]
    
    def edges(self):
        """
        Output the edges of the graph, generator
        """
        for fir in self.__graphDict:
            for sec in self.__graphDict[fir]:
                yield [fir, sec]
    
    def linkInfo(self, start, end):
        """
        get the link info between two linked nodes, pelink: [], dlink:[]
        """
        try:
            return self.__graphDict[start][end]
        except:
        #if start in self.__graphDict:
        #    if end in self.__graphDict[start]:
        #        return self.__graphDict[start][end]
        #    else:
        #        logger.error('There is no linkage between ' + str(start) + ' and ' + str(end),
        #                     exit_with_code=1
        #                     )
        #else:
            logger.error('node ' + str(start) + ' not on the graph',
                         exit_with_code=1
                         )
    
    def linkNum(self):
        """
        Return the number of edges of the graph
        """
        count = 0
        for key in self.__graphDict:
            for sec in self.__graphDict[key]:
                count += 1
        return count
    
    def degree(self, node=0):
        """
        Every time after deleting the node on the graph, the in and outdegree should be refreshed.
        """
        if not node:
            visited = set()
            ind = {}
            outd = {}
            for plus in self.__graphDict:
                if plus not in visited:
                    minus = reverseNode(plus, self.index)
                    o = len(self.__graphDict[plus])
                    if minus in self.__graphDict:
                        i = len(self.__graphDict[minus])
                    else:
                        i = 0
                    outd[plus], ind[plus], outd[minus], ind[minus] = o, i, i, o
                    visited.add(plus)
                    visited.add(minus)
            return ind, outd
        else:
            ind, oud = 0, 0
            if isinstance(node, int):
                if node in self.__graphDict:
                    oud = len(self.__graphDict[node])
                rev = reverseNode(node, self.index)
                if rev in self.__graphDict:
                    ind = len(self.__graphDict[rev])
            return (ind, oud)
    
    def indegree(self, node):
        rev = reverseNode(node, self.index)
        if node not in self.__graphDict and rev not in self.__graphDict:
            logger.error('node ' + str(node) + ' is not on the graph',
                         exit_with_code=1
                         )
        elif rev in self.__graphDict:
            return len(self.__graphDict[rev])
        elif node in self.__graphDict and rev not in self.__graphDict:
            return 0
    
    def outdegree(self, node):
        rev = reverseNode(node, self.index)
        if node in self.__graphDict:
            return len(self.__graphDict[node])
        if node not in self.__graphDict and rev not in self.__graphDict:
            logger.error('node ' + str(node) + ' is not on the graph',
                         exit_with_code=1
                         )
        elif node not in self.__graphDict and rev in self.__graphDict:
            return 0
    
    def __backTrace(self, rout, visitedNodes):
        """
        back trace the chain find branches
        """
        splitPaths = []
        restart = set()
        if len(rout) == 1:
            rev = reverseNode(rout[0], self.index)
            visitedNodes.update([rout[0], rev])
            if rout[0] in self.__graphDict:
                nextNodes = self.__graphDict[rout[0]].keys()
                restart.update(nextNodes)
            if rev in self.__graphDict:
                nextNodes = self.__graphDict[rev].keys()
                restart.update(nextNodes)
            return [rout], restart
        
        frag = []
        revRout = [reverseNode(i, self.index) for i in rout]
        while rout:
            node = rout.pop()
            revNode = reverseNode(node, self.index)
            frag.append(revNode)
            if revNode in self.__graphDict and len(self.__graphDict[revNode]) > 1:
                nextNodes = self.__graphDict[revNode].keys()
                splitPaths.append([reverseNode(i, self.index) for i in frag[::-1]])
                frag = []
                for ele in nextNodes:
                    if ele not in revRout:
                        restart.add(ele)
            if len(rout) == 0: # get the last node
                splitPaths.append([reverseNode(i, self.index) for i in frag[::-1]])
                break
        
        visitedNodes.update(rout)
        visitedNodes.update(revRout)
        return splitPaths, restart
    
    def compressGraph(self, linear=0):
        """
        Compress the graph by merging nodes that indegree = outdegree = 1
        output the simple chains at the same time
        """
        # find the loop
        loop = set()
        record = set()
        for node in self.__graphDict:
            if node not in record:
                record.add(node)
                for key in self.__graphDict[node]:
                    if key in self.__graphDict and node in self.__graphDict[key]:
                        loop.update([node, key])
        
        indegree, outdegree = self.degree()
        group = self.findSubgraph(out='gp')
        
        # traversal the each of the subraphs
        paths, circles = [], []
        for key in group:
            elements = group[key]
            visitedNodes = set()
            startNodes = collections.deque([i for i in elements if indegree[i] == 0])
            if not startNodes:
                # may be circle or may be loops exists
                startNodes = collections.deque(list(set(elements) & loop))
            if not startNodes:
                startNodes = collections.deque(elements)
            
            while startNodes:
                node = startNodes.pop()
                if node not in visitedNodes:
                    stack = [[node]]
                    while stack:
                        tmp = stack.pop()
                        last = tmp[-1]
                        if last in tmp[:-1]:
                            if len(tmp[:-1])*2 == len(group[key]): # a circle
                                splitPaths, restart = self.__backTrace(tmp[:-1], visitedNodes)
                                circles.extend(splitPaths)
                                for i in restart: startNodes.appendleft(i)
                                continue
                        if last in visitedNodes:
                            splitPaths, restart = self.__backTrace(tmp[:-1], visitedNodes)
                            paths.extend(splitPaths)
                            for i in restart: startNodes.appendleft(i)
                            continue
                        visitedNodes.update([last, reverseNode(last, self.index)])
                        if last in self.__graphDict:
                            nextNodes = self.__graphDict[last].keys()
                            if len(nextNodes) == 1:
                                longer = tmp + nextNodes
                                stack.append(longer)
                            else:
                                for ele in nextNodes:
                                    if ele not in visitedNodes:
                                        startNodes.appendleft(ele)
                                splitPaths, restart = self.__backTrace(tmp, visitedNodes)
                                paths.extend(splitPaths)
                                for i in restart: startNodes.appendleft(i)
                        else:
                            splitPaths, restart = self.__backTrace(tmp, visitedNodes)
                            paths.extend(splitPaths)
                            for i in restart: startNodes.appendleft(i)
            left = set(elements) - visitedNodes
            if left:
                print 'group ' + str(key)
                print left
                logger.error('There are nodes not visited',
                             exit_with_code=1
                             )
        
        simpleChain = []
        simpleChain.extend(circles)
        # check every path and delete the loop structure
        for path in copy.deepcopy(paths):
            sect = loop & set(path)
            if len(sect) > 1:
                if len(path)==2:
                    paths.remove(path)
                else:
                    pos = sorted([path.index[i] for i in sec])
                    paths.append(path[:pos[0]])
                    paths.append(path[pos[-1]+1:])
        
        graph = copy.deepcopy(self.__graphDict)
        # first, delete the circle nodes from the graph
        for cir in circles:
            for node in cir:
                revNode = reverseNode(node, self.index)
                if node in graph: del graph[node]
                if revNode in graph: del graph[revNode]
        
        # filter the paths, length, in and out degree of the two end nodes
        filteredPaths = set()
        for path in paths:
            if len(path) > 1: # check out the two end nodes
                flagF, flagT = 0, 0
                firstDeg = (indegree[path[0]], outdegree[path[0]])
                lastDeg = (indegree[path[-1]], outdegree[path[-1]])
                if firstDeg==(0,1) or firstDeg==(1,1):
                    flagF = 1
                if lastDeg==(1,0) or lastDeg==(1,1):
                    flagT = 1
                if not flagF:
                    path = path[1:]
                if not flagT:
                    path = path[:-1]
                if len(path) > 1:
                    filteredPaths.add(tuple(path))
        # first check if the paths containe fake connection
        for path in filteredPaths:    
            for pos in xrange(len(path)-1):
                if path[pos+1] not in graph[path[pos]]:
                    logger.error('Nodes: ' + str(path[pos]) + 'and ' + str(path[pos+1]) + ' are linked incorrectly.',
                                 exit_with_code=1
                                 )
        # check the path to merge the nodes
        combinedNodes = {}
        for path in filteredPaths:
            revPath = [reverseNode(i, self.index) for i in path[::-1]]
            if indegree[path[0]]==0 and outdegree[path[-1]]==0: # a very simple path, delete them directly
                simpleChain.append(path)
                for ele in path:
                    if ele in graph:
                        del graph[ele]
                for ele in revPath:
                    if ele in graph:
                        del graph[ele]
            elif indegree[path[0]]==1 and outdegree[path[-1]] == 1:
                try:
                    combinedNodes[path[0]] = path
                    combinedNodes[revPath[-1]] = revPath
                    graph[path[0]] = graph[path[-1]]
                    reverKey = reverseNode(graph[path[-1]].keys()[0], self.index) 
                    if len(graph[reverKey]) == 1:
                        graph[reverKey] = {revPath[-1]:1}
                    else:
                        del graph[reverKey][revPath[0]]
                        graph[reverKey][revPath[-1]] = 1
                    for i in path[1:]: del graph[i]
                    for i in revPath[:-1]: del graph[i]
                except:
                    pass
            elif indegree[path[0]]==1 and outdegree[path[-1]] == 0:
                combinedNodes[path[0]] = path
                for i in path[:-1]: del graph[i]
                combinedNodes[revPath[-1]] = revPath
                for i in revPath[:-1]:del graph[i]
            elif indegree[path[0]]==0 and outdegree[path[-1]]==1:
                combinedNodes[revPath[0]] = revPath
                for i in revPath[:-1]: del graph[i]
                combinedNodes[path[-1]] = path
                for i in path[:-1]: del graph[i]
        
        # check symetric
        links = set()
        for s in graph:
            for e in graph[s]:
                links.add((s,e))
        # check
        record = set()
        for ele in links:
            if ele not in record:
                revEle = tuple([reverseNode(i, self.index) for i in ele[::-1]])
                if revEle in links:
                    # ok, symetric
                    record.update([ele, revEle])
                else:
                    logger.error('no links between ' + str(revEle[0]) + ' ' + str(revEle[1]),
                                 exit_with_code=1
                                 )
        
        if linear == 0:
            return graph, combinedNodes
        elif linear == 1:
            return graph, combinedNodes, simpleChain
        else:
            logger.error('The parameter: linear of method compressGraph should be 1 or 0',
                         exit_with_code=1
                         )
    
    def shortestPath(self, start, end, step = 5):
        """
        Use bfs algorithm to get the true dirrection and order of two plink nodes
        Because we are using a PE lib of 180 bp. There is no reason the step was setted with a larger number.
        """
        visitedNodes = set()
        queue = collections.deque([[start]])
        searched = []
        path = []
        poll = set([end, reverseNode(end, self.index)])
        while len(queue) > 0:
            tmp = queue.pop()
            last = tmp[-1]
            if last in visitedNodes:
                continue
            visitedNodes.add(last)
            if last in self.__graphDict:
                for second in self.__graphDict[last]:
                    if second not in poll:
                        longerPath = tmp + [second]
                        if len(longerPath) < step + 1:
                            queue.appendleft(longerPath)
                    elif second in poll:
                        path = tmp + [second]
                        return path
        return path
    
    def allPaths(self, start, end, step = 10):
        """
        get all of the paths between two nodes
        """
        queue = collections.deque([[start]])
        pool = set([end, reverseNode(end, self.index)])
        paths = set()
        # there is no visitedNodes, no matter the redundant, beacuse if we set vistiedNodes
        # there will be a huge rate of unbalance link in two strand (caused by bubbule)
        while queue:
            tmp = queue.pop()
            last = tmp[-1]
            if last in pool:
                paths.add(tuple(tmp))
                continue
            if len(tmp) == step+1:
                continue
            if last in tmp[:-1]: # cycle
                continue
            if last in self.__graphDict:
                for second in self.__graphDict[last]:
                    queue.appendleft(tmp + [second])
        
        return paths
    
    def plinkNodes(self, node, pelink):
        """
        Get the pelinked nodes a the start node
        
        """
        if not node:
            logger.error('Give a node to the method: graph.plinkNodes()',
                         exit_with_code=1
                         )
        
        aroudNodes = set()
        rev = reverseNode(node, self.index)
        
        if node in pelink:
            aroudNodes.update(pelink[node].keys())
            aroudNodes.update([reverseNode(i, self.index) for i in pelink[node].keys()])
        
        if rev in pelink:
            aroudNodes.update(pelink[rev].keys())
            aroudNodes.update([reverseNode(i, self.index) for i in pelink[rev].keys()])
        
        return aroudNodes

    def plotGraph(self):
        """
        output the edges for plotting
        """
        with open('plotG','w') as af:
            for key in self.__graphDict:
                for second in self.__graphDict[key]:
                    af.write(str(key) + ' ' + str(second) + '\n')

def subGraph(graph, index, out='s'):
    """
    find subgraphs including sperating complementary strand
    """
    cl = {}
    group = {}
    count = 1
    tmp = []
    # undirected graph
    ug = {}
    for node in graph.iterkeys():
        cl[node] = 1
        if node not in ug:
            ug[node] = {}
        
        for ele in graph[node]:
            cl[ele] = 1
            if ele not in ug:
                ug[ele] = {}
            ug[node][ele] = 1
            ug[ele][node] = 1
    
    for node in cl:
        if cl[node] == 1:
            tmp.append(node)
            while tmp:
                last = tmp.pop()
                cl[last] = count + 1
                for second in ug[last]:
                    if cl[second] == 1:
                        tmp.append(second)
            count += 1
    
    for node in cl:
        if cl[node] not in group:
            group[cl[node]] = set()
        group[cl[node]].add(node)
    
    if out == 's':
        return group
    elif out == 'p':
        record = set()
        pairs = set()
        for key in group.iterkeys():
            if key not in record:
                rev = set([reverseNode(i,index) for i in group[key]])
                for cycle in group:
                    if group[cycle] == rev:
                        if key == cycle:
                            pairs.add(tuple([key]))
                            record.add(key)
                            break
                        else:
                            order = tuple(sorted([key, cycle]))
                            pairs.add(order)
                            record.update([key, cycle])
                            break
        
        count = len([k for i in pairs for k in i])
        if len(group) != count:
            logger.error('Group number is not included by the pairs, \
                         it may caused by nosymetric of the compressed function in the Class Graph',
                         exit_with_code=1
                        )
        
        return group, pairs

def checkFile(*f):
    """
    check if the parameters are exists
    """
    
    for tmp in f:
        if not os.path.exists(tmp) or not os.path.getsize(tmp):
            logger.error('file: ' + tmp + ' not exists or the file size is zero, please recheck it',
                         exit_with_code=1
                         )

def indexLinks(fn, index, kmersize, linkType='d'):
    """
    linkType: d -- arcfile (the Arc file generated by SOAPdenovo)
    linkType: p -- linkfile (the links file gnerated by SOAPdenovo)
    """
    if isinstance(kmersize, str):
        logger.error('Kmersize should be int',
                     exit_with_code=1
                     )
    elif kmersize <= 0:
        logger.error('Kmersize should be larger than 0 !',
                     exit_with_code=1
                     )
    
    graph = Graph(index)
    if linkType == 'd':
        logger.info('Indexing arc file ...')   
        # generate directed_dlink, undirectedlink and cl
        with open(fn,'r') as b_file:
            for line in b_file:
                arry = map(int, line.rstrip().split())
                # index directed graph
                for incre in xrange(1,len(arry),2):
                    graph.addEdge(arry[0], arry[incre], -kmersize)
    elif linkType == 'p':
        logger.info('Indexing link file ...')
        with open(fn,'r') as c_file:
            for line in c_file:
                info = map(int,line.rstrip().split())
                graph.addEdge(info[0], info[1], info[2:])
    
    return graph

def indexSmallGraph(index, *files):
    """
    indexing the graph of target genome
    """
    # the link info: pelink: [gap, suport number, libsize], dlink: number, newlink: [overlap numer]
    graph = Graph(index)    
    for filename in files:
        if os.path.exists(filename) and os.path.getsize(filename):
            logger.info('reading file: ' + filename )
            with open(filename, 'r') as a_file:
                for line in a_file:
                    info = []
                    for ele in line.rstrip().split():
                        try:
                            info.append(int(ele))
                        except:
                            info.append(float(ele))
                    if len(info) == 3:
                        graph.addEdge(info[0], info[1], int(info[2]))
                    else:
                        graph.addEdge(info[0], info[1], info[2:])
    try:
        graph.checkSymetric()
    except:
        logger.error('The link file if not symetric, please recheck it',
                     exit_with_code=1
                     )
    
    return graph

def reverseNode(node, index):
    """
    Get the node id of it's  reverse complementery strand node
    """
    if index[node] == 0:
        return node + 1
    elif index[node] == 2:
        return node - 1
    else:
        return node
