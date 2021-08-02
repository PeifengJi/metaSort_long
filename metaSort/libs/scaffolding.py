############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import config
import collections
import copy
import sys
import os
import math
import cPickle as pickle
import indexGraph
import config
from indexGraph import subGraph, Graph
import itertools
from indexGraph import indexSmallGraph
from statistic import contigInfo, scaffoldInfo
from log import get_logger

# set log file
logger = get_logger()

class corGraph:
    
    def __init__(self, seqs):
        
        self.seqs = seqs
        
    def combination(self, size):
        """
        Generate 256 types of permutation of tertra nuclotides
        """
        combination = []
        S = ['A','T','G','C']
        for seq in itertools.product(S, repeat = size):
                combination.append(''.join([i for i in seq]))
        combination.sort()
        return combination
    
    def standarFreq(self, freq, size):
        """
        Normalize the frequency of the frequency table
        """
        sumOccurences = sum(freq.values())
        coeff = math.pow(4, size)/sumOccurences
        for key in freq:
            normalized = freq[key]*coeff
            freq[key] = normalized
    
    def revCom(self, string):
        """
        reverse complement the given string
        """
        comp = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C', 'N':'N'}
        result = ''
        for i in string[::-1]:
            result += comp[i]
        return result
    
    def zScore(self, seq, size):
        """
        Caculate the z-scores of the sequence, scaf length cutoff > 5000
        both of the two strand should be caculated
        """
        # get the frequency of the tera, tri, di- oligo of the seq
        revSeq = self.revCom(seq)
        # tertra or tri
        arr1 = ( seq[i:i+size].upper() for i in xrange(len(seq)-size+1) )
        arr2 = ( revSeq[i:i+size].upper() for i in xrange(len(seq)-size+1) )
        oligoF = collections.defaultdict(int)
        oligoS = collections.defaultdict(int)
        oligoT = collections.defaultdict(int)
        for i in arr1: oligoF[i] += 1
        for i in arr2: oligoF[i] += 1
        # remove the 'N' contained kmers
        for key in itertools.ifilter(lambda x: set(x) - set('ATCG'), oligoF.keys()):
            del oligoF[key]
        # normalize the frequency
        self.standarFreq(oligoF, size)
        for key, value in oligoF.iteritems():
            oligoS[key[:size-1]] += value
            oligoT[key[:size-2]] += value
        # caculate the z-score: Teeling et al. BMC Bioinformatics 2004; 5:163.
        # Zxyzw = (Nxyzw - Exyzw)/Var(Nxyzw)**0.5
        # Exyzw = Nxyz*Nyzw/Nyz
        # Var(Nxyzw) = Exyzw*(Nyz-Nxyz)(Nyz-Nyzw)/Nyz**2
        zScores = []
        for XYZW in self.combination(size):
            try:
                XYZ, YZW, YZ = XYZW[:size-1], XYZW[1:], XYZW[1:size-1]
                exp = oligoS[XYZ]*oligoS[YZW]/oligoT[YZ]
                var = exp*(oligoT[YZ]-oligoS[XYZ])*(oligoT[YZ]-oligoS[YZW])/math.pow(oligoT[YZ], 2)
                z = (oligoF[XYZW] - exp)/math.sqrt(var)
                zScores.append(z)
            except:
                zScores.append(0)
        
        return zScores

    def pearsonCorr(self, zX, zY):
        """
        Use the pearson correlation coefficent when the seq length is long
        """
        vectLen = len(zX)
        sumX, sumXsq, sumY, sumYsq = 0, 0, 0, 0
        for ele in zX:
            sumX += ele
            sumXsq += math.pow(ele,2)
        for ele in zY:
            sumY += float(ele)
            sumYsq += pow(float(ele),2)
        # add up the products
        product = sum([zX[pos]*zY[pos] for pos in xrange(vectLen)])
        # caculate the pearson correlation score
        numerator = product - (sumX*sumY/vectLen)
        denominator = math.sqrt((sumXsq - pow(sumX, 2)/vectLen)*(sumYsq - pow(sumY, 2)/vectLen))
        # can not have division by 0
        if denominator == 0: return 0
        cor = numerator/denominator
        return cor

    def buildGraph(self):
        
        graph = collections.defaultdict(dict)
        zS4, zS3 = {}, {}
        for key in self.seqs:
            z4 = self.zScore(self.seqs[key], 4)
            z3 = self.zScore(self.seqs[key], 3)
            zS4[key] = z4
            zS3[key] = z3
        
        combs = itertools.combinations(self.seqs.keys(), 2)
        for eles in combs:
            dis4 = self.pearsonCorr(zS4[eles[0]], zS4[eles[1]])
            if dis4 >= 0.9:
                graph[eles[0]][eles[1]] = graph[eles[1]][eles[0]] = dis4
            elif 0.7< dis4 < 0.9 and (len(self.seqs[eles[0]])< 30000 or len(self.seqs[eles[1]]) < 30000):
                dis3 = self.pearsonCorr(zS3[eles[0]], zS3[eles[1]])
                if dis3 >= 0.85:
                    graph[eles[0]][eles[1]] = graph[eles[1]][eles[0]] = dis3
        return graph

class scafWalking:
    
    def __init__(self, index, spalinks, soapScafs):
        self.pairSpa = {}
        self.relate = {}
        self.spa = spalinks # spa newlinks
        self.index= index
        self.allSoap = soapScafs # soapdenov scaf links
    
    def reverseNodes(self, node):
        
        if self.index[node] == 0:
            return node + 1
        elif self.index[node] == 2:
            return node - 1
        else:
            return node
    
    def prepare(self):
        """
        Make some transformation for scaffolding walking
        """
        allSides, ref, num = set(), {}, -1
        for ele in self.spa:
            path = self.spa[ele]
            ref[path[0]] = ref[path[-1]] = num
            self.pairSpa[num] = path
            revPath = [self.reverseNodes(i) for i in path[::-1]]
            ref[revPath[0]] = ref[revPath[-1]] = num-1
            allSides.update([path[0], path[-1], revPath[0], revPath[-1]])
            self.pairSpa[num-1] = revPath
            self.relate[num] = self.relate[num-1] = ele
            num -= 2
        
        allNodes, soapScaf, num = set(), {}, 1
        backTrace = {}
        for path in self.allSoap:
            soapScaf[num] = list(path)
            revPath = [self.reverseNodes(i) for i in path[::-1]]
            soapScaf[num+1] = revPath
            allNodes.update(path)
            allNodes.update(revPath)
            for ele in path:
                backTrace[ele] = num
            for ele in revPath:
                backTrace[ele] = num + 1
            num += 2
        
        return allSides, ref, allNodes, soapScaf, backTrace
    
    def spaExtend(self, source, visited, allNodes, soapScaf, backTrace, direct):
        
        swith = 1
        # add to visited
        if source%2!=0:
            visited.update([source, source-1])
        else:
            visited.update([source, source+1])
        # extend
        if direct == 1:
            rnode = self.pairSpa[source][-1]
            if rnode in allNodes:
                soapid = backTrace[rnode]
                if soapid in visited:
                    swith = 0
                    return swith, [] 
                path = soapScaf[soapid]
                if path.index(rnode)!= len(path)-1:
                    truncatePath = path[path.index(rnode)+1:]
                    return swith, truncatePath
                else:
                    swith = 0
                    return swith, []
            else:
                swith = 0
                return swith, []
        
        if direct == -1:
            lnode = self.pairSpa[source][0]
            if lnode in allNodes:
                soapid = backTrace[lnode]
                if soapid in visited:
                    swith = 0
                    return swith, []
                path = soapScaf[soapid]
                if path.index(lnode) != 0:
                    truncatePath = path[:path.index(lnode)]
                    return swith, truncatePath
                else:
                    swith = 0
                    return swith, []
            else:
                swith = 0
                return swith, []
    
    def soapExtend(self, nodes, visited, allSides, ref, backTrace, direct):
        
        swith = 1
        if not isinstance(nodes, list):
            logger.error('List should be given',
                         exit_with_code=1
                         )
        soapid = backTrace[nodes[-1]]
        if soapid%2 != 0:
            visited.update([soapid, soapid + 1])
        else:
            visited.update([soapid, soapid - 1])
        # extend
        inter = set(nodes).intersection(allSides) # one math
        if len(inter)==1:
            match = list(inter)[0]
            spaid = ref[match]
            if spaid in visited:
                swith = 0
                return swith, []
            if direct == 1:
                if match == self.pairSpa[spaid][0]:
                    truncate = nodes[:nodes.index(match)]
                    truncate.append(spaid)
                    return swith, truncate
                else:
                    swith = 0
                    return swith, []
            else:
                if match == self.pairSpa[spaid][-1]:
                    truncate = nodes[nodes.index(match)+1:]
                    return swith, [spaid]+truncate
                else:
                    swith = 0
                    return swith, []
        else:
            swith = 0
            return swith, []
    
    def filter(self):
        # 100% contain
        toDel, delSpa = set(), set()
        for eleSpa in self.spa.keys():
            for pos, soapPath in enumerate(self.allSoap):
                revPath = [self.reverseNodes(node) for node in soapPath[::-1]]
                spaPath = self.spa[eleSpa]
                inter = set(spaPath).intersection(set(soapPath))
                interRev = set(spaPath).intersection(set(revPath))
                if len(inter) == len(soapPath):
                    toDel.add(pos)
                elif len(interRev) == len(soapPath):
                    toDel.add(pos)
        
        for pos in sorted(toDel, reverse=True):
            self.allSoap.pop(pos)
    
    def walking(self):
        """
        Merge the 2 types of scaffolds, select a seed and extend it
        """
        # filter the whole contained scaffolds
        #pre = len(self.allSoap)
        #self.filter()
        #logger.info('{0} are deleted'.format(pre-len(self.allSoap)))
        
        allSides, ref, allNodes, soapScaf, backTrace = self.prepare()
        # filter the spa scaffolds that can not be extended
        paths, visited = [], set()
        for spaid in self.pairSpa:
            if spaid in visited:
                continue
            lswith = rswith = 1
            stack = [[spaid]]
            while stack:
                tmp = stack.pop()
                if rswith:
                    last = tmp[-1]
                    if last < 0:
                        rswith, enlongth = self.spaExtend(last, visited, allNodes, soapScaf, backTrace, 1)
                        if rswith:
                            tmp.extend(enlongth)
                            stack.append(tmp)
                    else:
                        for ele in tmp[::-1]:
                            if ele <0:
                                pos = tmp.index(ele)
                                break
                        nodes = tmp[pos+1:]
                        rswith, enlongth = self.soapExtend(nodes, visited, allSides, ref, backTrace, 1)
                        if rswith:
                            tmp = tmp[:pos+1]
                            tmp.extend(enlongth)
                            stack.append(tmp)
                if lswith:
                    first = tmp[0]
                    if first < 0:
                        lswith, enlongth = self.spaExtend(first, visited, allNodes, soapScaf, backTrace, -1)
                        if lswith:
                            enlongth.extend(tmp)
                            stack.append(enlongth)
                    else:
                        for ele in tmp:
                            if ele <0:
                                pos = tmp.index(ele)
                                break
                        nodes = tmp[:pos]
                        lswith, enlongth = self.soapExtend(nodes, visited, allSides, ref, backTrace, -1)
                        if lswith:
                            tmp = tmp[pos:]
                            enlongth.extend(tmp)
                            stack.append(enlongth)
                if rswith == lswith == 0:
                    paths.append(tmp)
                    break
        # add other soap scaffolds
        plus = set()
        for ele in soapScaf.keys():
            if ele in visited:
                del soapScaf[ele]
            else:
                if ele %2 != 0:
                    plus.add(ele)
        
        toDel = set()
        for eleSpa in self.spa:
            for eleSoap in plus:
                soapPath = soapScaf[eleSoap]
                revPath = [self.reverseNodes(node) for node in soapPath[::-1]]
                spaPath = self.spa[eleSpa]
                inter = set(spaPath).intersection(set(soapPath))
                interRev = set(spaPath).intersection(set(revPath))
                if len(inter)*1.0/len(soapPath)>0.5 or len(interRev)*1.0/len(soapPath)>0.5:
                    toDel.add(eleSoap)
        
        for ele in plus:
            if ele not in toDel:
                paths.append(soapScaf[ele])
        
        return paths

class overlapGraph:
    
    def __init__(self, index):
        
        self.__soap = {}
        self.__ingap = {}
        self.graph = {}
        self.index = index
        
    def allNodes(self):
        
        return set(self.__soap.keys()) | set(self.__ingap.keys()) 
    
    def indexGraph(self, soapPaths, allPaths):
        
        # index ingap scafs
        count = 1
        for path in allPaths:
            revPath = [reverseNode(i, self.index) for i in path[::-1]]
            self.__ingap[count] = path
            self.__ingap[count+1] = revPath
            count += 2
        # index soap scafs
        count = -1
        for path in soapPaths:
            revPath = [reverseNode(i, self.index) for i in path[::-1]]
            self.__soap[count] = path
            self.__soap[count-1] = revPath
            count -= 2
    
    def revStrandNode(self, node):
        """
        Get the reverse complement strand of the nodes of different types
        """
        if node > 0:
            if node % 2 == 0:
                return node - 1
            else:
                return node + 1
        if node < 0:
            if node % 2 == 0:
                return node + 1
            else:
                return node - 1
    
    def degree(self):
        """
        Output two dict: in and out degree
        """
        indegree = {}
        outdegree = {}
        for key in self.graph:
            if key not in outdegree:
                outdegree[key] = 0
            outdegree[key] += len(self.graph[key].keys())
            for sec in self.graph[key]:
                if sec not in indegree:
                    indegree[sec] = 0
                indegree[sec] += 1
        for ele in indegree:
            if ele not in outdegree:
                outdegree[ele] = 0
        for ele in outdegree:
            if ele not in indegree:
                indegree[ele] = 0
        return indegree, outdegree
    
    def checkSymmertic(self):
        """
        Check if two starnds are reverse complementry
        """
        pairs = []
        for key in self.graph:
            for sec in self.graph[key]:
                pairs.append([key, sec])
        for ele in pairs:
            rev = [self.revStrandNode(node) for node in ele[::-1]]
            if rev not in pairs:
                logger.error('The overlap graph is not symmetric, check it',
                             exit_with_code=1
                             )
    
    def __contains__(self, key):
        
        if not isinstance(key, int):
            logger.error('Please provide a int to represent the node on the graph',
                         exit_with_code=1
                         )
        if key in self.graph:
            return 1
        else:
            return 0
    
    def __getitem__(self, key):
        """
        get the keys of dict of dict, key must be provided
        """
        if not key:
            logger.error('Please provide a key to graph, like: graph[key]',
                         exit_with_code=1
                         )
        
        if key in self.graph:
            return self.graph[key]
        else:
            logger.error('Key ' + str(key) + ' is not on the graph',
                         exit_with_code=1
                         )
    
    def buildEdges(self, soapPaths, allPaths):
        """
        build the graph
        """
        
        self.indexGraph(soapPaths, allPaths)
        
        # get the overlap
        overlaps =[]
        covered = set()
        for key1 in self.__soap:
            for key2 in self.__ingap:
                inter = set(self.__soap[key1]) & set(self.__ingap[key2])
                if inter: # there are overlaps, check the types of overlap
                    posSoap = sorted([self.__soap[key1].index(i) for i in inter])
                    posIngap =  sorted([self.__ingap[key2].index(i) for i in inter])
                    # get the left and right region of flank regions
                    flankSoap = (posSoap[0], len(self.__soap[key1])-1-posSoap[-1])
                    flankIngap = (posIngap[0], len(self.__ingap[key2])-1-posIngap[-1])
                    if flankSoap[0] >= flankIngap[0] and flankSoap[1] >= flankIngap[1]:
                        covered.add(key2)
                    elif flankSoap[0] <= flankIngap[0] and flankSoap[1] <= flankIngap[1]:
                        covered.add(key1)
                    elif flankSoap[0] >= flankIngap[0] and flankSoap[1] < flankIngap[1]:
                        # soap -> ingap
                        overlaps.append(tuple([key1, key2, posSoap[-1], posIngap[-1], posIngap[-1]-posIngap[0]+1]))
                    elif flankSoap[0] < flankIngap[0] and flankSoap[1] >= flankIngap[1]:
                        # ingap -> soap
                        overlaps.append(tuple([key2, key1, posIngap[-1], posSoap[-1], posIngap[-1]-posIngap[0]+1]))
        # build a overlap graph
        graph = {}
        for ele in overlaps:
            if ele[0] not in graph:
                graph[ele[0]] = {}
            graph[ele[0]][ele[1]] = ele[2:]
        # filter the multiple outdegree nodes
        toDel = set()
        pairs = set()
        for key in graph:
            for sec in graph[key]:
                pairs.add(tuple([key, sec, graph[key][sec]]))
            if len(graph[key]) >= 2:
                linkOrder = sorted(graph[key], key=lambda x:graph[key][x][2])
                for ele in linkOrder[:-1]:
                    toDel.add(ele)
                    toDel.add(self.revStrandNode(ele))
           
        for ele in pairs:
            if not toDel.intersection(ele[:2]):
                self.graph[ele[0]] = {}
                self.graph[ele[0]][ele[1]] = ele[2]
        
        self.checkSymmertic()
        return covered, toDel
    
    def iterkeys(self):
        """
        Only get the nodes the outdegree larger than one, return a generator
        """
        return self.graph.iterkeys()
    
    def seq(self, node):
        if node > 0:
            return self.__ingap[node]
        else:
            return self.__soap[node]

class serialization:
    
    def __init__(self, soap, paras):
        
        self.soap = soap
        self.paras = paras
        self.dvertices = self.soap.dDlink.vertices()
        self.pvertices = self.soap.dPlink.vertices()
        self.distance = {}
        
    def getSeqs(self, allGenomeNodes):
        
        # soap seqs
        soapSeqs = self.soap.fetchSeqs(allGenomeNodes)
        """
        visited = set()
        with open(self.soap.transedContig, 'r') as af:
            for line in af:
                info = line.split()
                if int(info[0]) in allGenomeNodes:
                    seq = info[-1].rstrip()
                    visited.add(int(info[0]))
                    soapSeqs[int(info[0])] = seq
                    rev = self.soap.revNode(int(info[0]))
                    soapSeqs[rev] = revCom(seq)
                    visited.add(rev)
        left = allGenomeNodes - visited
        """
        
        left = allGenomeNodes - set(soapSeqs.keys())
        if left:
            logger.notice(' '.join(map(str, left)) + ' do not get their seqs')
            for ele in left:
                soapSeqs[ele] = 'N'*self.soap.kmersize
            
        # spa seqs
        spaSeqs = {}
        with open(self.paras.refn, 'r') as af:
            for line in af:
                if  line[0] == '>':
                    seqid = line[1:-1]
                else:
                    seq = line.rstrip()
                    spaSeqs[seqid] = seq
        
        for ele in spaSeqs.keys():
            if 'NODE' not in ele:
                del spaSeqs[ele]
        
        return soapSeqs, spaSeqs
    
    def linkType(self, first, second):
        
        if first in self.dvertices  and second in self.dvertices:
            try:
                linkinfo = self.soap.dDlink[first][second]
                return linkinfo
            except:
                pass
        if first in self.pvertices and second in self.pvertices:
            try:
                linkinfo = self.soap.dPlink[first][second]
                return linkinfo
            except:
                pass
        if first in self.distance and second in self.distance[first]:
            lap = self.distance[first][second]
            if lap < 0:
                return lap
            else:
                return [lap, 0, 0]
        
        return [20, 0, 0]
    
    def scafComp(self, walks, finalPaths):
        
        com, tmpSpaScafId = [], set()
        for path in finalPaths:
            tmp = []
            for pos in xrange(len(path)):
                if path[pos] > 0:
                    tmp.append(path[pos])
                else:
                    tmp.extend(walks.pairSpa[path[pos]])
                    tmpSpaScafId.add(walks.relate[path[pos]])
            com.append(tmp)
        
        scafNodes, num = set(), 0
        with open(config.COMPOSITION, 'w') as af:
            for path in com:
                scafNodes.update(path)
                af.write('scaf' + str(num) + ' ' + ' '.join(map(str, path)) + '\n')
                num += 1
        
        spaScafId = set()
        for ele in tmpSpaScafId:
            if ':' not in ele:
                spaScafId.add(ele)
            else:
                for key in ele.split(':'):
                    if 'NODE' in key:
                        spaScafId.add(key)
        return scafNodes, spaScafId
    
    def outSeqs(self, walks, gNodes, recoverSwith):
        
        # get the seqs of mergered spa scaffolds
        newlinks = self.paras.mappedLinks(self.soap.index, gNodes)
        aftSeqs = self.paras.newSpaSeqs(newlinks, self.soap)
        
        # get a mixed scaffolds
        finalPaths = walks.walking()
        scafNodes, sapScafId = self.scafComp(walks, finalPaths)
        scafNodes.update([self.soap.revNode(i) for i in scafNodes])
        
        # get seqs
        gNodes.update(scafNodes)
        soapSeqs, spaSeqs = self.getSeqs(gNodes)
        
        # scaffolding
        linkedScafs, num = {}, 1
        for path in finalPaths:
            if path[0] > 0:
                seq = soapSeqs[path[0]]
            else:
                spaid = walks.relate[path[0]]
                seq = aftSeqs[spaid] if path[0] %2 != 0 else revCom(aftSeqs[spaid])
            for pos in xrange(1, len(path)):
                fir, sec = path[pos-1], path[pos]
                if fir > 0 and sec > 0:
                    linkinfo = self.linkType(fir, sec)
                    if isinstance(linkinfo, int): # dlink
                        if linkinfo > 0:
                            seq += 'N'*linkinfo
                            seq += soapSeqs[sec]
                        else:
                            seq += soapSeqs[sec][-linkinfo:]
                    else: # pelink
                        nu = overlaps(soapSeqs[fir], soapSeqs[sec])
                        if nu:
                            seq += soapSeqs[sec][nu:]
                        else:
                            seq += 'N'*linkinfo[0]
                            seq += soapSeqs[sec]
                elif fir < 0 and sec < 0: # note the spa directions in the finalPaths
                    source, sink = walks.pairSpa[fir][-1], walks.pairSpa[sec][0]
                    direcRightSeq = aftSeqs[walks.relate[sec]] if sec%2 !=0 else revCom(aftSeqs[walks.relate[sec]])
                    linkinfo = self.linkType(source, sink)
                    if isinstance(linkinfo, int):
                        if linkinfo > 0:
                            seq += 'N'*linkinfo
                            #seq += direcRightSeq[self.soap.kmersize:]
                            seq += direcRightSeq
                        else:
                            seq += direcRightSeq[-linkinfo:]
                    else: # pelink
                        nu = overlaps(soapSeqs[source], soapSeqs[sink])
                        if nu:
                            seq += direcRightSeq[nu:]
                        else:
                            seq += 'N'*linkinfo[0]
                            seq += direcRightSeq
                elif fir > 0 and sec < 0:
                    sink = walks.pairSpa[sec][0]
                    direcRightSeq = aftSeqs[walks.relate[sec]] if sec%2 !=0 else revCom(aftSeqs[walks.relate[sec]])
                    linkinfo = self.linkType(fir, sink)
                    if isinstance(linkinfo, int):
                        if linkinfo > 0:
                            seq += 'N'*linkinfo
                            #seq += direcRightSeq[self.soap.kmersize:]
                            seq += direcRightSeq
                        else:
                            seq += direcRightSeq[-linkinfo:]
                    else:
                        nu = overlaps(soapSeqs[fir], soapSeqs[sink])
                        if nu:
                            seq += direcRightSeq[nu:]
                        else:
                            seq += 'N'*linkinfo[0]
                            seq += direcRightSeq
                elif fir < 0 and sec > 0:
                    source = walks.pairSpa[fir][-1]
                    linkinfo = self.linkType(source, sec)
                    if isinstance(linkinfo, int):
                        if linkinfo > 0:
                            seq += 'N'*linkinfo
                            #seq += soapSeqs[sec][self.soap.kmersize:]
                            seq += soapSeqs[sec]
                        else:
                            seq += soapSeqs[sec][-linkinfo:]
                    else:
                        nu = overlaps(soapSeqs[source], soapSeqs[sec])
                        if nu:
                            seq += soapSeqs[sec][nu:]
                        else:
                            seq += 'N'*linkinfo[0]
                            seq += soapSeqs[sec]
            # out seqs
            #na = 'scaf{0}'.format(num)
            linkedScafs[num] = {}
            linkedScafs[num]['seq'] = seq
            linkedScafs[num]['order'] = path
            num += 1
        
         # delete the scaffolds that not belongs to the target genome
        mappedNodes = set([k for i in self.paras.nodes for k in i, self.soap.revNode(i)])
        if recoverSwith:
            errorNodes = correct(linkedScafs, mappedNodes, self.soap.length)
            errorNodes = set([k for i in errorNodes if i > 0 for k in (i, self.soap.revNode(i))])
        else:
            errorNodes = set()
        
        # generate singletons
        gNodes = gNodes - errorNodes
        spaSingle = {i:spaSeqs[i] for i in self.paras.spa - sapScafId}
        soapSinId = set([ i for i in gNodes-scafNodes if self.soap.index[i]==0])
        
        newSeqs = self.soap.fetchSeqs([i for i in soapSinId | self.paras.nodes if i not in soapSeqs])
        # generate seeds
        with open(config.SEED, 'w') as af:
            for ele in self.paras.nodes:
                outid = '>{}'.format(ele)
                try:
                    af.write(outid+'\n'+soapSeqs[ele]+'\n')
                except:
                    af.write(outid+'\n'+newSeqs[ele]+'\n')
        
        # genearte recovered
        gNodes = set([i for i in gNodes if self.soap.index[i]==0])
        with open(config.RECOVERED,'w') as af:
            for ele in gNodes:
                outid = '>{}'.format(ele)
                try:
                    af.write(outid+'\n'+soapSeqs[ele]+'\n')
                except:
                     af.write(outid+'\n'+newSeqs[ele]+'\n')
        
        # generate final results
        with open(config.FINALRES, 'w') as af:
            for ele in linkedScafs:
                outid = '>scaf{0}'.format(ele)
                af.write(outid + '\n' + linkedScafs[ele]['seq']+'\n')
            for ele in spaSingle:
                outid = '>{0}'.format(ele)
                af.write(outid + '\n' + spaSingle[ele]+'\n')
            for ele in soapSinId:
                outid = '>{0}'.format(ele)
                #if len(soapSingle[ele]) > 100:
                try:
                    af.write(outid + '\n' + soapSeqs[ele]+'\n')
                except:
                    af.write(outid+'\n'+newSeqs[ele]+'\n')
        
        #output MDA seeds
        flag = 0
        mdascafs =  open(config.MDASEEDS, 'w')
        with open(self.paras.refn, 'r') as af:
            for line in af:
                if '>' in line:
                    if 'NODE' in line:
                        mdascafs.write(line)
                        flag = 1
                    else :
                        flag = 0
                else:
                    if flag == 1:
                        mdascafs.write(line)
        mdascafs.close()
        
        logger.info('Outputing done: ')
        logger.info('The seed contigs are writed into file: {0}'.format(config.SEED))
        logger.info('The recovered contigs are writed into file: {0}'.format(config.RECOVERED))
        logger.info('The scaffolds are writed into file: {0}'.format(config.FINALRES))
        logger.info('done !')

class correlation:
    
    def __init__(self, x, y, size):
        
        self.x = x
        self.y = y
        self.size = size
        self.combs = self.combination()
        
    def combination(self):
        """
        Generate 256 types of permutation of tertra nuclotides
        """
        combination = []
        S = ['A','T','G','C']
        for seq in itertools.product(S, repeat = self.size):
                combination.append(''.join([i for i in seq]))
        combination.sort()
        return combination
    
    def standarFreq(self, freq):
        """
        Normalize the frequency of the frequency table
        """
        sumOccurences = sum(freq.values())
        coeff = math.pow(4, self.size)/sumOccurences
        for key in freq:
            normalized = freq[key]*coeff
            freq[key] = normalized
    
    def revCom(self, string):
        """
        reverse complement the given string
        """
        comp = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C', 'N':'N'}
        result = ''
        for i in string[::-1]:
            result += comp[i]
        return result
    
    def zScore(self, seq):
        """
        Caculate the z-scores of the sequence, scaf length cutoff > 5000
        both of the two strand should be caculated
        """
        # get the frequency of the tera, tri, di- oligo of the seq
        revSeq = self.revCom(seq)
        # tertra or tri
        arr1 = ( seq[i:i+self.size].upper() for i in xrange(len(seq)-self.size+1) )
        arr2 = ( revSeq[i:i+self.size].upper() for i in xrange(len(seq)-self.size+1) )
        oligoF = collections.defaultdict(int)
        oligoS = collections.defaultdict(int)
        oligoT = collections.defaultdict(int)
        for i in arr1: oligoF[i] += 1
        for i in arr2: oligoF[i] += 1
        # remove the 'N' contained kmers
        for key in itertools.ifilter(lambda x: set(x) - set('ATCG'), oligoF.keys()):
            del oligoF[key]
        # normalize the frequency
        self.standarFreq(oligoF)
        for key, value in oligoF.iteritems():
            oligoS[key[:self.size-1]] += value
            oligoT[key[:self.size-2]] += value
        # caculate the z-score: Teeling et al. BMC Bioinformatics 2004; 5:163.
        # Zxyzw = (Nxyzw - Exyzw)/Var(Nxyzw)**0.5
        # Exyzw = Nxyz*Nyzw/Nyz
        # Var(Nxyzw) = Exyzw*(Nyz-Nxyz)(Nyz-Nyzw)/Nyz**2
        zScores = []
        for XYZW in self.combs:
            try:
                XYZ, YZW, YZ = XYZW[:self.size-1], XYZW[1:], XYZW[1:self.size-1]
                exp = oligoS[XYZ]*oligoS[YZW]/oligoT[YZ]
                var = exp*(oligoT[YZ]-oligoS[XYZ])*(oligoT[YZ]-oligoS[YZW])/math.pow(oligoT[YZ], 2)
                z = (oligoF[XYZW] - exp)/math.sqrt(var)
                zScores.append(z)
            except:
                zScores.append(0)
        
        return zScores
    
    def pearsonCorr(self):
        """
        Use the pearson correlation coefficent when the seq length is long
        """
        zX = self.zScore(self.x)
        zY = self.zScore(self.y)
        vectLen = len(zX)
        sumX, sumXsq, sumY, sumYsq = 0, 0, 0, 0
        for ele in zX:
            sumX += ele
            sumXsq += math.pow(ele,2)
        for ele in zY:
            sumY += float(ele)
            sumYsq += pow(float(ele),2)
        # add up the products
        product = sum([zX[pos]*zY[pos] for pos in xrange(vectLen)])
        # caculate the pearson correlation score
        numerator = product - (sumX*sumY/vectLen)
        denominator = math.sqrt((sumXsq - pow(sumX, 2)/vectLen)*(sumYsq - pow(sumY, 2)/vectLen))
        # can not have division by 0
        if denominator == 0: return 0
        cor = numerator/denominator
        return cor

def correct(scaffolds, mappedNodes, length):
    """
    Use a reliable reference scaf as seed to binning the other scaffolds
    """
    print('Caculating distance among scaffolds ...')
    iden = refSeq = nu =  0
    for key in scaffolds:
        order =  scaffolds[key]['order']
        tmpNu = sum([length[i] for i in order if i in mappedNodes])
        if tmpNu > nu:
            nu = tmpNu
            refSeq = scaffolds[key]['seq']
            iden = key
    # 0.75 can be a cutoff of scafffolds that longer than 100k
    distance = {}
    for key in scaffolds:
        if set(scaffolds[key]['order']) & mappedNodes:
            continue
        seq = scaffolds[key]['seq']
        if len(seq) >= 15000:
            cor = correlation(refSeq, seq, 4)
            distance[key] = (len(seq), cor.pearsonCorr())
    
    # discard wrong scaffolds
    record, delKey = set(), []
    for key in distance:
        if distance[key][0] < 100000 and distance[key][1] <= 0.7:
            delKey.append(key)
        if distance[key][0] >= 100000 and distance[key][1] <= 0.7:
            delKey.append(key)
    
    logger.debug('{0} scaffolds are deleted.'.format(len(delKey)))
    record = set([k for i in delKey for k in scaffolds[i]['order']])
    
    for key in delKey:
        del scaffolds[key]
    
    return record

def delCands(soap, paths, mappedNodes, suspNodes, lc=10000, ic=0.7):
    """
    use corelation to delete candidates if the nodes if large
    """
    # only use seqs that larger than 15k
    suNodes = set([i for i in suspNodes if soap.index[i]==0 and soap.length[i]>=lc])
    if not suNodes:
        return set()
    
    nu = 0
    for path in paths:
        if set(path) & mappedNodes:
            tmpNu = sum([soap.length[i] for i in path])
            if tmpNu > nu:
                nu = tmpNu
                order = path
    
    # 0.75 can be a cutoff of scafffolds that longer than 100k
    seqs = soap.fetchSeqs(set(order) | suNodes)
    
    refSeq = seqs[order[0]]
    for i in xrange(1, len(order[1:])):
        try:
            overlap = soap.dDlink[order[i-1]][order[i]]
        except:
            overlap = soap.dPlink[order[i-1]][order[i]][0]
        if overlap < 0:
            refSeq += seqs[order[i]][abs(overlap):]
        else:
            refSeq += 'N'*overlap
            refSeq += seqs[order[i]]
    
    delNodes = set()
    for node in suNodes:
        cor = correlation(refSeq, seqs[node], 4)
        if cor.pearsonCorr() <= ic:
            delNodes.add(node)
    delNodes.update([i+1 for i in delNodes if soap.index[i]==0])
    logger.debug('{0} candidates should be deleted'.format(len(delNodes)))
    
    return delNodes

def plotSubgraph(graph, elements):
    """
    plot the subgraph for igraph visualization
    """
    if isinstance(elements, list):
        elements = set(elements)
    elif isinstance(elements, int):
        logger.error('Subgraph should contain more than one node, check it',
                     exit_with_code=1
                     )
    
    with open('subgraphPlot', 'w') as af:
        for i in graph:
            for k in graph[i]:
                edge = [i, k]
                if elements.intersection(edge):
                    af.write(' '.join(map(str, edge)) + '\n')

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

def revCom(string):
    """
    reverse complement the given string
    """
    comp = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C', 'N' : 'N',
            'a':'T', 't':'A', 'c':'G', 'g':'C', 'n':'N'
            }
    
    result = ''
    for i in string[::-1]:
        result += comp[i]
    
    return result

def comlinks(graph, plink, pool, mergedNodes):
    # find the pelinks among different nodes (combinedNodes) in one single elements (group)
    transferMerge = {}
    for ele in mergedNodes:
        for i in mergedNodes[ele]:
            transferMerge[i] = ele
    
    links = {}
    for ele in pool:
        arEle = set()
        if ele in mergedNodes:
            acEle = set(mergedNodes[ele])
            for i in acEle:
                arEle.update(graph.plinkNodes(i, plink))
        else:
            acEle = set([ele])
            arEle =  graph.plinkNodes(ele, plink)
        arEle = arEle - acEle
        if arEle:
            for sec in arEle:
                if sec in pool:
                    if ele not in links:
                        links[ele] = {}
                    links[ele][sec] = 1
                else:
                    if sec in transferMerge:
                        if ele not in links:
                            links[ele] = {}
                        links[ele][transferMerge[sec]] = 1
    
    return links

def bfs_search(start, graph, step = 5):
    
    queue = collections.deque([[start]])
    searched = []
    while queue:
        tmp = queue.pop()
        last = tmp[-1]
        if last in graph:
            for second in graph[last]:
                longerPath = tmp + [second]
                if len(longerPath) == step + 1:
                    searched.append(longerPath)
                else:
                    queue.appendleft(longerPath)
        else:
            searched.append(tmp)
    
    return searched   # return a structer link [[1,2,4],[3,5,6],[2,3,4]]

def indegree(node, graph, index):
    
    rev = reverseNode(node, index)
    if rev not in graph:
        return 0
    else:
        return len(graph[rev])

def backTrace(graph, index, links, chain, visitedNodes):
    """
    back trace the chain find branches
    """
    splitPaths = set()
    if len(chain) == 1:
        rev = reverseNode(chain[0], index)
        visitedNodes.update([chain[0], rev])
        if rev in graph:
            nextNodes = set(graph[rev].keys())
            return set([(chain[0],)]), nextNodes
        else:
            return set([(chain[0],)]), set()
    
    frag = []
    restart = set()
    revChain = [reverseNode(i, index) for i in chain]
    while chain:
        node = chain.pop()
        revNode = reverseNode(node, index)
        frag += [revNode]
        if revNode in graph and len(graph[revNode]) > 1:
            nextNodes = graph[revNode].keys()
            nextLink = set([k for i in frag if i in links for k in links[i]])
            if not nextLink: # can not guild the traversal
                splitPaths.add(tuple([reverseNode(i, index) for i in frag[::-1]]))
                frag = []
                for ele in nextNodes:
                    if ele not in revChain:
                        restart.add(ele)
            else:
                extPaths = bfs_search(revNode, graph, 2)
                linked = []
                for path in extPaths:
                    if nextLink.intersection(path[1:]):
                        linked.append(path[1])
                if len(linked) == 1 and linked[0] not in revChain:
                    splitPaths.add(tuple([reverseNode(i, index) for i in frag[::-1]]))
                    frag = []
                    for ele in nextNodes:
                        if ele not in revChain:
                            restart.add(ele)
                elif len(linked) == 1 and linked[0] in revChain:
                    for ele in nextNodes:
                        if ele not in revChain:
                            restart.add(ele)
                else: # len(linked) == 0 or len(linked) > 1
                    splitPaths.add(tuple([reverseNode(i, index) for i in frag[::-1]]))
                    frag = []
                    for ele in nextNodes:
                        if ele not in revChain:
                            restart.add(ele)
        if len(chain) == 0: # get the last node
            splitPaths.add(tuple([reverseNode(i, index) for i in frag[::-1]]))
            break
    
    doneNodes = set()
    for path in splitPaths:
        doneNodes.update(path)
        doneNodes.update([reverseNode(i,index) for i in path])
    visitedNodes.update(doneNodes)
    
    return splitPaths, restart

def plinkGuildPaths(oriGraph, comGraph, index, pelink, combinedNodes, elements):
    """
    Use pelink to supervise the traversal
    """
    links = comlinks(oriGraph, pelink, elements, combinedNodes)
    paths = set()
    startNodes = collections.deque([i for i in elements if indegree(i, comGraph, index) == 0])
    if not startNodes:
        #  lot of node that point to each others
        loop = set()
        record = set()
        for node in comGraph:
            if node not in record:
                record.add(node)
                for key in comGraph[node]:
                    if key in comGraph and node in comGraph[key]:
                        loop.update([node, key])
        startNodes = collections.deque(list(set(elements) & loop))
    # may be loop
    if not startNodes:
        startNodes = collections.deque([elements[0]])
    # go to tarversal
    visitedNodes = set()
    while startNodes:
        node = startNodes.pop()
        if node not in visitedNodes:
            stack = [[node]]
            while stack:
                tmp = stack.pop()
                last = tmp[-1]
                if last in visitedNodes or last in tmp[:-1]:
                    splitPaths, restart = backTrace(comGraph, index, links, tmp[:-1], visitedNodes)
                    paths.update(splitPaths)
                    for i in restart: startNodes.appendleft(i)
                    continue
                visitedNodes.update([last, reverseNode(last, index)])
                if last not in comGraph: # the end of chain
                    splitPaths, restart = backTrace(comGraph, index, links, tmp, visitedNodes)
                    paths.update(splitPaths)
                    for i in restart: startNodes.appendleft(i)
                else:
                    fullTmpNodes = tmp + [reverseNode(i, index) for i in tmp]
                    nextLink = set([k for i in tmp if i in links for k in links[i] if k not in fullTmpNodes])
                    nextNodes = comGraph[last].keys()
                    if len(nextNodes)==1: # may be weak connection (no pelink)
                        longer = tmp + nextNodes
                        stack.append(longer)
                    else: #len(nextNodes) > 1
                        extPaths = bfs_search(last, comGraph, 2)
                        linked = []
                        for path in extPaths:
                            if nextLink.intersection(path[1:]):
                                linked.append(path[1])
                        if len(linked) == 1:
                            longer = tmp + [linked[0]]
                            for ele in nextNodes:
                                if ele not in visitedNodes:
                                    startNodes.appendleft(ele)
                            stack.append(longer)
                        else: # end of the chain, reverse tarversal
                            splitPaths, restart = backTrace(comGraph, index, links, tmp, visitedNodes)
                            paths.update(splitPaths)
                            for i in restart: startNodes.appendleft(i)
                            for ele in nextNodes:
                                if ele not in visitedNodes:
                                    startNodes.appendleft(ele)
    
    if len(visitedNodes) == len(elements) :
        return paths
    else:
        visitedNodes.update([reverseNode(i, index) for i in visitedNodes])
        logger.error('Some of the nodes are not visited, please re-check your code',
                     exit_with_code=1
                     )

def deleteLoops(comGraph, index, combinedNodes):
    """
    delete the edeges that two nodes point to each other
    """
    # delete the loop from the group, because it will interfere the travesaling
    loop = set()
    record = set()
    for node in comGraph:
        if node not in record:
            record.add(node)
            for key in comGraph[node]:
                if key in comGraph and node in comGraph[key]:
                    loop.update([node, key])
    # delete the loop edges
    oldNodes = set()
    edges = set()
    for source in comGraph:
        for sink in comGraph[source]:
            edges.add(tuple([source, sink]))
            oldNodes.update([source, sink])
    # get a new graph
    newGraph = {}
    newNodes = set()
    for edge in edges:
        if not loop.intersection(set(edge)):
            if edge[0] not in newGraph:
                newGraph[edge[0]] = {}
            newGraph[edge[0]][edge[1]] = comGraph[edge[0]][edge[1]]
            newNodes.update(edge)
    
    deletedNodes = set([i for i in (oldNodes-newNodes) if index[i]==0])
    singStrand = []
    for ele in deletedNodes:
        if ele in combinedNodes:
            singStrand.append(combinedNodes[ele])
        else:
            singStrand.append([ele])
    
    return newGraph, singStrand

def dPaths(oriGraph, soap, fn=False):
    """
    Out put the final paths, use PElink to supervise the traversal
    """
    # handel scaffolds
    chains = []
    # first compress the graph, formate of simpleChain: [[1,2,3],[4,5,6] ...]
    comGraph, combinedNodes, simpleChain = oriGraph.compressGraph(linear=1)
    chains.extend(simpleChain)
    # handel the complex subgraph
    newGraph, loops = deleteLoops(comGraph, soap.index, combinedNodes)
    if loops: chains.extend(loops)
    # use pelink to supervise the traveral
    group, pairs = subGraph(newGraph, soap.index, out='p')
    graphPaths = set()
    for pair in pairs:
        elements = list(group[pair[0]])
        if len(pair) == 2:
            elements.extend(list(group[pair[1]]))
        paths = plinkGuildPaths(oriGraph, newGraph, soap.index, soap.dPlink, combinedNodes, elements)
        graphPaths.update(paths)
    
    if not graphPaths: logger.info('NOTE: no nodes can be linked by the graph')
    # output the path to chain
    for fragment in graphPaths:
        chains.append([])
        for ele in fragment:
            if ele in combinedNodes:
                chains[-1].extend(combinedNodes[ele])
            else:
                chains[-1].append(ele)
    #output the super contig compositions
    dlinkSingletons = set()
    dpath = []
    for path in chains:
        if len(path) == 1:
            dlinkSingletons.update([path[0], soap.revNode(path[0])])
        elif len(path) > 1:
            dpath.append(path)          
    if fn:
        with open(fn,'wb') as af:
            pickle.dump(dpath, af)
    
    return dlinkSingletons, dpath

def pPaths(oriGraph, index, fn=False):
    """
    output the pelink paths, for there is no addtional info to guild the
    cross site on the graph break when encounter a cross site
    """
    paths = []
    plinkSingltons = set()
    comGraph, combinedNodes, simpleChain = oriGraph.compressGraph(linear=1)
    paths.extend(simpleChain)
    record = set()
    for ele in combinedNodes:
        if ele not in record:
            revEle = reverseNode(ele, index)
            record.update([ele, revEle])
            paths.append(combinedNodes[ele])
    # output the composition of each scaffold
    pathNodes = set()
    for path in paths:
        revNodes = set([reverseNode(i, index) for i in path])
        pathNodes.update(path)
        pathNodes.update(revNodes)
    if fn:
        with open(fn, 'wb') as af:
            pickle.dump(paths, af)
    # output the singlets
    plinkSingltons.update([i for i in (oriGraph.vertices() - pathNodes) if index[i]==0])
    plinkSingltons.update([i+1 for i in plinkSingltons])
    return plinkSingltons, paths

def overlaps(x, y):
    """
    get the end maches between two seqs
    """
    if not isinstance(x, str) or not isinstance(y, str):
        logger.error('The overlaps function only accept two string objects',
                     exit_with_code=1
                     )
    
    maxlimit = 100
    minlimit = 5
    for m in xrange(maxlimit, minlimit, -1):
        if x[-m:] == y[:m]:
            return m
    else:
        return 0
 
def reviseScaf(scafs, index, mappedNodes, *graphs):
    """
    Revise and extend the scafs: first, error correction; second, extend the scaffold
    """
    genomeNodes = set()
    genomeNodes.update(mappedNodes)
    dgraph, pgraph = graphs
    if dgraph:
        genomeNodes.update(dgraph.vertices())
    if pgraph:
        genomeNodes.update(pgraph.vertices())
    # choose the scaf that belongs to the genome
    recovered = [] # a feasible links, just check relocation
    for ele in scafs:
        inter = genomeNodes.intersection(ele)
        if inter:
            rate = len(inter)*1.0/len(ele)
            if 0.4 <= rate:
                recovered.append(ele)
    
    if dgraph:
        dNodes= dgraph.vertices()
        step2Scafs = []
        for path in recovered:
            step2Scafs.append([])
            for pos in xrange(len(path)-1):
                if path[pos] in dNodes and path[pos+1] in dNodes:
                    connects = dgraph.allPaths(path[pos], path[pos+1], step=100)
                    if connects and len(connects)==1:
                        step2Scafs[-1].extend(connects[0][:-1])
                    else:
                        step2Scafs[-1].append(path[pos])
                else:
                    step2Scafs[-1].append(path[pos])
            else:
                step2Scafs[-1].append(path[-1])
        return step2Scafs
    else:
        return recovered

def extendScaf(scafs, soap, *roads):
    """
    Extend the scafs by the LinkedContigs
    """
    # first build a graph use dgraph and pgraph
    allPaths = []
    if len(roads) == 1:
        dpath = roads[0]
        allPaths.extend(dpath)
    elif len(roads) == 2:
        dpath, ppath = roads
        allPaths.extend(dpath)
        allPaths.extend(ppath)
    # build a overlap graph
    graph = overlapGraph(soap.index)
    covered, deleted = graph.buildEdges(scafs, allPaths)
    # traverl the graph and output the connected nodes
    paths = []
    visitedNodes = set()
    indegree, outdegree = graph.degree()
    startNodes = [ key for key in graph.iterkeys() if indegree[key]==0 and outdegree[key] > 0]
    for node in startNodes:
        if node not in visitedNodes:
            stack = [[node]]
            while stack:
                tmp = stack.pop()
                last = tmp[-1]
                if last in visitedNodes:
                    paths.append(tmp)
                    continue
                visitedNodes.add(last)
                if last > 0:
                    revNode = last-1 if last%2==0 else last+1
                else:
                    revNode = last+1 if last%2==0 else last-1
                visitedNodes.update([last, revNode])
                if last in graph:
                    for sec in graph[last]:
                        longer = tmp + [sec]
                        stack.append(longer)
                else:
                    paths.append(tmp)
    # check two end of scaffolds
    newPaths = []
    for path in paths:
        path = collections.deque(path)
        if len(path) == 2:
            shares = graph[path[0]][path[1]]
            seq1 = graph.seq(path[0])
            seq2 = graph.seq(path[1])
            if shares[0] < len(seq1) - 1:
                l1 = sum([soap.length[i] for i in seq1[shares[0]:]])
                l2 = sum([soap.length[i] for i in seq2[shares[1]:]])
                if l1 > l2:
                    delNode = path.pop()
                    deleted.update([delNode, graph.revStrandNode(delNode)])
                    newPaths.append(path)
                    continue
            if shares[1] != 0:
                l1 = sum([soap.length[i] for i in seq1[:shares[0]]])
                l2 = sum([soap.length[i] for i in seq2[:shares[1]]])
                if l1 < l2 :
                    delNode = path.popleft()
                    deleted.update([delNode, graph.revStrandNode(delNode)])
            newPaths.append(path)
        else: # len path > 2
            shares = graph[path[-2]][path[-1]]
            seq1 = graph.seq(path[-2])
            seq2 = graph.seq(path[-1])
            if shares[0] < len(seq1) - 1:
                l1 = sum([soap.length[i] for i in seq1[shares[0]:]])
                l2 = sum([soap.length[i] for i in seq2[shares[1]:]])
                if l1 > l2:
                    delNode = path.pop()
                    deleted.update([delNode, graph.revStrandNode(delNode)])
            shares = graph[path[0]][path[1]]
            seq1 = graph.seq(path[0])
            seq2 = graph.seq(path[1])
            if shares[1] != 0:
                l1 = sum([soap.length[i] for i in seq1[:shares[0]]])
                l2 = sum([soap.length[i] for i in seq2[:shares[1]]])
                if l1 < l2 :
                    delNode = path.popleft()
                    deleted.update([delNode, graph.revStrandNode(delNode)])
            newPaths.append(path)
    # generate final scaffolds
    finalScafs = {}
    count = 1
    for path in newPaths:
        if len(path) == 1:
            seq = graph.seq(path[0])
        else:
            truncate = []
            seq = []
            for pos in xrange(len(path)-1):
                infoUp = graph[path[pos]][path[pos+1]]
                if pos == 0: # compare which scaf is longer
                    truncate.append([0, infoUp[0]])
                else:
                    truncate[-1].append(infoUp[0])
                truncate.append([infoUp[1]])
            for pos, value in enumerate(truncate):
                if len(value) == 2:
                    node = path[pos]
                    seq.extend(graph.seq(node)[value[0]:value[1]])
                else:
                    node = path[pos]
                    seq.extend(graph.seq(node)[value[0]:])
        finalScafs[count] = seq
        count += 1
    # get the unlinked scaffold of ingap
    left = graph.allNodes() - (covered | deleted)
    scafNodes = set()
    for path in newPaths:
        scafNodes.update(path)
    scafNodes.update([graph.revStrandNode(node) for node in scafNodes])
    left = left - scafNodes    
    left = set([i for i in left if i%2==1])
    for node in left:
        seq = graph.seq(node)
        nu = len([node for node in seq if soap.length[node] > 2*soap.kmersize])
        if nu > 0:
            finalScafs[count] = seq
            count += 1
    
    return finalScafs

def report():
    """
    Report the assembly result and comparsion with previous 40% nodes 
    """
    # gather the info of every file
    # seqNu, TotalLen, maxLen, N50, averLen, GCrate
    preInfo = contigInfo(config.SEED)
    aftInfo = contigInfo(config.RECOVERED)
    if os.path.exists(config.MDASEEDS) and os.path.getsize(config.MDASEEDS):
        mdaInfo = contigInfo(config.MDASEEDS)
    else:
        mdaInfo = [0 for i in xrange(10)]
    # (seqNu, TotalLen, maxLen, NLenInfo[4][0], averLen, GCrate),\ NLenInfo,\ (Arate, Trate, Crate, Grate, Nrate),\
    #(scaflg100, ratelg100, scaflg500, ratelg500, scaflg1k, ratelg1k, scaflg100k, ratelg100k, scaflg1m, ratelg1m)
    normalInfo, NLenInfo, baseRate, lgRate =  scaffoldInfo(config.FINALRES)
    # make the final report
    try:
        logfile = open(config.STATISTIC,'w')
        oldStdout = sys.stdout
        sys.stdout = logfile
        
        # first output the recovered contig info
        if os.path.exists(config.MDASEEDS):
            print('<-- Information of sequecnes before assembly -->\n')
            print("%-25s %-15s %-15s" % ('sequences ', 'from_metaO', 'from_metaS'))
            print("%-25s %-15d %-15d" % ('Sequece Number:', preInfo[0], mdaInfo[0]))
            print("%-25s %-15d %-15d" % ('Total Length:', preInfo[1], mdaInfo[1]))
            print("%-25s %-15d %-15d" % ('Longest Length:', preInfo[2], mdaInfo[2]))
            print("%-25s %-15d %-15d" % ('N50 Length:', preInfo[3], mdaInfo[3]))
            print("%-25s %-15d %-15d" % ('Average Length:', preInfo[4], mdaInfo[4]))
            print("%-25s %.2f%% %-8s %.2f%%\n" % ('GC Content:', preInfo[5], ' ', mdaInfo[5]))
        else:
            print('<-- Information of sequecnes before assembly -->\n')
            print("%-25s %-15d" % ('Contig Number:', preInfo[0]))
            print("%-25s %-15d" % ('Total Length:', preInfo[1]))
            print("%-25s %-15d" % ('Longest Length:', preInfo[2]))
            print("%-25s %-15d" % ('N50 Length:', preInfo[3]))
            print("%-25s %-15d" % ('Average Length:', preInfo[4]))
            print("%-25s %.2f%%\n" % ('GC Content:', preInfo[5]))
        
        # then output the comparions of before and after
        print('<-- MGA assembly information -->\n')
        print("%-25s %-15s %-14s" % ('sequences ', 'Contig', 'Scaffold'))
        print("%-25s %-15d %-15d" % ('Number:', aftInfo[0], normalInfo[0]))
        print("%-25s %-15d %-15d" % ('Total Length:', aftInfo[1], normalInfo[1]))
        print("%-25s %-15d %-15d" % ('Longest Length:', aftInfo[2], normalInfo[2]))
        print("%-25s %-15d %-15d" % ('N50 Length:', aftInfo[3], normalInfo[3]))
        print("%-25s %-15d %-15d" % ('Average Length:', aftInfo[4], normalInfo[4]))
        print("%-25s %.2f%% %-8s %.2f%%" % ('GC Content:', aftInfo[5], ' ', normalInfo[5]))
        print("At last, addtional %d bases are recovered.\n" % (normalInfo[1] - preInfo[1]))
        
        # print the info of scaffolds   
        print('<-- Information for assembly Scaffold (cut_off_length < 100bp) -->')
        print("%-25s %-15d %.2f%%" % ('scaffolds>100', lgRate[0], lgRate[1]))
        print("%-25s %-15d %.2f%%" % ('scaffolds>500', lgRate[2], lgRate[3]))
        print("%-25s %-15d %.2f%%" % ('scaffolds>1K', lgRate[4], lgRate[5]))
        print("%-25s %-15d %.2f%%" % ('scaffolds>10K', lgRate[6], lgRate[7]))
        print("%-25s %-15d %.2f%%" % ('scaffolds>100K', lgRate[8], lgRate[9]))
        print("%-25s %-15d %.2f%%\n" % ('scaffolds>1M', lgRate[10], lgRate[11]))
        print("%-25s %-15d %.1f%%" % ('Nucleotide_A', baseRate[0], baseRate[1]))
        print("%-25s %-15d %.1f%%" % ('Nucleotide_T', baseRate[2], baseRate[3]))
        print("%-25s %-15d %.1f%%" % ('Nucleotide_C', baseRate[4], baseRate[5]))
        print("%-25s %-15d %.1f%%" % ('Nucleotide_G', baseRate[6], baseRate[7]))
        print("%-25s %-15d %.1f%%\n" % ('sNucleotide_N', baseRate[8], baseRate[9]))
        print("%-25s %-15d %-15d" % ('N10', NLenInfo[0][0], NLenInfo[0][1]))
        print("%-25s %-15d %-15d" % ('N20', NLenInfo[1][0], NLenInfo[1][1]))
        print("%-25s %-15d %-15d" % ('N30', NLenInfo[2][0], NLenInfo[2][1]))
        print("%-25s %-15d %-15d" % ('N40', NLenInfo[3][0], NLenInfo[3][1]))
        print("%-25s %-15d %-15d" % ('N50', NLenInfo[4][0], NLenInfo[4][1]))
        print("%-25s %-15d %-15d" % ('N60', NLenInfo[5][0], NLenInfo[5][1]))
        print("%-25s %-15d %-15d" % ('N70', NLenInfo[6][0], NLenInfo[6][1]))
        print("%-25s %-15d %-15d" % ('N80', NLenInfo[7][0], NLenInfo[7][1]))
        print("%-25s %-15d %-15d" % ('N90', NLenInfo[8][0], NLenInfo[8][1]))
    finally:
        if logfile:
            logfile.close()
        if oldStdout:
            sys.stdout = oldStdout
    
    logger.info('Statistic info has been writed into file: {}'.format(config.STATISTIC))
    logger.info('\nDone')

# for no MDA scaffolds scaffolding
def generateSeqs(scafs, soap, targetGraph, mappedNodes, distance):
    """
    Gernerate the seqs of scafs to binning
    """
    scafNodes = set([v for i in scafs.itervalues() for k in i for v in (k, soap.revNode(k))])
    
    # prepare data for scaffolding
    allNodes = targetGraph.vertices() | mappedNodes | scafNodes 
    
    # get seqs
    graphSeqs = soap.fetchSeqs(allNodes)
    
    # for discarding error scaffolds
    scaffolds = {}
    for key, value in scafs.iteritems():
        scaffolds[key] = {}
        scaffolds[key]['order'] = value
        # initialize the sequence
        seq = graphSeqs[value[0]]
        for pos in xrange(len(value)-1):
            source, sink = value[pos], value[pos+1]
            if source in targetGraph and sink in targetGraph[source]:
                lap = targetGraph[source][sink]
                if isinstance(lap, int):
                    if lap < 0:
                        seq += graphSeqs[sink][-lap:]
                    else:
                        seq += 'N'*lap + graphSeqs[sink]
                else:# list
                    lap = lap[0]
                    if lap < 0:
                        seq += graphSeqs[sink][-lap:]
                    else:
                        seq += 'N'*lap + graphSeqs[sink]
            elif source in distance and sink in distance[source]:
                lap = distance[source][sink]
                if lap < 0:
                    seq += graphSeqs[sink][-lap:]
                else:
                    seq += 'N'*lap + graphSeqs[sink]
            elif source in soap.dPlink and sink in soap.dPlink[source]:
                nu = overlaps(graphSeqs[source], graphSeqs[sink])
                if nu:
                    seq += graphSeqs[sink][nu:]
                else:
                    linkinfo = soap.dPlink.linkInfo(source, sink)
                    if isinstance(linkinfo, int):
                        seq += 'N'*10 + graphSeqs[sink]
                    elif isinstance(linkinfo, list):
                        seq += 'N'*linkinfo[0]
                        seq += graphSeqs[sink]
            else:
                seq += 'N'*10
        scaffolds[key]['seq'] = seq
    
    # delete the scaffolds that not belongs to the target genome
    #errorNodes = correct(scaffolds, mappedNodes, soap.length)
    #errorNodes = set([k for i in errorNodes if i > 0 for k in (i, soap.revNode(i))])
    # refresh allGenomeNodes
    #allNodes = allNodes - errorNodes
    
    # get error nodes
    newScafNodes = set([ele for key in scaffolds for ele in scaffolds[key]['order']])
    newScafNodes.update([soap.revNode(i) for i in newScafNodes])
    
    # get singletons
    soapSinId = set([ i for i in (allNodes-newScafNodes) if soap.index[i]==0])
    
    # output sequences of every part
    logger.info('Outputing sequneces ...')
    
    with open(config.SEED, 'w') as af:
        for ele in (i for i in mappedNodes if soap.index[i] == 0):
            outid = '>{}'.format(ele)
            af.write(outid+'\n'+graphSeqs[ele]+'\n')
    
    with open(config.RECOVERED, 'w') as af:
        for ele in (i for i in allNodes if soap.index[i]==0):
            af.write('>' + str(ele) + '\n' + graphSeqs[ele] + '\n')
    
    with open(config.FINALRES, 'w') as seqfile:
        for key in scaffolds:
            seqid = '>' + str(key)
            seq = scaffolds[key]['seq']
            seqfile.write(seqid + '\n' + seq + '\n')
        
        if soapSinId:
            for ele in soapSinId:
                if len(graphSeqs[ele]) > 100:
                    seqfile.write('>' + str(ele) + '\n' + graphSeqs[ele] + '\n')
    
    logger.info('The seed contigs are writed into file: {0}'.format(config.SEED))
    logger.info('The recovered contigs are writed into file: {0}'.format(config.RECOVERED))
    logger.info('The scaffolds are writed into file: {0}'.format(config.FINALRES))
    logger.info('Done: ')
