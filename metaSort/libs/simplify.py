############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import math
import copy
import config
import shutil
import tempfile
import itertools
import indexGraph
import platform
import subprocess
import collections
#import scipy.optimize
import cPickle as pickle
from log import get_logger
from contextlib import contextmanager
# set log file
logger = get_logger()

noNeedCheck = set()

class similarity:
    """
    Stores the parameters for an alignment scoring function
    """
    def __init__(self, kmersize, simCutOff = 0.8, match = 2, mismatch = -1, gap = -5):
        
        self.simCutOff = simCutOff
        self.gap = gap
        self.match = match
        self.mismatch = mismatch
        self.kmersize = kmersize
        self.blPath = self.checkBl2seq()
    
    def __matchchar(self, a, b):
        """Return the score for aligning character a with b"""
        assert len(a) == len(b) == 1
        if a == b:
            return self.match
        else:
            return self.mismatch
    
    def __localAlign(self, x, y):
        """
        Use the dynamic programing algorithm to get the identity between two sequences
        """ 
        # create a zero-filled matrix
        A = [[0]*(len(y)+1) for i in xrange(len(x)+1)]
        # record the best socre and position in the matrix
        best = 0
        pos = (0, 0)
        # fill A from up left
        for i in xrange(1, len(x)+1):
            for j in xrange(1, len(y)+1):
                # the local alignment recurrance rule:
                A[i][j] = max(
                    A[i-1][j] + self.gap,
                    A[i][j-1] + self.gap,
                    A[i-1][j-1] + self.__matchchar(x[i-1], y[j-1]),
                    0
                )
                
                # track the cell with the largest score and position
                if A[i][j] >= best:
                    best = A[i][j]
                    pos = (i, j)
        
        # traceback step, find with operation resulted in the value of the cell and proceed to corresponding cell
        matchx = collections.deque()
        matchy = collections.deque()
        i, j = pos
        while A[i][j] > 0:
            # check up, left, digonal
            if A[i][j] == A[i-1][j] + self.gap: # gap from left
                matchx.appendleft(x[i-1])
                matchy.appendleft('-')
                i -= 1
            elif A[i][j] == A[i][j-1] + self.gap: # grap from up
                matchx.appendleft('-')
                matchy.appendleft(y[j-1])
                j -= 1
            elif A[i][j] == A[i][j-1] + self.mismatch:
                matchx.appendleft(x[i-1])
                matchy.appendleft(y[j-1])
                i -= 1
                j -= 1
            else:
                matchx.appendleft(x[i-1])
                matchy.appendleft(y[j-1])
                i -= 1
                j -= 1
        
        matchNum = len([i [0] for i in itertools.izip(matchx, matchy) if i[0] == i[1] != '-'])
        matchxlen = len([i for i in matchx if i != '-'])
        matchylen = len([i for i in matchy if i != '-'])
        matchLen = matchxlen if matchxlen > matchylen else matchylen
        return matchNum, matchLen
    
    def __sharedKmer(self, x, y):
        """
        use the shared kmer number to represent the similarity between two seqs
        """
        len1, len2 = len(x), len(y)
        ar1 = set([x[i:i+self.kmersize] for i in xrange(len1-self.kmersize+1)])
        ar2 = set([y[i:i+self.kmersize] for i in xrange(len2-self.kmersize+1)])
        revSeq2 = revCom(y)
        ar3 = set([revSeq2[i:i+self.kmersize] for i in xrange(len2-self.kmersize+1)])
        shared1 = len(ar1 & ar2)
        shared2 = len(ar1 & ar3)
        shared = shared1 if shared1 > shared2 else shared2
        len11 = len(ar1)
        len22 = len(ar2)
        short = len11 if len11 < len22 else len22
        similarity = shared*1.0/short
        return similarity
    
    def checkBl2seq(self):
        """
        check if bl2seq is available, if yes use it algin seqs
        """
        # check if the bl2seq program exists in the system path
        #with open(os.devnull, "w") as f:
        #    proc = subprocess.Popen(['bl2seq'], stdout=f, stderr=f)
        #    if proc:
        #        bl2seqPath = 'bl2seq'
        # use the default directory of bl2seq
        if platform.system() == 'Darwin':
            bl2seqPath = os.path.join(config.libsLocation, 'align/macosx/bl2seq')
        elif platform.system() == 'Linux':
            bl2seqPath = os.path.join(config.libsLocation, 'align/linux/bl2seq')
        
        with open(os.devnull, "w") as f:
            proc = subprocess.Popen([bl2seqPath], stdout=f, stderr=f)
            if proc:
                return bl2seqPath
            else:
                return 0

    @contextmanager
    def named_pipes(self, count):
        dirname = tempfile.mkdtemp()
        try:
            paths = []
            for i in range(count):
                paths.append(os.path.join(dirname, 'named_pipe' + str(i)))
                os.mkfifo(paths[-1])
            yield paths
        finally:
            shutil.rmtree(dirname)

    def bl2seq(self, seq1, seq2):
        """
        use bl2seq in the blast suit to algin two seqs
        """
        with self.named_pipes(2) as paths:
            cmd = [self.blPath, '-p', 'blastn', '-F', 'F',  '-i', paths[0], '-j', paths[1], '-D', '1']
            proc = subprocess.Popen(cmd, stdout = subprocess.PIPE)
            with open(paths[0], 'w') as af:
                af.write('>seq1\n' + seq1)
            
            with open(paths[1], 'w') as af:
                af.write('>seq2\n' + seq2)
            out = proc.communicate()
            
            # extract identity fraction
            out = out[0].split('\n')
            for line in out:
                if not line:
                    idenRatio = 0
                    break
                if line[0] != '#':
                    ae = line.split()
                    idenRatio = round(1 - int(ae[4])*1.0/int(ae[3]), 4)
                    break

            
            return idenRatio
    
    def identity(self, x, y, lenRatio = 0.8):
        """
        Caculate the similarity between two seqs
        """
        seq1, seq2 = x, y
        # check if seq1 and seq2 are exists
        if not seq1 or not seq2:
            return 100
        
        if self.blPath:
            return self.bl2seq(x, y)
        else:
            if len(seq1) >= 300 and len(seq2) >= 300:
                sim1 = self.__sharedKmer(seq1, seq1)
                if sim1 < 0.5:
                    matchNum, lencut = self.__localAlign(seq1, seq2)
                    reflen = len(seq1) if len(seq1) < len(seq2) else len(seq2)
                    if matchNum >= reflen*lenRatio:
                        sim = matchNum*1.0 / lencut
                        return 1 if sim >= self.simCutOff else 0
                    else:
                        return 0
                else:
                    return 1
            else:
                matchNum, lencut = self.__localAlign(seq1, seq2)
                reflen = len(seq1) if len(seq1) < len(seq2) else len(seq2)
                if matchNum >= reflen*lenRatio:
                    sim = matchNum*1.0 / lencut
                    if sim >= self.simCutOff:
                        return 1 if sim >= self.simCutOff else 0
                else:
                    return 0

# bubbles feature calling

def indexBubs(soap, bubbles):
    """
    Using num 0 -- len(bubbles) to index the bubbles
    """
    indexedBubs, count = {}, 1
    # remove the bubles that included by other paths
    for key in bubbles:
        tmp = bubbles[key]
        exclude = set()
        for com in itertools.combinations(tmp, 2):
            if set(com[0]).issubset(set(com[1])):
                exclude.add(com[0])
            elif set(com[1]).issubset(set(com[0])):
                exclude.add(com[1])
        # select the longest no overlap paths
        newtmp, bubNodes = [], []
        for value in tmp:
            if value not in exclude:
                newtmp.append(value)
                bubNodes.extend(value)
        # check if the bubble can be decomposed
        counter = collections.Counter(bubNodes)
        sharedNodes = [i for i in counter if counter[i]==len(newtmp)]
        if len(sharedNodes) == 2:
            if count in indexedBubs:
                count += 1
            indexedBubs[count] = newtmp
        else:
            try:
                sharedNodes.sort(key=lambda x: newtmp[0].index(x))
                pairs = [(sharedNodes[i], sharedNodes[i+1]) for i in range(len(sharedNodes)-1)]
                for pair in pairs:
                    decomposed = set()
                    for path in newtmp:
                        s, e = [path.index(i) for i in pair]
                        decomposed.add(path[s: e+1])
                    decomposed = [i for i in decomposed if len(i)>2]
                    if len(decomposed) > 1:
                        if count in indexedBubs:
                            count += 1
                        indexedBubs[count] = decomposed
            except:
                # there is a loop on one path
                if count in indexedBubs:
                    count += 1
                indexedBubs[count] = newtmp
        
    # make the bubble one strand
    #for num in indexedBubs:
    #    tmp = []
    #    for path in indexedBubs[num]:
    #        tmp.append(tuple([i-1 if soap.index[i]==2 else i for i in path]))
    #    indexedBubs[num] = tmp

    # accept only dlink bubbles
    delNum = []
    for k, v in indexedBubs.iteritems():
        count = len(v)
        paths = []
        for path in v:
            for pos in xrange(len(path)-1):
                try:
                    soap.dDlink[path[pos]][path[pos+1]]
                except:
                    count -= 1
                    paths.append(path)
                    break
        if count <= 1:
            delNum.append(k)
        elif count < len(v):
            for path in paths:
                indexedBubs[k].remove(path)
    for ele in delNum:
        del indexedBubs[ele]

    for ele in indexedBubs.keys():
        if len(indexedBubs[ele]) > 3:
            del indexedBubs[ele]

    return indexedBubs

# subspecies numbers
def setXMax(xMax, binWidth):
    return int((math.floor(xMax / binWidth)) * binWidth)

def getFirstXMax(vector, binWidth):
    maxCov = subMaxCov = 0
    for i in vector:
        if i > maxCov:
            subMaxCov = maxCov
            maxCov = i
    
    xMax = setXMax(subMaxCov, binWidth) + binWidth * 5
    return xMax

def setWidthByXMax(xMax):
    listWidth = [0, 0]  # [binWidth, widthMovAve]
    if xMax > 300:
        listWidth = [6, 5]
    if xMax <= 300:
        listWidth = [4, 3]
    if xMax <= 120:
        listWidth = [2, 3]
    if xMax <= 100:
        listWidth = [1, 1]
    
    return listWidth

def histo(vector, xMin, xMax, binWidth):
    """
    Get the number into bins
    """
    bins = range(xMin, xMax, binWidth)
    dicHisto = { x : 0 for x in bins}
    for cov in vector:
        if cov < xMin or cov > xMax:
            continue
        
        for x in bins:
            if x <= cov < x+binWidth:
                dicHisto[x] += 1
                break
    
    return dicHisto

def smoothingHisto(dicHisto, xMin, xMax, binWidth, widthMovAve):
    
    dicSmoothHisto = {}
    listMovAve = []
    
    for x in xrange(xMin, xMax, binWidth):
        listMovAve.append(dicHisto[x])
        if len(listMovAve) < widthMovAve:
            continue
        dicSmoothHisto[x - binWidth*((widthMovAve - 1)/2)] = sum(listMovAve)/float(widthMovAve)
        listMovAve.pop(0)
    
    return dicSmoothHisto

def detectPeakPandS(dicHisto, xMin, xMax, binWidth, thresHeight, listPeakPandS):
    
    countIncrease = countDecrease = 0
    thresIncrease = thresDecrease = 3
    beforeHeight = -1
    flagPeakStart = False
    peakHeight = peakCov = 0
    
    for x in range(xMax - binWidth, xMin - binWidth, -1 * binWidth):
        if beforeHeight == -1:
            beforeHeight = dicHisto[x]
            continue
        
        if not flagPeakStart:
            if dicHisto[x] >= thresHeight:
                if dicHisto[x] >= beforeHeight:
                    countIncrease += 1
                    if countIncrease >= thresIncrease:
                        countIncrease = 0
                        flagPeakStart = True
            beforeHeight = dicHisto[x]
        
        if flagPeakStart:
            if dicHisto[x] >= peakHeight:
                peakHeight = dicHisto[x]
                peakCov = x
            else:
                countDecrease += 1
                if countDecrease >= thresDecrease:
                    for i in range(2):
                        if listPeakPandS[i] == -1:
                            tmpBias = float(binWidth) / 2
                            listPeakPandS[i] = peakCov + tmpBias
                            peakHeight = 0; peakCov = 0
                            break
                    if listPeakPandS[1] != -1:
                        return listPeakPandS
                    countDecrease = 0
                    flagPeakStart = False
    
    return listPeakPandS

def peakFinding(coverVector):
    """
    Finding the peaks in the vector, used for bubble classification
    """
    listPeak = []
    binWidth = 4
    xMin = 0
    xMax = 1000
    widthMovAve = 5
    listPeakPandS = [-1, -1]
    thresHeight = 0
    while True:
        if len(listPeak) == 0:
            xMax = getFirstXMax(coverVector, binWidth)
        # Set width and xMax
        listWidth = setWidthByXMax(xMax)
        binWidth, widthMovAve = listWidth
        xMax = setXMax(xMax, binWidth)
        # make smoothed histogram
        xMin = 0
        dicHisto = histo(coverVector, xMin, xMax, binWidth)
        dicSmoothHisto = smoothingHisto(dicHisto, xMin, xMax, binWidth, widthMovAve)
        xMin += binWidth*((widthMovAve - 1)/2)
        xMax -= binWidth*((widthMovAve - 1)/2)
        # get thresHeight
        if len(listPeak) == 0:
            thresHeight = dicSmoothHisto[xMax - binWidth]
        # detect (primary and ) secondary peak
        listPeakPandS = detectPeakPandS(dicSmoothHisto, xMin, xMax, binWidth, thresHeight, listPeakPandS)
        # record peak
        if len(listPeak) == 0:
            listPeak.append(listPeakPandS[0])
        listPeak.append(listPeakPandS[1])
        # when could not detect secondary peak, break
        if listPeakPandS[1] == -1:
            listPeak.pop(-1)
            return listPeak
        # prepare for enxt peak
        listPeakPandS[0] = listPeakPandS[1]
        listPeakPandS[1] = -1
        xMax = listPeakPandS[0]

def aveStd(vector):
    """
    Caculate the average and standard variation of a list
    """
    if not vector:
        raise Exception('Warning: the vector length is zero')
    for ele in vector:
        if isinstance(ele, str):
            raise Exception('Waring: there are strings in the vector')
    mean = sum(vector)*1.0/len(vector)
    var = sum([(i-mean)**2 for i in vector])/(len(vector)-1)
    sd = var**0.5
    return (mean, sd)

def getSubspBubbles(bubbles, soap, mappedNodes):
    """
    Distingish the number of subspeies from the bubble info
    use 3 criterias to filter the fake bubbules
    """
    leftNodes = set()
    cleaned = {}
    for key in bubbles:
        kcs = []
        for path in bubbles[key]:
            kcs.append(sum([soap.cov[i] for i in path[1:-1]]))
        kcs.append(soap.cov[path[-1]])
        kcs.insert(0,soap.cov[path[0]])
        # 1: the kmer coverage of every bubble node should larger than 0
        if 0 not in kcs:
            innerCov = sum([i for i in kcs[1:-1]])
            # 2: the kmer of the two end of nodes should larger than the sum of the two bubble nodes
            if innerCov <= kcs[0] + 5 and innerCov <= kcs[-1] + 5:
                cleaned[key] = kcs
    # 3: use the kmer coverage to get the subspecies level bubbles
    if len(cleaned) < 5:
        return 0, 0
    else:
        num, kcMax = spNumRan(config.STATUS, soap, mappedNodes)
        if num == 0: return 1, bubbles
        # filter the bubbles with kcmax
        for key in cleaned.keys():
            allcov = cleaned[key]
            if allcov[0] > kcMax or allcov[-1] > kcMax:
                del cleaned[key]
    # because the kmer coverage will not extrac when the length of node is small, so it impossible to distingish them
    for ele in bubbles.keys():
        if ele not in cleaned:
            del bubbles[ele]
    if len(bubbles) == 0:
        return 0, 0
    else:
        return num, bubbles
 
def spNumRan(statusfn, soap, mappedNodes):
    """
    Get the total kmer coverage of subspcies by using long contigs, the should be the shared region of subspecies
    """
    sharedRegion = set()

    with open(statusfn, 'r') as af:
        for line in af:
            info = line.rstrip().split()
            if int(info[1]) == 0:
                sharedRegion.add(int(info[0]))
    candi = [i for i in mappedNodes if i in sharedRegion and soap.index[i]==0]
    #candi.sort(key=lambda x: length[x])
    covPoll = [soap.cov[i] for i in candi if soap.length[i]> config.CLEANLENCUTOFF]
    # find the peaks
    peaks = peakFinding(covPoll)
    if len(peaks) == 1:
        # no subspecies or the kmer of subspecies are similar
        if len(covPoll) <= 3:
            return 0, 0
        else:
            mean, std = aveStd(covPoll)
            return 1, mean + 3*std
    else:
        bins = []
        peaks.sort()
        # the number of bins are the same as peaks
        for i in xrange(len(peaks)):
            bins.append([])
        for ele in covPoll:
            if ele <= peaks[0]:
                bins[0].append(ele)
            elif ele >= peaks[-1]:
                bins[-1].append(ele)
        for pos in xrange(0, len(peaks)-1):
            low, high = peaks[pos], peaks[pos+1]
            for ele in covPoll:
                if low < ele < high:
                    if abs(ele-low) < abs(ele-high):
                        bins[pos].append(ele)
                    else:
                        bins[pos+1].append(ele)
        # this step can only get the conserved subspecies, because of two subspeices have a similar kmer coverage
        conservedSpNum = len(peaks) # the last one is the sum of subspecies kmer coverage
        kcRan = 0
        for ele in bins:
            if len(ele) > 2:
                mean, std = aveStd(ele)
                kcMax = mean + 3*std
            if len(ele) == 2:
                kcMax = max(ele)
            elif len(ele) == 1:
                kcMax = ele[0] + 1
        return conservedSpNum, kcMax

# tips and bubbles
def revCom(string):
    """
    reverse complement the given string
    """
    comp = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C',
            'n':'N', 'N':'N',
            'a' : 'T', 't' : 'A', 'c' : 'G', 'g' : 'C'
            }
    result = ''
    for i in string[::-1]:
        result += comp[i]
    
    return result

def scoreMatrix(path, identity, candidate, cov, out = 'bubble'):
    """
    Used for filter the bubbles and tips, select the best nodes represent the real genome.
    output a score matrix which contains four columns: 1, number of nodes in mapped ndoes; 2 average identity
    5:number of nodes in the candidate; 4: kmer coverage distance with the two end of the nodes on the path
    1 to 4 representt the importance weight, 1 indicate the most important, if a path have a higher value of 1 than another
    path, it means it is the best one no matter the value of other 3 value.
    """
    if out == 'bubble':
        mappedNum = len([i for i in path[1:-1] if i in identity])
        candiNum = len([i for i in path[1:-1] if i in candidate])
        if mappedNum == 0:
            averIdentity = 0
        else:
            averIdentity = sum([identity[i] for i in path[1:-1] if i in identity])*1.0/mappedNum
        refCov = (cov[path[0]] + cov[path[-1]])/2
        queryCovList = [cov[i] for i in path[1:-1]]
        queryCov = refCov - sum(queryCovList)/len(queryCovList)
        return (mappedNum, averIdentity, candiNum, queryCov)
    elif out == 'tip':
        if path[1] in identity:
            iden = identity[path[1]]
        else:
            iden = 0
        if path[1] in candidate:
            candi = 1
        else:
            candi = 0
        queryCov = 0 - (cov[path[0]] - cov[path[1]] )
        return (iden, candi, queryCov)
    elif out == 'mirro':
        mappedNum = len([i for i in path if i in identity])
        candiNum = len([i for i in path if i in candidate])
        if mappedNum == 0:
            averIdentity = 0
        else:
            averIdentity = sum([identity[i] for i in path[1:-1] if i in identity])*1.0/mappedNum
        refCov = (cov[path[0]] + cov[path[-1]])/2
        queryCovList = [cov[i] for i in path[1:-1]]
        queryCov = refCov - sum(queryCovList)/len(queryCovList)
        return (mappedNum, averIdentity, candiNum, queryCov)

def removeTips(graph, soap, identity, candidate):
    """
    There will be common regions (Tips) from other genomes be searched in the eariler steps.
    definition: short dead end (and also have a different kmer cov in meta graph) On meta graph: almost tips are subspecies
                Strategies: if a node have a outdegree >=2, and have 1 follow node that its outdegree is 0
                           it will be taken as tip. Only one not more than one follow node, because we don't know
                           the result if we delet one more follow nodes (nodes are different in length)
    """
    candiTips, candiNodes = {}, set()
    # select the node randomly but the outdegree of this node should larger than 1
    for node in graph.keys():
        if graph.outdegree(node) > 1:
            # first find the tips
            deadEnd, aliveEnd = set(), set()
            for second in graph[node].keys():
                if graph.degree(second) == (1,0):
                    deadEnd.add(second)
                else:
                    aliveEnd.add(second)
            # get the types of the tip
            if deadEnd and not aliveEnd:
                candiTips[node] = [list(deadEnd)]
                candiNodes.update(deadEnd)
            elif deadEnd and aliveEnd:
                candiTips[node] = [list(deadEnd), list(aliveEnd)]
                candiNodes.update(deadEnd)
                candiNodes.update(aliveEnd)
    
    # find the true subspecies
    graphSeqs = soap.fetchSeqs(candiNodes)
    tips, discard = {}, set()
    #print('a')
    sim = similarity(soap.kmersize)
    #print('b')
    for key, value in candiTips.iteritems():
        trueSub = set()
        # all of the point are dead end
        if len(value) == 1:
            for com in itertools.combinations(value, 2):
                try:
                    seqIdentity = sim.identity(graphSeqs[com[p1]], graphSeqs[com[p2]])
                except:
                    seqIdentity = 1
                
                if seqIdentity > 0:
                    trueSub.update(com)
            
            # check which nodes will be removed
            if trueSub:
                tips[node] = list(trueSub)
                maxScore, leftId = 0, 0
                for ele in trueSub:
                    scoreList = scoreMatrix([key, ele], identity, candidate, soap.cov, out = 'tip')
                    if scoreList > maxScore:
                        maxScore, leftId = scoreList, ele
                discard.update([i for i in trueSub if i != leftId])
        elif len(value) == 2:
            deadEnd, aliveEnd = value
            # not all of the nodes are deadEnd, delete the dead end nodes
            for v1 in deadEnd:
                for v2 in aliveEnd:
                    try:
                        if sim.identity(graphSeqs[v1], graphSeqs[v2]) > 0:
                            trueSub.update([v1,v2])
                            discard.add(v1)
                    except:
                        discard.add(v1)
            
            if trueSub:
                 tips[node] = list(trueSub)
    
    # remove dead ends
    _ = graph.deleteNodes(discard)
    
    return tips

def branchCut(graph, soap, identity, candidate):
    """
    Detect long tips and remove them
    """
    logger.info('Detecting long tips ...')
    
    toDelete, deletedNodes, reservedSp = set(), set(), set()
    passedSp = {}
    
    # formate of subSpecies: {0: [ (1,2,3), (5) ], 7:[ (8,9,10),(11,22,45) ],  ...}
    subSpecies = detectLongSpBranches(graph, soap)
    
    # filter and select the represent branch, and remove other subspecies branches
    for sourceNode in subSpecies:
        # filter, all of the tips should be single chain (outdegree of the last node should be zero)
        elements = subSpecies[sourceNode]
        for ele in elements:
            if isinstance(ele, int):
                break
            elif len(ele) == 1 or graph.outdegree(ele[-1]) != 0:
                break
        else:
            reservedSp.add(sourceNode)
            passedSp[sourceNode] = elements
            representId, score = {}, []
            for pos, chain in enumerate(elements):
                mappedNum = len([i for i in chain if i in identity])
                candiNum = len([i for i in chain if i in candidate])
                averIdentity = 0 if mappedNum == 0 else sum([identity[i] for i in chain if i in identity])*1.0/mappedNum
                scorelist = (mappedNum, averIdentity, candiNum)
                score.append(scorelist)
                representId[scorelist] = pos
            # select the best candidate
            score.sort()
            best = representId[score[-1]]
            for pos, ele in enumerate(elements):
                if pos != best:
                    toDelete.update(ele)
    
    logger.info(str(len(toDelete)) + ' nodes are deleted from the graph')
    if toDelete:
        _ = graph.deleteNodes(toDelete)
        logger.info('Graph size: ' + str(len(graph.vertices())) + ', now')
    
    return passedSp

def detectLongSpBranches(oriGraph, soap):
    """
    Sometimes the variations from subspecies not represented as the formate of bubbles but tips
    Compare them using shared kmer ratio
    """
    # compress the graph
    comGraph, combinedNodes = oriGraph.compressGraph() 
    
    # find long tips
    candiTips = {}
    for node in comGraph.iterkeys():
        if node in noNeedCheck:
            continue
        
        if len(comGraph[node]) > 1:
            deadEnd, aliveEnd = set(), set()
            sp = [] # subspecies
            for second in comGraph[node].keys():
                if second not in comGraph:
                    rev = soap.revNode(second)
                    if len(comGraph[rev]) == 1:
                        deadEnd.add(second)
                else:
                    aliveEnd.add(second)
            
            if deadEnd:
                tmp = [deadEnd]
                if aliveEnd:
                   tmp.append(aliveEnd)
                candiTips[node] = tmp
    
    # fetch the nodes sequences
    nodes = set([e for i in candiTips.itervalues() for k in i for e in k])
    for key in copy.copy(nodes):
        if key in combinedNodes:
            nodes.update(combinedNodes[key])
    graphSeqs = soap.fetchSeqs(nodes)
    
    # extend the merged sequences
    for key in graphSeqs:
        if key in combinedNodes:
            chain = combinedNodes[key]
            seq = graphSeqs[chain[0]]
            # check the link types, note the newlink overlaps, like: -32
            for i in xrange(0, len(chain)-1):
                try: # dlink
                    if oriGraph[chain[i]][chain[i+1]] > 0:
                        try:
                            seq += 'N'*oriGraph[chain[i]][chain[i+1]]
                            seq += graphSeqs[chain[i+1]]
                        except:
                            seq += 'N'
                    else:
                        try:
                            overlap = abs(oriGraph[chain[i]][chain[i+1]])
                            seq += graphSeqs[chain[i+1]][overlap:]
                        except:
                            seq += 'N'
                except: # pelink
                    gap = oriGraph[chain[i]][chain[i+1]][0]
                    if gap > 0:
                        try:
                            seq += 'N'*gap
                            seq += graphSeqs[chain[i+1]]
                        except:
                            seq += 'N'
                    else:
                        try:
                            seq += graphSeqs[chain[i+1]][-gap:]
                        except:
                            seq += 'N'
                    
            graphSeqs[key] = seq
    
    sim = similarity(soap.kmersize)
    subSpecies, sp = {}, []
    for node, value in candiTips.iteritems():
        # value = [set(), set()] or [set()]
        deadEnd = [k for i in value for k in i]
        # get the similarity of the appendant nodes, a 2-d matrix
        for pos in xrange(len(deadEnd)-1):
            for pos2 in xrange(pos+1, len(deadEnd)):
                try:
                    seqIdentity = sim.identity(graphSeqs[deadEnd[pos]], graphSeqs[deadEnd[pos2]])
                except:
                    seqIdentity = 1
                if  seqIdentity > 0:
                    sp.append([deadEnd[pos], deadEnd[pos2]])
            
            if len(sp) == 1:
                subSpecies[node] = sp[0]
            elif len(sp) > 1:
                # handel the sp list by a tree structure, and get the subgraphs
                g, tmp, links, cl, groups = 1, [], {}, {}, {}
                for ele in sp:
                    if ele[0] not in links:
                        links[ele[0]] = {}
                    if ele[1] not in links:
                        links[ele[1]] = {}
                    cl[ele[0]] = cl[ele[1]] = 1
                    links[ele[0]][ele[1]] = 1
                    links[ele[1]][ele[0]] = 1
                
                for node in cl:
                    if cl[node] == 1:
                        tmp.append(node)
                        while tmp:
                            last = tmp.pop()
                            cl[last] = g + 1
                            if last in links:
                                tmp.extend([ x for x in links[last] if cl[x] == 1])
                        g += 1
                
                for ele in cl.iterkeys():
                    if cl[ele] not in groups:
                        groups[cl[ele]] = []
                    groups[cl[ele]].append(ele)
                
                for v in groups.itervalues():
                    if len(v) > 1:
                        subSpecies[node] = v
    
    for node in subSpecies:
        elements = subSpecies[node]
        backedEle = []
        for ele in elements:
            if ele in combinedNodes:
                backedEle.append(combinedNodes[ele])
            else:
                backedEle.append((ele))
        subSpecies[node] = backedEle
    
    # formate of subSpecies: {0: [ (1,2,3), (5) ], 7:[ (8,9,10),(11,22,45) ],  ...}
    return subSpecies

def allPathsIn2Nodes(start, end, graph):
    """
    Find all of the paths between two nodes, they are bubbles
    """
    paths = set()
    stack = [[start]]
    while stack:
        tmp = stack.pop()
        last = tmp[-1]
        if last == end:
            paths.add(tuple(tmp))
            continue
        
        if last in graph:
            for second in graph[last].keys():
                longer = tmp + [second]
                if len(longer) < 10:
                    stack.append(longer)
            
    return paths

def removeMirrorNodes(graph, soap, identity, candidate):
    """
    Remove the mirror structur (subspecies also)
    The structure likes: 1 -> 2; 1->3 ; 5 -> 2 ; 5 -> 3.
    1 and 5 are subspecies, remove the short path containing the 1 or 5
    """
    # get all of the nodes on the graph
    pairs = {}
    for node in graph.keys():
        if graph.outdegree(node) == 2:
            nodes = tuple(graph[node].keys())
            if nodes not in pairs:
                pairs[nodes] = []
            pairs[nodes].append(node)
    
    # choose the candidate mirror structure
    candPair = set()
    for k, v in pairs.iteritems():
        if len(v) == 2:
            if soap.length[v[0]] > 5*soap.length[v[1]] or soap.length[v[1]] > 5*soap.length[v[0]]:
                continue
            else:
                candPair.add(k)
    graphSeqs = soap.fetchSeqs([k for i in candPair for k in pairs[i]])
    
    # the format the counter is: (5,8):[node1, node3]
    sim = similarity(soap.kmersize)
    sp, perfect, half = [], set(), set()
    deleteNodes = set()
    for k in candPair:
        v = pairs[k]
        try:
            seqIdentity = sim.identity(graphSeqs[v[0]], graphSeqs[v[1]])
        except:
            seqIdentity = 1
        
        if seqIdentity > 0:
            sp.append(tuple(v))
            model = set([graph.indegree(i) for i in v])
            if model == set([0]):
                # delete the shorter node
                v.sort(key=lambda x: soap.length[x])
                perfect.add(v[0])
            elif 0 in model:
                half.update([i for i in v if graph.indegree(i)==0])
            else:
                # don't know which node to delete, delete the two nodes
                deleteNodes.update(k)
    
    # delete one of the pair nodes
    deleteNodes.update(perfect)
    deleteNodes.update(half)
    
    if deleteNodes:
        _ = graph.deleteNodes(deleteNodes)
    
    return sp

def backward(graph, soap, ref, bubbles):
    """
    Back trace the visited roads, and find the bubbles
    """
    tmpBub = set()
    border = 10
    visitedNodes = set(ref[:-1])
    sourceNode = soap.revNode(ref[-1])
    stack = [[sourceNode]] # source node
    while stack:
        tmp = stack.pop()
        last = tmp[-1]
        if soap.revNode(last) in visitedNodes:
            if len(tmp) > 2:
                # a bubbles is found
                revLast = soap.revNode(last)
                pos = ref.index(revLast)
                revStrand = [soap.revNode(i) for i in tmp[::-1]]
                if 2 < len(ref[pos:]) < 10: 
                    tmpBub.update([tuple(ref[pos:]), tuple(revStrand)])
            continue
        if last in graph:
            for sec in graph[last]:
                longer = tmp + [sec]
                if len(longer) > border:
                    continue
                stack.append(longer)
    if tmpBub:
        # in case of , find the most left [(495021, 281577, 192783, 548550, 216361, 529483, 150603, 163431),
                #(245339, 159056, 569792, 491008, 163431),(495021, 245339, 159056, 569792, 491008, 163431)]
        outcomes = itertools.permutations(tmpBub, 2)
        for ele in outcomes:
            if set(ele[0]).issubset(ele[1]):
                tmpBub.discard(ele[0])
        left = set([i[0] for i in tmpBub])
        tmpBub = list(tmpBub)
        if len(left) > 1:
            newBub = []
            for path in tmpBub:
                if left.issubset(path):
                    # the path can determine the order
                    order = sorted([path.index(i) for i in left])
                    leftNodes = path[order[0]:order[-1]+1]
                    break
            # complement the shorter path
            for path in tmpBub:
                if path[0] == leftNodes[0]:
                    newBub.append(path)
                else:
                    most = leftNodes.index(path[0])
                    newBub.append(leftNodes[:most] + path)
            if leftNodes[0] not in bubbles:
                bubbles[leftNodes[0]] = []
                bubbles[leftNodes[0]].append(newBub)
        else:
            fir = next(iter(left))
            if fir not in bubbles:
                bubbles[fir] = []
            bubbles[fir].append(tmpBub)

def findBubble(graph, soap):
    """
    Find the bubbles on the graph, first travel the deepest path, and trace back 
    """
    logger.info('Finding bubbles ')
    logger.info('Results ... ')
    bubbles = {}
    group = graph.findSubgraph(out='gp')
    for key, elements in group.iteritems():
        visitedNodes = set()
        startNodes = set([i for i in elements if graph.indegree(i)==0])
        if not startNodes:
            loop = set()
            for node in elements:
                if key in graph:
                    for key in graph[node]:
                        if key in graph and node in graph[key]:
                            loop.update([node, key])
            startNodes = set(elements) & loop
        if not startNodes: # loop
            startNodes = set(elements)
        while startNodes:
            node = startNodes.pop()
            if node not in visitedNodes:
                stack = [[node]]
                while stack:
                    tmp = stack.pop()
                    last = tmp[-1]
                    if last in tmp[:-1]:
                        # a loop
                        continue
                    if last in visitedNodes:
                        if len(tmp) > 2 :
                            bub = backward(graph, soap, tmp, bubbles)
                            if bub:
                                bubbles.append(bub)
                        continue
                    visitedNodes.update([last, soap.revNode(last)])
                    if last in graph:
                        for sec in graph[last]:
                            longer = tmp + [sec]
                            stack.append(longer)
    
    # formate of bubbles: {1: [[(1,2,3),(1,2,4)]],
    #                      10: [ [(12,15,16), (12,17,16)],      [(12,18,19), (12,110,19)] ],
    #                      31: [ [(31,32,33,334),(31,35,334), (31,91, 92,334)],  [(31,38,33), (31,39,33)] ]}
    
    record, flag = set(), 1
    keys = set(bubbles.keys())
    while flag:
        flag = 0
        for key in keys:
            if key not in record:
                internalNodes = set()
                for path in bubbles[key]:
                    for ele in path:
                        internalNodes.update(ele[1:-1])
                nest = internalNodes.intersection(keys) - record
                if nest:
                    flag = 1
                    for ele in nest:
                        bubbles[key].extend(bubbles[ele])
                        record.add(ele)
    
    for ele in record:
        del bubbles[ele]
    
    logger.info('#'*60)
    logger.info('{0:,} bubbles structure are found'.format(len(bubbles)))
    tanglyBubs = {}
    for key in bubbles.keys():
        if len(bubbles[key]) > 1:
            tmpG, revG = {}, {}
            graphNodes = set()
            for bub in bubbles[key]:
                for path in bub:
                    for pos in xrange(len(path)-1):
                        graphNodes.update([path[pos], path[pos+1]])
                        if path[pos] not in tmpG: tmpG[path[pos]] = {}
                        tmpG[path[pos]][path[pos+1]] = 1
                        if path[pos+1] not in revG: revG[path[pos+1]] = {}
                        revG[path[pos+1]][path[pos]] = 1
            source = [i for i in graphNodes if i not in revG]
            sink = [i for i in graphNodes if i not in tmpG]
            if len(source) == len(sink) == 1:
                # traversal the graph
                paths = []
                stack = [[source[0]]]
                while stack:
                    tmp = stack.pop()
                    last = tmp[-1]
                    if last in tmp[:-1]:
                        continue
                    if last in tmpG:
                        for sec in tmpG[last]:
                            longer = tmp + [sec]
                            stack.append(longer)
                    else:
                        paths.append(tuple(tmp))
                        # in the case of cycle
                        if len(paths) == 100:
                            paths = []
                            break
                if paths:
                    bubbles[key] = [paths]
                else:
                    tanglyBubs[key] = bubbles[key]
            else:
                tanglyBubs[key] = bubbles[key]

    logger.info('{0:,} tangly structure in it'.format(len(tanglyBubs)))
    for key in tanglyBubs:
        del bubbles[key]

    bubNodes = set([z for i in bubbles.itervalues() for k in i for v in k for z in v])
    bubNodes.update([soap.revNode(i) for i in bubNodes])
    # combine the paths in the complex structure
    delNodes, delPathNum, allPathNum = set(), 0, 0
    for key in tanglyBubs:
        local = {}
        for paths in tanglyBubs[key]:
            for path in paths:
                if path[0] not in local:
                    local[path[0]] = {'p': set(), 'n': set()}
                local[path[0]]['p'].add(tuple(path))
                allPathNum += 1
                local[path[0]]['n'].update([k for i in path for k in (i, soap.revNode(i))])
        # handle one by one, remove the paths that have no branches
        for sourceNode in local:
            purePath, dirtyPath = [], []
            for path in local[sourceNode]['p']:
                # check if the path has branches
                for node in path[1:-1]:
                    if node in graph:
                        if set(graph[node].keys()) - local[sourceNode]['n']:
                            dirtyPath.append(path)
                            break
                    if soap.revNode(node) in graph:
                        if set(graph[soap.revNode(node)].keys())- local[sourceNode]['n']:
                            dirtyPath.append(path)
                            break
                else:
                    purePath.append(path)
            if purePath:
                if dirtyPath:
                    reservedNodes = set([node for path in dirtyPath for node in path])
                    delPathNum += len(purePath)
                else:
                    reservedNodes = set([node for node in purePath[0]])
                    delPathNum += len(purePath)-1
                # remove path nodes
                reservedNodes.update([soap.revNode(i) for i in reservedNodes])
                delNodes.update(local[sourceNode]['n'] - reservedNodes)
    
    logger.info('{0:,} out of the {1:,} paths are deleted'.format(delPathNum, allPathNum))
    
    # may be some to be deleted nodes have it's own position in bubble nodes, exclude them
    _ = graph.deleteNodes(delNodes - bubNodes)

    # check if the bubble are cause by subspecies, using similarity
    for bub in bubbles.keys():
        if len(bubbles[bub]) == 1:
            path = bubbles[bub][0]
            bubbles[bub] = path
    
    graphSeqs = soap.fetchSeqs([x for v in bubbles.itervalues() for k in v for x in k])
    identity = []
    sim = similarity(soap.kmersize)
    for bub in bubbles.keys():
        # fetch rout sequences
        seqs = []
        bubblePaths = bubbles[bub]
        for path in bubblePaths:
            bPath = path[1:-1]
            if len(bPath) == 1:
                seq = graphSeqs[bPath[0]] if bPath[0] in graphSeqs else ''
            else:
                seq = graphSeqs[bPath[0]] if bPath[0] in graphSeqs else ''
                for p in xrange(1, len(bPath)):
                    pre, aft = bPath[p-1], bPath[p]
                    overlap = graph[pre][aft]
                    if isinstance(overlap, int):
                        try:
                            seq += graphSeqs[aft][abs(overlap):]
                        except:
                            seq += 'N'
                    elif isinstance(overlap, int):
                        overlap = overlap[0]
                        if overlap < 0:
                            try:
                                seq += graphSeqs[aft][abs(overlap):]
                            except:
                                seq += 'N'
                        else:
                            seq += 'N'*overlap
                            try:
                                seq += graphSeqs[aft]
                            except:
                                seq += 'N'
            seqs.append(seq)
        # compare the similarity among the sequences
        similar = set()
        for i in xrange(len(seqs)-1):
            for k in xrange(i+1, len(seqs)):
                try:
                    seqIdentity = sim.identity(seqs[i], seqs[k])
                except:
                    seqIdentity = 1
                
                if seqIdentity:
                    similar.update([i, k])
                    identity.append(seqIdentity)
        
        dele = set(xrange(len(seqs))) - similar 
        if dele:
            for pos in sorted(dele, reverse=True):
                bubblePaths.pop(pos)
            
            if bubblePaths:
                bubbles[bub]= bubblePaths
            else:
                del bubbles[bub]
        
    bubbles = cleanBubbles(bubbles, soap)
    logger.info('At last {0} bubbles are found'.format(len(bubbles)))
    logger.info('#'*60)

    return bubbles, identity

def cleanBubbles(bubbles, soap):
    # clean the bubbles like: 502290: [(502290, 610068, 502289), (502290, 610067, 502289)]}
    logger.info('Cleaning bubbles structures ...')
    # if there are reverse complement nodes in a bubbles, it is a fake bubble
    for key in bubbles.keys():
        paths = bubbles[key]
        allnodes = []
        for path in paths:
            allnodes.append(set(path[1:-1]))
        for pos in xrange(len(allnodes)-1):
            rev = set([soap.revNode(i) for i in allnodes[pos+1]])
            if rev.intersection(allnodes[pos]): # has intersection
                del bubbles[key]
                break
    
    return bubbles

def flateBubble(graph, soap, bubbles, identity, candidate):
    """
    Select the best bubbles sequence to represent the path
    only the node in the bubble have a indegree and outdegree of one will be taken as bubble path
    """
    nodesTobeDelete = set()
    for source in bubbles:
        bub = bubbles[source]
        bubbleNodes = [ node for path in bub for node in path[1:-1] ]
        counter = collections.defaultdict(int)
        for node in bubbleNodes:
            counter[node] += 1
        branchNodes = [i for i in bubbleNodes if counter[i] == 1]
        model = set([graph.degree(node=i) for i in branchNodes])
        if model == set([(1,1)]):
            score = []
            representId = {}
            for pos, path in enumerate(bub):
                scoreList = scoreMatrix(path, identity, candidate, soap.cov)
                representId[scoreList] = pos
                score.append(scoreList)
            score.sort()
            best = representId[score[-1]]
            ref = bub[best]
            # if the nodes on the bubble not present in the ref path, then delete them
            nodesTobeDelete.update([node for node in bubbleNodes if node not in ref ])
        else:
            # may be existing a perfect bubble, some of the inner node extend to another paths that not belong to the bubble paths
            allNodes = set([node for path in bub for node in path])
            dirtyNodes = set()
            for path in bub:
                for node in path[1:-1]:
                    rev = soap.revNode(node)
                    for i in graph[node]:
                        if i not in allNodes:
                            dirtyNodes.add(i)
                    for i in graph[rev]:
                        if i not in allNodes:
                            dirtyNodes.add(i)
            mark = set()
            for pos, path in enumerate(bub):
                if dirtyNodes.intersection(path):
                    mark.add(pos)
            
            if 0 < len(mark) < len(bub): # can delete bubbles
                allpos = set(range(len(bub)))
                reservedNodes = set()
                for i in (allpos - mark):
                    reservedNodes.update([bub[i]])
                nodesTobeDelete.update([node for node in bubbleNodes if node not in reservedNodes ])
    
    _ = graph.deleteNodes(nodesTobeDelete)

def recordNodeStatus(graph, soap, tips, mirros, bubbles, branches):
    """
    Record the status of source node of bubbles, tips, branches, for later classification
    0 cleaned node (no tips,no bubbles, no branches);
    1 tips; 2: bubbles 3 branches 4: mirror
    """
    # mark all of the nodes
    nodeStatus = {}
    for ele in graph.vertices():
        nodeStatus[ele] = 0
    # mark tips
    tipNodes = set()
    for key in tips:
        tipNodes.add(key)
        for node in tips[key]:
            tipNodes.add(node)
    for ele in tipNodes:
        rev = soap.revNode(ele)
        nodeStatus[ele] = 1
        nodeStatus[rev] = 1
    # mark bubbles
    bubbleNodes = set()
    for key in bubbles:
        for path in bubbles[key]:
            bubbleNodes.update(path)
    for ele in bubbleNodes:
        rev = soap.revNode(ele)
        nodeStatus[ele] = 2
        nodeStatus[rev] = 2
    # mark branches
    branchNodes = set()
    for key in branches:
        branchNodes.add(key)
        for path in branches[key]:
            branchNodes.update(path)
    for ele in branchNodes:
        rev = soap.revNode(ele)
        nodeStatus[ele] = 3
        nodeStatus[rev] = 3
    
    # mark mirror nodes
    mirNodes = set()
    for ele in mirros:
        mirNodes.update(ele)
    for ele in mirNodes:
        rev = soap.revNode(ele)
        nodeStatus[ele] = 4
        nodeStatus[rev] = 4
    return nodeStatus

def merge(graph, soap, identity, candidate):
    """
    find and distinguish the subspeice of the genome
    """
    # find tips
    logger.info('Trimming tips ...')
    tips = removeTips(graph, soap, identity, candidate)
    
    # remove mirror structure, they are subspecies
    logger.info('Removing mirror nodes')
    mirrors = removeMirrorNodes(graph, soap, identity, candidate)
    
    # after remove mirror nodes, tips will be genreated
    tips2 = removeTips(graph, soap, identity, candidate)
    for key, value in tips2.iteritems():
        if key not in tips:
            tips[key] = value
        else:
            tips[key].extend(value)
    
    # fourth, find bubbles
    bubbles, identity = findBubble(graph, soap)
    
    # merge bubbles
    flateBubble(graph, soap, bubbles, identity, candidate)
    
    # transform the bubbles
    indexedBubs = indexBubs(soap, bubbles)
    
    # find longtips, they are may be subspcies
    passedSP = branchCut(graph, soap, identity, candidate)
    
    # record the status of every node
    passedSP = set()
    nodeStatus = recordNodeStatus(graph, soap, tips, mirrors, bubbles, passedSP)
    with open(config.STATUS,'w') as af:
        for ele in nodeStatus:
            af.write(str(ele) + ' ' + str(nodeStatus[ele]) + '\n')
    
    # store the tips, bubbles
    with open(config.VARS, 'wb') as af:
        pickle.dump(mirrors, af) if mirrors else pickle.dump({}, af)
        pickle.dump(tips, af) if tips else pickle.dump({}, af)
        pickle.dump(passedSP, af) if passedSP else pickle.dump({}, af) 
        pickle.dump(indexedBubs, af) if bubbles else pickle.dump({}, af)

