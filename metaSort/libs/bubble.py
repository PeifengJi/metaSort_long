############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import re
import sys
import math
import config
import genemark
import indexGraph
import platform
import itertools
import subprocess
import collections
import numpy as np
#import scipy.optimize
import cPickle as pickle
from log import get_logger
# set log file
logger = get_logger()

class SNVs:
    """
    Stores the parameters for an alignment scoring function
    """
    def __init__(self, match = 2, mismatch = -1, gap = -5):
        
        self.gap = gap
        self.match = match
        self.mismatch = mismatch
    
    def __matchchar(self, a, b):
        """Return the score for aligning character a with b"""
        assert len(a) == len(b) == 1
        if a == b:
            return self.match
        else:
            return self.mismatch
    
    def localAlign(self, x, y):
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
        
        return matchx, matchy
    
    def call(self, x, y):
        """
        Caculate the similarity between two seqs
        """
        mx, my = self.localAlign(x, y)

        count, vars, visited = 0, {}, set()        
        for e1, e2 in itertools.izip(mx, my):
            if e1 != e2 and e1 != '-' and e2 != '-':
                vars[count] = '{0}/{1}'.format(e1, e2)
            elif e1 == '-':
                vars[count] = './{0}'.format(e2)
            elif e2 == '-':
                 vars[count] = '{0}/.'.format(e1)
            count += 1    

        # find the continuesly insertion
        keys = sorted([i for i, k in vars.iteritems() if k[0]=='.'])
        if not keys:
            return vars

        pre, tmp = keys[0], [[keys[0]]]
        for key in keys[1:]:
            if key - pre == 1:
                tmp[-1].append(key)
            else:
                tmp.append([key])
            pre = key

        for k in xrange(len(tmp)):
            delLen = len(tmp[k])
            if len(tmp[k]) == 1:
                vars[tmp[k][0]-1] = vars[tmp[k][0]]
                del vars[tmp[k][0]]
            else:
                vars[tmp[k][0]-1] = './'
                for ele in tmp[k]:
                    vars[tmp[k][0]-1] += vars[ele][-1]
                    del vars[ele]
            # change the index of other snv that larger than the region
            for key in sorted(vars.keys()):
                if key > tmp[k][-1]:
                    vars[key - delLen] = vars[key]
                    del vars[key]
            # renew the index list (tmp)
            newt = []
            for p, e in enumerate(tmp):
                if p >= k+ 1:
                    newt.append([i-delLen for i in e])
                else:
                    newt.append(e)
            tmp = newt

        # find the continusely deletion
        keys = sorted([i for i, k in vars.iteritems() if k[-1]=='.'])
        if not keys:
            return vars

        pre, tmp = keys[0], [[keys[0]]]
        for key in keys[1:]:
            if key - pre == 1:
                tmp[-1].append(key)
            else:
                tmp.append([key])
            pre = key
        tmp = [i for i in tmp if len(i) > 1]
        for pos in tmp:
            vars[pos[0]] = vars[pos[0]][0]
            for ele in pos[1:]:
                vars[pos[0]] += vars[ele][0]
                del vars[ele]
            vars[pos[0]] += '/.'

        return vars

class Pipeline():
    """
    bwa pipeline to map reads to reference to get the sam file
    """
    def __init__(self, ref):

        if platform.system() == 'Darwin':
            self.bwa = os.path.join(config.libsLocation, 'bwa/macosx/bwa-spades')
        else:
            self.bwa = os.path.join(config.libsLocation, 'bwa/linux/bwa-spades')
        
        self.ref = ref

    def isIndexed(self):
        """
        Checks to see if a given reference is indexed already
        """
        ref_ext = set(['.amb', '.ann', '.bwt', '.pac', '.sa'])
        ref_indexes = glob.glob(self.ref + '.*' )

        if not ref_indexes:
            return False

        ext = set([os.path.splitext(f)[1] for f in ref_indexes])
        intersec = ref_ext & ext
        # Return true only if the intersection of the found extensions is equal to
        # the expected ref_ext set
        return intersec == ref_ext

    def indexRef(self):
        """
        Index a given reference
        """
        if self.isIndexed():
            print('{0} is already indexed'.format(self.ref))
            return 1

        # index
        print('Indexing reference file ...')
        indexCMD = [self.bwa, 'index', self.ref]
        #don't print the stdout to the screen
        with open(os.devnull, "w") as f:
            returnCode = subprocess.call(indexCMD, stdout=f, stderr=f)

        if returnCode != 0:
            raise Exception('Bwa index reference file failed!')

    def algin(self, p1, p2, outPrefix, cpu):
        """
        Index reference and sai, sam files
        """
        # aln
        s1 = '{0}_1.sai'.format(outPrefix)
        s2 = '{0}_2.sai'.format(outPrefix)
        cmd1 = [self.bwa, 'aln', self.ref, '-t', cpu, '-f', s1, p1]
        cmd2 = [self.bwa, 'aln', self.ref, '-t', cpu, '-f', s2, p2]
        print('Creating sai files ...')

        # don't print the stdout to the screen
        with open(os.devnull, "w") as f:
            f1 = subprocess.Popen(cmd1, stdout=f, stderr=f)

        with open(os.devnull, "w") as f:
            f2 = subprocess.Popen(cmd2, stdout=f, stderr=f)

        f1.wait()
        f2.wait()
        outsam = outPrefix + '.sam'
        # generate sam file
        cmdSam = [self.bwa, 'sampe',
                  self.ref,
                  '-f', outsam,
                  s1, s2, p1, p2]

        with open(os.devnull, "w") as f:
            proc = subprocess.call(cmdSam)
            print('Sam file is generated!')

        self.__clean(s1, s2)

    def mem(self, p1, p2, outPrefix, cpu):

        cmdSam = [self.bwa, 'mem', '-t', str(cpu), self.ref, p1, p2]

        outfn = outPrefix + './sam'
        with open(outfn, 'w') as af:
            subprocess.call(cmdSam, stdout=af)

    def __clean(self, *files):
        """
        remove the temple files
        """
        for fn in files:
            os.remove(fn)

class ParseSAM:
    """
    This class provides fuctions to parsing/processing/transforming SAM format file
    """
    _reExpr1=re.compile(r'\s+') #>=1 spaces
    _reExpr2=re.compile(r'^\s*$') #blank line
    _splicedHit = re.compile(r'^(\d+)M(\d+)N(\d+)M$',re.IGNORECASE) #regular expression for spliced mapped reads
    _monoHit = re.compile(r'^(\d+)M$',re.IGNORECASE)                #regular expresion for Non-spliced mapped reads

    def __init__(self,samFile):

        self.fileName = os.path.basename(samFile)
        self.f=open(samFile,'r')

    def stat (self):
        """
        Calculate mapping statistics
        """

        _numTotalReads=0        #Total number of reads. both mapped and unmapped.
        _numMapReads=0          #Number of mapped reads
        _numUnmapReads=0        #Number of unmapped reads
        _numMultiReads=0        #Number of non-unique hit reads
        _numPlusReads=0         #Number of reads mapped to + strand. useful for strand specific sequencing.
        _numMinusReads=0        #Number of reads mapped to - strand. useful for
        _numOddReadsTotal=0     #number of reads that NOT properly paired
        _numMapReads1=0         #number of first reads mapped
        _numMapReads2=0         #number of second reads mapped
        #_numUnmapReads1=0      #number of first reads Unmapped
        #_numUnmapReads2=0      #number of second reads Unmapped
        _numPairedReads=0       #number of reads that properly paired
        _numPaired_PP=0         #both mapped to +
        _numPaired_PM=0         #each mapped to +/-
        _numPaired_MM=0         #both mapped to -
        _numSingle=0            #single hit. no mate
        _numMonoHit=0           #read not splicing mapped
        _numSplitHit=0          #read was splicing mapped
        
        for line in self.f:
            line=line.rstrip()
            if line[0] == '@':
                continue      #skip head lines
            elif ParseSAM._reExpr2.match(line):
                continue     #skip blank lines
            else:
                _numTotalReads = _numTotalReads+1
                field = line.split()
                flagCode=string.atoi(field[1])
                if (flagCode & 0x0004) != 0: #Query unmapped. count and skipped
                    _numUnmapReads += 1
                    continue
                else:  #Query mapped
                    if (flagCode & 0x0010 != 0):  #map to strand -
                        _numMinusReads += 1
                        _numMapReads +=1
                    if (flagCode & 0x0010 == 0):  #map to forward +
                        _numPlusReads += 1
                        _numMapReads +=1
                    if (flagCode & 0x0100)!= 0:
                        _numMultiReads += 1  #not primary hit
                    if (ParseSAM._splicedHit.match(field[5])):  #Splicing mappedreads
                        _numSplitHit += 1
                    else:  #mono mapped reads
                        _numMonoHit += 1
                    if (flagCode & 0x0001) != 0: #if this is paired sequencing
                        if (flagCode & 0x0040)!= 0:
                            _numMapReads1 += 1   #first read mapped
                        if (flagCode & 0x0080)!= 0:
                            _numMapReads2 += 1   #second read mapped
                        if (flagCode & 0x0002)!= 0:
                            _numPairedReads += 1  #mapped in proper pair
                        if (flagCode & 0x0008)!= 0: #onle one end could map
                            _numSingle +=1
                        else:     #both end could map
                            if(flagCode & 0x0010)!=0 and (flagCode & 0x0020)!=0 :  #both end mapped to -
                                _numPaired_MM +=1
                            if(flagCode & 0x0010)==0 and (flagCode & 0x0020)==0:   #both end mapped to +
                                _numPaired_PP +=1
                            if(flagCode & 0x0010)!=0 and (flagCode & 0x0020)==0 :  #first end -, second +
                                _numPaired_PM +=1
                            if(flagCode & 0x0010)==0 and (flagCode & 0x0020)!=0 :  #first end +, second -
                                _numPaired_PM +=1

        print >>sys.stderr,"\n#================================================="
        print >>sys.stderr,"#=================Report=========================="
        print >>sys.stderr, "Total Reads:\t",_numTotalReads
        print >>sys.stderr, "Total Mapped:\t",_numMapReads
        print >>sys.stderr, "Total UnMapped:\t",_numUnmapReads
        print >>sys.stderr, "Plus Strand Hits:\t", _numPlusReads
        print >>sys.stderr, "Minus Strand Hits:\t",_numMinusReads
        print >>sys.stderr, "Non-splicing Hits:\t",_numMonoHit
        print >>sys.stderr, "Splicing Hits:\t",_numSplitHit
        print >>sys.stderr, "Primary Hits:\t",_numTotalReads-_numMultiReads,"\n"
        print >>sys.stderr, "Mapped read-1:\t",_numMapReads1
        print >>sys.stderr, "Mapped read-2:\t",_numMapReads2
        print >>sys.stderr, "Proper Paired Hits:\t",_numPairedReads
        print >>sys.stderr, "Paired Hits (+/+):\t",_numPaired_PP
        print >>sys.stderr, "Paired Hits (-/-):\t",_numPaired_MM
        print >>sys.stderr, "Paired Hits (+/-):\t",_numPaired_PM
        print >>sys.stderr, "single end mapped:\t",_numSingle,"\n"
        print >>sys.stderr,"#================================================="
        self.f.seek(0)

    def samTobed(self,outfile=None):
        """
        Convert SAM file to BED file. BED file will be saved as xxx.sam.bed unless otherwise specified.
        reads spliced more than once were NOT recognized in current version
        """
        if outfile is None:
            outfile=self.fileName + ".bed"

        print >>sys.stderr,"Writing bed entries to\"",outfile,"\"...",
        FO=open(outfile,'w')
        for line in self.f:
            line=line.rstrip()
            field=line.split()
            if line[0] == '@':
                continue #skip head lines
            if ParseSAM._reExpr2.match(line):
                continue #skip blank lines
            if (string.atoi(field[1]) & 0x0004)!=0: 
                continue    #skip unmapped line
            flag = string.atoi(field[1])
            if (ParseSAM._splicedHit.match(field[5])):          #read is spliced mapped
                match1 = ParseSAM._splicedHit.match(field[5]).group(1)
                gap = ParseSAM._splicedHit.match(field[5]).group(2)
                match2 = ParseSAM._splicedHit.match(field[5]).group(3)
                chrom = field[2]
                chromStart = string.atoi(field[3])-1
                chromEnd = chromStart + string.atoi(match1) + string.atoi(gap) + string.atoi(match2)
                name = field[9]
                score = field[4]
                if(flag & 0x0010)==0:
                    strand = '+'
                else:
                    strand = '-'
                thickStart = chromStart
                thickEnd = chromEnd
                itemRgb = "0,255,0"
                blockCount = 2
                blockSizes = match1 + ',' + match2
                blockStarts = '0' +','+ str(string.atoi(match1)+string.atoi(gap))
            
            if(ParseSAM._monoHit.match(field[5])): #read mapped as whole
                match=ParseSAM._monoHit.match(field[5]).group(1)
                chrom = field[2]
                chromStart = string.atoi(field[3])-1
                chromEnd = chromStart + string.atoi(match)
                name = field[9]
                score = field[4]
                if(flag & 0x0010)==0:
                    strand = '+'
                else:
                    strand = '-'
                thickStart = chromStart
                thickEnd = chromEnd
                itemRgb = "0,0,0"
                blockCount = '1'
                blockSizes = match
                blockStarts = '0'
            print >>FO, string.join((str(i) for i in [chrom,chromStart,chromEnd,name,score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts]),sep="\t")

        print >>sys.stderr, "Done"
        self.f.seek(0)
        FO.close()

    def samTowig(self,outfile=None,log2scale=False,header=True):
        """
        Convert SAM file to wig file. WIG file will be saved as xxx.sam.wig unless otherwise specified.
        """
        if outfile is None:
            outfile = self.fileName + ".wig"
        
        FO=open(outfile,'w')
        print >>sys.stderr, "Writing wig file to\"",outfile,"\"..."

        _splicedHit = re.compile(r'^(\d+)M(\d+)N(\d+)M$',re.IGNORECASE)     #match CIGAR string if read is splicing mapped
        _monoHit = re.compile(r'^(\d+)M$',re.IGNORECASE)                    #match CIGAR string if read is Non-splicing mapped
        headline="track type=wiggle_0 name=" + outfile + " track_label description='' visibility=full color=255,0,0"
        wig=collections.defaultdict(dict)
        for line in self.f:
            hits=[]
            if line[0] == '@':
                continue #skip head lines
            if ParseSAM._reExpr2.match(line):
                continue       #skip blank lines
            field=line.rstrip().split()
            if (string.atoi(field[1]) & 0x0004)!=0: 
                continue    #skip unmapped line
            
            if(_monoHit.match(field[5])):                   #if this is a single match. i.e. no splicing match
                match=_monoHit.match(field[5]).group(1)
                hits=range(string.atoi(field[3]),string.atoi(field[3])+string.atoi(match))
            
            if (_splicedHit.match(field[5])):
                match1 = _splicedHit.match(field[5]).group(1)
                gap = _splicedHit.match(field[5]).group(2)
                match2 = _splicedHit.match(field[5]).group(3)
                hits=range(string.atoi(field[3]),string.atoi(field[3])+string.atoi(match1))
                hits.extend(range(string.atoi(field[3])+string.atoi(match1)+string.atoi(gap), string.atoi(field[3])+string.atoi(match1)+string.atoi(gap)+string.atoi(match2)))
            for i in hits:
                if wig[field[2]].has_key(i):
                    wig[field[2]][i] +=1
                else:
                    wig[field[2]][i]=1

        if header:
            FO.write(headline + "\n")

        for chr in sorted(wig.keys()):
            print >>sys.stderr, "Writing ",chr, " ..."
            FO.write('variableStep chrom='+chr+'\n')
            for coord in sorted(wig[chr]):
                if log2scale:
                    FO.write("%d\t%5.3f\n" % (coord,math.log(wig[chr][coord],2)))
                else:
                    FO.write("%d\t%d\n" % (coord,wig[chr][coord]))
        self.f.seek(0)
        FO.close()

    def getUnmap(self, outfile=None,fastq=True):
        """
        Extract unmapped reads from SAM file and write to fastq [when fastq=True] or fasta [when fastq=False] file
        """

        if outfile is None:
            if fastq: 
                outfile = self.fileName + ".unmap.fq"
            else: 
                outfile = self.fileName + ".unmap.fa"
        
        FO=open(outfile,'w')
        unmapCount=0
        print >>sys.stderr, "Writing unmapped reads to\"",outfile,"\"... ",

        for line in self.f:
            hits=[]
            if line[0] == '@':
                continue #skip head lines
            if ParseSAM._reExpr2.match(line):
                continue #skip blank lines
            field=line.rstrip().split()
            flagCode = string.atoi(field[1])
            seq=field[9]
            qual=field[10]
            if (flagCode & 0x0004) != 0: #read unmap
                unmapCount +=1
                if (flagCode & 0x0001) != 0: #paried in sequencing
                    if (flagCode & 0x0040)!=0:
                        seqID=field[0] + '/1' #first read
                    if (flagCode & 0x0080)!=0:
                        seqID=field[0] + '/2' #second read
                else: 
                    seqID=field[0]

                if fastq: 
                    FO.write('@' + seqID + '\n' + seq +'\n' + '+' +'\n' + qual+'\n')
                else: 
                    FO.write('>' + seqID + '\n' + seq +'\n')

        print >>sys.stderr, str(unmapCount) + " reads saved!\n"
        FO.close()
        self.f.seek(0)

    def getProperPair(self,outfile=None):
        """
        Extract proper paried mapped reads.
        """
        if outfile is None:
            outfile = self.fileName + ".PP.sam"
        
        FO=open(outfile,'w')
        PPcount=0
        print >>sys.stderr, "Writing proper paired reads to\"",outfile,"\"... ",
        
        for line in self.f:
            hits=[]
            if line[0] == '@':
                continue  #skip head lines
            if ParseSAM._reExpr2.match(line):
                continue  #skip blank lines
            field=line.rstrip('\n').split()
            flagCode=string.atoi(field[1])
            if ((flagCode & 0x0001) != 0) and ((flagCode & 0x0002)!=0):
                PPcount +=1
                FO.write(line)
        
        FO.close()
        print >>sys.stderr, str(PPcount) + " reads were saved!\n",
        self.f.seek(0)

    def collapseSAM(self, outfile=None,collevel=10):
        """
        At most collevel[default=10] identical reads will be retained in outputting SAM file
        The original SAM file must be sorted before hand. if not, using linux command like (sort -k3,3 -k4,4n myfile.sam >myfile.sorted.sam)
        """
        if outfile is None:
            outfile = self.fileName + ".collapsed.sam"
        
        print >>sys.stderr, "Writing collapsed SAM file to\"",outfile,"\"... "
        FO=open(outfile,'w')
        flag=""
        
        for line in self.f:
            if line[0] == '@':
                continue #skip head lines
            if ParseSAM._reExpr2.match(line):
                continue #skip blank lines
            field=line.rstrip('\n').split()
            if (string.atoi(field[1]) & 0x0004)!=0: 
                continue #skip unmapped line
            id=field[2] + field[3] + field[5]
            if (id != flag):
                FO.write(line)
                flag=id
                skipTrigger=0
            else:
                skipTrigger+=1
                if skipTrigger < collevel:
                    FO.write(line)
                else:
                    continue
        FO.close()
        self.f.seek(0)

def bubFeatureCal(graph, soap, indexedBubs, align = False):
    """
    Caculate the density of bubbles containing graphs, and average bubbles distance
    """
    # transfer the formate the bubbles into (start, end) formate
    bub1start, bub1end, cl = set(), set(), {}
    for k, v in indexedBubs.iteritems():
        bub1start.add(v[0][0])
        cl[v[0][0]] = k
        bub1end.add(v[0][-1])
    
    # get the subgraphs
    density, distance, lenRate = [], [], []
    group = graph.findSubgraph(out='gp')
    for key in group:
        inter = set(group[key]) & bub1start
        if inter:
            # get the inner nodes of the bubbles on this subgraph
            innerNodes = set([n for i in inter for p in indexedBubs[cl[i]] for n in p[1:-1]])
            density.append(round(len(innerNodes)*2.0/len(group[key]), 2))
            bubSumLen = sum([soap.length[i] for i in innerNodes])
            totalLen = sum([soap.length[i] for i in group[key]])
            lenRate.append(round(bubSumLen*1.0/totalLen, 3))
            # calculate the distance between bubbles
            if len(inter) > 1:
                # use a temple graph
                tmpGraph = indexGraph.Graph(soap.index)
                for ele in group[key]:
                    if ele in graph:
                        for sec in graph[ele]:
                            tmpGraph.addEdge(ele, sec, 1)
                # delete the bubble inner node and get the subgraphs
                singlts = tmpGraph.deleteNodes(innerNodes)
                if singlts:
                    for x in singlts:
                        # if x in the inter node of two continously bubbles
                        if x in bub1start and x in bub1end:
                            distance.append(soap.length[x])
                # transfromt the reduced graph into a bidirected formate
                tmpGroup = tmpGraph.findSubgraph(out='gp')
                biGraph = {}
                for ele in tmpGraph.iterkeys():
                    if ele not in biGraph: 
                        biGraph[ele] = {}
                    for sec in tmpGraph[ele]:
                        biGraph[ele][sec] = 1
                        if sec not in biGraph: biGraph[sec] = {}
                        biGraph[sec][ele] = 1
                # get the far distance between two nodes, there is no way that
                # a bub start locates on the same subgraph of the bub end of a same bub
                for value in tmpGroup.itervalues():
                    sPoll, ePoll = set(value) & bub1start, set(value) & bub1end
                    if sPoll and ePoll:
                        for node in ePoll:
                            queue = collections.deque([[node]])
                            visitedNodes, flag = set(), 0
                            while queue:
                                tmp = queue.pop()
                                last = tmp[-1]
                                if last in visitedNodes:
                                    continue
                                visitedNodes.add(last)
                                for sec in biGraph[last]:
                                    longer = tmp + [sec]
                                    if sec in sPoll or soap.revNode(sec) in sPoll:
                                        distance.append(sum([soap.length[i] for i in longer]))
                                        sPoll.discard(sec)
                                        sPoll.discard(soap.revNode(sec))
                                        if not sPoll:
                                            flag = 1
                                            break
                                        else:
                                            queue.appendleft(longer)
                                    else:
                                        queue.appendleft(longer)
                                if flag==1: 
                                    break
    
    if not align: return lenRate, density, distance

    # calculate the similarity
    bubNodes = set([k for p in indexedBubs.itervalues() for i in p for k in i[1:-1]])
    graphSeqs = soap.fetchSeqs(bubNodes)
    identity = []
    for value in indexedBubs.itervalues():
        tmpSeqs = []
        for path in value:
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
            tmpSeqs.append(seq)
        
        # compare the similarity among the sequences
        sim = similarity(soap.kmersize)
        bubidens, iden = [], 0
        for i in xrange(len(tmpSeqs)-1):
            for k in xrange(i+1, len(tmpSeqs)):
                try:
                    iden = sim.identity(tmpSeqs[i], tmpSeqs[k])
                except:
                    logger.error('Error happends when align sequences',
                                    to_stderr=True,
                                    exit_with_code=3
                        )
                if iden:
                    bubidens.append(iden)
        
        if bubidens:
            identity.extend(bubidens)

    return lenRate, density, distance, identity

def allelCov(soap, allelfn, f1, f2, cpu):
    """
    use bwa mapping reads to allels and calculate the coverage of the allels\
    """
    # mapping by using bwa mem
    bwa = Pipeline(allelfn)
    bwa.mem(allelfn, f1, f2, 'result', cpu)
    
    # convert sam to bam
    cmdConvert = ['samtools', 'import', 'result.sam', 'result.bam']
    subprocess.call(cmdConvert)
    os.remove('result.sam')
    
    # using samtools add the MD tag
    cmdMD = ['samtools', 'calmd', 'result.bam', allelfn]
    with open('md.sam', 'w') as af:
        subprocess.call(cmdMD, stdout=af)
    os.remove('result.bam')

    # filter the unmapped and multiple mapped reads
    counter = collectons.defaultdict(int)
    with open('md.dam') as af:
        for line in af:
            if line[0] != '@':
                break
        for line in af:
            ae = line.rstrip().split()
            if int(ae[5]) > 0:
                if 'A' not in ae[-1] and 'T' not in ae[-1] and 'C' not in ae[-1] and 'G' not in ae[-1] :
                    counter[ae[2]] += int(ae[-1].split(':')[-1])
    
    # calculate coverage of the allels
    length = {}
    with open(allelfn, r) as af:
        for line in af:
            if line[0] == '>':
                tid = line[1:-1]
            else:
                length[tid] = len(line[:-1])

    cov = {}
    for key, value in counter.iteritems():
        cov[key] = round(value*1.0/length(key),2)

    return cov

def bubSnvs(soap, indexedBubs, outfn):
    """
    Transfomr the format of bubbles using increasing number
    """
    #link the inner nodes of bubbles
    bubSeqs = soap.fetchSeqs([v for i in indexedBubs.itervalues() for k in i for v in k])
    out, sim = open(outfn, 'w'), SNVs()
    snvs = {}
    for k, v in indexedBubs.iteritems():
        seqs, snvs[k] = [], {}
        for p in v:
            seq = bubSeqs[p[1]]
            for ele in p[1:-1]:
                seq += bubSeqs[ele][soap.kmersize:]
            seqs.append(seq)
        for p, s in enumerate(seqs):
            out.write('>{0}_{1}\n{2}\n'.format(k, p, s))
        # take the first as reference
        for x in xrange(1, len(seqs)):
            snvs[k][x] = sim.call(seqs[0], seqs[x])

    return indexedBubs, snvs

def allelAbundance(soap, indexedBubs, depth=False):
    """
    use nnls to estimate the abundance of every path in the bubble
    """
    logger.info('Estimating abundance of the variations ...')

    abundance = {}
    for key in indexedBubs:
        paths = indexedBubs[key]
        uniqNodes = list(set([i for k in paths for i in k]))
        
        # initialize respone vector
        if depth:
            try:
                y =  np.array([depth[i] for i in uniqNodes])
            except:
                # there must be node that no sam mapper coverage
                continue

        else:
            y =  np.array([soap.cov[i] for i in uniqNodes])
        
        # initialize design matrix
        A =  np.zeros((len(uniqNodes), len(paths)))
        position = 0
        for node in uniqNodes:
            for p, v in enumerate(paths):
                if node in v:
                    A[position][p] = 1
            position += 1
        
        # initialize the weighted matrix
        weights = [1-math.exp(-soap.length[i]*0.5/100) for i in uniqNodes]
        # http://en.wikipedia.org/wiki/Linear_least_squares_%28mathematics%29#Weighted_linear_least_squares
        Aw = A*np.sqrt(weights)[:, None]
        Yw = y*np.sqrt(weights)
        # use non-negative linear least squares method to get the optimized parasmeters
        #  finds path weights that minimize the sum of the squared residuals
        # min||WAx - Wb||2 s.t.x > 0
        abundance[key] = scipy.optimize.nnls(Aw, Yw)[0]
        #abundance[key] = scipy.optimize.nnls(A, y)[0]

    logger.info('Done!')
    return abundance

def genePrediction(infile):
    """
    Use MetaGeneMark to predict the orf of every contig, and get the codon usage of every contig
    """
    # get the absolut path of the predict files
    genemarkPath = genemark.getPath()
    gmhmmp = genemarkPath + '/gmhmmp'
    modfile = '/'.join(genemarkPath.split('/')[:-1]) + '/MetaGeneMark_v1.mod'
    # predict
    protein = {}
    seq = ''
    ids = collections.deque()
    cmdMetaGeneMark = [gmhmmp, '-d', '-f', 'G', '-m', modfile, '-o', 'refGenes.lst', infile]
    logger.info('Predicting proteins of the fasta file: ' + infile)
    logger.print_command_line(cmdMetaGeneMark, wrap_after=None)
    
    child = subprocess.call(cmdMetaGeneMark)
    if child != 0:
        sys.stderr.write(out)
        logger.error('Error happened while using MetaGeneMark', exit_with_code=1)
    else:
        # gather the line info
        genes = {}
        with open('refGenes.lst','r') as af:
            for line in af:
                if line[0] == '>':
                    ae = line.rstrip().split('|')
                    gid, s = int(ae[-2]), ae[0][1:]
                    nextae = ae[-1].split()
                    t, refid = int(nextae[0]), nextae[1][1:]
                    if refid not in genes:
                        genes[refid] = {}
                    genes[refid][gid] = [s, t]
    return genes

def locateVars(soap, indexedBubs, bubfn, refn, gffn=False):
    """
    Give the positio of the bubbles on the reference
    """
    # algin the bubble nodes with reference
    if platform.system() == 'Darwin':
        mummerPath = os.path.join(config.libsLocation, 'MUMmer3.23-osx')
    else:
        mummerPath = os.path.join(config.libsLocation, 'MUMmer3.23-linux')
    nucPath = os.path.join(mummerPath, 'nucmer')
    snpPath = os.path.join(mummerPath, 'show-snp')
    
    cmdAlgn = [nucPath, '-c', '65', '-l', '65', '--maxmatch', '-p', 'var', refn, bubfn]
    subprocess.call(cmdAlgn)
    
    cmdSnp = []

    # extract the position of the bubble nodes
    position, varOnGenome = {}, {}
    with open('varNucmer.result','r') as af:
        for i in xrange(4):
            _ = af.readline()
        for line in af:
            ae = line.rstrip().split()
            if float(ae[10]) >= 95:
                if ae[-1] not in varOnGenome:
                    varOnGenome[int(ae[-1])] = set()
                if int(ae[-1]) not in position:
                    position[int(ae[-1])] = []
                position[int(ae[-1])].append([int(ae[0]), int(ae[1])])
                varOnGenome[int(ae[-1])].add(ae[-2])

    # identity the uniq mapping position of multiple mapping nodes
    passed, passedBub = set(), set()
    for key in status:
        uniqPos, multiNodes = [], set()
        for v in status[key]:
            if v in position and v not in passed and len(position[v]) > 1:
                multiNodes.add(v)
            elif v in position and len(position[v]) == 1:
                uniqPos.extend(position[v][0])
                position[v] = position[v][0]
                passed.add(v)
        
        if multiNodes:
            if uniqPos:
                uniqPos.sort()
                for node in multiNodes:
                    sp, sv = 1000, 10000000000
                    for p, v in enumerate(position[node]):
                        meanPos = sum(v)/2
                        meanUniq = sum(uniqPos)/2
                        dis = abs(meanUniq - meanPos)
                        if dis < sv:
                            sv, sp =  dis, p
                    position[node] = position[node][sp]
                    passed.add(node)
                passedBub.add(key)
        else:
            passedBub.add(key)

    # the left node bubbles must locate on multiple positon of the genome 
    # they may be repeatative sequences, do not use them

    # bubbles on genome
    bubOnGenome = {}
    for k in passedBub:
        tmp = collections.Counter([v for i in status[k] if i in varOnGenome for v in varOnGenome[i]])
        # take the mode number
        if tmp:
            bubOnGenome[k] = tmp.most_common(1)[0][0]

    # record the position and status of each node
    with open('location', 'w') as af:
        for key in sorted(passedBub):
            for p, ele in enumerate(status[key]):
                note = 'm'
                if p == 0:
                    note = 's'
                elif p == len(status[key]) -1:
                    note = 'e'
                if ele in position:
                    # node, cluster number, positon on the bubble, reference sequence id , start on the reference, end on the reference
                    af.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(ele, key, note, bubOnGenome[key], position[ele][0], position[ele][1]))

    # annotate the sequences to get the codon region
    if gffn:
        scafs = {}
        with open(gffn) as af:
            for line in af:
                if line[0] == '#':
                    continue
                ae = line.split('\t')
                if ae[2] == 'CDS' or ae[2] == 'gene':
                    s, e = int(ae[3]), int(ae[4])
                    if e-s >= 500:
                        if 'Genbank:' in ae[-1]:
                            nextae = ae[-1].split(';')
                            if ae[0] not in scafs:
                                scafs[ae[0]] = {}
                            for ele in nextae:
                                if 'Genbank:' in ele:
                                    gid = ele.split(':')[1]
                                    scafs[ae[0]][gid] = [s, e]
    else:
        scafs = genePrediction(refn)

    # statistic info (the bubble number on this gene and the id of the bubbles)
    identifer, parser = scafs.keys()[0], 0
    if '|' not in identifer: parser = 1

    bubinfo, bubNote = {}, {}
    with open('location') as af:
        for line in af:
            ae = line.rstrip().split()
            #  0           1                  2                     3                             4                    5   
            # node, cluster number, structure on the bubble, reference sequence id , start on the reference, end on the reference
            newScafId = ae[3].split('|')[-2] if parser else ae[3]
            bubNote[int(ae[0])] = ae[2]
            if newScafId not in bubinfo:
                bubinfo[newScafId] = {}
            if int(ae[1]) not in bubinfo[newScafId]:
                bubinfo[newScafId][int(ae[1])] = {}
            if int(ae[0]) not in bubinfo[newScafId][int(ae[1])]:
                bubinfo[newScafId][int(ae[1])][int(ae[0])] = [int(i) for i in ae[4:]]

    geneBubs, visited = {}, set()
    for scaf in scafs:
        if scaf in bubinfo:
            for gene in scafs[scaf]:
                geneRange = scafs[scaf][gene]
                for cl in bubinfo[scaf]:
                    if cl in visited:
                        continue
                    for n, v in bubinfo[scaf][cl].iteritems():
                        if geneRange[0] <= v[0] and v[1] <= geneRange[1]:
                            if scaf not in geneBubs:
                                geneBubs[scaf] = {}
                            if gene not in geneBubs[scaf]:
                                geneBubs[scaf][gene] = {}
                            if cl not in geneBubs[scaf][gene]:
                                geneBubs[scaf][gene][cl] = []
                            geneBubs[scaf][gene][cl].append(n)
                            visited.add(cl)
    # delete the clsuter that noly one node mapped on the gene
    for scaf in geneBubs.keys():
        for gene in geneBubs[scaf].keys():
            for cl, nodes in geneBubs[scaf][gene].items():
                if len(nodes) == 1:
                    del geneBubs[scaf][gene][cl]
            if len(geneBubs[scaf][gene].keys()) == 0:
                del geneBubs[scaf][gene]

    # bubble position on the gene
    bubPosition = {}
    for scaf in geneBubs:
        if scaf not in bubPosition:
            bubPosition[scaf] = {}
        for gene in geneBubs[scaf]:
            if gene not in bubPosition[scaf]:
                bubPosition[scaf][gene] = {}
            for cl in geneBubs[scaf][gene]:
                tmp = []
                for node in geneBubs[scaf][gene][cl]:
                    #if bubNote[node] == 'm':
                    tmp.extend(bubinfo[scaf][cl][node])
                tmp.sort()
                bubPosition[scaf][gene][cl] = [tmp[0], tmp[-1]]

    # write the bubble position file 
    with open('bubPosition','w') as af:
        for scaf in bubPosition:
            for gene in bubPosition[scaf]:
                for node, v in bubPosition[scaf][gene].iteritems():
                    af.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(scaf, gene, node, v[0], v[1]))

    # store variables
    with open('geneBubs','wb') as af:
        pickle.dump(indexedBubs, af)
        pickle.dump(cls, af)
        pickle.dump(bubinfo, af)
        pickle.dump(bubPosition, af)
        pickle.dump(geneBubs, af)

def alphaDiversity(cov, outfn=False):
    """
    Calculate the alpha diverity of the bubbles
    formula: H = - sum( (Pi * ln(Pi)))
    """
    # transform the coverage
    cl = {}
    for key, value in cov.iteritems():
        ae = key.split('_')
        if int(ae[0]) not in cl:
            cl[int(ae[0])] = []
        cl[int(ae[0])].append(value)

    alpha = {}
    for key, value in cl.iteritems():
        # add filter, no zero and there must one allel coverage larger than 3
        if len(value) > 1 and 0 not in value:
            for ele in value:
                if ele >= 5:
                    break
            else:
                continue
            H,sumAbun = 0, sum(value)
            for i in value:
                pi = i*1.0/sumAbun
                if i > 0 and sumAbun > 0:
                    H += (pi) * math.log(pi)
            alpha[key] = round(-H, 2)

    if outfn:
        with open(outfn, 'w') as af:
            for ele in alpha:
                af.write('{0}\t{1:.2}\n'.format(ele, alpha[ele]))

    return alpha

