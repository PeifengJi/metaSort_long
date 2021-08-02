############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import collections
import subprocess
import config
import platform
import copy
import re
import glob
import operator
import tempfile
import shutil
import itertools
import createPelinks
from contextlib import contextmanager
import indexGraph
import cPickle as pickle
from log import get_logger
# set log file
logger = get_logger()

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
            logger.info('{0} is already indexed'.format(self.ref))
            return 1
        
        # index
        logger.info('Indexing reference file ...')
        indexCMD = [self.bwa, 'index', self.ref]
        #don't print the stdout to the screen
        with open(os.devnull, "w") as f:
            returnCode = subprocess.call(indexCMD, stdout=f, stderr=f)
        
        if returnCode != 0:
            logger.error('Bwa index reference file failed!',
                    to_stderr=True,
                    exit_with_code=3
                    )
    
    def algin(self, p1, p2, outPrefix, cpu):
        """
        Index reference and sai, sam files
        """
        # aln
        s1 = '{0}_1.sai'.format(outPrefix)
        s2 = '{0}_2.sai'.format(outPrefix)
        cmd1 = [self.bwa, 'aln', self.ref, '-t', cpu, '-f', s1, p1]
        cmd2 = [self.bwa, 'aln', self.ref, '-t', cpu, '-f', s2, p2]
        logger.info(' '.join(cmd1))
        logger.info(' '.join(cmd2))
        logger.info('Creating sai files ...')
        
        # don't print the stdout to the screen
        with open(os.devnull, "w") as f:
            f1 = subprocess.Popen(cmd1, stdout=f, stderr=f)
        
        with open(os.devnull, "w") as f:
            f2 = subprocess.Popen(cmd2, stdout=f, stderr=f)
        
        f1.wait()
        f2.wait()
        # generate sam file
        logger.info('Creating sam files ...')
        cmdSam = [self.bwa, 'sampe',
                  '-f', '{0}.sam'.format(outPrefix),
                  self.ref, 
                  s1, s2, p1, p2]
        
        with open(os.devnull, "w") as f:
            proc = subprocess.call(cmdSam, stdout=f, stderr=f)
        
        logger.info('Sam file is generated!')
        
        self.__clean(s1, s2)
        self.__clean(self.ref+'.amb', self.ref+'.ann', self.ref+'.bwt', self.ref+'.pac', self.ref+'.sa')
    
    def __clean(self, *files):
        """
        remove the temple files
        """
        for fn in files:
            os.remove(fn)
    
    def bam(self):
        """
        generate sorted bam file
        """
        logger.info('Generating sorted bam file ...')
        cmd = ['samtools', 'faidx', self.ref]
        p1 = subprocess.call(cmd)
        samfn = outPrefix + '.sam'
        b1 = outPrefix + '.bam'
        cmd = ['samtools', 'import', self.ref, samfn, b1]
        p2 = subprocess.call(cmd)
        b2 = 'sorted_' + outPrefix 
        cmd = ['samtools', 'sort', b1, b2]
        p3 = subprocess.call(cmd)
        self.__clean(b1)
        return b2

class similarity2:
    """
    Stores the parameters for an alignment scoring function
    """
    def __init__(self, kmersize = 17, simCutOff = 0.8, match = 2, mismatch = -1, gap = -5):
        
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
        check if exisits bl2seq, if yes use it algin seqs
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

class similarity:
    """
    Stores the parameters for an alignment scoring function
    """
    def __init__(self, simCutOff = 0.8, match = 2, mismatch = -1, gap = -5):
        self.simCutOff = simCutOff
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
    
    def identity(self, x, y):
        """
        Caculate the similarity between two seqs
        """
        matchNum, lencut = self.__localAlign(x, y)
        return lencut, round(matchNum*1.0/lencut, 2)

class SOAPtransformer():
    
    def __init__(self, metaO, library, threads):
        
        self.metaO = metaO
        self.library = library
        self.tmpContig = 'tmp.fa'
        self.transed_metaO = 'g.contig'
        self.BINFILE = 'Graphs'
        self.cpu = threads
        self.absPath = os.getcwd()
        self.libfn = 'libSetting'
    
    def revCom(self, string):
        """
        reverse complement the given string
        """
        comp = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C',
                'n': 'N', 'N' : 'n',
                'a' : 'T', 't' : 'A', 'c' : 'G', 'g' : 'C'
                }
        result = ''
        for i in string[::-1]:
            result += comp[i]
        
        return result
    
    def palindromic(self, seq, clen):
        """
        Check if the seq is a palindromic seq, criterion:
        1: length if multiple of 2;
        2: the half of the seq is same as the reverse complementary of the other half seq
        """
        # criterion 1
        if clen%2:
            return 0
        else:
            # criterion 2
            fir10bp, last10bp = seq[:10], self.revCom(seq[-10:])
            half = clen/2
            midFir10bp, midLast10bp = seq[half-10: half], self.revCom(seq[half:half+10])
            if fir10bp == last10bp and midFir10bp == midLast10bp:
                return 1
            else:
                return 0
    
    def transFasta(self, infn):
        """
        transform the splited lines of seq into one line and build index
        """
        output = open(self.tmpContig, 'w')
        header, seq, index, length = 1, '', {}, {}
        with open(infn, 'r') as af:
            for line in af:
                if line[0] == '>':
                    if seq:
                        seqid, clen = header, len(seq)
                        flag = self.palindromic(seq, clen)
                        if not flag: # not self reverse complementary
                            index[seqid], length[seqid]  = 0, clen
                            header += 2
                        else:        # self reverse complementary
                            index[seqid], length[seqid] = 1, clen
                            header += 1
                        output.write('>{0}\n{1}\n'.format(seqid, seq))
                        seq = ''
                else:
                    seq += line.rstrip()
            else:
                seqid, clen = header, len(seq)
                flag = self.palindromic(seq, clen)
                if not flag: 
                    index[seqid], length[seqid]  = 0, clen
                else:
                    index[seqid], length[seqid] = 1, clen
                output.write('>{0}\n{1}\n'.format(seqid, seq))

        output.close()
        return index, length
    
    def readLib(self):
        """
        Read the libConfig file and output the lib info
        """
        libs = {}
        with open(self.library, 'r') as af:
            for line in af:
                if line[0] == '#':
                    continue
                if line[:3] == 'LIB':
                    libID = int(line.rstrip().split('=')[-1])
                    libs[libID] = {'ins': 0,
                                        'len' : 0,
                                        'reads' : []
                                        }
                    counter = 1
                elif line[:3] == 'rea':
                    libs[libID]['len'] = int(line.split('=')[-1].strip())
                elif line[:3] == 'ins':
                    libs[libID]['ins'] = int(line.split('=')[-1].strip())
                elif line[:3] == 'pai':
                    # must be the lib file
                    if counter%2:
                        libs[libID]['reads'].append([])
                    fn = line.rstrip().split('=')[1].strip()
                    libs[libID]['reads'][-1].append(fn)
                    counter += 1
        
        logger.info('There are {0} libraries'.format(len(libs)))
        for key in libs:
            logger.info('  Lib id: {0}, read length: {1}, ins_len: {2}, files:'.format(key, libs[key]['len'], libs[key]['ins']))
            for fn in libs[key]['reads']:
                logger.info('  {0}, {1}'.format(fn[0], fn[1]))
        
        return libs
    
    def checkFiles(self, *files):
        """
        check if the file is ready for using
        """
        for fn in files:
            if not os.path.getsize(fn) or not os.path.exists(fn):
                logger.error('File: {0} error, please check it'.format(fn),
                             to_stderr=True,
                             exit_with_code=1
                             )
    
    def bwaMapping(self):
        """
        Get the coverage of the contigs from sam file
        """
        sam = Pipeline(self.transed_metaO )
        sam.indexRef()
        for key in self.libs:
            prefix = 'lib{0}'.format(key)
            for pair in self.libs[key]['reads']:
                logger.info('Mapping files {0} and {1}'.format(pair[0], pair[1]))
                sam.algin(pair[0], pair[1], prefix, self.cpu)
        logger.info('Done!')
    
    def dlinkCreation(self):
        """
        Finding dlinks by finding the mode of the overlap
        """
        # store the 100bp on the two end of the sequence
        ids = collections.deque()
        seqs, seq = {}, ''
        with open(self.tmpContig , 'r') as af:
            for line in af:
                if line[0] == '>':
                    ids.append(int(line.rstrip()[1:]))
                    if seq:
                        cid = ids.popleft()
                        if len(seq) < 200:
                            seqs[cid] = seq
                            if self.index[cid] != 1:
                                seqs[cid+1] = self.revCom(seq)
                        else:
                            seqs[cid] = seq[:100] + seq[-100:]
                            if self.index[cid] != 1:
                                revSeq = self.revCom(seq)
                                seqs[cid+1] = revSeq[:100] + revSeq[-100:]
                        seq = ''
                else:
                    seq = line.rstrip()
            else:
                cid = ids.popleft()
                if len(seq) < 200:
                    seqs[cid] = seq
                    if self.index[cid] != 1:
                        seqs[cid+1] = self.revCom(seq)
                else:
                    seqs[cid] = seq[:100] + seq[-100:]
                    if self.index[cid] != 1:
                        revSeq = self.revCom(seq)
                        seqs[cid+1] = revSeq[:100] + revSeq[-100:]
        
        # pair algin to find the best kmersize and construct dlink graph
        links, ks = [], []
        nodes = set(seqs.iterkeys())
        for size in xrange(100, 40, -1):
            kmerLeft, kmerRight = {}, {}
            for node in nodes:
                if len(seqs[node]) > size:
                    if seqs[node][:size] not in kmerLeft:
                        kmerLeft[seqs[node][:size]] = []
                    kmerLeft[seqs[node][:size]].append(node)
                    
                    if seqs[node][-size:] not in kmerRight:
                        kmerRight[seqs[node][-size:]] = []
                    kmerRight[seqs[node][-size:]].append(node)
            
            # key is kmer string
            for key in kmerRight:
                if key in kmerLeft:
                    for up in kmerRight[key]:
                        for dow in kmerLeft[key]:
                            links.append([up, dow, -size])
                            ks.append(size)
                            nodes.discard(dow)
                        nodes.discard(up)
        
        # the mode should be the kmersize
        kmersize = collections.Counter(ks).most_common()[0][0]
        logger.info('Kmer size is: {0}'.format(kmersize))
        logger.info('Building dlink graph ...')
        
        dDlink = indexGraph.Graph(self.index)
        for link in links:
            if link[-1] == -kmersize:
                dDlink.addEdge(link[0], link[1], link[2])
        
        return kmersize, dDlink
    
    def revNode(self, node):
        """
        Get the reverse complementary node
        """
        if self.index[node] == 0:
            return node + 1
        elif self.index[node] == 2:
            return node - 1
        else:
            return node
    
    def scanKmeroverlap(self, seqs, link, compare):
        """
        Scan two seqs, find if they have a perfect overlap
        """
        # link: (205027, -1, 271959, 1)
        kmer = self.kmersize if self.kmersize > 0 else -kmersize
        direct =  (link[1], link[3])
        # find best link
        overlap = [[]]
        s1 = seqs[link[0]] if direct[0] == 1 else self.revCom(seqs[link[0]])
        s2 = seqs[link[2]] if direct[1] == -1 else self.revCom(seqs[link[2]])
        overlap[-1].extend(compare.identity(s1[-kmer:], s2[:kmer]))
        overlap[-1].extend(direct)
        
        revDir = (-link[1], -link[3])
        s1 = seqs[link[0]] if revDir[0] == 1 else self.revCom(seqs[link[0]])
        s2 = seqs[link[2]] if revDir[1] == -1 else self.revCom(seqs[link[2]])
        overlap.append([])
        overlap[-1].extend(compare.identity(s1[-kmer:], s2[:kmer]))
        overlap[-1].extend(revDir)
        # sort the result, and find the best one
        overlap.sort(key= lambda x: operator.itemgetter(1, 0))
        
        return overlap[0][2:]
    
    def plink(self):
        """
        Generate pelinks based on read mapping result
        """
        logger.info('Generating pelink and estimating gap size ...')
        pelink, contigCov = {}, {}
        for key in self.libs:
            logger.info('Parsing libs {} ...'.format(key))
            libSPlink, miu, sigma, cov = createPelinks.pelinksFromSam(self.samFiles[key], key, self.length,
                                    self.kmersize, self.libs[key]['ins'], self.libs[key]['len'])
            self.libs[key]['miu'], self.libs[key]['sigma'], contigCov = miu, sigma, cov
            # there are situations that more than one peaks in pair linked two links, select the best one
            # also there are reduant in the pelink, redundancy, like 13:{1:15:-1:[xxx], -1:{15:1:[yyy]}}, one is error
            dis = set()
            for c1 in libSPlink:
                newarr = set()
                for d1 in libSPlink[c1]:
                    for c2 in libSPlink[c1][d1]:
                        for d2 in libSPlink[c1][d1][c2]:
                            arr = libSPlink[c1][d1][c2][d2]
                            newarr.add((c1, d1, c2, d2))
                            # find multiple mapping peaks
                            if len(arr) > 1: 
                                # use dlink for selection
                                paths = self.dDlink.allPaths(c1, c2)
                                if paths: # there are dlink path
                                    gaps = []
                                    for path in paths:
                                        if len(path) == 2:
                                            gaps.append(self.dDlink[path[0]][path[1]])
                                        else:
                                            gap = self.dDlink[path[0]][path[1]]
                                            for pos in xrange(1, len(path)-1):
                                                gap += self.length[path[pos]] + self.dDlink[path[pos]][path[pos+1]]
                                            gaps.append(gap)
                                    usedGap = min(gaps)
                                    minDis, minPos = 0, 0
                                    for pos, value in enumerate(arr):
                                        relaLen = abs(usedGap - value[1])
                                        if minDis > relaLen:
                                            minDis = relaLen
                                            minPos = pos
                                    libSPlink[c1][d1][c2][d2] = [arr[minPos]]
                                else:
                                    # no dlink path can connect two contigs
                                    # if they have overlaped gaps merge them
                                    arr.sort(key=lambda x:x[1])
                                    libSPlink[c1][d1][c2][d2] = [arr[0]]
                # find cross pair mapping results
                counter = collections.Counter([i[2] for i in newarr])
                out = [i for i in counter if counter[i] == 2] # candidate
                for ele in out:
                    o, t = [i for i in newarr if i[2] == ele]
                    if o[1] == -t[1] and o[3] == -t[3]:
                        dis.add((c1, o[1], o[2], o[3]))
            
            # store seqs
            seqs, nodes = {}, set([str(i) for k in dis for i in k])
            with open(self.tmpContig , 'r') as af:
                for line in af:
                    if line[0] == '>':
                        cid = line.rstrip()[1:]
                        flag = 1 if cid in nodes else 0
                    else:
                        if flag:
                            if len(line) > 200:
                                seqs[int(cid)] = line[:100] + line.rstrip()[-100:]
                            else:
                                seqs[int(cid)] = line.rstrip()
            
            # find the best overlaps
            sim = similarity()
            visited = set()
            for link in dis:
                if (link[0], link[2]) not in visited:
                    rigOder = self.scanKmeroverlap(seqs, link, sim)
                    del libSPlink[link[0]][-rigOder[0]][link[2]][-rigOder[1]]
                    del libSPlink[link[2]][-rigOder[1]][link[0]][-rigOder[0]]
                    visited.add((link[0], link[2]))
                    visited.add((link[2], link[0]))
            
            # transfer the formate of pelinks
            transedPlinks = {}
            for c1 in libSPlink:
                for d1 in libSPlink[c1]:
                    for c2 in libSPlink[c1][d1]:
                        for d2 in libSPlink[c1][d1][c2]:
                            arr = libSPlink[c1][d1][c2][d2]
                            fir = c1 if d1 == 1 else self.revNode(c1)
                            sec = c2 if d2 == 1 else self.revNode(c2)
                            tmp = [sec, arr[0][0], round(arr[0][1], 1), round(arr[0][2],1), key]
                            if fir not in transedPlinks:
                                transedPlinks[fir] = []
                            transedPlinks[fir].append(tmp)
            
            # combline different libs
            if not pelink:    
                pelink = copy.deepcopy(transedPlinks)
            else:
                for c1, v1 in transedPlinks.iteritems():
                    if c1 not in pelink:
                        pelink[c1] = v1
                    else:
                        pelink[c1].extend(v1)
        
        # remove the reads mapped to same contigs and the redundancy that same link shared by different libs
        # 3193493: [[3175737, -57, 3.0, 1], [3175737, -87, 20, 2]] for example
        linkScore = {}
        for ele in pelink:
            rev = self.revNode(ele)
            position = {}
            if ele not in linkScore:
                linkScore[ele] = {}
            if len(self.libs) > 1:
                for pos, end in enumerate(pelink[ele]):
                    if end[0] != rev and end[0] != ele:
                        if end[0] not in position:
                            position[end[0]] = []
                        position[end[0]].append(pos)
                for key, pos in position.iteritems():
                    if len(pos) > 1: # several libs support the link
                        tmp = [pelink[ele][i] for i in pos]
                        tmp.sort(key = lambda x: x[-1]) # c2, support number, mean gap, sigma, libid
                        linkScore[ele][key] = [tmp[0][2], tmp[0][3], 3, tmp[0][-1]]
                    else:
                        inf = pelink[ele][pos[0]]
                        if inf[3] >= 4:
                            linkScore[ele][key] =[inf[2], inf[3], 2, inf[-1]]
                        else:
                            linkScore[ele][key] =[inf[2], inf[3], 1, inf[-1]]
            else:
                for inf in pelink[ele]:
                    if inf[1] >= 4:
                        linkScore[ele][inf[0]] =[inf[2], inf[3], 3, inf[-1]]
                    else:
                        linkScore[ele][inf[0]] =[inf[2], inf[3], 2, inf[-1]]
            # if link contains only self mapping
            if len(linkScore[ele]) == 0:
                del linkScore[ele]
        
        # make two nodes on the same strand
        dPlink = indexGraph.Graph(self.index)
        for key in linkScore:
            for sec, value in linkScore[key].iteritems():
                dPlink.addEdge(key, self.revNode(ele), value)
                dPlink.addEdge(ele, self.revNode(key), value)
        
        return dPlink, contigCov
    
    def trans(self):
        """
        Check the files and generate the files that not exists
        """
        #ref_ext = set(['.contig', '.ContigIndex', '.links', '.Arc'])
        #ext = set([os.path.splitext(f)[1] for f in glob.glob('g' + '.*' )])
        #intersec = ref_ext - ext
        #if intersec:
        logger.info('Indexing contig and caculating contig length ...')
        # index and contig length
        self.index, self.length = self.transFasta(self.metaO)
        
        # get the libraries info
        self.libs = self.readLib()
        
        # check the fq files
        self.samFiles = {}
        for key in self.libs:
            prefix = 'lib{0}'.format(key)
            self.samFiles[key] = prefix + '.sam'
            for pair in self.libs[key]['reads']:
                for fn in pair:
                    self.checkFiles(fn)
        
        # map reads and calculate contig sequencing depth
        logger.info('Mapping reads to reference ...')
        sam = Pipeline(self.tmpContig)
        if not sam.isIndexed():
            sam.indexRef()
        for key in self.libs:
            if os.path.exists(self.samFiles[key]) and os.path.getsize(self.samFiles[key]):
                logger.info('Sam file {0} exists, skip mapping!'.format(self.samFiles[key]))
                continue
            for pair in self.libs[key]['reads']:
                prefix = 'lib{0}'.format(key)
                logger.info('Mapping {0} and {1}'.format(pair[0], pair[1]))
                sam.algin(pair[0], pair[1], prefix, self.cpu)
        
        # construct dlink graph
        for key in self.index.keys():
            if self.index[key] == 0:
                self.index[key+1] = 2
                self.length[key+1] = self.length[key]

        self.kmersize, self.dDlink = self.dlinkCreation()
        
        # finding pelinks and also calculating coverage info using the sam file
        self.dPlink, self.cov = self.plink()
        
        # create the g.contig file
        tmpf = open(self.tmpContig, 'r')
        seq, seqids = '', collections.deque()
        with open(self.transed_metaO, 'w') as af:
            for line in tmpf:
                if line[0] == '>':
                    seqids.append(int(line[1:-1]))
                    if seq:
                        tid = seqids.popleft()
                        af.write('>{0} length {1} cvg_{2}_tip_0\n{3}\n'.format(tid, self.length[tid], self.cov[tid], seq))
                        seq = ''
                else:
                    seq = line.rstrip()
            else:
                tid = seqids.popleft()
                af.write('>{0} length {1} cvg_{2}_tip_0\n{3}\n'.format(tid, self.length[tid], self.cov[tid], seq))
                seq = ''
        tmpf.close()
        os.remove(self.tmpContig)
        
        # output the files
        update = 0 # for finding lost dlinks
        logger.info('Outputing graphs ...')
        with open(self.BINFILE,'wb') as af:
            pickle.dump(self.kmersize, af, True)
            pickle.dump(self.index, af, True)
            pickle.dump(self.dPlink, af, True)
            pickle.dump(self.dDlink, af, True)
            pickle.dump(update, af, True)
            pickle.dump(self.libs, af, True)
        
        # output lib setting
        with open(self.libfn, 'wb') as af:
            pickle.dump(self.libs, af)
