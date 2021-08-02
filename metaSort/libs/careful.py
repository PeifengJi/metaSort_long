#!/usr/bin/python
############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import re
import subprocess
import itertools
import collections
import math
import numpy as np

import config
from log import get_logger
# set log file
logger = get_logger()

try:
    from sklearn.mixture import GMM
except:
    logger.warning('Package sklearn is not imported. GMM step will be passed.')

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
        logger.info('Creating sai files ...')
        
        # don't print the stdout to the screen
        with open(os.devnull, "w") as f:
            f1 = subprocess.Popen(cmd1, stdout=f, stderr=f)
        
        with open(os.devnull, "w") as f:
            f2 = subprocess.Popen(cmd2, stdout=f, stderr=f)
        
        f1.wait()
        f2.wait()
        # generate sam file
        cmdSam = [self.bwa, 'sampe',
                  self.ref, 
                  s1, s2, p1, p2]
        
        with open(os.devnull, "w") as f:
            proc = subprocess.Popen(cmdSam, stdout=subprocess.PIPE, stderr=f)
        
        with open('{0}.sam'.format(outPrefix), 'aw') as af:
            while 1:
                line = proc.stdout.readline()
                if line[:3] != '@SQ' and line[:3] != '@PG':
                    break
            af.write(line)
            while 1:
                line = proc.stdout.readline()
                if not line:
                    break
                af.write(line)
        
        logger.info('Sam file is generated!')
        
        self.__clean(s1, s2)
    
    def __clean(self, *files):
        """
        remove the temple files
        """
        for fn in files:
            os.remove(fn)

class BreakOut:
    """
    Detect the break points of bwa mapping results
    """
    def __init__(self, overlap):
        
        self.range = {}
        self.error = {}
        # dn: deleteion and link use 'N'
        # dl: deletion and link use two reads
        # b : break the error site
        self.overlap = overlap if overlap < 0 else 0
    
    def parse(self, lines):
        
        # parse the cigar: readID, refID, MappedPos, ciga, seq
        lines.sort(key=lambda x:x[2])
        cigars = (re.sub('\d','', lines[0][3]), re.sub('\d','',lines[1][3]))
        if lines[0][1] not in self.range:
            self.range[lines[0][1]] = []
        if lines[0][1] not in self.error:
            self.error[lines[0][1]] = []
        
        if cigars==('M', 'M'):
            self.case1(lines)
        elif cigars==('MS','M') or cigars==('M','SM'):
            self.case2(lines, cigars)
        elif cigars==('SM','M') or cigars==('M','MS'):
            self.case3(lines, cigars)
        elif cigars==('MS','MS') or cigars==('SM','SM'):
            self.case4(lines, cigars)
        elif cigars==('MS','SM') or cigars==('SM','MS'):
            self.case5(lines, cigars)
    
    def case1(self, tmp):
        
        lenUp = int(re.search(r'(\d+)M',tmp[0][3]).group(1))
        lenDown = int(re.search(r'(\d+)M', tmp[1][3]).group(1))
        distance = tmp[1][2] - tmp[0][2] - lenUp
        self.range[tmp[0][1]].append([tmp[0][2], tmp[1][2]+lenDown])
    
    def case2(self, tmp, cigar):
        # no way two reads overlapped
        lenUp = int(re.search(r'(\d+)M',tmp[0][3]).group(1))
        lenDown = int(re.search(r'(\d+)M',tmp[1][3]).group(1))
        if cigar == ('MS','M'):
            self.error[tmp[0][1]].append([tmp[0][2]+lenUp, tmp[1][2], 'dn'])
            #self.range[tmp[0][1]].append([tmp[0][2]+lenUp, tmp[1][2]+lenDown])
        else: # M SM
            self.error[tmp[0][1]].append([tmp[0][2]+lenUp, tmp[1][2], 'dn'])
            #self.range[tmp[0][1]].append([tmp[0][2]+lenUp, tmp[1][2]+lenDown])
    
    def case3(self, tmp, cigar):
        lenDown = int(re.search(r'(\d+)M',tmp[1][3]).group(1))
        if cigar == ('SM', 'M'):
            self.error[tmp[0][1]].append([tmp[0][2]])
            #self.range[tmp[0][1]].append([tmp[0][2], tmp[1][2]+lenDown])
        else: # M MS
            self.error[tmp[0][1]].append([tmp[1][2]+lenDown])
            #self.range[tmp[0][1]].append([tmp[0][2], tmp[1][2]+lenDown])
    
    def case4(self, tmp, cigar):
        # no way two reads overlapped
        lenUp = int(re.search(r'(\d+)M',tmp[0][3]).group(1))
        lenDown = int(re.search(r'(\d+)M',tmp[1][3]).group(1))
        if cigar == ('MS','MS'):
            self.error[tmp[0][1]].append([tmp[0][2]+lenUp, tmp[1][2]+lenDown, 'b'])
            #self.range[tmp[0][1]].append([tmp[0][2], tmp[0][2]+lenUp])
        elif cigar==('SM','SM'):
            self.error[tmp[0][1]].append([tmp[0][2], tmp[1][2], 'b'])
            #self.range[tmp[0][1]].append([tmp[1][2], tmp[1][2]+lenDown])
    
    def case5(self, tmp, cigar):
        lenUp = int(re.search(r'(\d+)M',tmp[0][3]).group(1))
        lenDown = int(re.search(r'(\d+)M',tmp[1][3]).group(1))
        if cigar==('SM','MS'): # no matter they are overlapped
            self.error[tmp[0][1]].append([tmp[0][2], tmp[1][2]+lenDown, 'b'])
        else: # MS','SM'
            if self.overlap < 0:
                stackLen = stack(tmp[0][-1], tmp[1][-1])
                self.error[tmp[0][1]].append([tmp[0][2]+lenUp, tmp[1][2], 'dl', stackLen, tmp[0][-1], tmp[1][-1]])
            else:
                stackLen == 0
                self.error[tmp[0][1]].append([tmp[0][2]+lenUp, tmp[1][2], 'dn'])
            
            #self.range[tmp[0][1]].append([tmp[0][2], tmp[0][2]+lenUp])
            #self.range[tmp[0][1]].append([tmp[1][2], tmp[1][2]+lenDown])
    
    def stack(self, seq1, seq2):
        total = 0
        score = {
            'TT':1, 'AA':1, 'CC':1, 'GG':1,
            'AT':-1, 'AC':-1, 'AG':-1,
            'CT':-1, 'GT':-1, 'CG':-1,
        }
        
        length = len(seq1) if len(seq1) <= len(seq2) else len(seq2)
        for i in xrange(3, length):
            subseq1 = seq1[-i:]
            subseq2 = seq2[:i+1]
            for comb in itertools.izip(subseq1, subseq2):
                string = ''.join(comb)
                total += score[string]
            
            if total == i:
                return 0 - total
            else:
                total = 0
        else:
            return 0
    
    def mergeRange(self):
        """
        use overlaps to detect the break out point, and merge the range
        """
        coverage = {}
        for key in self.range:
            if len(self.range[key]) == 1:
                coverage[key] = self.range[key]
                continue
            new = sorted(breaks.range[key])
            merged, record = [], [0]
            for pos in xrange(len(new)-1):
                pre = new[pos]
                aft = new[pos+1]
                if aft[0] < pre[1]: # overlaped
                    if pos+1 == len(new)-1:
                        leftmost = new[record[-1]][0]
                        rightmost = sorted([x[1] for x in new[record[-1]:]])[-1]
                        merged.append([leftmost, rightmost])
                    else:
                        continue
                else:
                    mergeLen = pos - record[-1]
                    if mergeLen == 0:
                        merged.append(pre)
                    else:
                        leftmost = new[record[-1]][0]
                        rightmost = sorted([x[1] for x in new[record[-1]:pos+1]])[-1]
                        merged.append([leftmost, rightmost])
                    record.append(pos+1)
                    if pos+1 == len(new)-1:
                        merged.append(new[pos+1])
            # store the range
            coverage[key] = merged
        return coverage
    
    def getErrors(self):
        
        # delete the unerrored contigs
        for key in self.error.keys():
            if len(self.error[key]) == 0:
                del self.error[key]
        
        newError = {}
        coverage = self.mergeRange()
        # to ensure the error site no well mapped reads
        for refid in self.error:
            if refid in coverage:
                for site in self.error[refid]:
                    if len(site) == 2:
                        for flate in coverage[refid]:
                            if flate[0] < site[0] < flate[1]:
                                break
                        else:
                            if refid not in newError:
                                newError[refid] = []
                            newError[refid].append(site)
                    else:
                        for flate in coverage[refid]:
                            if flate[0] < site[0] and site[1] < flate[1]:
                                break
                        else:
                            if refid not in newError:
                                newError[refid] = []
                            newError[refid].append(site) 
        
        return newError

def univeralSysCall(cmd, outfn=None, errfn=None):
    """
    Runs cmd and redirects stdout to out_filename (if specified),
    stderr to err_filename (if specified), or to log otherwise
    """
    if not isinstance(cmd, list):
        raise Exception('A list format of cmd should be given')
    if outfn:
        stdout = open(outfn, 'w')
    else:
        stdout = subprocess.PIPE
    if errfn:
        stderr = open(errfn, 'w')
    else:
        stderr = subprocess.PIPE
    
    proc = subprocess.call(cmd, stdout=stdout, stderr=stderr)
    
    if outfn:
        stdout.close()
    if errfn:
        stderr.close()

def runBwa(setting):
    """
    bwa mapping
    """
    bwaPath = os.path.join(config.libsLocation, 'bwa')
    # check the binary file
    if os.path.isfile(bwaPath) and os.access(bwaPath, os.X_OK):
        # bwa comand lines
        print('Running bwa ...')
        univeralSysCall( [bwaPath, 'index', '-a', 'is', setting['scaf']])
        univeralSysCall( [bwaPath, 'aln', setting['scaf'], setting['lfn'], '-t', setting['t']], 'tmp1.sai')
        univeralSysCall( [bwaPath, 'aln', setting['scaf'], setting['rfn'], '-t', setting['t']], 'tmp2.sai')
        univeralSysCall( [bwaPath, 'sampe', setting['scaf'], 'tmp1.sai', 'tmp2.sai', setting['lfn'], setting['rfn']], 'tmp.sam')
        logger.info('Mapping done')
    else:
        logger.error('Could not find bwa path, check it',
                     exit_with_code=1
                     )

def pilelup(samfile, breaks):
    """
    Caculate the depth of every single base of the scaffolds
    """
    minECLen = 1000
    #minECLen = config.minECLen
    acceptCigar = ('M', 'MS','SM', 'SMS')
    refs = {}
    pairs, tmp = [], []
    with open(samfile, 'r') as af:
        # header
        for line in af:
            if '@SQ' in line and 'SN:' in line and 'LN:' in line:
                continue
            elif '@PG' in line and 'ID:' in line and 'PN:' in line:
                continue
            else:
                break
        # mapping info
        for line in af:
            info = line.split()
            parser = bin(int(info[1]))
            if parser[-4:] == '0011' and info[6] == '=':
               pairs.append(parser[-6:])
               # readID, refID, MappedPos, ciga, seq
               tmp.append([info[0], info[2], int(info[3]), info[5], info[9]])
               #tmp.append([info[2], int(info[3]), info[5]])
               if len(tmp) == 2:
                    # check if we get a pair of reads
                    if tmp[0][0] != tmp[1][0]:
                        #print('One of two reads is not paired: {0} and {1}, pop the former'.format(tmp[0][0],tmp[1][0]))
                        tmp = tmp[1:]
                        pairs = pairs[1:]
                        continue
                    if set(pairs) == set(['100011', '010011']):
                        breaks.parse(tmp)
                        # parse the mapping results
                        refid = tmp[0][1]
                        reflen = int(re.search(r'length_(\d+)_cov', refid).group(1))
                        if reflen >= minECLen:
                            if refid not in refs: refs[refid] = {}
                            c1, c2 = re.sub('\d','', tmp[0][3]), re.sub('\d','', tmp[1][3])
                            if c1 in acceptCigar and c2 in acceptCigar:
                                for onePair in tmp:
                                    cigar = re.sub('\d','', onePair[3])
                                    lenDw = int(re.search(r'(\d+)M', onePair[3]).group(1))
                                    for pos in xrange(onePair[2]-1, onePair[2]+lenDw-1):
                                        if pos not in refs[refid]:
                                            refs[refid][pos] = 0
                                        refs[refid][pos] += 1
                        tmp = []
                        pairs = []    
    return refs

def isize(samfn):
    """
    Estimate the inner length between the paired end reads
    """
    maxNum, count = 200000, 0
    insertLen = []
    readLen = []
    with open(samfn, 'r') as af:
        for line in af:
            info = line.split()
            try:
                parser = bin(int(info[1]))
                if parser[-4:] == '0011' and info[6] == '=':
                   insertLen.append(abs(int(info[8])))
                   group = re.findall(r'(\d+)', info[5])
                   rlen = sum(map(int, group))
                   readLen.append(rlen)
                   count += 1
                   if count >= maxNum:
                    break
            except:
                continue
    averInnerSpa = int(sum(insertLen)*1.0/count - sum(readLen)*2.0/count)
    return averInnerSpa

def cluster(args):
    # read parameters
    c, data = args
    
    # use EM algorithm to fit the module
    gmm = GMM(n_components=c, covariance_type='full', n_init=5, n_iter=20000).fit(data)
    
    # get the bic criterion score
    bic = gmm.bic(data)
    
    if not gmm.converged_:
        print(("Clustering into {0} clusters did not "
                         "converge, consider increasing the number "
                         "of iterations.").format(c))
        print >> sys.stderr, "Cluster {0} did not converge".format(c)
    
    return bic,c, gmm.converged_

def estimate(points):
    
    clusterArgs = [[i, points] for i in range(1,6)]
    # parallels to run gmm
    results = []
    for arg in clusterArgs:
        results.append(cluster(arg))
    
    # determine the optimal cluster number by using BIC
    min_bic, optimalC, converged = min(results,key=lambda x: x[0])
    gmm = GMM(n_components=optimalC,
              covariance_type='full',
              n_init=5,
              n_iter=20000,
              ).fit(points)
    
    # use the largest peak
    clnum = sorted([(k, v) for k, v in zip(gmm.weights_, xrange(optimalC))])[-1][-1]
    # use the posterior probabilities of each data point to cluster
    clusters = gmm.predict(points)
    peak = []
    for k, v in itertools.izip(clusters, points):
        if k == clnum:
            peak.append(v)
    
    mean = np.array(peak).mean()
    std = np.array(peak).std()
    return mean - std, mean + std

def peakCluster(pileupfn, graphLim=0.8, mod = 1):
    """
    Estimated the coverage distribution of the scaffolds
    using a slding window, only used for scaffolds > 5k
    """
winSize, lenLim = 100, 5000
stat, paras,  = {}, collections.deque()
with open(pileupfn, 'r') as af:
    for line in af:
        ae = line.rstrip().split()
        if ae[0] not in stat:
            if paras and len(paras)*100 >= lenLim:
                _ = paras.popleft()
                _ = paras.pop()
                stat[iden] = np.array(paras)
            # initialization
            stat[ae[0]] = []
            counter, paras, tmp = 1, collections.deque(), []
        else:
            iden = ae[0]
            if winSize*(counter-1) < int(ae[1]) < winSize*counter+1:
                tmp.append(float(ae[2]))
            else:
                counter += 1
                if tmp:
                    paras.append(np.array(tmp).mean())
                else:
                    paras.append(0)
                tmp = []
    else:
        if paras and len(paras)*100 >= lenLim:
            _ = paras.popleft()
            _ = paras.pop()
            stat[iden] = np.array(paras)
    
    # for visualization
    with open('scaffold.cov','w') as af:
        for key, value in stat.iteritems():
            counter = 1
            for ele in value:
                af.write('{0}\t{1}\t{2}\n'.format(key, counter, ele))
                counter += 1
        
        if mod == 1:
            # compare the coverage distribution of each scaffolds and clustering
            # model the peak using em algorithm
            stat_ran = {}
            for k, v in stat.iteritems():
                ar = np.array([[i] for i in v])
                stat_ran[k] = estimate(ar)
            
            # cluster
            backup, graph = {}, {}
            for com in itertools.combinations(stat_ran.keys(), 2):
                s, t = sorted([stat_ran[com[0]], stat_ran[com[1]]])
                if t[0] < s[1]: # overlap
                    mint = sorted([s[1], t[1]])[0]
                    overlap = mint - t[0]
                    smallRan = sorted([s[1]-s[0], t[1]-t[0]])[0]
                    overlapRatio = round(overlap*1.0/smallRan, 2)
                    if overlapRatio > graphLim:
                        if com[0] not in graph:
                            graph[com[0]] = {}
                        if com[1] not in graph:
                            graph[com[1]] = {}
                        graph[com[0]][com[1]] = graph[com[1]][com[0]] = overlapRatio
                    else:
                        if com[0] not in backup:
                            backup[com[0]] = {}
                        if com[1] not in backup:
                            backup[com[1]] = {}
                        backup[com[0]][com[1]] = backup[com[1]][com[0]] = overlapRatio
                    
                    # find subgraph to cluster the coverage profiles
            cl = {i:1 for i in graph}
            group = 1
            stack = []
            for key in cl.keys():
                if cl[key] == 1:
                    stack.append(key)
                    while stack:
                        last = stack.pop()
                        cl[last] = group + 1
                        for ele in graph[last]:
                            if cl[ele] == 1:
                                stack.append(ele)
                    group += 1
            
            clusters = {}
            for k, v in cl.iteritems():
                if v not in clusters:
                    clusters[v] = []
                clusters[v].append(k)
            
            if len(clusters) > 1:
                print('mod 1')
                # use the contamination to absort miss contaminations
                for k in set(backup.keys()) - set([v for i in clusters.itervalues() for v in i]):
                    # find the maximum cluster
                    v = backup[k]
                    fin = sorted(v.keys(), key=lambda x: v[x], reverse=True)[0]
                    clusters[cl[fin]].append(k)
                
            return clusters
        else:
            # too close the contamination and postive data
            mean = {k: np.mean(v) for k, v in stat.iteritems()}
            sample_mean = np.array(mean.values()).mean()
            sample_std = np.array(mean.values()).std()
            contaminations = set()
            for i, v in mean.iteritems():
                if v <= sample_mean-3*sample_std or v >= sample_mean+3*sample_std:
                    contaminations.add(i)
    elif mod == 2:
        print('mod 2')
        # too close the contamination and postive data
        mean = {k: np.mean(v) for k, v in stat.iteritems()}
        sample_mean = np.array(mean.values()).mean()
        sample_std = np.array(mean.values()).std()
        contaminations = set()
        for i, v in mean.iteritems():
            if v <= sample_mean-2*sample_std or v >= sample_mean+2*sample_std:
                contaminations.add(i)
    
    return contaminations

def filterCon(scafn, pileupfn):
    """
    Remove the contamination by coverage
    """
    res = peakCluster(pileupfn)
    
seqs = {}
with open(scafn,'r') as af:
    for line in af:
        if line[0] == '>':
            iden = line.rstrip()[1:]
        else:
            seqs[iden] = line.rstrip()

if isinstance(res, set):
    cleanfn = open('cleaned.fa','w')
    confn = open('conta.fa','w')
    for k, v in seqs.iteritems():
        if k not in res:
            cleanfn.write('>{0}\n{1}\n'.format(k, v))
        else:
            confn.write('>{0}\n{1}\n'.format(k, v))
    
    cleanfn.close()
    confn.close()
else:
    for key in res:
        with open('cluster{0}.fa'.format(key), 'w') as af:
            for node in res[key]:
                af.write('>{0}\n{1}\n'.format(node, seqs[node]))
            
        














