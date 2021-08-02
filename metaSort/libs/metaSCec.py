############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import re
import subprocess
import collections
import config
import math
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

class gcBreak:
    
    """
    Return the percent GC between a pair of bases if the GC content of the
    range is more than 2.5 standard deviations from the average GC content
    of the sequence. It returns the average GC content of the sequence
    otherwise.
    """
    def __init__(self, seq, winSize):
        
        self.seq = seq.upper()
        self.seqLen = len(seq)
        self.winSize = winSize
        
    def averGC(self):
        
        numG = self.seq.count('G')
        numC = self.seq.count('C')
        rate = (numC+numG)*1.0/self.seqLen
        return rate
    
    def caculateSD(self):
        
        GC = frozenset(['G','C'])
        winCount = self.seqLen - self.winSize
        sum_so_far, gc_so_far, current_gc_count = 0, 0, 0
        
        for i in xrange(-self.winSize, winCount):
            if i > 0:
                preChar = self.seq[i-1]
                if preChar in GC:
                    current_gc_count -= 1
            
            newChar = self.seq[i+self.winSize]
            if newChar in GC:
                current_gc_count += 1
            
            if i >= 0:
                this_value = current_gc_count*1.0/self.winSize
                sum_so_far += this_value*this_value
                gc_so_far += this_value
        
        gc_average = gc_so_far/winCount
        sd = math.sqrt(sum_so_far*1.0/winCount - gc_average*gc_average)
        
        return sd
    
    def GScan(self):
        
        """
        Scan the genome and count the sites that larger than 2.5 times of sd
        """
        
        sd = self.caculateSD()
        gc_average = self.averGC()
        
        bias = []
        
        for i in xrange(self.seqLen-self.winSize+1):
            gc_count = 0
            subSeq = self.seq[i:i+self.winSize+1]
            gc_count += subSeq.count('G')
            gc_count += subSeq.count('C')
            gc_content = gc_count*1.0/self.winSize
            if abs(gc_content-gc_average) < sd*2:
                bias.append(0)
            else:
                bias.append(1)
        
        # get the continous news that have gc bias
        region = []
        for pos, value in enumerate(bias):
            if value == 1:
                if len(region) == 0:
                    region.append([pos])
                else:
                    if pos-region[-1][-1]==1:
                        region[-1].append(pos)
                    else:
                        region.append([pos])
        # merge range
        new = []
        for ele in region:
            if len(ele) == 1:
                newEle = [ele[0], ele[0]+self.winSize+1]
            else:
                newEle = [ele[0], ele[-1]+self.winSize+1]
            new.append(newEle)
        
        # merge them
        new.sort()
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
    
        return merged

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
    
def relocation(lfn, rfn, refn, threads=1):
    """
    This script was designed to find and revise the misassemble assembly results
    of mini-meta assembly by using big-meta reads
    """
    setting = {'lfn':lfn,
              'rfn':rfn,
              'scaf':refn,
              't':str(threads)
              }
    
    # bwa maps the reads to reference
    runBwa(setting)
    
    # only find the relocation of on the scaffolds
    sites = {}
    with open('tmp.sam', 'r') as af:
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
               # readID, refID, MappedPos, ciga, distance
               tmp.append([info[0], info[2], int(info[3]), info[5], int(info[8])])
               if len(tmp) == 2:
                    # check if we get a pair of reads
                    if tmp[0][0] != tmp[1][0]:
                        tmp = tmp[1:]
                        pairs = pairs[1:]
                        continue
                    if set(pairs) == set(['100011', '010011']):
                        if re.findall(r'\D',tmp[0][3]) == ['M'] and re.findall(r'\D',tmp[1][3]) == ['M']:
                                if tmp[0][-1] > 1000: # find the relocation
                                    if tmp[0][1] not in sites:
                                        sites[tmp[0][1]] = []
                                    sites[tmp[0][1]].append([tmp[0][0], tmp[0][2], tmp[1][2]])
    
def errorScan(spafn, ecfn):
    """
    Break the large scaffolds based on the error sites
    """    
    if not os.path.exists(spafn):
        logger.error('File: {0} not exists, please recheck it'.format(spafn),
                     exit_with_code=1
                     )
    
    # store the scaffold error sites
    logger.info('Error sites are writed into file: erroSites')
    output = open(ecfn, 'w')
    segments, scaffolds = [], {}
    with open(spafn,'r') as af:
        for line in af:
            if not line:
                continue
            if line[0] == '>':
                sid = line[:-1]
                scaffolds[sid] = ''
            else:
                 scaffolds[sid] += line[:-1]
    
    # rename the scaffolds into the format: NODE_xx_length_xx_cov_xx_ID_xx
    mgadir = os.path.dirname(ecfn)
    errout = open(os.path.join(mgadir, 'erroSites'),'w')
    idnumber = 0
    for key in scaffolds:
        seq = scaffolds[key]
        idnumber += 1
        if len(seq) < config.minECLen:
            seqid = 'NODE_{0}_length_{1}_cov_NULL_ID_1'.format(idnumber, len(seq))
            errout.write('{0}\t{1}\n'.format(key, seqid))
            output.write('>{0}\n{1}\n'.format(seqid, seq))
        else:  
            breaks = gcBreak(seq, config.winSize)
            breakPoints = breaks.GScan()
            if len(breakPoints) > 0:
                # output the sites
                seqid = 'NODE_{0}_length_{1}_cov_NULL'.format(idnumber, len(seq))
                errout.write('{0}\t{1}\n'.format(key, seqid))
                # get seqs
                for pair in breakPoints:
                    errout.write(' '.join(map(str,pair)) + '\n')
                # extract the partial sequences
                for pos, value in enumerate(breakPoints):
                    if pos == 0:
                        segments.append(seq[:value[0]])
                        segments.append(seq[value[0]:value[1]])
                    else:
                        segments.append(seq[breakPoints[pos-1][1]:value[0]])
                        segments.append(seq[value[0]:value[1]])
                # rename them and output
                count = 1
                for ele in segments:
                    if len(ele) >= 1000:
                        fragid = '>NODE_{0}_length_{1}_cov_NULL_ID_{2}'.format(idnumber, len(ele), count)
                        count += 1
                        output.write(fragid + '\n' + ele + '\n')
                segments = []
            else:
                seqid = 'NODE_{0}_length_{1}_cov_NULL_ID_1'.format(idnumber, len(seq))
                errout.write('{0}\t{1}\n'.format(key, seqid))
                output.write('>{0}\n{1}\n'.format(seqid, seq))
    
    output.close()
    errout.close()