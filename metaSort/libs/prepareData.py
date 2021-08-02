############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import genemark
import cPickle as pickle
import glob
import multiprocessing
import collections
import subprocess
import config
import shutil
import itertools
import genemark
import gapEstimate
import indexGraph
import metaSCec
import platform
import lostLinks
from pybloom import ScalableBloomFilter
from log import get_logger

# set log file
logger = get_logger()

fivebase = ['A','T','C','G','N']
codon = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA',
         'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC',
         'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG',
         'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT',
         'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']

def unwrap_kmerOrder(worker):
    return worker.kmerOrder()

def unwrap_kmerFreq(worker):
    return worker.kmerFreq()

def unwrap_kmerComp(worker):
    return worker.kmerComposition()

def unwrap_codon(worker, genemarkPath):
    return worker.proteinFormate(genemarkPath)

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

class metaBAT:
    """
    metaBAT binning
    """
    def __init__(self, reference, output, threads=1):
        
        self.refn = os.path.abspath(reference)
        self.output = os.path.abspath(output)
        self.clusters = 'metaO_cluster'
        self.threads = str(threads)
        
    def checkPlatform(self):
        """
        Check if metabat is available
        """
        if platform.system() == 'Darwin':
            metabatPath = os.path.join(config.libsLocation, 'metabat', 'macosx')
        elif platform.system() == 'Linux':
            metabatPath = os.path.join(config.libsLocation, 'metabat', 'linux')
        
        metabatExc = os.path.join(metabatPath, 'metabat')
        with open(os.devnull, "w") as f:
            proc = subprocess.Popen([metabatExc], stdout=f, stderr=f)
            if proc:
                return metabatPath
            else:
                return 0
    
    def cleanUp(self, *files):
        """
        remoe unused files
        """
        for fn in files:
            if os.path.exists(fn):
                os.remove(fn)
        
    def checkSamtools(self):
        """
        Check if samtools is available
        """
        with open(os.devnull, "w") as f:
            proc = subprocess.Popen(['samtools'], stdout=f, stderr=f)
            if proc:
                return 1
            else:
                return 0
    
    def sortSAM(self, samfn):
        """
        Generate sorted bam file
        """
        cmd = ['samtools', 'faidx', self.refn]
        retunCode = subprocess.call(cmd)
        
        cmd = ['samtools', 'import', self.refn, samfn, 'result.bam']
        retunCode = subprocess.call(cmd)
        
        cmd = ['samtools', 'sort', 'result.bam', 'sorted']
        retunCode = subprocess.call(cmd)
        
        os.remove('result.bam')
        for fn in glob.glob(self.refn + '.*'):
            self.cleanUp(fn)
    
    def binning(self, samfn):
        """
        Use metaBAT bins the scaffolds
        """
        binPath = self.checkPlatform()
        
        if self.checkSamtools():
            # go on
            logger.info('Detect samtools, go on ...')
        else:
            logger.warning('Does not detect samtools, metaBAT could not be used !')
            return 0
        
        if os.path.exists('sorted.bam') and os.path.getsize('sorted.bam'):
            logger.info('File: sorted.bam is found, skip sorting sam step ...')
        else:
            logger.info('Generating bam files, please wait ...')
            self.sortSAM(samfn)
        
        curDir = os.getcwd()
        if not os.path.isdir(self.output):
            os.makedirs(self.output)
        else:
            shutil.rmtree(self.output)
            os.makedirs(self.output)
        os.chdir(self.output)
        
        logger.info('Using metaBAT binning contigs, please wait ...')
        # calculating contig depth
        calDepth = os.path.join(binPath, 'jgi_summarize_bam_contig_depths')
        cmd = [calDepth, '--outputDepth', 'contig_cov.txt', os.path.join(curDir, 'sorted.bam')]
        retunCode = subprocess.call(cmd)
        
        # binning
        prefix = 'bin'
        getBins = os.path.join(binPath, 'metabat')
        cmd = [getBins, '-o', prefix, '-i', self.refn, '-a', 'contig_cov.txt', '-t',
               self.threads, '--saveCls', '--superspecific', '-m', str(config.metaO_minLen)
               ]
        retunCode = subprocess.call(cmd)
        
        shutil.move('contig_cov.txt', curDir)
        shutil.move(prefix, curDir)
        os.chdir(curDir)
        
        cluster = {}
        with open('bin','r') as af:
            for line in af:
                ae = [int(i) for i in line.rstrip().split()]
                if ae[1] > 0:
                    key = 'o' + str(ae[1])
                    if key not in cluster:
                        cluster[key] = set()
                    cluster[key].add(ae[0])
        
        # store the clusters for contaminants removal
        with open(self.clusters, 'wb') as af:
            pickle.dump(cluster, af)
        
        os.remove('sorted.bam')
        os.remove(samfn)
        logger.info('binning done !')

class SOAP:
    """
    Store and transform the format of soap files
    """
    def __init__(self, database, threads=1):
        
        #self.soapdir = soapdir
        self.threads = threads
        self.database = database
        self.update = 0
        self.BINFILE = os.path.join(self.database,'Graphs')
        self.transedContig = os.path.join(self.database, 'contigfile')
        self.PROFREQ = os.path.join(self.database, 'profreq')
        self.TETRAFREQ = os.path.join(self.database, 'tetrafreq')
        self.tetraOrder = os.path.join(self.database, 'tetraOrder')
        self.nufile = '/'.join(os.path.abspath(self.transedContig).split('/')[:-1]) + '/proSVMseqs.fa'
        self.profile = '/'.join(os.path.abspath(self.transedContig).split('/')[:-1]) + '/pro.lst'
        self.bins = os.path.join(self.database, config.LMETACL)
        self.clusters = os.path.join(self.database, 'metaO_cluster')
        self.libStat = os.path.join(self.database, 'libStat.txt')
        self.kmerFile = os.path.join(self.database, 'fiftyKmers')
        self.idIndexFile = os.path.join(self.database, 'idIndex')
        self.arcfile, self.readOnContigFn, self.linkfile, self.contigfile, \
        self.indexfile, self.scafile, self.scafSeqFile = ['' for i in xrange(7)]
    
    def checkFile(self, *fns):
        """
        Check if the file exists
        """
        for fn in fns:
            if not os.path.exists(fn) or not os.path.getsize(fn):
                logger.error('File: {0} not exists or size of the file is zero.'.format(fn),
                             to_stderr=True,
                             exit_with_code=3
                             )
    
    def getCovLen(self):
        """
        Get the kmer coverage info from the contig file
        """
        print('Getting kmer coverage info of every nodes on the graph ...')
        self.cov = collections.defaultdict(float)  # for storing kmer coverage info of every nodes
        self.length = collections.defaultdict(int)
        
        with open(self.transedContig, 'r') as contig_file:
            for line in contig_file:
                contigId, contigLen, value, _ = line.rstrip().split()
                contigId, contigLen, value = int(contigId), int(contigLen), float(value)
                rev = self.revNode(contigId)
                self.cov[contigId] = self.cov[rev] = value
                self.length[contigId] = self.length[rev] = contigLen
    
    def oneStrandCovLen(self):
        """
        Get the coverage and length info generated by SOAPdenovo
        """
        print('Getting kmer coverage info of every nodes on the graph ...')
        with open(self.BINFILE,'rb') as af:
            _ = pickle.load(af)
            self.index = pickle.load(af)
        
        self.cov = collections.defaultdict(float)  # for storing kmer coverage info of every nodes
        self.length = collections.defaultdict(int)
        
        with open(self.transedContig, 'r') as contig_file:
            for line in contig_file:
                contigId, contigLen, value, _ = line.rstrip().split()
                contigId, contigLen, value = int(contigId), int(contigLen), float(value)
                rev = self.revNode(contigId)
                self.cov[contigId] = self.cov[rev] = value
                self.length[contigId] = self.length[rev] = contigLen
    
    def transContigfile(self, contigFile, target):
        """
        Transform the contig file into a file including four columns:
        1:id; 2:length; 3 kmer coverage; 4: sequence
        Used for fast reading and convience, less code content
        """
        if os.path.exists(target) and os.path.getsize(target):
            return 1
        
        if not (os.path.exists(target) and os.path.getsize(target)) \
            or not (os.path.exists(self.idIndexFile) and os.path.getsize(self.idIndexFile)):
            out = open(target,'w')
            flag = collections.deque()
            seq, kmersize, idIndex = '', '', {}
            
            with open(contigFile, 'r') as af:
                for line in af:
                    if line[0] == '>':
                        flag.append(line[1:])
                        info = flag.popleft().split()
                        kmerCov = info[3].split(r'_')[1]
                        idIndex[int(info[0])] = out.tell()
                        kmersize = int(line.split()[2]) - 1
                        out.write(' '.join([info[0], info[2], kmerCov, seq]) + '\n')
                        seq = ''
                    else:
                        seq += line[:-1]
                else:
                    info = flag.popleft().split()
                    kmerCov = info[3].split(r'_')[1]
                    idIndex[int(info[0])] = out.tell()
                    out.write(' '.join([info[0], info[2], kmerCov, seq]) + '\n')
            out.close()
            
            update = 0
            with open(self.idIndexFile, 'wb') as af:
                pickle.dump(idIndex, af)
                pickle.dump(update, af)
            
            return kmersize
    
    def proteinFormate(self, genemarkPath):
        """
        Predict the protein coding region of the contigs for SVM
        """
        protein, seq, ids = {}, '', collections.deque()
        gmhmmp = genemarkPath + '/gmhmmp'
        modfile = '/'.join(genemarkPath.split('/')[:-1]) + '/MetaGeneMark_v1.mod'
        
        svmfile = open(self.nufile, 'w')
        with open(self.transedContig, 'r') as c_f:
            for line in c_f:
                identifier, _, _, seq = line[:-1].split()
                svmfile.write('>' + identifier + '\n' + seq + '\n')
        svmfile.close()
        # using genemark to predict genes
        cmdMetaGeneMark = [gmhmmp, '-d', '-f', 'G', '-m', modfile, '-o', self.profile, self.nufile]
        logger.info('Predicting proteins from the contigs using gmhmmp')
        logger.print_command_line(cmdMetaGeneMark, wrap_after=None)
        
        child = subprocess.call(cmdMetaGeneMark)
        if child != 0:
            logger.error('Error happened while using MetaGeneMark',
                         to_stderr=True,
                         exit_with_code=1
                         )
        else:
            with open(self.profile,'r') as lfile:
                for i in xrange(6): _ = lfile.readline()
                for line in lfile:
                    if line[0] != '\n' and line[0] != '#':
                        if not ids:
                            ids.append(line.split()[0])
                            continue
                        
                        info = line.split()
                        if ids[-1] != info[0]:
                            ids.append(info[0])
                        if seq:
                            usage = collections.defaultdict(int)
                            for po in xrange(0,len(seq)-2,3):
                                usage[seq[po:po+3]] += 1
                            protein[ids.popleft()] = [usage.get(ele,0) for ele in codon]
                            seq = ''
                    elif line[0] != '\n' and line[0:2] == '##' and line[2] in fivebase:
                        seq += line[2:-1]
                else:
                    usage = collections.defaultdict(int)
                    for po in xrange(0,len(seq)-2,3):
                        usage[seq[po:po+3]] += 1
                    protein[ids.popleft()] = [usage.get(ele,0) for ele in codon]
                    seq = ''
        
        # dumple it 
        with open(self.PROFREQ, 'wb') as af:
            pickle.dump(protein, af)
    
    def cleanUp(self):
        """
        Delete the unused tmple files
        """
        if os.path.exists(self.nufile):
            os.remove(self.nufile)
        
        if os.path.exists(self.profile):
            os.remove(self.profile)
    
    def permutation(self, size = 4):
        """
        Generate the combination of a,t,c,g
        """
        combination = []
        S = ['A','T','G','C']
        for seq in itertools.product(S, repeat = size):
            combination.append(''.join([i for i in seq]))
        
        merged = []
        for ele in combination:
            revcom = self.revCom(ele)
            if revcom not in merged:
                merged.append(ele)
        
        merged.sort()
        return merged
    
    def calFreq(self, seq, kmerSpectrum, size = 4, out = 'both'):
        """
        caculate the kmer frequency and sort the list by the dict value
        """
        if size % 2 == 0:
            number = (4 ** size + 4 ** (size / 2)) / 2
        else:
            number = (4 ** size) / 2
        
        numberindex = xrange(number)
        # get the kmer frequency
        arr = ( seq[i:i+size] for i in xrange(len(seq)-size+1) )
        table = collections.defaultdict(int)
        for i in arr:
            table[i] += 1
        
        # remove the 'N' contained kmers
        for key in itertools.ifilter(lambda x: set(x) - set('ATCG'), table.keys()):
            del table[key]
        
        for i in (set(table) - set(kmerSpectrum)):
            table[self.revCom(i)] += table[i]
            del table[i]
        
        get = table.get
        # use the total count to normalize
        allCount = sum(table.values())
        if out == 'freq':
            #get the number of each type of kmer in the kmerSpectrum to ready for SVM
            kmerNumber = ( get(i, 0) for i in kmerSpectrum )
            return kmerNumber
        elif out == 'order':
            # This step is to make sure the order of kmers that number is zero is the same among differenct contigs
            orders = sorted(table, key = get) # from low to high (value), returns the oreded keys
            if len(orders) < number:
                newKmers = [i for i in kmerSpectrum if i not in orders]
                newKmers.extend(orders)
                return dict(zip(newKmers, numberindex))
            else:
                return dict(zip(orders, numberindex))
        elif out == 'both':
            kmerNumber = ( get(i, 0) for i in kmerSpectrum )
            orders = sorted(table, key = get)
            if len(orders) < number:
                newKmers = [i for i in kmerSpectrum if i not in orders]
                newKmers.extend(orders)
                return dict(zip(newKmers, numberindex)), kmerNumber
            else:
                return dict(zip(orders, numberindex)), kmerNumber
    
    def revCom(self, string):
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
            try:
                result += comp[i]
            except:
                 result += i
        
        return result
    
    def kmerOrder(self):
        # generate the tetranuclotied order of the contigs
        kmerSpectrum = self.permutation()
        # transform the store method
        outfile = open(self.tetraOrder, 'w')
        with open(self.transedContig, 'r') as af:
            for line in af:
                info = line.split()
                if int(info[1]) >= 300:
                    seqOrder = self.calFreq(info[-1][:-1], kmerSpectrum, out='order')
                    outfile.write(info[0] + ' ' + ' '.join(map(str, [seqOrder[i] for i in kmerSpectrum])) + '\n')
        outfile.close()
    
    def kmerFreq(self):
        # generate the tetranucleotide frequency of the contigs
        kmerSpectrum = self.permutation()
        # transform the store method
        outfile = open(self.TETRAFREQ, 'w')
        with open(self.transedContig, 'r') as af:
            for line in af:
                info = line.split()
                if int(info[1]) >= 100:
                    seqFreq = self.calFreq(info[-1][:-1], kmerSpectrum, out='freq')
                    outfile.write(info[0] + ' ' + ' '.join(map(str, seqFreq)) + '\n')
        outfile.close()
    
    def loadData(self, samFn = False):
        """
        Check if the files in soapdir exists and transform the file formate
        also predict the code usage of each contig
        """
        self.contigfile = os.path.join(self.database, 'g.contig')
        
        # less code content also get the assembly kmersize for later use
        self.kmersize = self.transContigfile(self.contigfile, self.transedContig)
        
        # handling index file, first check if the files were already generated
        if os.path.exists(self.BINFILE) and os.path.getsize(self.BINFILE):
            logger.info('Loading graphs from file: ' + self.BINFILE)
            with open(self.BINFILE,'rb') as af:
                self.kmersize = pickle.load(af)
                self.index = pickle.load(af)
                self.dPlink = pickle.load(af)
                self.dDlink = pickle.load(af)
                self.update = pickle.load(af)
                self.libs = pickle.load(af)
            # load coverage and length info
            self.getCovLen() 
            if not self.update:
                # update other assemblers, find lost dlinks
                self.dDlink = lostLinks.findLostDlinks(self)
                self.update = 1
                with open(self.BINFILE,'wb') as af:
                    pickle.dump(self.kmersize, af, True)
                    pickle.dump(self.index, af, True)
                    pickle.dump(self.dPlink, af, True)
                    pickle.dump(self.dDlink, af, True)
                    pickle.dump(self.update, af, True)
                    pickle.dump(self.libs, af, True)
        
        # handeling protein coding file, first check if the files were already generated
        if not os.path.exists(self.PROFREQ) or not os.path.getsize(self.PROFREQ):
            # install genemark
            genemarkPath = genemark.getPath()
            p1 = multiprocessing.Process(target=unwrap_codon, args=(self, genemarkPath))
        else:
            logger.info('\nCode usage file existis, skip!')
            p1 = 0
        
        # generate the tetranuclotied order of the contigs
        if not os.path.exists(self.tetraOrder) or not os.path.getsize(self.tetraOrder):
            logger.info('\nCaculating tetra-nucleotide orders ...')
            p2 = multiprocessing.Process(target=unwrap_kmerOrder, args=(self,))
        else:
            logger.info('Tetra oder file existis, skip!')
            p2 = 0
        
        # generate the tetranucleotide frequency of the contigs
        if not os.path.exists(self.TETRAFREQ) or not os.path.getsize(self.TETRAFREQ):
            logger.info('Caculating tetra-nucleotide frequency ...')
            p3 = multiprocessing.Process(target=unwrap_kmerFreq, args=(self,))
        else:
            logger.info('Tetra freq file existis, skip!')
            p3= 0
        
        # binning the big scaffolds used for mda binning
        if os.path.exists(self.bins):
            logger.info('The binning files of meta-O exists, skip!')
            p4 = 0
        else:
            logger.info('\nBinning the meta-O sequences ...')
            metabat = metaBAT(self.contigfile, self.bins, self.threads)
            p4 = multiprocessing.Process(target=metabat.binning, args=(samFn,))
        
        # get kmer compositon file
        if not os.path.exists(self.kmerFile) or not os.path.getsize(self.kmerFile):
            p5 = multiprocessing.Process(target=unwrap_kmerComp, args=(self,))
        else:
            logger.info('Kmer file existis, skip!')
            p5 = 0
        
        if p1: p1.start()
        if p2: p2.start()
        if p3: p3.start()
        if p4: p4.start()
        if p5: p5.start()
        
        # join
        if p1: p1.join()
        if p2: p2.join()
        if p3: p3.join()
        if p4: p4.join()
        if p5: p5.join()
        
        # delete tmple files
        self.cleanUp()
    
    def revNode(self, node):
        
        if self.index[node] == 0:
            return node + 1
        elif self.index[node] == 2:
            return node - 1
        else:
            return node
    
    def loadScaf(self, genomeNodes=[]):
        """
        Load scaffolds file
        """
        distance = {} # the last nuclotide of the uper node with the first nuclotide of the down node
        scafs, links, ids = [], [], []
        with open(self.scafile, 'r') as af:
            for line in af:
                if '>' in line:
                    ids.append(line[1:-1].split()[0])
                    if links:
                        scafs.append([])
                        for pos in xrange(len(links)-1):
                            up = links[pos]
                            dow = links[pos+1]
                            overlap = 0 - (up[1] - (dow[2] - up[2]))
                            if up[0] not in distance:
                                distance[up[0]] = {}
                            distance[up[0]][dow[0]] = overlap
                            scafs[-1].append(links[pos][0])
                        else:
                            scafs[-1].append(links[-1][0])
                        links = []
                else:
                    info = line.rstrip().split()
                    # unify the direction
                    node, start, direction, length = int(info[0]), int(info[1]), info[2], int(info[3])
                    scafNode = node if direction == '+' else self.revNode(node)
                    links.append([scafNode, length, start])
            else:
                scafs.append([])
                for pos in xrange(len(links)-1):
                    up = links[pos]
                    dow = links[pos+1]
                    overlap = 0 - (up[1] - (dow[2] - up[2]))
                    if up[0] not in distance:
                        distance[up[0]] = {}
                    distance[up[0]][dow[0]] = overlap
                    scafs[-1].append(links[pos][0])
                else:
                    scafs[-1].append(links[-1][0])
        # recover the target genome scaffolds of soapdenovo 
        if genomeNodes:
            genomeScaf = []
            for ele in scafs:
                inter = genomeNodes.intersection(ele)
                if inter:
                    rate = sum([self.length[i] for i in inter])*1.0/sum([self.length[i] for i in ele])
                    if rate >=0.6:
                        genomeScaf.append(ele)
            return genomeScaf, distance
        else:
            return {i[0]:i[1] for i in itertools.izip(ids, scafs)}
    
    def fetchSeqs(self, ids):
        """
        Get the sequence of id
        """
        identifer = set()
        # check the formate of ids
        if isinstance(ids, int):
            identifer.add(ids)
        elif isinstance(ids, set) or isinstance(ids, list) \
                                  or isinstance(ids, tuple):
            try:
                identifer = set([int(i) for i in ids])
            except:
                logger.error('Find erros when fetching the sequences, please check your node ids',
                             exit_with_code = 1
                             )
        
        with open(self.idIndexFile, 'rb') as af:
            idIndex = pickle.load(af)
            update =  pickle.load(af)
        
        if not update:
            for key in idIndex.keys():
                idIndex[self.revNode(key)] = idIndex[key]
            
            with open(self.idIndexFile, 'wb') as af:
                pickle.dump(idIndex, af)
                pickle.dump(1, af)
        
        position = [idIndex[i] for i in identifer if i in idIndex]
        
        # in case of the contigs that is kmer size length and the seq is not given by soapdenovo
        noSeqId = [i for i in identifer if i not in idIndex]
        if noSeqId:
            logger.debug('There are {0} nodes that no seq'.format(len(noSeqId)))
        
        # try to get the seq by link info
        pairs = set()
        for node in noSeqId:
            if node in self.dDlink:
                for nextNode, value in self.dDlink[node].iteritems():
                    if -value == self.kmersize:
                        break
                pairs.add((node, nextNode, 1))
            else:
                rev = self.revNode(node)
                for nextNode, value in self.dDlink[node].iteritems():
                    if -value == self.kmersize:
                        break 
                pairs.add((node, nextNode, -1))
        
        identifer.update([i[1] for i in pairs])
        position.extend([idIndex[i[1]] for i in pairs if i[1] in idIndex])
        position.sort()
        
        sequences = {}
        with open(self.transedContig,'r') as af:
            for pos in position:
                af.seek(pos)
                line = af.readline()
                ae = line.rstrip().split()
                if int(ae[0]) in identifer:
                    sequences[int(ae[0])] = ae[-1].rstrip()
                rev = self.revNode(int(ae[0]))
                if rev in identifer:
                    sequences[rev] = self.revCom(ae[-1].rstrip())
        
        for pair in pairs:
            try:
                if pair[2] == 1:
                    sequences[pair[0]] = sequences[pair[1]][:self.kmersize]
                else:
                    sequences[pair[0]] = self.revCom(sequences[pair[1]][:self.kmersize])
            except:
                continue
        
        return sequences
    
    def kmerComposition(self):
        """
        split the contigs into kmers which used for checking
        if the target genome used later is included by the soap contigs
        """
        logger.info('Splitting seqs into kmers, it will take a while on the first time of running ...')
        ksize = 50
        sbf = ScalableBloomFilter(initial_capacity=500000,
                                  mode=ScalableBloomFilter.LARGE_SET_GROWTH
                                  )
        
        # record the kmers
        with open(self.transedContig, 'r') as af:
            for line in af:
                _, _, _, seq = line.rstrip().split()
                arr = ( seq[i:i+ksize] for i in xrange(len(seq)-ksize+1) )
                for ele in arr:
                    _ = sbf.add(ele)
        
        logger.info('{0:,} types of kmer are generated, done!'.format(len(sbf)))
        with open(self.kmerFile, 'w') as tmpf:
            sbf.tofile(tmpf)
    
    def kmerCMP(self, paras):
        """
        Caculate the fraction of the shared kmers between target genome and soap contigs
        ksize = 50, species: 63.2%, genus: 7.3%, family: 2.3, order: 0.06, class: 0.02, phylum:0.01
        use a conserved species level fraction : 30%
        """
        
        # check if the kmer file exists
        if not os.path.exists(self.kmerFile) or not os.path.getsize(self.kmerFile):
            logger.warning('Kmer file: {0} not exists, re-generated it!'.format(self.kmerFile))
            self.kmerComposition()
        
        if not paras.spa:
            # the seed seqs are SOAPdenovo scaffolds, pass
            logger.info('Becausing SOAPdenovo scaffolds are used, no need to compare kmers, skip ...')
            return 1
        
        # first load kmer file
        logger.info('Loading kmer file: {0}'.format(self.kmerFile))
        with open(self.kmerFile, 'r') as tmpf:
            sbf = ScalableBloomFilter.fromfile(tmpf)
        
        logger.info('Caculating shared kmer fraction ...')
        # read the target genome file
        seq, ksize = '', 50
        matches, total = 0, 0
        
        with open(paras.refn, 'r') as af:
            for line in af:
                if line[0] == '>':
                    if seq:
                        arr = ( seq[i:i+ksize] for i in xrange(len(seq)-ksize+1) )
                        for ele in arr:
                            total += 1
                            if ele in sbf:
                                matches += 1
                            else:
                                rev = self.revCom(ele)
                                if rev in sbf:
                                    matches += 1
                        seq = ''
                else:
                    seq += line.rstrip()
        
        # return the fraction of the shared kmers to judge if it is necessary go on the next steps
        fraction = matches*1.0/total
        logger.info('Shared kmer fraction is {0:.3}.'.format(fraction))
        return fraction

class BAFchecker:
    
    """
    Store the mapping and binning files of current project and also
    direct the Paras(a class) instance find these files
    """
    def __init__(self, projectDir):
        
        self.metaSdir = projectDir
        if not os.path.exists(self.metaSdir):
            os.makedirs(self.metaSdir)
        
        self.metaSdir = os.path.join(self.metaSdir)
        self.erfn = os.path.join(self.metaSdir, config.errFreefn)
    
    def checkFns(self):
        
        # libpaths, identity, fishing file, exclfn
        logger.info('Checking mapping analysis files:')
        
        guildlinksFn =  os.path.join(self.metaSdir, config.GUILDLINKS)
        if os.path.exists(guildlinksFn) and os.path.getsize(guildlinksFn):
            logger.info('Position file: {0} ........................ passed'.format(config.GUILDLINKS))
        else:
            logger.error('File: {0} not exists or size of this file is zero, please check it.'.format(guildlinksFn),
                        to_stderr=True,
                        exit_with_code=3
                         )
        
        fishFn = os.path.join(self.metaSdir, config.FISHCLUSTERS)
        if os.path.exists(fishFn) and os.path.getsize(fishFn):
            logger.info('Position file: {0} ........................... passed'.format(config.FISHCLUSTERS))
        else:
            logger.error('File: {0} not exists or size of this file is zero, please check it.'.format(fishFn),
                        to_stderr=True,
                        exit_with_code=3
                         )
        
        idenFn = os.path.join(self.metaSdir, config.IDENTITY)
        if os.path.exists(idenFn) and os.path.getsize(idenFn):
            logger.info('Position file: {0} .......................... passed'.format(config.IDENTITY))
        else:
            logger.error('File: {0} not exists or size of this file is zero, please check it.'.format(idenFn),
                        to_stderr=True,
                        exit_with_code=3
                         )
    
    def writeFns(self, guilLinks, fishing, identity):
        """
        Write mapping analysis files into MGA directory for later uage
        """
        logger.warning('Caustion: new mapping analysis files are writen and the following files will be overwriten!')
        guildlinksFn =  os.path.join(self.metaSdir, config.GUILDLINKS)
        fishFn = os.path.join(self.metaSdir, config.FISHCLUSTERS)
        idenFn = os.path.join(self.metaSdir, config.IDENTITY)
        
        logger.info('    {0}'.format(guildlinksFn))
        logger.info('    {0}'.format(fishFn))
        logger.info('    {0}'.format(idenFn))
        
        with open(guildlinksFn, 'wb') as af:
            pickle.dump(guilLinks, af)
        
        with open(fishFn, 'wb') as af:
            pickle.dump(fishing, af)
        
        with open(idenFn, 'wb') as af:
            pickle.dump(identity, af)
    
    def ec(self, spadesfn):
        """
        Perform error correction to the MDA data
        """
        if not os.path.exists(self.erfn) or not os.path.getsize(self.erfn):
            logger.info('Detecting error sites on the scaffolds')
            metaSCec.errorScan(spadesfn, self.erfn)
        else:
            logger.warning('Breaked file exists, skipping!')

class Paras:
    
    def __init__(self, ref, metaSdir):
        
        logger.info('Loading FAB files ...')
        if not ref:
            logger.error('Please provide the fasta file of the partial genome which you want to recover',
                         to_stderr=True,
                         exit_with_code=1
                         )
        
        if not os.path.exists(ref):
            logger.error('Reference file: {0} not exist!'.format(ref),
                            to_stderr=True,
                            exit_with_code=1
                            )
        elif not os.path.getsize(ref):
            logger.error('Size of the reference file: {0} is zero!'.format(ref),
                        to_stderr=True,
                        exit_with_code=1
                        )
        self.refn = os.path.abspath(ref)
        self.erfn = os.path.join(metaSdir, config.errFreefn)
        self.metaSdir = metaSdir
        
        # check if the ref are MDA assembly result
        self.renew = 0
        with open(ref, 'r') as af:
            for line in af:
                if line[0] == '>':
                    if 'NODE' in line: # spa marker
                        continue
                    else:
                        try:
                            int(line[1:-1])
                        except:
                            # other character in the header of fasta, not from MDA
                            self.renew = 1
                            break
    
    def loadData(self, soap):
        """
        Load the mapping data
        """
        logger.info('Loading data from {0}'.format(self.metaSdir))
        # read paras file
        guildlinksFn = os.path.join(self.metaSdir, config.GUILDLINKS)
        fishFn = os.path.join(self.metaSdir, config.FISHCLUSTERS)
        idenFn = os.path.join(self.metaSdir, config.IDENTITY)
        
        disapearCounter = 0
        if os.path.exists(guildlinksFn):
            with open(guildlinksFn, 'rb') as af:
                self.guildLinks = pickle.load(af)
        else:
            self.guildLinks = {}
            disapearCounter += 1
        
        if os.path.exists(fishFn):
            with open(fishFn,'rb') as af:
                self.fishCls = pickle.load(af)
        else:
            self.fishCls = {}
            disapearCounter += 1
        
        if os.path.exists(idenFn):
            with open(idenFn,'rb') as af:
                self.identity = pickle.load(af)
        else:
            self.identity = {}
            disapearCounter += 1
        
        self.mislib = {}
        
        if disapearCounter == 4:
            logger.warning('Mapping files are not found!')
    
        if not self.renew:
            self.spa, self.nodes, self.length = set(), set(), {}
            tmp = collections.deque()
            with open(self.refn, 'r') as af:
                for line in af:
                    if '>' in line:
                        tmp.append(line[1:-1])
                        tmp.append('')
                        if len(tmp) == 4:
                            if 'NODE' in tmp[0]:
                                self.spa.add(tmp[0])
                            else:
                                self.length[int(tmp[0])] = len(tmp[1])
                            tmp.popleft()
                            tmp.popleft()
                    else:
                        tmp[-1] += line[:-1]
                
                if 'NODE' in tmp[0]:
                    self.spa.add(tmp[0])
                else:
                    self.length[int(tmp[0])] = len(tmp[1])
            
            self.nodes = set(self.length.keys())
            # get the mean length of the soap nodes
            self.meanLen = int(sum(self.length.values())*1.0/len(self.length))
            
            for ele in self.spa:
                seqlen = int(ele.split('_')[3])
                self.length[ele] = seqlen
            
            self.marker = 'm' # parsing MDA scaf data
            logger.info('Loading done!')
        else:
            # accept only the soapdenovo scaffolds binning results
            self.refids, self.nodes = set(), set()
            with open(self.refn) as af:
                for line in af:
                    if line[0] == '>':
                        if 'scaffold' in line:
                            self.refids.add(line[1:-1])
                        else:
                            self.nodes.add(int(line[2:-1]))
            soapScaf = soap.loadScaf()
            self.nodes.update([k-1 if soap.index[k]==2 else k for i in self.refids for k in soapScaf[i]])
            self.length= {i:soap.length[i] for i in self.nodes}
            self.spa = {}
            self.meanLen = int(sum(self.length.values())*1.0/len(self.length))
            self.marker = 's' # parsing SOAPdenovo scaf data
            logger.info('Loading done!')
    
    def mappedLinks(self, index, gNodes):
        """
        Output the links of the node that mapped to reference.
        """
        newlinks = {}
        for ele, value in self.guildLinks.iteritems():
            order = [self.revNodes(i[-1], index) if i[-2] == "-1" else i[-1] for i in value]
            if ele in self.spa:
                newlinks[ele] = order
            else:
                # alos use the recovered nodes to add more spaids
                if len(set(order) & gNodes)*1.0/len(set(order)) > 0.5:
                    newlinks[ele] = order
        
        # direct link the mapping units by overlap
        occerence = collections.defaultdict(int)
        for key in newlinks:
            s, e = newlinks[key][0], newlinks[key][1]
            rs = self.revNodes(s, index)
            re = self.revNodes(e, index)
            occerence[s] += 1
            occerence[rs] += 1
            occerence[e] += 1
            occerence[re] += 1
        # find overlap betweeen two links
        multId = [x for x in occerence if occerence[x]>1]
        multId.sort()
        overlap, used = [], set()
        for pos in xrange(0,len(multId)-1,2):
            ids = []
            pair = set([multId[pos], multId[pos+1]])
            for ele in newlinks:
                if newlinks[ele][0] in pair:
                    ids.extend([ele, 0])
                elif newlinks[ele][-1] in pair:
                    ids.extend([ele, -1])
            if len(ids) != 4:
                continue
            if ids[1]==ids[3] == 0 and ids[0]!=ids[2]:
                overlap.append([-1, ids[2], 1, ids[0]])
            elif ids[1]==ids[3] == -1 and ids[0]!=ids[2]:
                overlap.append([1, ids[0], -1, ids[2]])
            elif ids[1] == -1 and ids[3] == 0 and ids[0]==ids[2]:
                overlap.append([1, ids[0], 1, ids[2]])
            elif ids[1] == 0 and ids[3] == -1 and ids[0]==ids[2]:
                overlap.append([1, ids[2], 1, ids[0]])
            used.update([ids[0],ids[2]])
        # out
        finalLinks = {}
        for ele in overlap:
            newStrain = []
            newId = ':'.join(map(str, ele))
            if ele[0] == 1:
                newStrain.extend(newlinks[ele[1]][:-1])
            elif ele[0] == -1:
                newStrain.extend([self.revNodes(i, index) for i in newlinks[ele[1]][::-1][:-1]])
            if ele[2] == 1:
                newStrain.extend(newlinks[ele[3]])
            elif ele[2] == -1:
                newStrain.extend([self.revNodes(i, index) for i in newlinks[ele[3]][::-1]])
            finalLinks[newId] = newStrain
        
        for ele in newlinks:
            if ele not in used:
                finalLinks[ele] = newlinks[ele]
        
        return finalLinks
    
    def update(self, soap, recovery):
        """
        add the recovered nodes to paras.nodes and update the files
        """
        # update paras variables
        self.nodes.update(recovery)
        for i in recovery:
            self.length[i] = soap.length[i]
        
        self.meanLen = int(sum(self.length.values())*1.0/len(self.length))
        logger.info('{0} new nodes are recovered.'.format(len(recovery)))
        
        # updating seed seqs
        seqs = soap.fetchSeqs(recovery)
        with open(self.refn, 'aw') as af:
            for key in seqs:
                af.write('>{0}\n{1}\n'.format(key, seqs[key]))
        
        return recovery
    
    def revCom(self, string):
        """
        reverse complement the given string
        """
        comp = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C',
                 'a' : 'T', 't' : 'A', 'c' : 'G', 'g' : 'C',
                 'n' : 'N', 'N' : 'N'
                }
        result = ''
        for i in string[::-1]:
            result += comp[i]
        
        return result
    
    def revNodes(self, node, index):
        
        if index[node] == 0:
            return node + 1
        elif index[node] == 2:
            return node - 1
        else:
            return node
    
    def newSpaSeqs(self, newlinks, soap):
        """
        Refresh the sequences of linked mapping cluster
        and also output singletons
        """
        preSeqs, flag = {}, 0
        with open(self.refn,'r') as af:
            for line in af:
                if 'NODE' in line:
                    seqid = line[1:-1]
                    flag = 1
                else:
                    if flag:
                        preSeqs[seqid] = line.rstrip()
                        flag = 0
        
        addtional = set()
        for key in set(newlinks.keys()) - self.spa:
            if ':' not in key:
                addtional.add(key)
            else:
                ae = key.split(':')
                for pos in xrange(1, len(ae), 2):
                    addtional.add(ae[pos])
        
        flag = 0
        with open(self.erfn) as af:
            for line in af:
                if line[0] == '>':
                    tmpid = line[1:-1]
                    if tmpid in addtional:
                        flag = 1
                else:
                    if flag:
                        preSeqs[tmpid] = line.rstrip()
                        flag = 0
        
        # get plus direction links
        genGuildLinks, tmpids = {}, set()
        for spaId in preSeqs:
            if spaId in self.guildLinks:
                order = []
                for ele in self.guildLinks[spaId]:
                    tmpids.add(ele[-1])
                    if ele[-2] == '1':
                        newEle = []
                        newEle.extend(ele[2:4])
                        newEle.extend([ele[5], ele[-1]])
                        order.append(newEle)
                    else:
                        rev = soap.revNode(ele[-1])
                        newEle = []
                        newStart = ele[5]-ele[2]+1
                        newEnd = ele[5]-ele[3]+1
                        newEle.extend([newStart, newEnd, ele[5]])
                        newEle.append(rev)
                        order.append(newEle)
                genGuildLinks[spaId] = order
        
        # get soap contig seqs
        tmpids.update([soap.revNode(i) for i in tmpids])
        tmpseq = soap.fetchSeqs(tmpids)
        for key, value in tmpseq.iteritems():
            preSeqs[key] = value
        
        # generate new sequences
        aftSeqs = {}
        for ele in newlinks:
            if ':' not in ele:
                aftSeqs[ele] = preSeqs[ele]
            else:
                info, allOrder, linkedId = ele.split(':'), [], []
                for pos in xrange(0, len(info),2):
                    linkedId.append(info[pos+1])
                    if info[pos] == '1':
                        allOrder.append(genGuildLinks[info[pos+1]])
                    else:
                        contigOrder = collections.deque()
                        for key in genGuildLinks[info[pos+1]]:
                            revOrder = []
                            newStart = key[2]-key[0]+1
                            newEnd = key[2]-key[1]+1
                            revOrder.extend([newEnd, newStart, key[2]])
                            revId = soap.revNode(key[-1])
                            revOrder.append(revId)
                            contigOrder.appendleft(revOrder)
                        allOrder.append(contigOrder)
                # fill the gap
                finalOrder= []
                for pos in xrange(1, len(allOrder)):
                    bef = allOrder[pos-1][-1]
                    aft = allOrder[pos][0]
                    if bef[-1] != aft[-1]:
                        logger.notice('Conflict when merging spades scafs')
                    gap = aft[0] - bef[1]
                    if gap > 0:
                        overlap = preSeqs[bef[-1]][bef[1]:aft[0]+1]
                        finalOrder.append(overlap)
                    else:
                        finalOrder.append(gap)
                # generate seqs
                seq = preSeqs[linkedId[0]]
                for pos in xrange(1,len(linkedId)):
                    if isinstance(finalOrder[pos-1], int):
                        seq += preSeqs[linkedId[pos]][abs(finalOrder[pos-1]):]
                    else:
                        seq += finalOrder[pos-1]
                        seq += preSeqs[linkedId[pos]]
                aftSeqs[ele] = seq
        
        return aftSeqs
