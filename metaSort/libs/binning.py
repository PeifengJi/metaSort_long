############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import config
import collections
import subprocess
import copy
import itertools
import cPickle as pickle
import shutil
import platform
import initialize
from log import get_logger

# set log file
logger = get_logger()

def subgraphs(graph):
    """
    split the graph into several small subgraphs
    """
    count = 1
    tmp = []
    cl = {i:1 for i in graph.iterkeys()}
    for node in cl:
        if cl[node] is 1:
            tmp.append(node)
            while tmp:
                last = tmp.pop()
                cl[last] = count + 1
                for second in graph[last]:
                    if cl[second] is 1:
                        tmp.append(second)
            count += 1
    group = {}
    for key in cl:
        if cl[key] not in group:
            group[cl[key]] = set()
        group[cl[key]].add(key)
    return group 

def reverseNodes(node, index):
    
    if index[node] is 0:
        return node + 1
    elif index[node] is 2:
        return node - 1
    else:
        return node

def revCom(string):
    """
    reverse complement the given string
    """
    comp = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'}
    result = ''
    for i in string[::-1]:
        result += comp[i]
    
    return result

def tetracalc(fastafn, threads):
    """
    Use the published tools to binning the contigs
    """
    # Vaughn iverson et al. 2012 SCIENCE SEAStAR software
    logger.info('Binning with parameters: merge_tart = {0}'.format(config.MERGELEVEL))
    tool = 'tetracalc'
    tooldir = 'tetra'
    if config.platform_name in ['Darwin', 'macosx']:
        toolPath = os.path.join(config.libsLocation, tooldir, 'macosx', tool)
    else:
        toolPath = os.path.join(config.libsLocation, tooldir, 'linux', tool)
    
    if not os.path.exists(toolPath):
        logger.error('Can not find script: tetracalc, please check it.',
                     exit_with_code=1
                     )
    
    cmd = [toolPath,
           '--merge_tar', str(config.MERGELEVEL),
           '-m', str(config.MINLEN),
           '--num_threads', str(threads),
           fastafn
           ]
    # output cmd
    logger.print_command_line(cmd, wrap_after=None)
    
    child = subprocess.Popen(cmd, stdout = subprocess.PIPE)
    out, err = child.communicate()
    if err:
        sys.stderr.write(out)
        logger.error('Error happened while binning.',
                     exit_with_code=1
                     )
    else:
        line = out.replace('"','')
        line = line.replace('{clusters : [[','')
        line = line.replace(']]}\n','')
        info = line.split('],[')
    
    return [ele.split(',') for ele in info]

class metaBAT:
    """
    metaBAT binning
    """
    def __init__(self, reference, threads=1):
        
        self.refn = os.path.abspath(reference)
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
    
    def binning(self, samfn = False, covFile = False):
        """
        Use metaBAT bins the scaffolds
        """
        binPath = self.checkPlatform()
        
        if samfn and not covFile:
            if self.checkSamtools():
                # go on
                logger.info('samtools is detected, go on ...')
            else:
                logger.warning('Does not find samtools, metaBAT could not be used !')
                return 0
        
        curDir = os.getcwd()
        output = 'temple_dir'
        if not os.path.isdir(output):
            os.makedirs(output)
        else:
            shutil.rmtree(output)
            os.makedirs(output)
        
        if os.path.exists(covFile) and os.path.getsize(covFile):
            logger.info('File: {0} is found, jump to binning step directly ...'.format(covFile))
            shutil.copy(covFile, output)
            os.chdir(output)
        else:
            if os.path.exists('sorted.bam') and os.path.getsize('sorted.bam'):
                logger.info('File: sorted.bam is found, skip sorting sam step ...')
            else:
                logger.info('Generating bam files, please wait ...')
                self.sortSAM(samfn)
            if os.path.exists(samfn) and os.path.getsize(samfn):
                os.remove(samfn)
        
            shutil.move('sorted.bam', output)
            os.chdir(output)
            
            logger.info('Caculating contig coverage, please wait ...')
            # calculating contig depth
            calDepth = os.path.join(binPath, 'jgi_summarize_bam_contig_depths')
            cmd = [calDepth, '--outputDepth', 'contig_cov.txt', 'sorted.bam']
            retunCode = subprocess.call(cmd)
        
        # binning
        logger.info('Running metaBAT, please wait ...')
        prefix = 'cluster'
        getBins = os.path.join(binPath, 'metabat')
        cmd = [getBins, '-o', prefix, '-i', self.refn, '-a', 'contig_cov.txt', '-t',
               self.threads, '--saveCls', '--superspecific', '-m', str(config.metaS_minLen)
               ]
        retunCode = subprocess.call(cmd)
        
        cluster = {}
        with open(prefix,'r') as af:
            for line in af:
                ae = [i for i in line.rstrip().split()]
                if ae[1] not in cluster and int(ae[1]) > 0:
                    cluster[ae[1]] = set()
                elif ae[1] in cluster:
                    cluster[ae[1]].add(ae[0])
        
        logger.info('binning done !')
        os.chdir(curDir)
        shutil.rmtree(output)
        
        return cluster

def nucmerPaser(mapfile):
    """
    Parse the nucmer mapping file, get overlap, link info
    """
    print('importing marked nodes from nucmer mapping result ...')
    #       0        1        2         3           4           5          6         
    #   refstart, refend, querystart, queryend, refmatches, querymathces, identity,  
    #       7        8        9        10       11           12        13      14
    #    reflen, querylen, refcov, querycov, refdirect, querydirect, refid, queryid
    # whole mapping
    identity = {}
    occrence, wholeAlign = collections.defaultdict(list), []
    with open(mapfile,'r') as b_file:
        for i in xrange(4):
            b_file.readline()
        for line in b_file:
            info = line.rstrip().split()
            if float(info[6]) >= 95 and int(info[7]) >= config.metaS_minLen:
                identity[int(info[-1])] = float(info[6])
                if float(info[9]) >= 99:
                    wholeAlign.append(line)
                elif float(info[10]) >= 99 and int(info[8]) >= 2000:
                    wholeAlign.append(line)
                    occrence[info[-1]].append(line)
    
    # filter multiple mppping seq, the seq will not be included in graph
    multIds = set([x for x in occrence if len(occrence[x])>1])
    logger.info()
    logger.info('Mapping statistics:')
    logger.info('#'*60)
    logger.info('Number of multiple mapping:                    {0:,}'.format(len(multIds)))
    
    overlaps = {}
    for line in wholeAlign:
        ae = line.rstrip().split()
        if float(ae[10]) >= 99:
            pos = [int(ae[0]), int(ae[1])] if int(ae[0]) < int(ae[1]) else [int(ae[1]), int(ae[0])]
            info = [pos[0], pos[1], int(ae[-1])]
            if ae[-2] not in overlaps:
                overlaps[ae[-2]] = []
            overlaps[ae[-2]].append(info)
    
    # build the graph
    graph = {}
    for line in copy.deepcopy(wholeAlign):
        link = line.split()
        if link[-1] not in multIds:
            if link[-1] not in graph:
                graph[link[-1]] = {}
            if link[-2] not in graph:
                graph[link[-2]] = {}
            graph[link[-1]][link[-2]] = graph[link[-2]][link[-1]] = 1
        else:
            wholeAlign.remove(line)
    
    clusters = subgraphs(graph)
    
    return clusters, identity, {}

def MDAbinning(fastafn, soap, threads = '1'):
    """
    Use metaBAT to binning the sequences
    """
    prefixDir = os.path.split(os.path.abspath(fastafn))[0]
    metabat = metaBAT(fastafn, threads)
    if os.path.exists('contig_cov.txt') and os.path.getsize('contig_cov.txt'):
        info = metabat.binning(covFile = 'contig_cov.txt')
    else:
        # using metaBAT for binning
        sam = initialize.Pipeline(fastafn)
        sam.indexRef()
        for key in soap.libs:
            prefix = '{1}/lib{0}'.format(key, prefixDir)
            for pair in soap.libs[key]['reads']:
                logger.info('Mapping files {0} and {1}'.format(pair[0], pair[1]))
                sam.algin(pair[0], pair[1], prefix, str(threads))
        logger.info('Done!')
        info = metabat.binning(samfn = '{0}.sam'.format(prefix))
    
    return info

def fishing(extendCls, mdaCls, mdaSeqfn, soap):
    """
    Use the fishing and clustering combined method to binning the sequences
    """
    # use cls info to link the different bins
    directGraph = {}
    for line in mdaCls.itervalues():
        for pair in itertools.combinations(line, 2):
            if pair[0] not in directGraph:
                directGraph[pair[0]] = {}
            if pair[1] not in directGraph:
                directGraph[pair[1]] = {}
            directGraph[pair[0]][pair[1]] = directGraph[pair[1]][pair[0]] = 1
    
    # add the cls infomation
    for line in extendCls.itervalues():
        for pair in itertools.combinations(line, 2):
            if pair[0] not in directGraph:
                directGraph[pair[0]] = {}
            if pair[1] not in directGraph:
                directGraph[pair[1]] = {}
            directGraph[pair[0]][pair[1]] = directGraph[pair[1]][pair[0]] = 1
    
    groups = subgraphs(directGraph)
    for i in groups.keys():
        groups['s' + str(i)] = groups[i]
        del groups[i]
    
    # merge the metaO binning results, by calculating the overlapped soap contigs
    with open(soap.clusters, 'rb') as af:
        soapClusters = pickle.load(af)
    
    # discard large clusters
    lenCutoff = 4000000 # 8 M
    for key in soapClusters.keys():
        accLen = sum([soap.length[i] for i in soapClusters[key]])
        if accLen >= lenCutoff:
            del soapClusters[key]
    
    passed = set()
    for key, value in groups.iteritems():
        size1 = sum([soap.length[int(i)] for i in value if 'N' not in i and soap.length[int(i)] >= config.metaO_minLen])
        size2 = sum([int(i.split('_')[3]) for i in value if 'N' in i and int(i.split('_')[3]) >= config.metaS_minLen])
        if size1 > lenCutoff or size2 > lenCutoff:
            passed.add(key)
        elif size2 <= 600000:
            passed.add(key)
    
    link = {}
    for pair in itertools.product(groups.keys(), soapClusters.keys()):
        if pair[0] not in passed:
            trans = set([int(i) for i in groups[pair[0]] if 'N' not in i and soap.length[int(i)] >= config.metaO_minLen])
            overlap = len(trans & soapClusters[pair[1]])
            if overlap == len(trans):
                if pair[0] not in link:
                    link[pair[0]] = {}
                if pair[1] not in link:
                    link[pair[1]] = {}
                link[pair[0]][pair[1]] = link[pair[1]][pair[0]]= overlap
    subLinks = subgraphs(link)
    
    # output the final clusters
    merged, count = {}, 1
    for ele in subLinks.itervalues():
        new = set()
        for i in ele:
            if i[0] == 's':
                new.update([int(k) if 'N' not in k else k for k in groups[i]])
            else:
                new.update([k for k in soapClusters[i]])
        # count size
        size = 0
        for i in new:
            try:
                size += soap.length[int(i)]
            except:
                size += int(i.split('_')[3])
        if  size >= config.SIZECUF:
            key = 'm' + str(count)
            merged[key] = new
            count += 1
    
    for ele in set(groups.keys()) - set(link.keys()):
        new = set([int(k) if 'N' not in k else k for k in groups[ele]])
        # count size
        size = 0
        for i in new:
            try:
                size += soap.length[int(i)]
            except:
                size += int(i.split('_')[3])
        if size >= config.SIZECUF:
            key = 'm' + str(count)
            merged[key] = new
            count += 1
    
    # output the final sequences
    prefix = config.CLUSTERSEQDIR
    if os.path.isdir(prefix):
        shutil.rmtree(prefix)
    else:
        os.mkdir(prefix)
    
    # fetch sequences
    soapIds, mdaIds = set(), set()
    for ele in merged.itervalues():
        for i in ele:
            try:
                soapIds.add(int(i))
            except:
                mdaIds.add(i)
    
    soapSeqs = soap.fetchSeqs(soapIds)
    flag, mdaSeqs = 0, {}
    with open(mdaSeqfn, 'r') as af:
        for line in af:
            if line[0] == '>':
                if line.rstrip()[1:] in mdaIds:
                    mid = line.rstrip()[1:]
                    flag = 1
            else:
                if flag:
                    mdaSeqs[mid] = line.rstrip()
                flag = 0
    
    for key, value in merged.iteritems():
        with open(os.path.join(prefix, key + '.fa'), 'w') as af:
            for ele in value:
                try:
                    int(ele)
                    af.write('>{0}\n{1}\n'.format(ele, soapSeqs[ele]))
                except:
                    af.write('>{0}\n{1}\n'.format(ele, mdaSeqs[ele]))
    return 1
