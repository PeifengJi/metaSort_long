############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import config
import mapper
import platform
import itertools
import collections
#import scipy.optimize
import cPickle as pickle
from log import get_logger
# set log file
logger = get_logger()
import prepareData
import simplify
import bubble
import pathRecover

def covRange(arr):
    """
    Merge the coved range of query seq on the reference
    """
    new = collections.deque(sorted(arr))
    uniq = { i[0]:i[1] for i in new}
    mappedRegion = []
    while new:
        first = new.popleft()
        if new:
            for pos, ele in enumerate(new):
                if ele[0] <= uniq[first[0]] <= uniq[ele[0]]:
                    uniq[first[0]] = uniq[ele[0]]
                    if pos == len(new) - 1:
                        mappedRegion.append([first[0], uniq[first[0]]])
                        new = []
                        break
                elif uniq[first[0]] < ele[0]:
                    mappedRegion.append([first[0], uniq[first[0]]])
                    for i in range(pos):
                        new.popleft()
                    break
        else:
            mappedRegion.append([first[0], uniq[first[0]]])
            break
        
    return mappedRegion

def shareRef(ref, c1, c2):
    """
    Extract the shared ref region of two sets of contigs
    """
    logger.info('Mapping contigs to reference ...')
    
    p1 = os.path.basename(c1).split('.')[0]
    p2 = os.path.basename(c2).split('.')[0]
    if not os.path.exists(p1+'.result'):
        mapper.pipeline(ref, c1, p1)
    else:
        logger.info('Mapping file: {0} exists, skip ...'.format(p1+'.result'))
    
    if not os.path.exists(p2+'.result'):
        mapper.pipeline(ref, c2, p2)
    else:
        logger.info('Mapping file: {0} exists, skip ...'.format(p2+'.result'))
    # fetch the mapped region on the ref of c1
    c1ref = {}
    with open(p1+'.result', 'r') as af:
        for i in xrange(4): _ = af.readline()
        for line in af:
            info = line.split()
            start, end, identity, cov, refid = int(info[0]), int(info[1]), float(info[6]), float(info[10]), info[-2]
            if identity >= 95 and cov >= 95:
                if refid not in c1ref:
                    c1ref[refid] = []
                c1ref[refid].append([start, end])
    c1map = {k: covRange(v) for k, v in c1ref.iteritems()}
    
    # fetch the mapped region on the ref of c2
    c2ref = {}
    with open(p2+'.result', 'r') as af:
        for i in xrange(4): _ = af.readline()
        for line in af:
            info = line.split()
            start, end, identity, cov, refid = int(info[0]), int(info[1]), float(info[6]), float(info[10]), info[-2]
            if identity >= 95 and cov >= 95:
                if refid not in c2ref:
                    c2ref[refid] = []
                c2ref[refid].append([start, end])
    c2map = {k: covRange(v) for k, v in c2ref.iteritems()}
    
    # fetch the shared regions
    sharedMap = {}
    for k in c1map:
        if k in c2map:
            tmp = collections.defaultdict(int)
            for reg in c1map[k]:
                for x in xrange(reg[0], reg[1] + 1):
                    tmp[x] += 1
            for reg in c2map[k]:
                for x in xrange(reg[0], reg[1] + 1):
                    tmp[x] += 1
            nums = sorted([i for i, v in tmp.iteritems() if v == 2])
            # dispart the number into continously linked regions
            conRegion = [[nums[0]]]
            for pos in xrange(1, len(nums)):
                if nums[pos] - nums[pos-1] == 1:
                    if pos == len(nums) - 1:
                        conRegion[-1].append(nums[pos])
                    continue
                else:
                    if nums[pos-1] - nums[pos-2] == 1:
                        conRegion[-1].append(nums[pos-1])
                    conRegion.append([nums[pos]])
            # filter the regions that small than 200 bp
            sharedMap[k] = [i for i in conRegion if len(i) == 2 and i[1]-i[0]>=200]
    
    #fetch the seqs of the regions
    logger.info('Outputing shared regions')
    outfn = open('sharedRegions.fa', 'w')
    flag, seq = 0, ''
    with open(ref, 'r') as af:
        for line in af:
            if line[0] == '>':
                ae = line.rstrip().split()[0]
                if ae[1:] in sharedMap:
                    if seq:
                        for reg in sharedMap[tid]:
                            nid = '>{0}_{1}_{2}'.format(tid, reg[0], reg[1])
                            outfn.write(nid + '\n')
                            outfn.write(seq[reg[0]: reg[1]+1] + '\n')
                            seq = ''
                    tid = ae[1:]
                    flag = 1
                else:
                    flag = 0
            else:
                if flag == 1:
                    seq += line.rstrip()
        else:
            if seq:
                for reg in sharedMap[tid]:
                    nid = '>{0}_{1}_{2}'.format(tid, reg[0], reg[1])
                    outfn.write(nid + '\n')
                    outfn.write(seq[reg[0]: reg[1]+1] + '\n')
    outfn.close()

def findVars(nref, contig, outfn):
    """
    Mapping the contigs to shared reference, and find bubbles, respectively
    """
    # new mapping
    logger.info('Parsing project: {0}'.format(contig))
    nref = 'sharedRegions.fa'
    prefix = os.path.basename(contig).split('.')[0] + '_2'
    if not os.path.exists(prefix+'.result'):
        mapper.pipeline(nref, contig, prefix)
    else:
        logger.info('Mapping file: {0} exists, skip ...'.format(prefix+'.result'))
    
    # fetch the mapped nodes
    c1Nodes = set()
    with open(prefix+'.result', 'r') as af:
        for i in xrange(4):
            _ = af.readline()
        for line in af:
            info = line.rstrip().split()
            identity, cov, gid = float(info[6]), float(info[10]), int(info[-1])
            if identity >= 95 and cov >= 95:
                c1Nodes.add(gid)
    # get the soap path
    soapdir = os.path.dirname(contig)
    soap = prepareData.SOAP(soapdir)
    soap.loadData()
    c1Nodes.update([i+1 for i in c1Nodes if soap.index[i]==0])
    targetGraph = pathRecover.genomeGraph(soap, c1Nodes, set(), step=1)
    simplify.merge(targetGraph, soap, {}, set())
    with open('{0}Graph'.format(prefix), 'wb') as af:
        pickle.dump(targetGraph, af)
    
    with open('variations', 'rb') as af:
        _ = pickle.load(af)
        _ = pickle.load(af)
        _ = pickle.load(af)
        indexedBubs = pickle.load(af)
    
    lenRate, density, distance = bubble.bubFeatureCal(targetGraph, soap, indexedBubs)
    with open(outfn, 'wb') as af:
        pickle.dump(lenRate, af)
        pickle.dump(density, af)
        pickle.dump(distance, af)
    
    logger.info('Analysis done and data is writen into file: {0}'.format(outfn))

