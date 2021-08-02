############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import math
import os
import itertools
import config
import collections
from log import get_logger
# set log file
logger = get_logger()

def likelihood(freq, clNum, clWeight, clMean, clSD):
    """
    freq: frequency of every number
    clNum: number of clusters
    clWeight: (k size) weight of every cluster
    clMean: (k size) mean of every cluster
    clSD: (k size) standard deviation of every cluster
    maxLikehood: maximum likelihood estimation
    """
    pi, maxLikehood = 3.1416, 0
    epsilon = 10**-10
    for ele in freq:
        tmpl = 0
        for num in xrange(clNum):
            if math.exp(-(ele-clMean[num])**2/2.0/clSD[num]**2) < epsilon:
                tmpl += clWeight[num]*epsilon/math.sqrt(2*pi*clSD[num]**2)
            else:
                tmpl += clWeight[num]*math.exp(-(ele-clMean[num])**2/2.0/clSD[num]**2)/math.sqrt(2*pi*clSD[num]**2)
        maxLikehood += freq[ele]*math.log(tmpl)
    
    return maxLikehood

def GMM(freq, clNum, ltol, maxiter):
    """
    freq: frequency of every number
    clNum: number of clusters
    ltol: ltol: the boundary of the ratio between the MLE of two steps
    maxiter: the boundary of the iteration
    
    clWeight: (k size) weight of every cluster
    clMean: (k size) mean of every cluster
    clSD: (k size) standard deviation of every cl
    Ll: maximum likelihood estimation
    numiter: the number of iteration
    """
    # init
    clMean, clSD, clWeight = [], [], []
    sumAll, sumForM, sumForV, sumForW = 0, 0, 0, 0
    pi, epsilon = 3.1416, 10**-10
    sumOfX = sum(freq.itervalues())
    
    # set init value
    sumtemp, counter = 0, 1
    for ele in sorted(freq.keys()):
        sumAll += freq[ele]
        sumForM += ele*freq[ele]
        sumForV += ele**2*freq[ele]
        sumForW += freq[ele]
        if sumAll >= counter*1.0/clNum*sumOfX:
            clMean.append(sumForM*1.0/sumForW)
            clSD.append(math.sqrt(sumForV*1.0/sumForW-(sumForM*1.0/sumForW)**2))
            clWeight.append(sumForW*1.0/sumOfX)
            counter += 1
            # init next
            sumForM = 0
            sumForV = 0
            sumForW = 0
    
    # main
    Ll = likelihood(freq, clNum, clWeight, clMean, clSD)
    Ln = Ll * 2
    numiter = 1
    expect = {}
    while abs(Ll*1.0/Ln-1)>ltol and numiter<=maxiter:
        # E-step
        for num in xrange(clNum):
            if abs(clSD[num]) < epsilon:
                clSD[num] = epsilon
            for ele in freq:
                if ele not in expect:
                    expect[ele] = [0 for i in xrange(clNum)]
                if math.exp(-(ele-clMean[num])**2/2.0/clSD[num]**2) < epsilon:
                    expect[ele][num] = clWeight[num]*epsilon/math.sqrt(2*pi*clSD[num]**2)
                else:
                    expect[ele][num] = clWeight[num]*math.exp(-(ele-clMean[num])**2/2.0/clSD[num]**2)/math.sqrt(2*pi*clSD[num]**2)
        
        eSum = collections.defaultdict(float)
        for num in xrange(clNum):
            for ele in freq:
                eSum[ele] += expect[ele][num]
        
        for num in xrange(clNum):
            for ele in freq:
                expect[ele][num] = expect[ele][num]*1.0/eSum[ele]
        
        # M step
        clMean = [0 for i in xrange(clNum)]
        clSD = [0 for i in xrange(clNum)]
        clWeight = [0 for i in xrange(clNum)]
        for num in xrange(clNum):
            for ele in freq:
                clWeight[num] += expect[ele][num]*freq[ele]
        
        for num in xrange(clNum):
            for ele in freq:
                clMean[num] += expect[ele][num]*ele*freq[ele]
        
        for num in xrange(clNum):
            clMean[num] = clMean[num]*1.0/clWeight[num]
        
        for num in xrange(clNum):
            dxM = {ele: ele-clMean[num] for ele in freq}
            dXsum = sum([freq[ele]*dxM[ele]**2*expect[ele][num] for ele in freq])
            clSD[num] = math.sqrt(dXsum*1.0/clWeight[num])
        
        for ele in xrange(clNum):
            clWeight[ele] = clWeight[ele]*1.0/sumOfX
        
        # cal loglikelihood
        Ln = Ll
        Ll = likelihood(freq, clNum, clWeight, clMean, clSD)
        numiter += 1
        
        return clWeight, clMean, clSD, Ll, numiter

def BIC(Ll, k, numSample):
    """
    Bayesian Information Criterion (BIC)
    The BIC is defined as:
        BIC = 2 * ln(L) - k * ln(n)
        where:
           * ln(L) is the log-likelihood of the estimated model
            * k is the number of component
            * n is the number of sample
    """
    return 2*math.log(-Ll, math.e) - k*math.log(numSample, math.e)

def pelinkCreat(soap, totalReads, sizeOfLib):
    """
    First read the mapping info from readOnContigs file, estimate the gap size by using the pairs that mapped
    to single contig, and delete the falase mapped reads, then use the estimated parameters to evaluate the gap
    size, between two linked contigs by one pair
    """
    # store parameters
    rangeOfLib = config.RANGEOFLIB
    rateForCuting = config.RATE4CUT
    pairsForDeleting = config.PAIRS4DEL
    pairsForIgnoring = config.PAIRS4IGN
    rateForSpliting = config.RATE4SPLIT
    rateForRemoving = config.RATE4REM
    rateForDivision = config.RATE4DIV
    
    fout = open(soap.libStat,'w')
    fout.write('Lib_id\tReads\tMean\tSD\tHigher_noise_ratio\tLower_noise_ratio\n')
    
    # check the readOnContigFn file format
    if soap.readOnContigFn[-2:] == 'gz': # gz compressed file
        import gzip
        fin = gzip.open(soap.readOnContigFn, 'rt')
    else:
        fin = open(soap.readOnContigFn, 'r')
    _ = fin.readline()
    
    #init libray
    totReadPairs, sumOfIns, sumOfInsSquare = 0, 0, 0
    insLength = collections.defaultdict(int)
    rec, temp = {}, [0]
    miuOfLib, sigmaOfLib = {}, {}
    rateOfHigherNoise, rateOfLowerNoise = {}, {}
    for lib in totalReads:
        miuOfLib[lib] = 0
        sigmaOfLib[lib] = 0
        rateOfHigherNoise[lib] = 0
        rateOfLowerNoise[lib] = 0
    
    # scan the file
    numOfLib = 1
    prelink = {}
    for line in fin:
        # array: ae -> 1 read id; 2 mapped contig id; 3 5'_position on contig of this read
        ae = map(int, line.rstrip().split()[:3])
        # record the paired reads mapped to same contig or different contigs, do not handel and record single mapped reads
        if ae[0] % 2:
            temp = ae[:]
        else:
            if ae[0] == temp[0] + 1:
                # paired reads
                if (soap.index[ae[1]], ae[1]) in [(0, temp[1]-1), (2, temp[1]+1), (1, temp[1])] \
                        and soap.length[ae[1]]-ae[2] > sizeOfLib[numOfLib]*rangeOfLib \
                        and soap.length[temp[1]]-temp[2] > sizeOfLib[numOfLib]*rangeOfLib:
                    # two reads on the same contig
                    value = soap.length[ae[1]] - ae[2] - temp[2] # the insert length of paired reads
                    insLength[value] += 1
                    totReadPairs += 1
                    sumOfIns += value
                    sumOfInsSquare += value**2
                else:
                    # paired reads mapped to different contig. The gap will be estimated after evaluating the library insert length
                    if ae[1] >= temp[1]:
                        tmpId = (ae[1], temp[1])
                        #tmpId = str(ae[1]) + '-' + str(temp[1])
                        if tmpId not in rec:
                            rec[tmpId] = []
                        rec[tmpId].append([ae[2], temp[2], ae[2]+temp[2]])
                    else:
                        tmpId = (temp[1], ae[1])
                        #tmpId = str(temp[1]) + '-' + str(ae[1])
                        if tmpId not in rec:
                            rec[tmpId] = []
                        rec[tmpId].append([temp[2], ae[2], ae[2]+temp[2]])
            temp = [0]
        
        # handle all of the reads in one library
        if ae[0] >= totalReads[numOfLib]:
            logger.info('Handling library: {0}, insert length: {1}'.format(numOfLib, sizeOfLib[numOfLib]))
            # get the mean and sigma of the sample of the population
            aboveInsPairs, belowInsPairs= 0, 0
            miuOfInitIns = sumOfIns*1.0/totReadPairs
            sigmaOfInitIns = math.sqrt(sumOfInsSquare*1.0/totReadPairs-miuOfInitIns**2)
            
            # remove outliers
            outliers = set()
            for key, value in insLength.iteritems():
                if key > miuOfInitIns + 4*sigmaOfInitIns:
                    aboveInsPairs += value
                    outliers.add(key)
                elif key < miuOfInitIns - 4*sigmaOfInitIns:
                    belowInsPairs += value
                    outliers.add(key)
            
            for ele in outliers:
                del insLength[ele]
            
            # find best cluster number
            minL = 0
            for ele in [1, 2, 3]:
                clWeight, clMean, clSD, Ll, numiter = GMM(insLength, ele, 0.001, 100)
                min2 = -Ll + ele*math.log(totReadPairs - aboveInsPairs - belowInsPairs)
                if not minL or minL> min2:
                    minL = min2
                    bestclMean, bestclSD = clMean, clSD # both of them are list
            
            # find the best gmm cluster that can present the library
            # the best gmm cluster is the cluster with the minmum distance with idear lib size
            sortedCombs = sorted(itertools.izip(bestclMean, bestclSD),  key=lambda x: abs(x[0]-sizeOfLib[numOfLib]))
            miuOfLib[numOfLib], sigmaOfLib[numOfLib] = sortedCombs[0]
            
            # other cluster is noise and get the proportion of the noise
            for ele in insLength:
                if ele > miuOfLib[numOfLib] + 3*sigmaOfLib[numOfLib]:
                    aboveInsPairs += insLength[ele]
                if ele < miuOfLib[numOfLib] - 3*sigmaOfLib[numOfLib]:
                    belowInsPairs += insLength[ele]
            
            rateOfHigherNoise[numOfLib] = aboveInsPairs*1.0/totReadPairs
            rateOfLowerNoise[numOfLib] = belowInsPairs*1.0/totReadPairs
            
            # wirte the parameters of each the lib
            fout.write('{0}\t{1}\t{2:.2}\t{3:.2}\t{4:.4}\t{5:.4}\n'.format(
                                numOfLib,
                                totalReads[numOfLib],
                                miuOfLib[numOfLib],
                                sigmaOfLib[numOfLib],
                                rateOfLowerNoise[numOfLib],
                                rateOfHigherNoise[numOfLib]
                                                        ))
            
            # use the estimated mean and sigma to partion and denoise the reads that mapped to different contigs
            for contigs in rec:
                c1, c2 = contigs
                newarr = sorted(rec[contigs], key=lambda x:x[2])
                cutLower = int(len(newarr)*rateOfHigherNoise[numOfLib]*rateForCuting)
                cutHigher = int(len(newarr)*rateOfLowerNoise[numOfLib]*rateForCuting)
                if cutLower >= 1:
                    newarr = newarr[cutLower:]
                if cutHigher >= 1:
                    newarr = newarr[:len(newarr)-cutHigher]
                
                sumOfnewarr, sumOfnewarrSquare = 0, 0
                for ele in newarr:
                    sumOfnewarr += ele[2]
                    sumOfnewarrSquare += ele[2]**2
                miuOfnewarr = sumOfnewarr*1.0/len(newarr)
                sigmaOfnewarr = math.sqrt(sumOfnewarrSquare*1.0/len(newarr)-miuOfnewarr**2)
                
                if (soap.index[c1], c1) in [(0, c2-1), (1, c2), (2, c2+1)] \
                          and soap.length[c1]-miuOfnewarr-2 < miuOfLib[numOfLib]+sigmaOfnewarr \
                          and soap.length[c1]-miuOfnewarr+2 > miuOfLib[numOfLib]-sigmaOfnewarr:
                    # wipe the pairs on the same contig with the right insert length
                    pass
               # support number if enough to estimate the gap
                elif len(newarr) > pairsForDeleting:
                     # compar the sample sd with population sigma, they should be similar, if not noise ratio is large, remove noise
                    if sigmaOfnewarr <= rateForSpliting*sigmaOfLib[numOfLib]:
                        # situation 1: the sample noise is low
                        # get the both ends (inner) mapping position
                        max1, max2 = newarr[0][0:2]
                        for pos in newarr[1:]:
                            if max1 < pos[0]:
                                max1 = pos[0]
                            if max2 < pos[1]:
                                max2 = pos[1]
                        
                        if soap.length[c1]-max1 > miuOfLib[numOfLib]+rateForRemoving*sigmaOfLib[numOfLib] or \
                             soap.length[c2]-max2 > miuOfLib[numOfLib]+rateForRemoving*sigmaOfLib[numOfLib]:
                            pass
                        else:
                            if c1 not in prelink:
                                prelink[c1] = []
                            prelink[c1].append([c2,
                                                soap.length[c1]-max1,
                                                soap.length[c2]-max2,
                                                len(newarr),
                                                len(newarr)*miuOfnewarr,
                                                miuOfLib[numOfLib]-soap.length[c1]-soap.length[c2]+miuOfnewarr,
                                                sigmaOfnewarr,
                                                numOfLib
                                                ])
                    elif sigmaOfnewarr > rateForSpliting*sigmaOfLib[numOfLib] and len(newarr) > pairsForIgnoring:
                        # situation 2: the sample noise is high
                        # group the insert length into serveral groups, by distance
                        cl = [-1]
                        for pos in xrange(len(newarr)-1):
                            if newarr[pos+1][2]-newarr[pos][2] > rateForDivision*sigmaOfLib[numOfLib]:
                                cl.append(pos)
                        cl.append(len(newarr)-1)
                        
                        removeall = 0
                        for pos in xrange(1, len(cl)):
                            bestclNum, templist, sumOfEachGroup, sumOfEachGroupSquare = 0, 0, 0, 0
                            numOfEachGroup = cl[pos] - cl[pos-1] # get the number of read pairs in each group
                            for np in xrange(cl[pos-1]+1, cl[pos]+1):
                                sumOfEachGroup += newarr[np][2]
                                sumOfEachGroupSquare += newarr[np][2]**2
                            # mean and sd of the insert length group
                            miuOfEachGroup = sumOfEachGroup*1.0/numOfEachGroup
                            sigmaOfEachGroup = math.sqrt(sumOfEachGroupSquare*1.0/numOfEachGroup-miuOfEachGroup**2)
                            
                            if numOfEachGroup > pairsForDeleting:
                                if (soap.index[c1], c1) in [(0, c2-1), (1, c2), (2, c2+1)] \
                                          and soap.length[c1]-miuOfEachGroup-2 < miuOfLib[numOfLib]+rateForDivision*sigmaOfLib[numOfLib] \
                                          and soap.length[c1]-miuOfEachGroup+2 > miuOfLib[numOfLib]-rateForDivision*sigmaOfLib[numOfLib]:
                                        bestclNum = sumOfEachGroup
                                        #removeall = 1
                                        #templist = 0
                                else:
                                    for np in xrange(cl[pos-1]+1, cl[pos]+1): # same as the step above
                                        if np == cl[pos-1] + 1:
                                            max1 = newarr[np][0]
                                            max2 = newarr[np][1]
                                        elif newarr[np][2] <= miuOfEachGroup+rateForDivision*sigmaOfLib[numOfLib]\
                                                and newarr[np][2] <= miuOfEachGroup-rateForDivision*sigmaOfLib[numOfLib]:
                                            if max1 < newarr[np][0]:
                                                max1 = newarr[np][0]
                                            if max2 < newarr[np][1]:
                                                max2 =  newarr[np][1]
                                    if soap.length[c1]-max1 > miuOfLib[numOfLib]+rateForRemoving*sigmaOfLib[numOfLib] \
                                           or soap.length[c2]-max2 > miuOfLib[numOfLib]+rateForRemoving*sigmaOfLib[numOfLib]:
                                        pass
                                    else:
                                        if bestclNum < sumOfEachGroup:
                                            bestclNum = sumOfEachGroup
                                            #templist = [c2,
                                            #            soap.length[c1]-max1,
                                            #            soap.length[c2]-max2,
                                            #            numOfEachGroup,
                                            #            sumOfEachGroup,
                                            #            miuOfLib[numOfLib]-soap.length[c1]-soap.length[c2]+miuOfEachGroup,
                                            #            sigmaOfEachGroup,
                                            #            numOfLib]
                                            if c1 not in prelink:
                                                prelink[c1] = []
                                            prelink[c1].append([c2,
                                                                soap.length[c1]-max1,
                                                                soap.length[c2]-max2,
                                                                numOfEachGroup,
                                                                sumOfEachGroup,
                                                                miuOfLib[numOfLib]-soap.length[c1]-soap.length[c2]+miuOfEachGroup,
                                                                sigmaOfEachGroup,
                                                                numOfLib
                                                                ])
            # init for next library
            numOfLib += 1
            totReadPairs, sumOfIns, sumOfInsSquare = 0, 0, 0
            insLength = collections.defaultdict(int)
            rec, temp = {}, [0]
    fin.close()
    fout.close()
    
    logger.info('PElinks are created !')
    return prelink, miuOfLib, sigmaOfLib
