############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import re
import math
import itertools
import collections
import config
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

def checkCigar(cigarValue):
    """
    Check if the reads mapped correctly
    """
    charNum = re.split('\D+', cigarValue)[:-1]
    charType = re.split('\d+', cigarValue)[1:]
    
    freqOfS = cigarValue.count('S')
    mismatchNum = sum([int(i[1]) for i in itertools.izip(charType, charNum) if i[0] != 'M'])
    
    if mismatchNum == 1 and (charType[0] == 'S' or charType[-1] == 'S'):
        return 1
    elif mismatchNum < 5:
        return 1
    else:
        return 0

def calMappingLength(cigarValue):
    """
    calculate the mapping length of the read
    """
    CHAR = frozenset(['M', 'S', 'D', 'I', 'H', 'P'])
    sumM, num = 0, '' 
    for ele in cigarValue:
        if ele not in CHAR:
            num += ele
        else:
            if ele == 'M':
                sumM += int(num)
            num = ''
    return sumM

def pelinksFromSam(samfn, libId, length, kmersize, sizeOfLib, readLen):
    """
    parser sam file and extrat pair mapped reads, estimate the insert length of the lib
    then use the estimated insert length to evaluate the gap size for finding lost dlinks
    also record the support number, mapping quality
    PS: also calculate the contig coverage
        Caculate the coverage of the contigs by the following formula:
        (total matches) / (contig length)
    """
    CHAR = frozenset(['M', 'S', 'D', 'I', 'H', 'P', 'N', 'X', '='])
    matches = collections.defaultdict(int)
    
    mappingQualityCut = config.QULIATY
    rateForDivision = config.RATE4DIV
    rateForRemoving = config.RATE4REM
    rateForSplitting = config.RATE4SPLIT
    lowerBoundOfSupport = config.SUPPORTNUM
    rateForCutting = config.RATE4CUT
    rateOfLoose = config.RATEOFLOOSE
    rateOfStrict = config.RATEOFSTRICT
    lenOfLargeContig = 2*kmersize + 1
    
    # parse the sam file and record the pair end mapped reads
    insLengthHist, pemapping = collections.defaultdict(int), {}
    readId = samFlag = contigId = mappingPos = mappingQuality = mappingCigar = insLength = 0
    cigarSucc = cigarFail = 0
    with open(samfn) as af:
        for line in af:
            # [0] readId [1] flag [2] contigId [3] pos [4] mappingQuality [5] ciagr [6] mrnm [7] mpos [8] insLen
            ae = line.rstrip().split()
            try:
                ae[1], ae[2], ae[3], ae[4], ae[8] = int(ae[1]), int(ae[2]), int(ae[3]), int(ae[4]), int(ae[8])
            except:
                continue
            ##################  contig coverage step begin ###############
            if int(ae[1]) & 0x0004 == 0 and int(ae[4]) > 5:
                sumM, num = 0, ''
                for ele in ae[5]:
                    if ele not in CHAR:
                        num += ele
                    else:
                        if ele == 'M':
                            sumM += int(num)
                        num = ''
                matches[int(ae[2])] += sumM
            ##################  contig coverage step end #############3
            
            # conditions take pair reads as mapped pair: flag and mapping quality
            if ae[0] == readId and not int(ae[1]%2**3/2**2) and not int(ae[1]%2**4/2**3) \
                    and ae[4] > mappingQualityCut and mappingQuality > mappingQualityCut:
                 # paired mapping reads
                if checkCigar(ae[5]) and checkCigar(mappingCigar):
                    cigarSucc += 1
                    if ae[2] == contigId:
                        # one same contig, for evaluating insert length
                        if int(ae[1]%2**5/2**4) + int(samFlag%2**5/2**4) == 1:
                            # check the 5th binary position, two reads on diff strands
                            if int(ae[1]%2**5/2**4):
                                if insLength > 0:
                                    insLengthHist[insLength] += 1
                                elif insLength == 0:
                                    newInsLength = ae[3] - mappingPos + calMappingLength(ae[5])
                                    if newInsLength > 0:
                                        insLengthHist[newInsLength] += 1
                            else:
                                # on the forward strand
                                if ae[8] > 0:
                                    insLengthHist[ae[8]] += 1
                                elif insLength == 0:
                                    newInsLength = mappingPos - ae[3] + calMappingLength(mappingCigar)
                                    if newInsLength > 0:
                                        insLengthHist[newInsLength] += 1
                    else:
                        # on different contigs, for estmating gap size
                        if ae[2] > contigId:
                            if int(ae[1]%2**5/2**4):
                                # mapping neg, Dir tag: -1, tail length: contig length - mapping pos - read length
                                firDir, firPos = -1, length[ae[2]] - ae[3] - readLen
                            else:
                                firDir, firPos = 1, ae[3]
                            
                            if int(samFlag%2**5/2**4):
                                secDir, secPos = -1, length[contigId] - mappingPos - readLen
                            else:
                                secDir, secPos = 1, mappingPos
                        else:
                            if int(ae[1]%2**5/2**4):
                                secDir, secPos = -1, length[ae[2]] - ae[3] - readLen
                            else:
                                secDir, secPos = 1, ae[3]
                            
                            if int(samFlag%2**5/2**4):
                                firDir, firPos = -1, length[contigId] - mappingPos - readLen
                            else:
                                firDir, firPos = 1, mappingPos
                        
                        tmp = (firPos, secPos, firPos+secPos, (mappingQuality+ae[4])/2)
                        
                        #initialize the mapping structure
                        if ae[2] not in pemapping:
                            pemapping[ae[2]] = {}
                        if firDir not in pemapping[ae[2]] :
                            pemapping[ae[2]][firDir] = {}
                        if contigId not in pemapping[ae[2]][firDir]:
                            pemapping[ae[2]][firDir][contigId] = {}
                        if secDir not in pemapping[ae[2]][firDir][contigId]:
                             pemapping[ae[2]][firDir][contigId][secDir] = []
                        pemapping[ae[2]][firDir][contigId][secDir].append(tmp)
                else:
                    cigarFail += 1
            # record read1 in the pair
            readId, samFlag, contigId, mappingPos = ae[0], ae[1], ae[2], ae[3]
            mappingQuality, mappingCigar, insLength = ae[4], ae[5], ae[8]
    
    logger.info('Happy pairs {0:,} out of {1:,} paired reads'.format(cigarSucc, cigarSucc+cigarFail))
    
    # pair end lib estimation by using EM algorithm
    logger.info('Insert size estimating ...')
    totPairs = sum(insLengthHist.values())
    
    # minL: the minimum likelihood score
    minL = 0
    for ele in [1, 2, 3]:
        clWeight, clMean, clSD, Ll, numiter = GMM(insLengthHist, ele, 0.001, 100)
        min2 = -Ll + ele*math.log(totPairs)
        if not minL or minL > min2:
            minL = min2
            bestclMean, bestclSD = clMean, clSD # both of them are list
    
    # the best gmm cluster is the cluster with the minmum distance with idear lib size
    minDeltaOfM = -1
    for ele in xrange(len(bestclMean)):
        if minDeltaOfM == -1 or minDeltaOfM > abs(bestclMean[ele] - sizeOfLib):
            minDeltaOfM = abs(bestclMean[ele] - sizeOfLib)
            miuOfLib = bestclMean[ele] 
            sigmaOfLib = bestclSD[ele]
    
    # count ratio of above insert length read pairs and that of below pairs, they are noise
    aboveInsPairs = belowInsPairs = 0
    for key, value in insLengthHist.iteritems():
        if key > miuOfLib + 3*sigmaOfLib:
            aboveInsPairs += value
        elif key < miuOfLib - 3*sigmaOfLib:
            belowInsPairs += value
    rateOfHigherNoise = aboveInsPairs*1.0/totPairs
    rateOfLowerNoise = belowInsPairs*1.0/totPairs
    
    # pair end mapping to pelink
    logger.info('Transforming pair mapping to pelinks ...')
    # record max support
    supportRec = {}
    for c1 in pemapping:
        for d1 in pemapping[c1]:
            for c2 in pemapping[c1][d1]:
                for d2 in pemapping[c1][d1][c2]:
                    if length[c1] > lenOfLargeContig and length[c2] > lenOfLargeContig:
                        if c1 not in supportRec:
                            supportRec[c1] = {}
                        if d1 not in supportRec[c1]:
                            supportRec[c1][d1] = 0
                        if c2 not in supportRec:
                            supportRec[c2] = {}
                        if d2 not in supportRec[c2]:
                            supportRec[c2][d2] = 0
                        if supportRec[c1][d1] < len(pemapping[c1][d1][c2][d2]):
                            supportRec[c1][d1] = len(pemapping[c1][d1][c2][d2])
                        if supportRec[c2][d2] < len(pemapping[c1][d1][c2][d2]):
                            supportRec[c2][d2] = len(pemapping[c1][d1][c2][d2])
    
    supportMax = {}
    sumOfLargeContig = 0
    for c in supportRec:
        for d in supportRec[c]:
            if length[c] > lenOfLargeContig:
                if supportRec[c][d] not in supportMax:
                    supportMax[supportRec[c][d]] = {}
                if c not in supportMax[supportRec[c][d]]:
                    supportMax[supportRec[c][d]][c] = {}
                supportMax[supportRec[c][d]][c][d] = 1
                if length[c] > lenOfLargeContig:
                    sumOfLargeContig += length[c]
    
    # confirm support pairs cutoff
    accumOfLargeContig = looseCutoffOfSupport = strictCutoffOfSupport = 0
    tmp = sorted(supportMax.keys())
    for sup in tmp:
        # get accumOfLargeContig
        if not (looseCutoffOfSupport and strictCutoffOfSupport):
            for c in supportMax[sup]:
                if not (looseCutoffOfSupport and strictCutoffOfSupport):
                    for d in supportMax[sup][c]:
                        accumOfLargeContig += length[c]
        
        if accumOfLargeContig*1.0/sumOfLargeContig > rateOfLoose and not looseCutoffOfSupport:
            looseCutoffOfSupport = sup + 1
        if strictCutoffOfSupport*1.0/sumOfLargeContig > rateOfStrict and not strictCutoffOfSupport:
            strictCutoffOfSupport = sup + 1
    
    if looseCutoffOfSupport < lowerBoundOfSupport:
        looseCutoffOfSupport = lowerBoundOfSupport
    
    if strictCutoffOfSupport < lowerBoundOfSupport:
        strictCutoffOfSupport = lowerBoundOfSupport
    #supportMax = {}
    #supportRec = {}
    
    # mapping to pelinks
    prePelink = {}
    for c1 in pemapping:
        for d1 in pemapping[c1]:
            for c2 in pemapping[c1][d1]:
                for d2 in pemapping[c1][d1][c2]:
                    arr = pemapping[c1][d1][c2][d2]
                    arr.sort(key=lambda x: x[2])
                    # remove the outliers
                    cutLowerPairs = int(len(arr)*rateOfHigherNoise*rateForCutting)
                    cutHigherPairs = int(len(arr)*rateOfLowerNoise*rateForCutting)
                    if cutLowerPairs >= 1:
                        arr = arr[cutLowerPairs:]
                    if cutHigherPairs >= 1:
                        arr = arr[:len(arr)-cutHigherPairs]
                    # get the miu of sigma of the arr
                    sumOfarr = sumOfarrSquare = sumOfQual = 0
                    qualList = []
                    for ele in arr:
                        sumOfarr += ele[2]
                        sumOfarrSquare += ele[2]**2
                        #sumOfQual += ele[3]
                        qualList.append(ele[3])
                    miuOfarr = sumOfarr*1.0/len(arr)
                    sigmaOfarr = math.sqrt(sumOfarrSquare*1.0/len(arr) - miuOfarr**2)
                    #miuOfQual = sumOfQual*1.0/len(arr)
                    
                    if sigmaOfarr <= rateForSplitting*sigmaOfLib and len(arr) >= min([looseCutoffOfSupport, lowerBoundOfSupport]):
                        # simga fits the lib
                        max1, max2 = arr[0][0:2]
                        for ele in xrange(1, len(arr)):
                            if max1 < arr[ele][0]:
                                max1 = arr[ele][0]
                            if max2 < arr[ele][1]:
                                max2 = arr[ele][1]
                        
                        if length[c1]-max1 > miuOfLib+rateForRemoving*sigmaOfLib or length[c2]-max2 > miuOfLib+rateForRemoving*sigmaOfLib:
                            # error in this pelink: mapping tail is larger than $miuOfLib+$rateForRemoving*$sigmaOfLib
                            pass
                        else:
                            tmp = [c2, d2, length[c1]-max1, length[c2]-max2, len(arr),
                                   len(arr)*miuOfarr, miuOfLib-length[c1]-length[c2]+miuOfarr,
                                   sigmaOfarr, qualList, sizeOfLib] 
                            
                            if c1 not in prePelink:
                                prePelink[c1] = {}
                            if d1 not in prePelink[c1]:
                                prePelink[c1][d1] = []
                            prePelink[c1][d1].append(tmp)
                    elif sigmaOfarr > rateForSplitting*sigmaOfLib and len(arr) >= strictCutoffOfSupport:
                        # sigma does not fit the lib, and the support is large enough
                        # find the trip point and put it into @cl(cluster)
                        cl = [-1]
                        for pos in xrange(len(arr)-1):
                            if arr[pos+1][2]-arr[pos][2] >= rateForDivision*sigmaOfLib:
                                cl.append(pos)
                        cl.append(len(arr)-1)
                        
                        # check every group
                        for pos in xrange(1, len(cl)):
                            qualList = []
                            bestclNum, templist, sumOfEachGroup, sumOfEachGroupSquare, sumOfEachGroupQual = 0, 0, 0, 0, 0
                            numOfEachGroup = cl[pos] - cl[pos-1]
                            for np in xrange(cl[pos-1]+1, cl[pos]+1):
                                sumOfEachGroup += arr[np][2]
                                sumOfEachGroupSquare += arr[np][2]**2
                                #sumOfEachGroupQual += arr[np][3]
                                qualList.append(arr[np][3])
                            miuOfEachGroup = sumOfEachGroup*1.0/numOfEachGroup
                            sigmaOfEachGroup = math.sqrt(sumOfEachGroupSquare*1.0/numOfEachGroup-miuOfEachGroup**2)
                            #miuOfEachGroupQual = sumOfEachGroupQual*1.0/numOfEachGroup
                            
                            if numOfEachGroup > min([looseCutoffOfSupport, lowerBoundOfSupport]):
                                for np in xrange(cl[pos-1]+1, cl[pos]+1): # same as the step above
                                    if np == cl[pos-1] + 1:
                                        max1 = arr[np][0]
                                        max2 = arr[np][1]
                                    elif arr[np][2] <= miuOfEachGroup+rateForDivision*sigmaOfLib\
                                            and arr[np][2] <= miuOfEachGroup-rateForDivision*sigmaOfLib:
                                        if max1 < arr[np][0]:
                                            max1 = arr[np][0]
                                        if max2 < arr[np][1]:
                                            max2 =  arr[np][1]
                                
                                if length[c1]-max1 > miuOfLib+rateForRemoving*sigmaOfLib \
                                    or length[c2]-max2 > miuOfLib+rateForRemoving*sigmaOfLib:
                                    pass
                                else:
                                    if bestclNum < sumOfEachGroup:
                                         bestclNum = sumOfEachGroup
                                         templist = [c2, d2, length[c1]-max1, length[c2]-max2, numOfEachGroup,
                                                     sumOfEachGroup, miuOfLib-length[c1]-length[c2]+miuOfEachGroup,
                                                     sigmaOfEachGroup, qualList, sizeOfLib]
                            if templist:
                                if c1 not in prePelink:
                                    prePelink[c1] = {}
                                if d1 not in prePelink[c1]:
                                    prePelink[c1][d1] = []
                                prePelink[c1][d1].append(templist)
                    elif sigmaOfarr > rateForSplitting*sigmaOfLib and len(arr) >= min([looseCutoffOfSupport, lowerBoundOfSupport]):
                        max1, max2 = arr[0][0], arr[0][1]
                        for np in xrange(len(arr)):
                            if max1 < arr[np][0]:
                                max1 = arr[np][0]
                            if max2 < arr[np][1]:
                                max2 = arr[np][1]
                            
                            if length[c1]-max1 > miuOfLib+rateForRemoving*sigmaOfLib or length[c2]-max2 > miuOfLib+rateForRemoving*sigmaOfLib:
                                # error in this pelink: mapping tail is larger than $miuOfLib+$rateForRemoving*$sigmaOfLib
                                pass
                            else:
                                tmp = [c2, d2, length[c1]-max1, length[c2]-max2, len(arr),
                                   len(arr)*miuOfarr, miuOfLib-length[c1]-length[c2]+miuOfarr,
                                   sigmaOfarr, qualList, sizeOfLib]
                                if c1 not in prePelink:
                                    prePelink[c1] = {}
                                if d1 not in prePelink[c1]:
                                    prePelink[c1][d1] = []
                                prePelink[c1][d1].append(tmp)
    
    # pelink filter
    logger.info('Filtering pelinks ...')
    pelink = {}
    tmp = sorted(prePelink.keys())
    for c1 in tmp:
        for d1 in prePelink[c1]:
            arr = prePelink[c1][d1]
            # prelink: {c1: d1: [ [c2,d2,left1,left2,num,sum, estimated_gap,pe_sigma,qual,library size] ]
            arr.sort(key=lambda x: x[0])
            # check the adjacent same contigs to detect the noise
            if len(arr) > 1:
                record = []
                for pos in xrange(1, len(arr)):
                    if arr[pos][0] == arr[pos-1][0]:
                        if arr[pos][4] <= looseCutoffOfSupport:
                            #arr.pop(pos)
                            record.append(pos)
                        elif arr[pos-1][4] <= looseCutoffOfSupport:
                            #arr.pop(pos - 1)
                            record.append(pos-1)
                        else:
                            if abs(arr[pos][6]-arr[pos-1][6]) < length[arr[pos][0]]-kmersize:
                                if sigmaOfLib*1.0/math.sqrt(arr[pos][4]) > sigmaOfLib*1.0/math.sqrt(arr[pos-1][4]):
                                    #arr.pop(pos)
                                    record.append(pos)
                                else:
                                    #arr.pop(pos-1)
                                    record.append(pos-1)
                record.sort(reverse=True)
                for ele in record:
                    arr.pop(ele)
            
            # check current library
            for ele in arr:
                if ele[6] < miuOfLib + rateForRemoving*sigmaOfLib:
                    lengthScore1 = int(length[c1]*1.0/lenOfLargeContig*100)*1.0/100
                    lengthScore2 = int(length[ele[0]]*1.0/lenOfLargeContig*100)*1.0/100
                if lengthScore1 > 1: lengthScore1 = 1
                if lengthScore2 > 1: lengthScore2 = 1
                if ele[7] < 1.5*sigmaOfLib:
                    # [0] num [1] est_gap [2] pe_sigma [3] avequal [4] est_sigma [5] lib_id
                    # [6] pe_gap [7] length_score [8] support_score
                    tmp = [ele[4], ele[6], ele[7], ele[8], sigmaOfLib*1.0/math.sqrt(ele[4]),
                           libId, ele[6], lengthScore1*lengthScore2]
                    if c1 not in pelink:
                        pelink[c1] = {}
                    if d1 not in pelink[c1]:
                        pelink[c1][d1] = {}
                    if ele[0] not in pelink[c1][d1]:
                        pelink[c1][d1][ele[0]] = {}
                    if ele[1] not in pelink[c1][d1][ele[0]]:
                        pelink[c1][d1][ele[0]][ele[1]] = []
                    
                    if ele[0] not in pelink:
                        pelink[ele[0]] = {}
                    if ele[1] not in pelink[ele[0]]:
                        pelink[ele[0]][ele[1]] = {}
                    if c1 not in pelink[ele[0]][ele[1]]:
                        pelink[ele[0]][ele[1]][c1] = {}
                    if d1 not in pelink[ele[0]][ele[1]][c1]:
                        pelink[ele[0]][ele[1]][c1][d1] = []
                    
                    pelink[c1][d1][ele[0]][ele[1]].append(tmp)
                    pelink[ele[0]][ele[1]][c1][d1].append(tmp)
    
    logger.info('Libaray: {0}'.format(libId))
    logger.info('  Miu of Lib: {0}'.format(int(miuOfLib*100)*1.0/100))
    logger.info('  Sigma of Lib: {0}'.format(int(sigmaOfLib*100)*1.0/100))
    logger.info('  Ratio of lower noise: {0}'.format(int(rateOfLowerNoise*10000)*1.0/100))
    logger.info('  Ratio of higher noise: {0}'.format(int(rateOfHigherNoise*10000)*1.0/100))
    logger.info('  Base cut off of support: {0}'.format(lowerBoundOfSupport))
    logger.info('  Loose of support: {0}'.format(looseCutoffOfSupport))
    logger.info('  Strict cut off of support: {0}'.format(strictCutoffOfSupport))
    
    # calculate coverage
    cov = {i:round(matches[i]*1.0/length[i], 2) for i in matches}
    for i in set(length.keys()) - set(cov.keys()):
        cov[i] = 0.0
    
    
    return pelink, miuOfLib, sigmaOfLib, cov