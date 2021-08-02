############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import math
import decompose
import config
import itertools
import copy
import collections
import subprocess
import cPickle as pickle
import numpy as np
from log import get_logger
# set log file
logger = get_logger()

try:
    from sklearn.mixture import GMM
except:
    logger.warning('Package sklearn is not imported. GMM step will be passed.')

def installSVM():
    """
    Check if libsvm was installed, etherwise, install it
    """
    logger.info('Check if libSVM is installed ...')
    currentDir = os.getcwd()
    svmPath = os.path.join(config.libsLocation, 'libsvm-master')
    binaries = ['svm-train', 'svm-predict', 'svm-scale']
    for requiredf in binaries:
        if not os.path.isfile(os.path.join(svmPath, requiredf)):
            logger.info('libSVM not installed, compiling it...')
            break
        else:
            logger.info('OK, the software has been installed, now continue analysising ...')
            return 1
    
    os.chdir(svmPath)
    returnCode1 = subprocess.call(['make'])
    os.chdir(os.path.join(svmPath, 'python'))
    returnCode2 = subprocess.call(['make'])
    if returnCode1 != 0 or returnCode2 != 0:
        os.chdir(currentDir)
        logger.error('Failed to compile libSVM ' +  svmPath +' Try to compile it manually.',
            to_stderr=True,
            exit_with_code=3
            )
    else:
        logger.info('OK, the software is installed, now continue analysising ...')
        os.chdir(currentDir)

def average(x):
    
    assert len(x) > 0
    return float(sum(x))/len(x)

def pearson_cor(x, y):
    
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avgX = average(x)
    avgY = average(y)
    diffpord, ydiff2, xdiff2 = 0, 0, 0
    for pos in xrange(n):
        xdiff = x[pos] - avgX
        ydiff = y[pos] - avgY
        diffpord += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff
    
    return diffpord/math.sqrt(xdiff2 * ydiff2)

def cluster(args):
    # read parameters
    c, covariance, inits, iters, data = args
    
    logger.debug('Doing GMM with a cluster number of {}'.format(c))
    # use EM algorithm to fit the module
    gmm = GMM(n_components=c,
              covariance_type=covariance,
              n_init=inits,
              n_iter=iters).fit(data)
    
    # get the bic criterion score
    bic = gmm.bic(data)
    
    if gmm.converged_:
        logger.debug("Clustering into {0} clusters converged.".format(c))
    else:
        logger.warning(("Clustering into {0} clusters did not "
                         "converge, consider increasing the number "
                         "of iterations.").format(c))
        print >> sys.stderr, "Cluster {0} did not converge".format(c)
    
    return bic,c, gmm.converged_

# single threads
class peaksSingle():
    
    def __init__(self, soap, lenCut=0):
        
        # read and set variables
        self.soap = soap
        
        # default gmm paras:
        self.lenCut = lenCut
        self.componentsNu = config.COMPONENTNUM
        self.iters = config.ITERS
        self.inits = config.INITS
        self.covariance = config.COVARIANCE
    
    def fit(self, points):
        
        # note the number of points, if it small than the number of components
        numComp = self.componentsNu if len(points)>self.componentsNu[-1] else range(len(points))
        
        # set paramerters of gmm module
        clusterArgs = []
        for com in numComp:
            clusterArgs.append([com,
                                self.covariance,
                                self.inits,
                                self.iters,
                                points])
        
        # parallels to run gmm
        results = []
        for arg in clusterArgs:
            results.append(cluster(arg))
        
        # determine the optimal cluster number by using BIC
        min_bic, optimalC, converged = min(results,key=lambda x: x[0])
        if not converged:
            logging.error(("Optimal bic score was reached for non convergent "
                           "cluster number {0}, exiting without clustering "
                           "output").format(optimalC),
                          exit_with_code=1
                          )
        
        logger.info('{0} clusters are generated'.format(optimalC))
        gmm = GMM(n_components=optimalC,
                  covariance_type=self.covariance,
                  n_init=self.inits,
                  n_iter=self.iters,
                  ).fit(points)
        
        # use the posterior probabilities of each data point to cluster
        clusters = gmm.predict(points)
        
        return clusters
    
    def label(self, data):
        
        # coverage gorups
        logger.info('Fitting the coverages ...')
        
        # the data should be a list of list like format
        for lis in data:
            if isinstance(lis, list) or isinstance(lis, tuple):
                break
            else:
                logger.error('The format of the data should like list of list',
                            exit_with_code=1
                             )
        
        ids, points = [], []
        # transform the format of data
        discard = 0
        for na, point in data:
            if self.soap.length[na] < self.lenCut:
                discard += 1
                continue
            ids.append(na)
            points.append([point])
        
        logger.warning('{0} of the contigs are filtered, because of short length'.format(discard))
        
        # label the points
        points = np.array(points)
        group = self.fit(points)
        cls = {}
        for ele, value in itertools.izip(ids, group):
            cls[ele] = value
        
        return cls

def EMseedFilter(soap, paras, recovered):
    """
    The kmer coverage of the target genome follows a normal distribution
    use gmm model to fit the distribution and use the mean and std to
    filter the false positive recovered nodes
    """
    if not recovered:
        return set()
    
    points = np.array([[soap.cov.get(i, 0)] for i in paras.nodes if soap.length[i] >= 300])
    clusterArgs = []
    for com in [1, 2, 3, 4, 5]:
        clusterArgs.append([com,
                            config.COVARIANCE,
                            config.INITS,
                            config.ITERS,
                            points])
    
    results = []
    for arg in clusterArgs:
        results.append(cluster(arg))
    
    # determine the optimal cluster number by using BIC
    min_bic, optimalC, converged = min(results,key=lambda x: x[0])
    if not converged:
        logging.error(("Optimal bic score was reached for non convergent "
                       "cluster number {0}, exiting without clustering "
                       "output").format(optimalC),
                      exit_with_code=1
                      )
    
    logger.info('{0} clusters are generated'.format(optimalC))
    gmm = GMM(n_components=optimalC,
              covariance_type = config.COVARIANCE,
              n_init = config.INITS,
              n_iter = config.ITERS,
              ).fit(points)
    
    # use the posterior probabilities of each data point to cluster
    if optimalC == 1:
        mean = gmm.means_[0]
        std = np.sqrt(gmm.covars_).ravel()[0]
        refMin, refMax = mean - 2*std, mean + 2*std
        accept = [i for i in recovered if refMin <= soap.cov.get(i, 0) <= refMax]
    else:
        means = gmm.means_.ravel()
        weights = gmm.weights_.ravel()
        stds = np.sqrt(gmm.covars_).ravel()
        res = zip(weights, means, stds)
        res.sort(key=lambda x: x[0], reverse=True)
        weight = 0
        for pos, value in enumerate(res):
            weight += value[0]
            if weight >= 0.7:
                break
        refRange = []
        for i in xrange(pos+1):
            refRange.append([ res[i][1]-2*res[i][2], res[i][1]+2*res[i][2] ])
        
        accept = []
        for ele in recovered:
            if ele not in soap.cov:
                continue
            for r in refRange:
                if r[0] <= soap.cov[ele] <= r[1]:
                    accept.append(ele)
                    break
    
    return set(accept)

def EMratioFilter(soap, paras, recovered, intra, inter, sizeCut):
    """
    use coverag to filter false positive
    """
    numcf = 500
    if len(recovered) < numcf:
        return set()
    
    gmm = peaksSingle(soap, lenCut=sizeCut)
    data = [(i, soap.cov[i]) for i in (recovered|paras.nodes) if i in soap.cov]
    cls = gmm.label(data)
    
    if set(cls.values())==1:
        # only one cluster, pass
        passed = set(cls.keys()) 
    else:
        transedCLS = {}
        for i in cls:
            if cls[i] not in transedCLS:
                transedCLS[cls[i]] = set()
            transedCLS[cls[i]].add(i)
        
        # narrow the select region of mapped nodes
        refcls = [cls[i] for i in cls if i in paras.nodes]
        refdis = collections.Counter(refcls)
        reliableCls, totalCount = set(), sum(refdis.values())
        for key, refnum in refdis.iteritems():
            if refnum*1.0/totalCount >= intra\
                or refnum*1.0/len(transedCLS[key])>= inter:
                reliableCls.add(key)
        
        # in case of no reliableCls
        while not reliableCls:
            inter -= 0.05
            if inter == 0:
                break
            for ele in refdis:
                if refdis[ele]*1.0/refnum >= inter:
                    reliableCls.add(ele)
        
        # choose the reliable cls
        passed = set()
        for i in transedCLS:
            if i in reliableCls:
                passed.update([ele for ele in transedCLS[i] if ele not in paras.nodes])
    
    return passed

def getGNodes(soap, paras, threads):
    """
    combined the decompose info and coverage cluster info to remove the
    fasle positive nodes generated by svm
    """
    installSVM()
    import increase
    # get candidates
    preL, candidateL = increase.classification(soap, paras, int(threads))
    
    # EM algorithm (kmer coverage)
    candL = EMratioFilter(soap, paras, candidateL, config.INTRACFL, config.INTERCFL, 300)
    
    cleanL = EMseedFilter(soap, paras, candL)
    
    if config.DEBUG:
        with open('cleanL', 'wb') as af:
            pickle.dump(cleanL, af)
    
    return cleanL
