############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import collections
import subprocess
import copy
import glob
import config
import itertools
import random
import shutil
import genemark
import cPickle as pickle
import scaffolding
import numpy as np
import multiprocessing
from log import get_logger

# set log file
logger = get_logger()
libs_dirpath = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.join(libs_dirpath, 'libsvm-master', 'tools'))
sys.path.append(os.path.join(libs_dirpath, 'libsvm-master', 'python'))
from grid import *
from svmutil import *
from scipy.stats import ks_2samp
# meachine learning algorithm: random forest
randomForest = 0
try:
    import numpy as np
    from sklearn.externals import joblib
    #from sklearn.ensemble import RandomForestClassifier
    from sklearn.grid_search import GridSearchCV
    from sklearn.cross_validation import KFold
    from sklearn import feature_selection
    from sklearn.pipeline import Pipeline
except:
    randomForest = 0
    logger.warning('Warning: python package scikit-learn is not found, then no random forest step!')

_penalty, _fold = config.PENALTY, config.FOLD

fivebase = ['A','T','C','G','N']
codon = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA',
        'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC',
        'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG',
        'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT',
        'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']

def cleanUp(*files):
    """
    remoe unused files
    """
    for fn in files:
        if os.path.exists(fn):
            os.remove(fn)

def permutation(size = 4):
    """
    Generate the combination of a,t,c,g
    """
    combination = []
    S = ['A','T','G','C']
    for seq in itertools.product(S, repeat = size):
        combination.append(''.join([i for i in seq]))
    
    merged = []
    for ele in combination:
        revcom = revCom(ele)
        if revcom not in merged:
            merged.append(ele)
    
    merged.sort()
    return merged

def calFreq(seq, kmerSpectrum, size = 4, out = 'both'):
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
        table[revCom(i)] += table[i]
        del table[i]
    
    get = table.get
    # use the total count to normalize
    allCount = sum(table.values())
    if out is 'freq':
        #get the number of each type of kmer in the kmerSpectrum to ready for SVM
        kmerNumber = ( get(i, 0) for i in kmerSpectrum )
        return kmerNumber
    elif out is 'order':
        # This step is to make sure the order of kmers that number is zero is the same among differenct contigs
        orders = sorted(table, key = get) # from low to high (value), returns the oreded keys
        if len(orders) < number:
            newKmers = [i for i in kmerSpectrum if i not in orders]
            newKmers.extend(orders)
            #return dict(zip(newKmers, numberindex))
            return {i[0]:i[1] for i in itertools.izip(newKmers, numberindex)}
        else:
            #return dict(zip(orders, numberindex))
            return {i[0]:i[1] for i in itertools.izip(orders, numberindex)}
    elif out is 'both':
        kmerNumber = ( get(i, 0) for i in kmerSpectrum )
        orders = sorted(table, key = get)
        if len(orders) < number:
            newKmers = [i for i in kmerSpectrum if i not in orders]
            newKmers.extend(orders)
            #return dict(zip(newKmers, numberindex)), kmerNumber
            return {i[0]:i[1] for i in itertools.izip(newKmers, numberindex)}, kmerNumber
        else:
            #return dict(zip(orders, numberindex)), kmerNumber
             return {i[0]:i[1] for i in itertools.izip(orders, numberindex)}, kmerNumber

def revCom(string):
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

def normalization(freq,bmin=0,bmax=1):
    """
    normalized to [0, 1]
    (bmax-bmin)*(x-xmin)/(xmax-xmin) + lower
    """
    tmp = list(freq)
    allCount = sum(tmp)
    # first normailize by total number
    tmp = [i*1.0/allCount for i in tmp]
    # then normalized to [0, 1]
    maxv, minv = max(tmp), min(tmp)
    constant = maxv - minv
    normaled = [(ele-minv)*1.0/constant for ele in tmp]
    return normaled

def spearMan(order1, order2):
    
    distance = 0
    for kmer in order1:
        if order1[kmer] > order2[kmer]:
            distance += order1[kmer] - order2[kmer]
        else:
            distance += order2[kmer] - order1[kmer]
    return distance

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
    logger.info('Predicting proteins of the contigs in the file: ' + infile)
    logger.print_command_line(cmdMetaGeneMark, wrap_after=None)
    
    child = subprocess.call(cmdMetaGeneMark)
    if child != 0:
        sys.stderr.write(out)
        logger.error('Error happened while using MetaGeneMark', exit_with_code=1)
    else:
        with open('refGenes.lst','r') as lfile:
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
    
    with open('refPro', 'wb') as af:
        pickle.dump(protein, af)
    
    # return the hash, key: contig id; value: the codon usage of the predicted protein
    return protein 

def intraGenomeDis(soap, paras, threads):
    """
    Caculate the spear man footrule distance of intra-genome
    """
    logger.info('Caculating intra-genome distance')
    #logger.info('Distance recovery range: {0}% and {1}%'.format(config.RANGEL*100, config.RANGES*100))
    
    kmerSpectrum = permutation()
    orderL, count, i = [], 0, 0
    with open(paras.refn, 'r') as af:
        for line in af:
            if line[0] != '>':
                seq = line[:-1]
                if len(line) >= 10000:
                    # split into 5k
                    count += 1
                    subseq = seq[i*5000:(i+1)*5000]
                    orderL.append(calFreq(subseq, kmerSpectrum, out='order'))
                    if count >= 1000:
                        break
                elif len(line) >= 1000:
                    count += 1
                    orderL.append(calFreq(seq, kmerSpectrum, out='order'))
                    if count >= 1000:
                        break
    
    # combination to get the distribution of the intra-genome distance
    result, count = [], 0
    for pair in itertools.combinations(orderL, 2):
        count += 1
        result.append(spearMan(pair[0], pair[1]))
        if count >= 3000:
            break
    return result

def child_initialize(refOrder, refDis):
    """
    in case of large dataset of graph is passed to the function (which will result only one cpu occupied)
    http://stackoverflow.com/questions/25825995/python-multiprocessing-only-one-process-is-running#
    the huge volume varible should be shared by child threads
    """
    global referenceOrders, referenceDistribution, gSoap, gParas
    referenceOrders = refOrder
    referenceDistribution = refDis

def poolClass(seqid, seqOrder):
    """
    Multiple process of caculate the spearman distance
    """
    try:
        dis = [spearMan(seqOrder, value) for key, value in referenceOrders.iteritems()]
        stat = ks_2samp(referenceDistribution, dis)
        if stat.pvalue > 0.1: # same distribution, accept it as candidate
            return seqid, -1 # candidate
        elif stat.pvalue < 0.0001:
            return seqid, np.mean(dis) # unused
    except:
        return seqid, 0

def getCandiNodes(soap, paras, threads):
    """
    use spearman distance to get the candidates
    """
    if config.MDADEBUG:
        if os.path.exists('spearman') and os.path.getsize('spearman'):
            logger.info('Find existing spearman file, it will be used !')
            with open('spearman', 'rb') as af:
                candL = pickle.load(af)
                truncateL = pickle.load(af)
            return candL, truncateL
    
    # get negative and candidate distance boundary 
    refDis = intraGenomeDis(soap, paras, threads)
    
    # reference kmer order which used for comparing
    kmerSpectrum, length, refOrder, refPosIds = permutation(), {}, {}, set()
    with open(paras.refn, 'r') as af:
        for line in af:
            if line[0] is '>':
                seqid = line[1:-1]
            else:
                seqlen = len(line[:-1])
                refPosIds.add(seqid) 
                if seqlen >= 5000:
                    refOrder[seqid] = calFreq(line[:-1][:5000], kmerSpectrum, out='order')
                    length[seqid] = seqlen
    
    # recover large contigs, using KS test to determin the candidates
    logger.info('Recovering contigs ...')
    pool = multiprocessing.Pool(processes = threads, initializer = child_initialize, initargs = (refOrder, refDis,))
    result = []
    with open(soap.tetraOrder, 'r') as af:
        for line in af:
            info = line.rstrip().split() # info[0] is contig id
            if str(info[0]) not in refPosIds:
                seqOrder = {x[0]:x[1] for x in itertools.izip(kmerSpectrum, map(int, info[1:]))}
                result.append(pool.apply_async(poolClass, (info[0], seqOrder,)))
    pool.close()
    pool.join()
    
    candL, negDis = set(), {}
    for key in result:
        ele = key.get()
        if ele:
            if ele[1] == -1:
                candL.add(ele[0])
            elif ele[1] > 0:
                negDis[ele[0]] = ele[1]
    
    # select the negatie data
    logger.info('Get {0:,} contigs'.format(len(candL)))
    orderedNegId = sorted(negDis, key=negDis.get)
    
    lenPos = len(refOrder) if paras.marker == 'm' else len(paras.refids) + len(set([str(i) for i in paras.nodes if soap.length[i] >= 300]))
    try:
        truncateL = set(random.sample(orderedNegId, _fold*lenPos))
    except:
        # increase of sample larger than population
        truncateL = set(orderedNegId)
    
    if config.MDADEBUG:
        with open('spearman', 'wb') as af:
            pickle.dump(candL, af, True)
            pickle.dump(truncateL, af, True)
    
    return candL, truncateL

def nuclData(soap, paras, cand, neg):
    """
    output the svm files
    """
    kmerSpectrum = permutation()
    logger.info('Recovering tetranuclotide frequency')
    if paras.marker == 'm':
        posids = set([str(i) for i, k in paras.length.iteritems() if k >= 300])
    else:
        posids = paras.refids | set([str(i) for i in paras.nodes if soap.length[i] >= 300])
    
    # positive data for random forest
    sampeSize = len(posids|neg)
    label = np.empty((sampeSize))
    train = np.empty((sampeSize, len(kmerSpectrum)))
    
    flag, counter1 = 0, 0
    trainfile = open(config.NUCTRAIN, 'w') # for svm
    with open(paras.refn,'r') as af:
        for line in af:
            if line[0] == '>':
                flag = 1 if line[1:-1] in posids else 0 
            else:
                if flag:
                    seq = line[:-1]
                    seqFreq = calFreq(seq, kmerSpectrum, out='freq')
                    vect = normalization(seqFreq)
                    train[counter1] = vect
                    flag, label[counter1] = 0, 1
                    counter1 += 1
                    # for svm pos data
                    trainfile.write('1 ')
                    for pos, ele in enumerate(vect):
                        if ele == 0: continue
                        trainfile.write('{0}:{1:.3f} '.format(pos+1, ele))
                    trainfile.write('\n')
    
    # candi and negtive data
    counter2 = 0
    canfile = open(config.NUCANDI, 'w') # for svm candidates
    test = np.empty((len(cand), len(kmerSpectrum)))
    testID = np.empty((len(cand)))
    with open(soap.TETRAFREQ, 'r') as af:
        for line in af:
            inf = line.split()
            if inf[0] in cand:
                seqFreq = map(float, inf[1:])
                vect = normalization(seqFreq)
                test[counter2], testID[counter2] = vect, inf[0]
                counter2 += 1
                # for svm candidates
                canfile.write(inf[0] + ' ')
                for pos, ele in enumerate(vect):
                    if ele == 0: continue
                    canfile.write('{0}:{1:.3f} '.format(pos+1, ele))
                canfile.write('\n')
            elif inf[0] in neg:
                seqFreq = map(float, inf[1:])
                vect = normalization(seqFreq)
                train[counter1], label[counter1] = vect, -1
                counter1 += 1
                trainfile.write('-1 ')
                for pos, ele in enumerate(vect):
                    if ele == 0: continue
                    trainfile.write('{0}:{1:.3f} '.format(pos+1, ele))
                trainfile.write('\n')
            else:
                if paras.marker == 's':
                    if inf[0] in posids:
                        seqFreq = map(float, inf[1:])
                        vect = normalization(seqFreq)
                        train[counter1], label[counter1] = vect, 1
                        counter1 += 1
                        trainfile.write('1 ')
                        for pos, ele in enumerate(vect):
                            if ele == 0: continue
                            trainfile.write('{0}:{1:.3f} '.format(pos+1, ele))
                        trainfile.write('\n')
    
    trainfile.close()
    canfile.close()
    
    return train, label, test, testID

def proData(paras, proteinFreq, cand, neg):
    """
    get protein freq svm files
    """
    logger.info('Recovering codon usage')
    
    if os.path.exists('refPro') and os.path.getsize('refPro'):
        with open('refPro', 'rb') as af:
            protein = pickle.load(af)
    else:
        protein = genePrediction(paras.refn)
    
    # for random forest
    posids = set([str(i) for i in paras.nodes if paras.length[i] >= 300])
    sampleSize = len(protein.keys()) + len(posids)  + len(neg)
    label = np.empty((sampleSize))
    pro_trainX = np.empty((sampleSize, len(codon)))
    
    counter1 = 0
    with open(config.PROTRAIN,'w') as pf:
        for ele in protein:
            vect = normalization(protein[ele])
            pro_trainX[counter1] = vect
            label[counter1] = 1
            counter1 += 1
            pf.write('1 ')
            for pos, ele in enumerate(vect):
                if ele is 0: continue
                pf.write('{0}:{1:.3f} '.format(pos+1, ele))
            pf.write('\n')
        
        for ele in posids:
            if ele in proteinFreq:
                vect = normalization(proteinFreq[ele])
                pro_trainX[counter1], label[counter1] = vect, 1
                counter1 += 1
                pf.write('1 ')
                for pos, ele in enumerate(vect):
                    if ele is 0: continue
                    pf.write('{0}:{1:.3f} '.format(pos+1, ele))
                pf.write('\n')
        
        for ele in neg:
            if ele in proteinFreq:
                vect = normalization(proteinFreq[ele])
                pro_trainX[counter1], label[counter1] = vect, -1
                counter1 += 1
                pf.write('-1 ')
                for pos, ele in enumerate(vect):
                    if ele is 0: continue
                    pf.write('{0}:{1:.3f} '.format(pos+1, ele))
                pf.write('\n')
    
    # NOTE: there are nodes that have no codon usage info
    pro_trainX = pro_trainX[:counter1]
    label = label[:counter1]
    
    counter2 = 0
    test = np.empty((len(cand), len(codon)))
    testID = np.empty((len(cand)))
    
    with open(config.PROCANDI, 'w') as candifile:
        for seqid in cand:                   
            if seqid in proteinFreq:
                vect = normalization(proteinFreq[seqid])
                test[counter2], testID[counter2] = vect, int(seqid)
                counter2 += 1
                candifile.write(seqid + ' ')
                for pos, ele in enumerate(vect):
                    if ele is 0: continue
                    candifile.write('{0}:{1:.3f} '.format(pos+1, ele))
                candifile.write('\n')
    test = test[:counter2]
    testID = testID[:counter2]
    
    return pro_trainX, label, test, testID

def graidRCF(trainX, label, test, threads, model=False):
    """
    random forest to predict the label
    """
    if os.path.exists(model) and os.path.getsize(model):
        logger.info('Optimal file exists, use it directly!')
        
        with open(model, 'r') as af:
            tunedParameters = pickle.load(af)
        
        P = tunedParameters['features__percentile']
        D = tunedParameters['rfc__max_depth']
        N = tunedParameters['rfc__n_estimators']
        W = tunedParameters['rfc__class_weight']
        
        fs = feature_selection.SelectPercentile(feature_selection.chi2, percentile = P)
        X_train_fs = fs.fit_transform(trainX, label)
        newTest = fs.transform(test)
        
        optimalEstimator = RandomForestClassifier(n_estimators = N,
                                                  max_depth = D,
                                                  class_weight = W
                                                  )
        optimalEstimator.fit(X_train_fs, label)
        
        # predicting
        y_pred = optimalEstimator.predict(newTest)
        return y_pred
    
    logger.warning('model directory {0} not exists, training a new one'.format(model))
    # tune the parameters: k, d, n and also do cross-validation
    rfc = RandomForestClassifier()
    selection = feature_selection.SelectPercentile(feature_selection.chi2)
    max_depth = np.linspace(5, 200, 10)
    n_estimators = [100, 500, 1000]
    CV = KFold(n = trainX.shape[0], n_folds = 10)
    
    pipeline = Pipeline([("features", selection), ("rfc", rfc)])
    param_grid = dict(features__percentile = range(20, 100, 20),
                       rfc__n_estimators = n_estimators,
                       rfc__max_depth = max_depth,
                       rfc__class_weight = ['subsample']
                       )
    
    classifier = GridSearchCV(pipeline,
                              param_grid=param_grid,
                              cv = CV,
                              n_jobs = threads,
                              scoring = 'accuracy')
        
    classifier.fit(trainX, label)
    
    P = classifier.best_params_['features__percentile']
    D = classifier.best_params_['rfc__max_depth']
    N = classifier.best_params_['rfc__n_estimators']
    logger.info('Optimal percentage of features:{0}'.format(P))
    logger.info('Optimal number of trees is {0}'.format(N))
    logger.info('Optimal max depth is {0}'.format(D))
    
    # reduce number of features
    fs = feature_selection.SelectPercentile(feature_selection.chi2, percentile = P)
    X_train_fs = fs.fit_transform(trainX, label)
    newTest = fs.transform(test)
    
    optimalEstimator = RandomForestClassifier(n_estimators = N,
                                              max_depth = D,
                                              class_weight = 'subsample'
                                              )
    
    optimalEstimator.fit(X_train_fs, label)
    
    # dump the model
    with open(model, 'wb') as af:
        pickle.dump(classifier.best_params_, af)
    
    # predicting
    y_pred = optimalEstimator.predict(newTest)
    
    return y_pred

def rcf(trainX, label, test, threads, model = False):
    """
    random forest to predict the label
    """
    prefix = 'rbf.pkl'
    if model:
        # no need to train a new model, load the old model to predict
        if os.path.isdir(model):
            logger.info('Random forest model exists, load it')
            curDir = os.getcwd()
            try:
                os.chdir(model)
                optimalEstimator = joblib.load(prefix)
                os.chdir(curDir)
                # predicting
                y_pred = optimalEstimator.predict(test)
                return y_pred
            except:
                logger.waring('no model found in the directory: {0}, trainning a new one'.format(model))
                os.chdir(curDir)
                shutil.rmtree(model)
    
    # de novo predict
    logger.warning('model directory {0} not exists, training a new one'.format(model))
    max_depth = np.linspace(5, 200, 20)
    n_estimators = [100, 500, 1000]
    
    mlt = RandomForestClassifier() # meachine learning tool
    # first grid search to find the optimal parameters
    parameters = dict(n_estimators = n_estimators,
                      max_depth = max_depth,
                      class_weight = ['subsample']
                      )
    
    # cross validation: shuffle data
    #CV = ShuffleSplit(trainX.shape[0],
    #                  n_iter=10,
    #                  test_size=0.2
    #                  )
    
    # cross validation: k-fold
    CV = KFold(n = trainX.shape[0],
               n_folds = 10
               )
    
    classifier = GridSearchCV(estimator = mlt,
                              cv = CV,
                              param_grid = parameters,
                              n_jobs = threads,
                              scoring = 'accuracy'
                              )
    
    classifier.fit(trainX, label)
    D = classifier.best_estimator_.max_depth
    N = classifier.best_estimator_.n_estimators
    
    logger.info('Optimal max depth is {0}, number of trees is {1}'.format(D, N))
    
    # weight = np.array([0.1 if i==1 else 1 for i in label])
    optimalEstimator = RandomForestClassifier(n_estimators = N,
                                              max_depth = D,
                                              class_weight = 'subsample'
                                              )

    # training
    optimalEstimator.fit(trainX, label)
    
    # dump the model
    curDir = os.getcwd()
    if os.path.isdir(model):
        shutil.rmtree(model)
    os.makedirs(model)
    os.chdir(model)
    joblib.dump(optimalEstimator, prefix)
    os.chdir(curDir)
    
    # predicting
    y_pred = optimalEstimator.predict(test)
    
    return y_pred

def svm(trainfn, candifn, penalty, loadSVM):
    
    logger.print_timestamp()
    if os.path.exists(loadSVM) and os.path.getsize(loadSVM):
        logger.info('Loading SVM parameters from: {}'.format(loadSVM))
        with open(loadSVM, 'r') as af:
            param = pickle.load(af)
    else:
        logger.notice('No SVM parameters file is found, trainning a new one!')
        rate, param = find_parameters(trainfn, '-gnuplot null')
        with open(loadSVM, 'w') as af:
            pickle.dump(param, af)
    
    # training
    logger.info('Running SVM with penalty = {}'.format(penalty))
    option = '-c {0} -g {1} -w1 {2} -w-1 1 -h 0'.format(param['c'], param['g'], penalty)
    y, x = svm_read_problem(trainfn)
    model = svm_train(y, x, option)
    
    # predict
    labels, candi = svm_read_problem(candifn)
    if not candi:
        logger.warning('SVM failed, no contig is recovered.')
        return set()
   
    #labels= set([int(x) for x in labels])
    pLabel, pAcc, pVal = svm_predict(labels, candi, model)
    
    cand = set([int(ele[0]) for ele in itertools.izip(labels, pLabel) if ele[1]>0])
    logger.info('Get {0:,} candidates\n'.format(len(cand)))
    
    return cand

def classification(soap, paras, threads):
    """
    two seperate classfication methonds: svm and random forest, and take the shared results
    """
     # determine the parameters
    logger.info('Average query length: {0:,}'.format(paras.meanLen))
    global _penalty
    if paras.meanLen >= 3000: _penalty = 1
    
    # output parameters for increaseing step
    if config.MDADEBUG:
        logger.debug('Increasing parameters:\n  intra miu1:  {0}\n  intra miu2: {1}\n  svm penalty: {2}\n'.\
                    format(config.INTRAmiu1, config.INTRAmiu2, _penalty))
        
    # get the coding frequencies of the contigs
    if os.path.exists(soap.PROFREQ) and os.path.getsize(soap.PROFREQ):
        with open(soap.PROFREQ, 'rb') as af:
            proteinFreq = pickle.load(af)
    else:
        logger.error('The file: ' + soap.PROFREQ + ' is not exists, please recheck it',
                     exit_with_code=1)
    
    # first recover the spearman candidates
    candL, negL = getCandiNodes(soap, paras, threads)
    
    # two machine learning methods
    if candL and negL:
        # get tetra candidates
        nuclTrain, label, cand, candID = nuclData(soap, paras, candL, negL)
        
        # svm codon usage
        nuclCand = svm(config.NUCTRAIN, config.NUCANDI, _penalty, 'paras_nucl')
        nuclCand = set([str(int(i)) for i in nuclCand])
        _, _, _, _ = proData(paras, proteinFreq, nuclCand, negL)
        svmL = svm(config.PROTRAIN, config.PROCANDI, _penalty, 'paras_pro')
        
        # random forest
        if randomForest:
            logger.info('Running random forest ...')
            nuclRes = rcf(nuclTrain, label, cand, threads, model='rbfNucl')
            #nuclRes2 = graidRCF(nuclTrain, label, cand, threads, model = 'RFCnucl')
            nuclPassed = set([str(int(ele[0])) for ele in itertools.izip(candID, nuclRes) if ele[1] == 1])
            # get the codon usage frequency
            proTrain, label, cand, candID = proData(paras, proteinFreq, nuclPassed, negL)
            proL = rcf(proTrain, label, cand, threads, model='rbfPro')
            #proL = graidRCF(proTrain, label, cand, threads, model = 'RFCpro')
            passedL = set([int(ele[0]) for ele in itertools.izip(candID, proL) if ele[1] == 1])
            logger.info('{0:,} candidates are got\n'.format(len(passedL)))
        else:
            passedL = []
        
        # combine the two sets of results or take the smaller dataset
        if randomForest:
            candidateL = svmL & passedL if randomForest else svmL
        else:
            candidateL = svmL
    else:
        candidateL = set()
    
    # remove files
    cleanUp('refGenes.lst', 'refPro', config.NUCTRAIN, config.NUCANDI, config.PROTRAIN, config.PROCANDI)
    candL =  set(map(int, candL))
    
    # store the file for later checking
    if config.MDADEBUG:
        with open('increaseFns','wb') as af:
            pickle.dump(candL, af)
            pickle.dump(candidateL, af)
    
    return candL, candidateL

