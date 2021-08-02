############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################
import config
from log import get_logger

# set log file
logger = get_logger()

def contigInfo(fn):
    """
    Get the info of fasta file including n50, average length et al.
    """
    seqLen = []
    gc = 0
    try:
        with open(fn,'r') as a_file:
            for content in a_file:
                if '>' in content:
                    seqLen.append(0)
                else:
                    seq = content.rstrip().upper()
                    seqLen[-1] += len(seq)
                    gc += len([i for i in seq if i=='G' or i == 'C'])
    except:
        logger.error('Can not open the file: ' + fn,
                     exit_with_code=1
                     )
    
    seqLen.sort(reverse = True)
    TotalLen = sum(seqLen)
    maxLen = seqLen[0]    
    averLen = TotalLen/len(seqLen)    
    GCrate = round(gc*100.0/TotalLen, 2)
    
    len_add = maxLen
    for value in seqLen:
        if len_add < TotalLen/2:
            len_add += value
        else:
            N50 = value
            break
    
    return len(seqLen), TotalLen, maxLen, N50, averLen, GCrate

def scaffoldInfo(fn):
    seqLen = []
    a = t = c = g = n = 0
    try:
        with open(fn,'r') as a_file:
            for content in a_file:
                if '>' in content:
                    seqLen.append(0)
                else:
                    seq = content.rstrip().upper()
                    seqLen[-1] += len(seq)
                    for i in seq:
                        if i == 'A':
                            a += 1
                        elif i == 'T':
                            t += 1
                        elif i == 'C':
                            c += 1
                        elif i == 'G':
                            g += 1
                        else:
                            n += 1
    except:
        logger.error('Can not open the file: ' + fn,
                     exit_with_code=1
                     )
    
    seqNu = len(seqLen)
    seqLen.sort(reverse = True)
    TotalLen = sum(seqLen)
    maxLen = seqLen[0]    
    averLen = TotalLen/seqNu    
    Arate = round(a*100.0/TotalLen, 2)
    Trate = round(t*100.0/TotalLen, 2)
    Crate = round(c*100.0/TotalLen, 2)
    Grate = round(g*100.0/TotalLen, 2)
    Nrate = round(n*100.0/TotalLen, 2)
    GCrate = round((g+c)*100.0/TotalLen, 2)
    scaflg100 = len([i for i in seqLen if i >= 100])
    ratelg100 = round(scaflg100*1.0/seqNu, 2)
    scaflg500 = len([i for i in seqLen if i >= 500])
    ratelg500 = round(scaflg500*1.0/seqNu, 2)
    scaflg1k = len([i for i in seqLen if i >= 1000])
    ratelg1k = round(scaflg1k*1.0/seqNu, 2)
    scaflg10k = len([i for i in seqLen if i >= 10000])
    ratelg10k = round(scaflg10k*1.0/seqNu, 2)
    scaflg100k = len([i for i in seqLen if i >= 100000])
    ratelg100k = round(scaflg100k*1.0/seqNu, 2)
    scaflg1m = len([i for i in seqLen if i >= 1000000])
    ratelg1m = round(scaflg1m*1.0/seqNu, 2)
    # n10, 20, 30, 40, 50, 60, 70, 80, 90
    Ncount = [[]]
    sumLen = 0
    Ndata = 10
    for ele in seqLen:
        sumLen += ele
        if sumLen >= (TotalLen*Ndata)/100.0:
            Ncount[-1].append(ele)
            Ndata += 10
            if Ndata == 100:
                break
            Ncount.append([])
        else:
            Ncount[-1].append(ele)
    try:
        NLenInfo = []
        for ele in Ncount:
            NLenInfo.append((ele[-1], len(ele)))
    except:
        Ncount = [[]]
    
    return (seqNu, TotalLen, maxLen, NLenInfo[4][0], averLen, GCrate),\
           NLenInfo,\
           (a, Arate, t, Trate, c, Crate, g, Grate, n, Nrate),\
           (scaflg100, ratelg100, scaflg500, ratelg500, scaflg1k, ratelg1k, scaflg10k, ratelg10k,  scaflg100k, ratelg100k, scaflg1m, ratelg1m)

