############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import datetime
import os
import platform
import sys

def check_python_version():
    if sys.version[0:3] not in supportedPythonVersions:
        sys.stderr.write("ERROR! Python version " + sys.version[0:3] + " is not supported!\n" +\
                         "Supported versions are " + ", ".join(supportedPythonVersions) + "\n")
        sys.exit(1)

mainLocation = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
libsLocation = os.path.join(mainLocation, 'libs')

# the versions of python are supported
supportedPythonVersions = ['2.7', '3.0']

"""
CONSTANT
"""

# logger names
logBUILD = 'build'
logFAB = 'fab'
logMGA = 'mga'

def softVersion():
    
    version_fpath = os.path.join(libsLocation, '..', 'VERSION')
    version = "unknown"
    build = "unknown"
    if os.path.isfile(version_fpath):
        version_file = open(version_fpath)
        version = version_file.readline()
        if version:
            version = version.strip()
        else:
            version = "unknown"
        build = version_file.readline()
        if build:
            build = build.strip().lower()
        else:
            build = "unknown"
    return version, build


#########################################################
#                    BUILD                                #
#########################################################

# decomposition paras
SMALLSIZE = 50  # decompose to smaller subgraphs
CRLEN = 200     # the determination length of crbranches

# choose the binning tools: tetracalc and metabat, default tetracal
# Note that if the meta-O reads were assembled by soapdenovo, only tetracal
# could be used. But if you want to use metabat, you can choose the -b parameter
# to tell metaSort, the assembly was obtained usingn other assemblers.

# Also note that if metabat was chose, the samtools should be installed
# on you computer and in you system path
BINMODE = 'metabat'

### tetrcala binning parameters
# binning parameter of tetracala softeware
# the smaller the value is, the higher of binning feasiblity is, but also
# the smaller bin size is
MERGELEVEL = 0.89 # for MDA
MERGELEVEL2 = 0.9 # for SOAPdenovo Contigs

# minmum length for contig binning of tetracala
# the larger the number, the higher of the accuracy, but also
# the smaller bin size
metaS_minLen = 5000 # for MDA, range from 2500 to 500000
metaO_minLen = 5000 # for SOAPdenovo Contigs, range from 2500 to 500000

#########################################################
#                    BAF                                #
#########################################################
DEBUG = 1

# min scaffolds length for gc breaks
# the smaller window size is, the more fragmented are

winSize = 300  # ranging from 100 to 500

# the minmum scaffolds length for gc break
minECLen = 50000000 # ranging from 50,000 to maximum contig length

# the fasta file storing the breaked scaffolds
errFreefn = 'EC.fa'

### Nucmer alignment files
# prefix of the generated nucmer mapping result
NUCMER = 'align'

# FAB script output tmple files
LMETACL = 'lmetaCls'
# A bin, which has a size smaller than this value, will be discarded
SIZECUF = 1000000 # ranging from 1 to 5,000,000
SOFTDIR = 'soapClseqs'

CLUSTERSEQDIR = 'clSeqs'
IDENTITY = 'identity'
PARTIALCUT = 200

# mapping info between soapdenovo contigs and SPAdes scaffolds
GUILDLINKS = 'guildLinks'

# mapping cluster name
FISHCLUSTERS = 'fishCls'
MIXEDBIN = 'mixedBin'


# paths library name
LIBPATHS = 'libpaths'

mincluster = 65
#########################################################
#                    MGA                                #
#########################################################

# for debuging
MDADEBUG = 0

# The shared kmer in a species level cutoff
# The database and the query are split into 51-mers, compared and then used for calculating shared kmers
SPECIES = 0.6 # range 0.5 to 0.8

# control the number of svm recovered positive data nubmer
FOLD = 20 # ranging from 1 to 50

# penalty value if imbanlance datasets are used
PENALTY = 0.45 # range from 0.1 to 1

# output files name
NUCTRAIN = 'nuclTraining'
NUCANDI = 'nuclCandi'
PROTRAIN = 'proTraining'
PROCANDI = 'proCandi'
PROTRAINS = 'prots'
PROCANDIS = 'procs'

# GMM paras
COMPONENTNUM = range(2,20,2)
COVARIANCE = 'full'
ITERS = 20000
INITS = 5
INTRACFL = 0.2
INTERCFL = 0.1

# max step size for dfs searching canidate nodes
SIZELEVEL1 = 3
SIZELEVEL2 = 2
SIZELEVEL3 = 1
FILTERSIZE = 6
conserPEstep = [1]
greedyPEstep = [1]

# step size for searching triangle structure
TRIANGELIMIT = 3

# contigs that compose the scaffold
DPATH = 'dpath'
PPATH = 'ppath'

# output the final file
MDAlenCf = 5000000
COMPOSITION = 'linkedContigs'
RECOVERED = 'recover.fa'
SEED = 'seeds.fa'
MDASEEDS  = 'mda.fa'
STATISTIC = 'scafStatistics'
FINALRES = 'targetGenome.fa'

# subspecies file names
STATUS = 'node_status'
VARS = 'variations'
CLEANLENCUTOFF = 5000

# detect the operating system
if platform.system() == 'Darwin':
    platform_name = 'macosx'
else:
    if platform.architecture()[0] == '64bit':
        platform_name = 'linux_64'
    else:
        platform_name = 'linux_32'


#########################################################
#                    LOST DLINK                         #
#########################################################

QULIATY = 20

RATEOFLOOSE = 0.05

RATEOFSTRICT = 0.1

SUPPORTNUM = 2

RANGEOFLIB = 1.5

RATE4CUT = 0.8

PAIRS4DEL = 1

PAIRS4IGN = 2

RATE4SPLIT = 2

RATE4REM = 3

RATE4DIV = 2

ALIGNCUTOFF = 0.9



