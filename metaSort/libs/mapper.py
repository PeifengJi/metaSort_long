############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################
import os
import platform
import subprocess
import config
from log import get_logger

# set log file
logger = get_logger()

# exce used 
binaries = ['nucmer', 'delta-filter', 'show-coords']
if platform.system() == 'Darwin':
    mummerPath = os.path.join(config.libsLocation, 'MUMmer3.23-osx')
else:
    mummerPath = os.path.join(config.libsLocation, 'MUMmer3.23-linux')

def binPath(fname):
    return os.path.join(mummerPath, fname)

def run_nucmer(refPath, contigPath, prefix):
    cmdNuc = [binPath('nucmer'),
         '-c', '65',
         '-l', '65',
         '--maxmatch',
         '-p', prefix,
         refPath,
         contigPath]
    
    # print cmd
    logger.print_command_line(cmdNuc, wrap_after=None)
    
    returnCode = subprocess.call(cmdNuc, stdout = subprocess.PIPE)
    if returnCode != 0:
        logger.error('Nucmer failed',
                to_stderr=True,
                exit_with_code=3
                )
    
    deltaf = prefix + '.delta'
    cmdCor = [binPath('show-coords'),
              '-c',
              '-l',
              '-r',
              '-T',
              '-d',
              deltaf]
    
    # print cmd
    logger.print_command_line(cmdCor, wrap_after=None)
    
    resultName = prefix + '.result'
    returnCode = subprocess.call(
                    cmdCor,
                    stdout=open(resultName, 'w'))
    
    if returnCode != 0:
        logger.error('Show-coords failed',
                to_stderr=True,
                exit_with_code=3
                )
    # dele the delta file
    os.remove(deltaf)

def binariesExists(mummerPath):
    for requiredf in binaries:
        if not os.path.isfile(os.path.join(mummerPath, requiredf)):
            return 0
    return 1

def pipeline(refPath, contigPath, prefix):
    
    logger.info('Check if MUMmer is installed ...')
    currentDir = os.getcwd()
    if not binariesExists(mummerPath):
        # making
        logger.info('MUMmeris not installed, compiling it...')
        
        os.chdir(mummerPath)
        returnCode = subprocess.call(['make', 'check'])
        returnCode = subprocess.call(['make', 'install'])
        
        if returnCode != 0 or not binariesExists(mummerPath):
            logger.error('Failed to compile MUMmer ' +  mummerPath +' Try to compile it manually.',
                to_stderr=True,
                exit_with_code=3
                )
    else:
        logger.info('OK, the software is installed, now mapping ...')
    
    os.chdir(currentDir)
    run_nucmer(refPath, contigPath, prefix)

