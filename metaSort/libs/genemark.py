############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import shutil
import config
from log import get_logger

# set log file
logger = get_logger()

def install_genemark(tool_dirpath):
    """Installation instructions for GeneMark Suite.

    Please, copy key "gm_key" into users home directory as:
    cp gm_key ~/.gm_key
    (genemark_suite_linux_XX/gmsuite/INSTALL)
    """
    import subprocess
    import filecmp
    gm_key_fpath = os.path.join(tool_dirpath, 'gm_key')
    gm_key_dst = os.path.expanduser('~/.gm_key')
    if not os.path.isfile(gm_key_dst) or \
        (not filecmp.cmp(gm_key_dst, gm_key_fpath) and os.path.getmtime(gm_key_dst) < os.path.getmtime(gm_key_fpath)):
        shutil.copyfile(gm_key_fpath, gm_key_dst)
    # checking the installation
    tool_exec_fpath = os.path.join(tool_dirpath, 'gmhmmp')
    proc = subprocess.Popen([tool_exec_fpath], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while not proc.poll():
        line = proc.stdout.readline()
        if line.find('license period has ended') != -1:
            logger.warning('License period for GeneMark has ended! \n' \
                           'To update license, please visit http://topaz.gatech.edu/license_download.cgi page and fill in the form.\n' \
                           'You should choose GeneMarkS tool and your operating system (note that GeneMark is free for non-commercial use).\n' \
                           'Download the license key and replace your ~/.gm_key with the updated version. After that you can restart using this softeware.\n')
            return False
    return tool_exec_fpath

def getPath():
    
    tool_name = 'MetaGeneMark'
    tool_dirname = 'metagenemark'
    tool_dirpath = os.path.join(config.libsLocation, tool_dirname, config.platform_name)
    if not os.path.exists(tool_dirpath):
        raise Exception('Sorry, can\'t use %s on this platform, skipping gene prediction.' % tool_name)
    else:
        successful = install_genemark(tool_dirpath)
        if not successful:
            logger.error('Sorry, genemark is not installed successfully.',
                         exit_with_code=1
                         )
        else:
            return tool_dirpath
            
