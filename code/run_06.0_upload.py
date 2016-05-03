#!/usr/bin/env python
import os
import yaml
import subprocess
import time
import glob

import numpy as np
import upload
from utils import mkdir

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-s','--sleep',default=1,type=float)
    parser.add_argument('-q','--queue',default='local')
    opts = parser.parse_args()

    force = '-f' if opts.force else ''

    config = yaml.load(open(opts.config))
    catdir = config['catdir']
    keydir = config['keydir']
    section = config['db']
    table = config['dbtable']
    index = config['dbindex']

    logdir = mkdir('upload')

    ### # Upload the catalog
    infiles = sorted(glob.glob(catdir+'/*.fits'))
    params = (force,table,section,' '.join(infiles))
    cmd = "upload.py %s -t %s -s %s %s"%params
    logfile = os.path.join(logdir,'catalog_upload.log')
    submit = 'csub -o %s -u %s "%s"'%(logfile,opts.queue,cmd)
    subprocess.call(submit,shell=True)
    time.sleep(1)

    # Upload the keys
    infiles = sorted(glob.glob(keydir+'/*.fits'))
    params = (force,index,section,' '.join(infiles))
    cmd = "upload.py %s -t %s -s %s %s"%params
    logfile = os.path.join(logdir,'index_upload.log')
    submit = 'csub -o %s -u %s "%s"'%(logfile,opts.queue,cmd)
    subprocess.call(submit,shell=True)
    time.sleep(1)
