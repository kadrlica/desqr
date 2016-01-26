#!/usr/bin/env python
import os
import yaml
import subprocess
import time
import glob

import numpy as np
import upload


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
    table = config['table']
    index = config['index']
    logdir = os.path.join('upload')

    #infiles = ['unique/hpx_01313.fits','unique/hpx_01314.fits','unique/hpx_01315.fits']
    #infiles = ['unique/hpx_02957.fits','unique/hpx_02958.fits','unique/hpx_02959.fits']
    #infiles = sorted(glob.glob(unidir+'/hpx_0295*.fits'))

    ### # Upload the catalog
    ### infiles = sorted(glob.glob(catdir+'/*.fits'))
    ### params = (table,force,' '.join(infiles))
    ### cmd = "upload.py -t %s %s %s"%params
    ### logfile = os.path.join(logdir,'catalog_upload.log')
    ### if opts.queue == 'local':
    ###     submit = "%s > %s"%(cmd,logfile)
    ### else:
    ###     submit = "csub -o %s %s"%(logfile,cmd)
    ### subprocess.call(submit,shell=True)
    ### time.sleep(1)

    # Upload the keys
    infiles = sorted(glob.glob(keydir+'/*.fits'))
    params = (index,force,' '.join(infiles))
    cmd = "upload.py -t %s %s %s"%params
    logfile = os.path.join(logdir,'index_upload.log')
    if opts.queue == 'local':
        submit = "%s >  %s"%(cmd,logfile)
    else:
        submit = "csub -o %s %s"%(logfile,cmd)
    subprocess.call(submit,shell=True)
    time.sleep(1)
