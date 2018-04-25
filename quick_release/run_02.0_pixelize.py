#!/usr/bin/env python
import os
import yaml
import subprocess
import time
import glob

from ugali.utils.shell import mkdir
from pixelize import pixelize

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-s','--sleep',default=1,type=float)
    parser.add_argument('-n','--njobs',default=None,type=int)
    parser.add_argument('-q','--queue',default='condor')
    args = parser.parse_args()

    config = yaml.load(open(args.config))
    rawdir = config['rawdir']
    hpxdir = config['hpxdir']

    for band in config['bands']:
        indir = os.path.join(rawdir,band)
        outdir = mkdir(os.path.join(hpxdir,band))
        logdir = mkdir(os.path.join(outdir,'log'))
        logfile = os.path.join(logdir,'pixelize_%s.log'%band)

        outbase = 'hpx_%s'%band+'_%05d.fits'
        cmd = 'pixelize.py %s %s -o %s -n %i'%(indir,outdir,outbase,config['nside'])

        if args.queue == 'local':
            submit = cmd
        else:
            submit = 'csub -o %s '%(logfile)
            if args.njobs: submit += '-n %s '%args.njobs
            submit += cmd
        subprocess.call(submit,shell=True)
        time.sleep(args.sleep)
        #break
