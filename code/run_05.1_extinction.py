#!/usr/bin/env python
"""
Apply zeropoints and extinction correction to catalog files.
"""
import os
from os.path import basename
import yaml
import subprocess
import time
import glob

import fitsio

from ugali.utils.shell import mkdir

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-s','--sleep',default=1,type=float)
    parser.add_argument('-n','--njobs',default=10,type=int)
    parser.add_argument('-q','--queue',default='condor')
    args = parser.parse_args()

    config = yaml.load(open(args.config))
    catdir = config['catdir']
    ebv    = config.get('ebv',None)
    bands  = config['bands']

    infiles = sorted(glob.glob(catdir+'/*.fits'))
    logdir = mkdir(os.path.join(catdir,'log'))

    for infile in infiles:
        if not args.force:
            fits = fitsio.FITS(infile)
            if 'EBV' in fits[1].get_colnames():
                print "Found column 'EBV'; skipping %s..."%(basename(infile))
                continue

        logbase = ('ebv_'+os.path.basename(infile)).replace('.fits','.log')     
        logfile = os.path.join(logdir,logbase)
        params = dict(infile=infile,force='-f' if args.force else '',
                      extinction = '-e %s'%ebv if ebv else '',
                      bands = '-b '+' -b '.join(bands))

        cmd = "extinction.py %(force)s %(extinction)s %(bands)s %(infile)s"%params

        if args.queue == 'local':
            print cmd
            submit = cmd
        else:
            submit = 'csub -o %s -n %s "%s"'%(logfile,args.njobs,cmd)
        subprocess.call(submit,shell=True)
        time.sleep(args.sleep)
