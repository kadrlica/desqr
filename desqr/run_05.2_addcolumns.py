#!/usr/bin/env python
"""
Add de-reddened columns
"""
import os
from os.path import basename
import yaml
import subprocess
import time
import glob

import fitsio

from ugali.utils.shell import mkdir

COLUMNS = ['MAG_PSF','WAVG_MAG_PSF','MAG_AUTO','WAVG_MAG_AUTO']

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-s','--sleep',default=1,type=float)
    parser.add_argument('-n','--njobs',default=10,type=int)
    parser.add_argument('-q','--queue',default='local')
    args = parser.parse_args()

    config = yaml.load(open(args.config))
    catdir = config['catdir']
    bands  = config['bands']

    infiles = sorted(glob.glob(catdir+'/*.fits'))
    logdir = mkdir(os.path.join(catdir,'log'))
    for col in COLUMNS:
        for b in bands:
            logbase = 'add_%s_sfd_%s.log'%(col.lower(),b.lower())
            logfile = os.path.join(logdir,logbase)
            params = dict(band=b.upper(),infiles=' '.join(infiles),
                          force='-f' if args.force else '',logfile=logfile,
                          column = col)
            params['outcol'] = '%(column)s_SFD_%(band)s'%params
            print("Adding %(outcol)s..."%params)

            # The nested use of single and double quotes is really
            # nasty with csub from the shell, and is (nearly?) to
            # impossible from a python call to subprocess
            cmd = """add_column.py %(force)s --column %(outcol)s --formula "data['%(column)s_%(band)s']-data['EXTINCTION_%(band)s']" %(infiles)s | tee %(logfile)s"""%params
            #print(cmd)
            subprocess.call(cmd,shell=True)
            print("\n")
