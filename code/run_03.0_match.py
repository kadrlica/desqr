#!/usr/bin/env python
"""
Match catalog objects within and across input files.
"""

import os
import time
import yaml
import subprocess
import time
import glob

import fitsio
import numpy as np
import healpy

from ugali.utils.shell import mkdir
from const import OBJECT_ID

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config')
    parser.add_argument('-f','--force',action='store_true')
    # Sometimes sleep needs to be 30s or more...
    parser.add_argument('-s','--sleep',default=5,type=float)
    parser.add_argument('-n','--njobs',default=15,type=int)
    parser.add_argument('-q','--queue',default='condor')
    args = parser.parse_args()

    config = yaml.load(open(args.config))
    hpxdir = config['hpxdir']
    logdir = mkdir(os.path.join(hpxdir,'log'))
    radius = config['radius']
    #OBJECT_ID = config['objid']
    force = '-f' if args.force else ''
    
    for pix in np.arange(healpy.nside2npix(config['nside'])):
        infiles = glob.glob(hpxdir+'/*/*%05d*.fits'%pix)
        if len(infiles) == 0: continue

        fits = fitsio.FITS(infiles[0])
        if OBJECT_ID in fits[1].get_colnames() and not args.force:
            print "Found column '%s'; skipping..."%OBJECT_ID
            fits.close()
            continue
        fits.close()

        logfile = os.path.join(logdir,'hpx_%05d.log'%pix)

        params=(radius,' '.join(infiles),force)
        cmd = 'match.py -r %s %s %s'%params
        if args.queue == 'local':
            submit = cmd
        else:
            submit = 'csub -o %s -n %s %s'%(logfile,args.njobs,cmd)
        subprocess.call(submit,shell=True)
        if args.queue != 'local': time.sleep(args.sleep)
