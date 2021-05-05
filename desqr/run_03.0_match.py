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
import healpy as hp

from ugali.utils.shell import mkdir
from const import OBJECT_ID

if __name__ == "__main__":
    from parser import Parser
    parser = Parser(description=__doc__)
    parser.set_defaults(njobs=10)
    parser.add_argument('-m','--mlimit',default=40,type=float,
                        help='memory limit (GB)')
    parser.add_argument('-p','--pix',action='append',type=int)
    args = parser.parse_args()

    config = yaml.load(open(args.config))
    hpxdir = config['hpxdir']
    logdir = mkdir(os.path.join(hpxdir,'log'))
    radius = config['radius']
    force = '-f' if args.force else ''
    mlimit = '-m %s'%args.mlimit if args.mlimit else ''

    if args.pix: pixels = args.pix
    else: pixels = np.arange(hp.nside2npix(config['nside']))
        
    for pix in pixels:
        infiles = glob.glob(hpxdir+'/*/*%05d*.fits'%pix)
        if len(infiles) == 0: continue

        done = (not args.force)
        for f in infiles:
            if not done: break
            fits = fitsio.FITS(f)
            try:
                done = (OBJECT_ID in fits[1].get_colnames())
            except IOError as e: 
                print(f)
                raise(e)
            fits.close()

        # Not really necessary to compare with force again...
        if done and not args.force: 
            print("Found column '%s'; skipping %05d..."%(OBJECT_ID,pix))
            continue

        logfile = os.path.join(logdir,'hpx_%05d.log'%pix)

        params=(radius,force,mlimit,' '.join(infiles))
        cmd = 'match.py -r %s %s %s %s'%params
        submit = 'csub -q %s -o %s -n %s %s'%(args.queue,logfile,args.njobs,cmd)
        if args.verbose: print(submit)
        if not args.dryrun: subprocess.call(submit,shell=True)
        time.sleep(args.sleep)
