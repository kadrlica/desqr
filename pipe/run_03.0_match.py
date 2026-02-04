#!/usr/bin/env python
"""
Match catalog objects within healpix.
"""
import os
import time
import subprocess
import glob
import yaml

import fitsio
import numpy as np
import healpy as hp

from utils import mkdir
from const import OBJECT_ID

if __name__ == "__main__":
    from parser import Parser
    parser = Parser(description=__doc__)
    parser.set_defaults(njobs=12) 
    parser.add_argument('-m','--mlimit',default=40,type=int,
                        help='memory limit (GB)')
    args = parser.parse_args()

    config = args.config
    hpxdir = config['hpxdir']
    logdir = mkdir(os.path.join(hpxdir,'log'))
    radius = config['radius']
    force = '-f' if args.force else ''
    mlimit = '-m %s'%args.mlimit if args.mlimit else ''

    if args.pix: pixels = args.pix
    else: pixels = np.arange(hp.nside2npix(config['nside']))

    print("Starting matching...")
    
    for i,pix in enumerate(pixels):
        infiles = glob.glob(hpxdir+'/*/*%05d.fits'%pix)
        if len(infiles) == 0: continue
        print("(%s/%s): hpx %05d"%(i+1,len(pixels), pix))

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
