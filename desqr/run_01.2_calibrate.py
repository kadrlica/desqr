#!/usr/bin/env python
"""
Run photometric calibration against ATLAS-REFCAT2.
"""
import os
import yaml
import subprocess
import time
import glob

import numpy as np
from utils import mkdir
import download

if __name__ == "__main__":
    import parser
    parser = parser.Parser(description=__doc__)
    parser.add_argument('-k','--chunk',default=None,type=int,
                        help='number of exposures to run in a chunk')
    parser.add_argument('-b','--band',default=None,type=str,action='append',
                        help='band to submit')
    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))
    zpsdir = mkdir(config['zpsdir'])
    explist = config['explist']
    tags = config.get('tags')
    section = config.get('db','bliss')

    chunk = 1 if args.chunk is None else args.chunk

    exposures = np.recfromcsv(explist)
    bands = config['bands'] if args.band is None else args.band

    for band in bands:
        expband = exposures[exposures['band'] == band]
        
        outdir = mkdir(os.path.join(zpsdir,band))
        logdir = mkdir(os.path.join(outdir,'log'))

        for exp in np.array_split(expband,len(expband)//chunk):
            edict = dict(zip(exp[0].dtype.names,exp[0]))

            outbase = config['zpsbase']
            outfile = os.path.join(outdir,outbase)

            # If only one exposure, should skip if done
            if len(exp)==1: outfile = outfile.format(**edict)
            if os.path.exists(outfile) and not args.force:
                print('Found %s; skipping...'%outfile)
                continue

            # Create a logfile named after the first exposure
            logfile = os.path.join(logdir,os.path.splitext(outbase)[0]+'.log')
            logfile = logfile.format(**edict)

            expnums = ' '.join(str(e) for e in exp['expnum'])
            params = dict(expnum=expnums,
                          force='-f' if args.force else '',
                          verbose='-v' if args.verbose else '',
                          outfile=outfile,explist=explist)

            cmd = "calibrate.py %(force)s %(verbose)s %(explist)s %(outfile)s -e %(expnum)s "%params
            submit = "csub -q %s -o %s "%(args.queue, logfile)
            if args.njobs: submit += "-n %s "%args.njobs
            submit += cmd
            if args.verbose: print(submit)
            subprocess.call(submit,shell=True)
            time.sleep(args.sleep)

            #break
        #break
        #time.sleep(300)

