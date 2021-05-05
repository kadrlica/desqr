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
import pandas as pd
from utils import mkdir

if __name__ == "__main__":
    import parser
    parser = parser.Parser(description=__doc__)
    parser.add_argument('-k','--chunk',default=None,type=int,
                        help='number of exposures to run in a chunk')
    parser.add_argument('-s','--survey',default='refcat2',
                        help='calibration survey')
    parser.add_argument('-b','--band',default=None,type=str,action='append',
                        help='band to submit')
    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))
    zpsdir = mkdir(config['zpsdir'])
    explist = config['explist']
    tags = config.get('tags')
    section = config.get('db','bliss')

    chunk = 1 if args.chunk is None else args.chunk

    exposures = pd.read_csv(explist).to_records(index=False)
    exposures.dtype.names = [str(n).lower() for n in exposures.dtype.names]
    bands = config['bands'] if args.band is None else args.band

    for band in bands:
        expband = exposures[exposures['band'] == band]

        if len(expband) == 0:
            print("WARNING: No exposures found for %s-band; skipping..."%band)
            continue

        outdir = mkdir(os.path.join(zpsdir,band))
        logdir = mkdir(os.path.join(outdir,'log'))

        for exp in np.array_split(expband,len(expband)//chunk):
            #edict = dict(zip(exp[0].dtype.names,exp[0]))

            outbase = config['zpsbase']
            outfile = os.path.join(outdir,outbase)

            expnums = []
            for e in exp:
                edict = dict(zip(e.dtype.names,e))
                filename = outfile.format(**edict)
                if os.path.exists(filename) and not args.force:
                    print('Found %s; skipping...'%filename)
                else:
                    expnums.append(e['expnum']) 

            if not len(expnums): continue

            # Create a logfile named after the first exposure
            logfile = os.path.join(logdir,os.path.splitext(outbase)[0]+'.log')
            logfile = logfile.format(**edict)

            # If only one exposure, should skip if done
            #if len(exp)==1: outfile = outfile.format(**edict)
            #if os.path.exists(outfile) and not args.force:
            #    print('Found %s; skipping...'%outfile)
            #    continue

            expnums = ' '.join(str(e) for e in expnums)
            params = dict(expnum=expnums,
                          force='-f' if args.force else '',
                          verbose='-v' if args.verbose else '',
                          survey=args.survey,
                          outfile=outfile,explist=explist)

            cmd = "calibrate.py %(force)s %(verbose)s %(explist)s %(outfile)s -e %(expnum)s -s %(survey)s "%params
            submit = "csub -q %s -o %s "%(args.queue, logfile)
            if args.njobs: submit += "-n %s "%args.njobs
            submit += cmd
            if args.verbose: print(submit)
            if not args.dryrun: subprocess.call(submit,shell=True)
            time.sleep(args.sleep)

            #break
        #break
        #time.sleep(300)

