#!/usr/bin/env python
"""
Run photometric calibration against reference catalog.
"""
import os
from os.path import getmtime
import yaml
import subprocess
import time
import glob

import numpy as np
import pandas as pd
from utils import mkdir, is_found, found

if __name__ == "__main__":
    import parser
    parser = parser.Parser(description=__doc__)
    parser.add_argument('-k','--chunk',default=None,type=int,
                        help='number of exposures to run in a chunk')
    parser.add_argument('-s','--survey',default='refcat2',
                        help='calibration survey')
    parser.add_argument('-b','--band',default=None,type=str,action='append',
                        help='band to submit')
    parser.add_argument('--expmin',default=0,type=int,
                        help='minimum exposure number')
    parser.add_argument('--expmax',default=np.inf,type=int,
                        help='maximum exposure number')
    parser.add_argument('--transform',default='linear',choices=['linear','interp'],
                        help='transformation type')
    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))
    zpsdir = mkdir(config['zpsdir'])
    explist = config['explist']

    chunk = 1 if args.chunk is None else args.chunk

    exposures = pd.read_csv(explist).to_records(index=False)
    exposures.dtype.names = [str(n).lower() for n in exposures.dtype.names]
    
    #Hack to subselect exposures
    sel = np.ones(len(exposures),dtype=bool)
    #sel &= (exposures['teldec'] > -40) & (exposures['teldec'] < -20)
    #sel &= (exposures['propid'] == '2012B-0001')
    exposures = exposures[sel]

    bands = config['bands'] if args.band is None else args.band
    outbase = config['zpsbase']
    #procbase = '{expnum}.log'
    procbase = 'D{expnum:08d}_{band}_25_r1p1_immask.fits.fz'

    for band in bands:
        expband = exposures[exposures['band'] == band]

        if len(expband) == 0:
            print("WARNING: No exposures found for %s-band; skipping..."%band)
            continue

        outdir = mkdir(os.path.join(zpsdir,band))
        logdir = mkdir(os.path.join(outdir,'log'))

        for exp in np.array_split(expband,len(expband)//chunk):
            #edict = dict(zip(exp[0].dtype.names,exp[0]))

            outfile = os.path.join(outdir,outbase)

            expnums = []
            for e in exp:
                if (e['expnum'] > args.expmax) or (e['expnum'] < args.expmin):
                    continue

                edict = dict(zip(e.dtype.names,e))
                filename = outfile.format(**edict)
                # If processing log modified more recently than zp
                # Use try/except to minimize filesystem calls
                # NOTE: Remove first os.path.exists(filename) for speed
                if not args.force and os.path.exists(filename):
                    procfile = os.path.join(e['path'],procbase.format(**edict))
                    try:
                        if getmtime(procfile) < getmtime(filename): 
                            # Both files exist and processing is older
                            found(filename)
                            continue
                        else:
                            # Both files exist and processing is newer
                            print("WARNING: New processing for %s..."%e['expnum'])
                    except OSError:
                        if not os.path.exists(procfile):
                            # No processing file
                            print("WARNING: %s not found; skipping...."%procfile)
                            continue
                        else:
                            # No zeropoint file
                            pass

                print("Adding %s..."%e['expnum'])
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
                          force='-f', # controlled above
                          verbose='-v' if args.verbose else '',
                          survey=args.survey, transform=args.transform,
                          outfile=outfile,explist=explist)

            cmd = "calibrate.py %(force)s %(verbose)s %(explist)s %(outfile)s -e %(expnum)s -s %(survey)s --transform %(transform)s"%params
            submit = "csub -q %s -o %s "%(args.queue, logfile)
            if args.njobs: submit += "-n %s "%args.njobs
            submit += cmd
            if args.verbose: print(submit)
            if not args.dryrun: subprocess.call(submit,shell=True)
            time.sleep(args.sleep)

            #break
        #break
        #time.sleep(300)

