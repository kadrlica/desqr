#!/usr/bin/env python
""" Downlink exposures that are already skimmed """
import os, shutil
from os.path import getmtime
import yaml
import subprocess
import time
import glob

import numpy as np
import pandas as pd

from utils import mkdir, found, is_found

if __name__ == "__main__":
    from parser import Parser
    parser = Parser(description=__doc__)
    parser.add_argument('--copy', action='store_true')
    args = parser.parse_args()

    config = args.config
    rawdir = config['rawdir']
    explist = config['explist']
    tags = config.get('tags')
    basename = config.get('rawbase','D{expnum:08d}_{band:s}_cat.fits')

    exposures = pd.read_csv(explist).to_records(index=False)
    for band in config['bands']:
        if args.bands and (band not in args.bands): continue

        outdir = mkdir(os.path.join(rawdir,band))
        logdir = mkdir(os.path.join(outdir,'log'))

        sel = exposures['band'] == band
        for i, exp in enumerate(exposures[sel]):
            print("(%s/%s): %s"%(i+1,sel.sum(), exp['expnum']))

            if (tags is not None) and (exp['tag'] not in tags):
                print("Skipping exposure with tag: '%s'"%exp['tag'])
                continue

            try:
                float(exp['path'])
                print("Missing path for: %(expnum)d"%exp)
                continue
            except:
                pass
            
            edict = dict(zip(exp.dtype.names,exp))
            outbase = basename.format(**edict)
            outfile = os.path.join(outdir,outbase)
            procfile = os.path.join(exp['path'])

            if os.path.exists(outfile):
                if args.force:
                    os.remove(outfile)
                else:
                    print("Found %s; skipping..."%outfile)
                    continue

            if args.copy:
                if args.verbose: print("Copying %(path)s"%exp)
                shutil.copy(exp['path'], outfile)
            else:
                if args.verbose: print("Linking %(path)s"%exp)
                os.symlink(exp['path'], outfile)

            time.sleep(args.sleep)

