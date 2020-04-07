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
    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))
    zpsdir = mkdir(config['zpsdir'])
    explist = config['explist']
    tags = config.get('tags')
    section = config.get('db','bliss')

    exposures = np.recfromcsv(explist)
    for band in config['bands']:
        outdir = mkdir(os.path.join(zpsdir,band))
        logdir = mkdir(os.path.join(outdir,'log'))
        for exp in exposures[exposures['band'] == band]:
            print("WARNING: selecting on expnum!")
            if exp['expnum'] < 600000: continue
            if (tags is not None) and (exp['tag'] not in tags):
                print("No exposures with tag '%s'"%tags)
                continue
            outbase = config['zpsbase']%exp
            outfile = os.path.join(outdir,outbase)
            if os.path.exists(outfile) and not args.force:
                print 'Found %s; skipping...'%outfile
                continue

            logfile = os.path.join(logdir,os.path.splitext(outbase)[0]+'.log')
            params = dict(expnum=exp['expnum'],tag=exp['tag'],
                          force='-f' if args.force else '',
                          verbose='-v' if args.verbose else '',
                          outfile=outfile,explist=explist)

            cmd = "calibrate.py %(force)s %(verbose)s -t %(tag)s -e %(expnum)s %(explist)s %(outfile)s"%params
            submit = 'csub -q %s -o %s '%(args.queue, logfile)
            if args.njobs: submit += '-n %s '%args.njobs
            submit += cmd
            if args.verbose: print submit
            subprocess.call(submit,shell=True)
            time.sleep(args.sleep)

            #break
        #break
        #time.sleep(300)

