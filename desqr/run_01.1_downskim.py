#!/usr/bin/env python
""" Downskim exposures to disk """
import os
import yaml
import subprocess
import time
import glob

import numpy as np
from utils import mkdir
import download

if __name__ == "__main__":
    from parser import Parser
    parser = Parser(description=__doc__)
    args = parser.parse_args()

    config = yaml.load(open(args.config))
    rawdir = config['rawdir']
    explist = config['explist']
    tags = config.get('tags')
    section = config.get('db','bliss')
    basename = config.get('rawbase','D{expnum:08d}_{band:s}_cat.fits')
    
    if os.path.exists(explist) and not args.force:
        print("Found %s; skipping download..."%explist)
    else:
        print("Didn't find %s; downloading..."%explist)
        query = download.exposure_query(tags)
        print query
        sqlfile = os.path.splitext(explist)[0]+'.sql'
        download.download(explist,query,sqlfile=sqlfile,section=section,force=args.force)
    
    exposures = np.recfromcsv(explist)
    for band in config['bands']:
        outdir = mkdir(os.path.join(rawdir,band))
        logdir = mkdir(os.path.join(outdir,'log'))
        for exp in exposures[exposures['band'] == band]:
            if (tags is not None) and (exp['tag'] not in tags):
                print("No exposures with tag '%s'"%tags)
                continue

            #outbase = basename%(exp)
            edict = dict(zip(exp.dtype.names,exp))
            outbase = basename.format(**edict)
            outfile = os.path.join(outdir,outbase)
            if os.path.exists(outfile) and not args.force:
                print('Found %s; skipping...'%outfile)
                continue

            logfile = os.path.join(logdir,outbase.replace('.fits','.log'))
            params = dict(expnum=exp['expnum'],tag=exp['tag'],
                          force='-f' if args.force else '',
                          verbose='-v' if args.verbose else '',
                          outfile=outfile,explist=explist)
                      
            cmd = "downskim.py %(force)s %(verbose)s -t %(tag)s -e %(expnum)s %(explist)s %(outfile)s"%params
            submit = 'csub -q %s -o %s '%(args.queue, logfile)
            if args.njobs: submit += '-n %s '%args.njobs
            submit += cmd
            if args.verbose: print(submit)
            if not args.dryrun: subprocess.call(submit,shell=True)
            time.sleep(args.sleep)
        
