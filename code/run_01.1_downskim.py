#!/usr/bin/env python
import os
import yaml
import subprocess
import time
import glob

import numpy as np
from utils import mkdir

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-s','--sleep',default=2,type=float)
    parser.add_argument('-n','--njobs',default=24,type=int)
    parser.add_argument('-q','--queue',default='condor')
    parser.add_argument('-v','--verbose',action='store_true')
    args = parser.parse_args()

    config = yaml.load(open(args.config))
    rawdir = config['rawdir']
    explist = config['explist']
    tags = config.get('tags')
    
    if os.path.exists(explist) and not args.force:
        print "Found %s; skipping download..."%explist
    else:
        msg = "Missing exposure list: %s"%explist
        raise Exception(msg)

    exposures = np.recfromcsv(explist)
    for band in config['bands']:
        outdir = mkdir(os.path.join(rawdir,band))
        logdir = mkdir(os.path.join(outdir,'log'))
        for exp in exposures[exposures['band'] == band]:
            if (tags is not None) and (exp['tag'] not in tags):
                continue
            unitname = 'D%08d'%(exp['expnum'])
            outbase = '%s_%s_cat.fits'%(unitname,band)
            outfile = os.path.join(outdir,outbase)
            if os.path.exists(outfile) and not args.force:
                print 'Found %s; skipping...'%outfile
                continue

            logfile = os.path.join(logdir,outbase.replace('.fits','.log'))
            params = dict(expnum=exp['expnum'],tag=exp['tag'],
                          force='-f' if args.force else '',
                          verbose='-v' if args.verbose else '',
                          outfile=outfile,explist=explist)
                      
            cmd = "downskim.py %(force)s %(verbose)s -t %(tag)s -e %(expnum)s %(explist)s %(outfile)s"%params
            if args.queue == 'local':
                submit = cmd
            else:
                submit = 'csub -o %s '%(logfile)
                if args.njobs: submit += '-n %s '%args.njobs
                submit += cmd
            if args.verbose: print submit
            subprocess.call(submit,shell=True)
            time.sleep(args.sleep)

            #break
        #break
        #time.sleep(300)
        
