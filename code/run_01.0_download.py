#!/usr/bin/env python
import os
import yaml
import subprocess
import time
import glob

import numpy as np
import download
from utils import mkdir

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-s','--sleep',default=2,type=float)
    parser.add_argument('-q','--queue',default='condor')
    parser.add_argument('-v','--verbose',action='store_true')
    args = parser.parse_args()

    config = yaml.load(open(args.config))
    rawdir = config['rawdir']
    explist = config['explist']
    tags = config.get('tags')
    
    if os.path.exists(explist) and not args.force:
        print "Found %s; skipping ..."%explist
    else:
        query = download.exposure_query(tags,program='survey')
        print query
        sqlfile = os.path.splitext(explist)[0]+'.sql'
        download.download(explist,query,sqlfile=sqlfile,force=args.force)

    exposures = np.recfromcsv(explist)
    for band in config['bands']:
        outdir = os.path.join(rawdir,band)
        logdir = mkdir(os.path.join(outdir,'log'))
        sqldir = mkdir(os.path.join(outdir,'sql'))
        for exp in exposures[exposures['band'] == band]:
            if (tags is not None) and (exp['tag'] not in tags):
                continue
            unitname = exp['unitname']
            outbase = '%s_%s_cat.fits'%(unitname,band)
            outfile = os.path.join(outdir,outbase)
            if os.path.exists(outfile) and not args.force:
                print 'Found %s; skipping...'%outfile
                continue
            sqlfile = os.path.join(sqldir,outbase.replace('.fits','.sql'))
            logfile = os.path.join(logdir,outbase.replace('.fits','.log'))
            params = (exp['tag'],exp['expnum'],exp['reqnum'],exp['attnum'],
                      '-f' if args.force else '',sqlfile,outfile)
            cmd = "download.py -t %s -e %s -r %s -a %s %s -l %s %s"%params
            if args.queue == 'local':
                submit = cmd
            else:
                submit = "csub -o %s %s"%(logfile,cmd)
            if args.verbose: print cmd
            subprocess.call(submit,shell=True)
            time.sleep(args.sleep)

            #break
        #break
        #time.sleep(300)
        
