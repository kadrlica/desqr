#!/usr/bin/env python
"""
Append MJD and EXPTIME to SE hpx files.

This is a temporary script and eventually these values should be included with the original download.
"""
import os
import yaml
import subprocess
import time
import glob

import fitsio

from ugali.utils.shell import mkdir

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-s','--sleep',default=0,type=float)
    parser.add_argument('-n','--njobs',default=15,type=int)
    parser.add_argument('-q','--queue',default='condor')
    args = parser.parse_args()

    config = yaml.load(open(args.config))
    hpxdir = config['hpxdir']
    mjdfile = config.get('mjdfile','data/desdm_mjd.fits')

    for band in config['bands']:
        indir = os.path.join(hpxdir,band)
        infiles = sorted(glob.glob(indir+'/*.fits'))
        logdir = mkdir(os.path.join(indir,'log'))
        for infile in infiles:
            if not args.force:
                fits = fitsio.FITS(infile)
                colname = 'MJD_OBS'
                if colname in fits[1].get_colnames():
                    print "Found column '%s'; skipping %s..."%(colname,os.path.basename(infile))
                    continue

            logbase = ('mjd_'+os.path.basename(infile)).replace('.fits','.log')
            logfile = os.path.join(logdir,logbase)

            params = dict(infile=infile,mjdfile='-m %s'%mjdfile,
                          force = '-f' if args.force else '')
            
            cmd = 'mjd.py %(force)s %(mjdfile)s %(infile)s '%params

            if args.queue == 'local':
                print cmd
                submit = cmd
            else:
                #submit = 'csub -o %s "%s"'%(logfile,cmd)
                submit = 'csub -o %s -n %s "%s"'%(logfile,args.njobs,cmd)
            subprocess.call(submit,shell=True)
            time.sleep(args.sleep)
        if args.queue != 'local': time.sleep(5)
            
