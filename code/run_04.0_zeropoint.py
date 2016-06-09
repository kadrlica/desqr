#!/usr/bin/env python
"""
Apply zeropoints and extinction correction to catalog files.
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
    parser.add_argument('-s','--sleep',default=5,type=float)
    parser.add_argument('-q','--queue',default='condor')
    args = parser.parse_args()

    config = yaml.load(open(args.config))
    hpxdir = config['hpxdir']
    zpfile = config['zpfile']
    blfile = config['blfile']
    ebv    = config.get('ebv',None)

    for band in config['bands']:
        indir = os.path.join(hpxdir,band)
        infiles = sorted(glob.glob(indir+'/*.fits'))
        logdir = mkdir(os.path.join(indir,'log'))
        for infile in infiles:

            if not args.force:
                fits = fitsio.FITS(infile)
                colname = 'MAG_ZERO'
                if colname in fits[1].get_colnames():
                    print "Found column '%s'; skipping %s..."%(colname,os.path.basename(infile))
                    continue

            logbase = ('zp_'+os.path.basename(infile)).replace('.fits','.log')     
            logfile = os.path.join(logdir,logbase)
            params = dict(infile=infile,zpfile=zpfile,
                          blacklist = '-b %s'%blfile if blfile else '',
                          extinction = '-e %s'%ebv if ebv else '',
                          force = '-f' if args.force else '')
            
            cmd = 'zeropoint.py %(force)s %(infile)s %(zpfile)s %(blacklist)s'%params
            if ebv:
                cmd += "; extinction.py %(force)s %(infile)s %(extinction)s"%params

            if args.queue == 'local':
                print cmd
                submit = cmd
            else:
                submit = 'csub -o %s "%s"'%(logfile,cmd)
            subprocess.call(submit,shell=True)
            time.sleep(args.sleep)
        if args.queue != 'local': time.sleep(5)
            