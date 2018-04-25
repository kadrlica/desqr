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
import download

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
    zpfile = config['zpfile']
    blfile = config['blfile']
    section= config['db']
    tags   = config['tags']
    ebv    = config.get('ebv',None)
    
    if os.path.exists(zpfile) and not args.force:
        print "Found %s; skipping download..."%zpfile
    else:
        query = download.zeropoint_query(tags)
        print query
        sqlfile = os.path.splitext(zpfile)[0]+'.sql'
        download.download(zpfile,query,sqlfile=sqlfile,section=section,force=args.force)

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

            if not blfile:
                blacklist = ''
            elif isinstance(blfile,basestring):
                blacklist = '-b %s'%blfile
            elif isinstance(blfile,list):
                blacklist = '-b '+' -b '.join(blfile)
            
            params = dict(infile=infile,zpfile=zpfile,
                          blacklist = blacklist,
                          extinction = '-e %s'%ebv if ebv else '',
                          force = '-f' if args.force else '')
            
            cmd = 'zeropoint.py %(force)s %(infile)s %(zpfile)s %(blacklist)s'%params
            if False: #if ebv:
                cmd += "; extinction.py %(force)s %(infile)s %(extinction)s"%params

            if args.queue == 'local':
                print cmd
                submit = cmd
            else:
                #submit = 'csub -o %s "%s"'%(logfile,cmd)
                submit = 'csub -o %s -n %s "%s"'%(logfile,args.njobs,cmd)
            subprocess.call(submit,shell=True)
            time.sleep(args.sleep)
        if args.queue != 'local': time.sleep(5)
            
