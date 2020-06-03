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
    from parser import Parser
    parser = Parser(description=__doc__)
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
        dirname = os.path.join(hpxdir,band)
        filenames = sorted(glob.glob(dirname+'/*.fits'))
        logdir = mkdir(os.path.join(dirname,'log'))
        for filename in filenames:
            basename = os.path.basename(filename)

            if not args.force:
                fits = fitsio.FITS(filename)
                colname = 'MAG_ZERO'
                if colname in fits[1].get_colnames():
                    print "Found column '%s'; skipping %s..."%(colname,basename)
                    continue

            logbase = ('zp_'+basename).replace('.fits','.log')     
            logfile = os.path.join(logdir,logbase)

            if not blfile:
                blacklist = ''
            elif isinstance(blfile,basestring):
                blacklist = '-b %s'%blfile
            elif isinstance(blfile,list):
                blacklist = '-b '+' -b '.join(blfile)
            
            params = dict(infile=filename,zpfile=zpfile,
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
            
