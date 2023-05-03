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

import download
from utils import isstring, mkdir

if __name__ == "__main__":
    from parser import Parser
    parser = Parser(description=__doc__)
    args = parser.parse_args()

    config  = yaml.safe_load(open(args.config))
    hpxdir  = config['hpxdir']
    zpfile  = config['zpfile']
    blfile  = config['blfile']
    ebv     = config.get('ebv',None)
    
    if not os.path.exists(zpfile):
        msg = "Couldn't find ZP file: %s"%zpfile
        raise Exception(msg)

    #query = download.zeropoint_query(config['tags'])
    #print(query)
    #sqlfile = os.path.splitext(zpfile)[0]+'.sql'
    #download.download(zpfile,query,sqlfile=sqlfile,section=config['db'],force=args.force)

    for band in config['bands']:
        dirname = os.path.join(hpxdir,band)
        logdir = mkdir(os.path.join(dirname,'log'))

        filenames = sorted(glob.glob(dirname+'/*.fits'))

        if args.pix:
            print("Running pixels: %s"%args.pix)
            pathname = os.path.join(dirname,config['hpxbase'].format(band=band))
            filenames = [pathname%p for p in args.pix]
            
        for filename in filenames:
            if not os.path.exists(filename):
                print("WARNING: File does not exist; skipping %s..."%filename)
                continue

            basename = os.path.basename(filename)

            if not args.force:
                fits = fitsio.FITS(filename)
                colname = 'MAG_ZERO'
                if colname in fits[1].get_colnames():
                    print("Found column '%s'; skipping %s..."%(colname,basename))
                    continue

            logbase = ('zp_'+basename).replace('.fits','.log')     
            logfile = os.path.join(logdir,logbase)

            if not blfile:
                blacklist = ''
            elif isstring(blfile):
                blacklist = '-b %s'%blfile
            elif isinstance(blfile,list):
                blacklist = '-b '+' -b '.join(blfile)
            
            params = dict(infile=filename,zpfile=zpfile,
                          blacklist = blacklist,
                          force = '-f' if args.force else '')
            
            cmd = 'zeropoint.py %(force)s %(infile)s %(zpfile)s %(blacklist)s'%params

            # Apply extinction to hpx files
            #if False: 
            #    params['extinction'] = '-e %s'%ebv if ebv else ''
            #    cmd += "; extinction.py %(force)s %(infile)s %(extinction)s"%params

            if args.queue == 'local': submit = cmd
            else: submit = 'csub -o %s -n %s "%s"'%(logfile,args.njobs,cmd)
                

            if args.verbose: print(submit)
                
            subprocess.call(submit,shell=True)
            time.sleep(args.sleep)
        if args.queue != 'local': time.sleep(5)
            
