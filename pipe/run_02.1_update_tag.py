#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import os
import time
import subprocess
import yaml

from utils import mkdir

COLUMN='TAG'

if __name__ == "__main__":
    from parser import Parser
    parser = Parser(description=__doc__)
    parser.set_defaults(njobs=12, queue='local')
    args = parser.parse_args()

    config = args.config
    #hpxdir = config['hpxdir']
    dirname = config['rawdir']
    #dirname = 'tmp'
    explist = config['explist']
    
    for band in config['bands']:
        if args.bands and (band not in args.bands): continue
            
        logdir = mkdir(os.path.join(dirname, 'log'))
        logfile = os.path.join(logdir,'update_tag_%s.log'%band)

        # Argument list too long so use xargs...
        #cmd  = ' merge_column.py %s/%s/*.fits -m %s'%(dirname, band, explist)
        cmd  = 'ls %s/%s/*.fits | xargs -n 1000 merge_column.py -m %s'%(dirname, band, explist)
        cmd += ' -c %s -f --nproc %s'%(COLUMN, args.njobs)
        cmd += ' -v' if args.verbose else ''

        #submit = 'csub -q %s -o %s %s'%(args.queue,logfile,cmd)
        submit = cmd
        if args.verbose: print(submit)
        if not args.dryrun: subprocess.call(submit, shell=True)
        time.sleep(args.sleep)
