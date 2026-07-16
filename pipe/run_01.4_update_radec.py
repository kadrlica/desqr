#!/usr/bin/env python
"""
Update RA, DEC from updated WCS.
"""
__author__ = "Alex Drlica-Wagner"
import os
import time
import subprocess
import yaml

from desqr.utils import mkdir

if __name__ == "__main__":
    from desqr.parser import Parser
    parser = Parser(description=__doc__)
    args = parser.parse_args()

    config = args.config
    dirname = config['rawdir']
    #dirname = 'tmp'

    for band in config['bands']:
        if args.bands and (band not in args.bands): continue
        print("Running %s-band..."%band)
        
        logdir = mkdir(os.path.join(dirname, 'log'))
        logfile = os.path.join(logdir,'update_radec_%s.log'%band)

        # Argument list too long so use xargs...
        cmd  = 'ls %s/%s/D*.fits | xargs -n 1000 update_radec.py -f'%(dirname, band)
        cmd += ' -v' if args.verbose else ''

        #submit = 'csub -q %s -o %s %s'%(args.queue,logfile,cmd)
        submit = cmd
        if args.verbose: print(submit)
        if not args.dryrun: subprocess.call(submit, shell=True)
        time.sleep(args.sleep)
    print("Done.")
