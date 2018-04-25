#!/usr/bin/env python
import os
import yaml
import subprocess
import time
import glob

import numpy as np
from ugali.utils.shell import mkdir
from utils import found

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-s','--sleep',default=1,type=float)
    parser.add_argument('-q','--queue',default='condor')
    parser.add_argument('-v','--verbose',action='store_true')
    args = parser.parse_args()

    force = '-f' if args.force else ''
    config = yaml.load(open(args.config))
    catdir = config['catdir']

    skmdir = mkdir(config['skmdir'])
    columns = config['skimcol'] if 'skimcol' in config else config['skim']
    select = config.get('skimsel')

    logdir = mkdir(os.path.join(skmdir,'log'))

    infiles = sorted(glob.glob(catdir+'/*.fits'))

    for infile in infiles:
        basename = os.path.basename(infile)
        outfile = os.path.join(skmdir,basename)
        logfile = os.path.join(logdir,basename.replace('.fits','.log'))
        if os.path.exists(outfile) and not args.force:
            found(outfile)
            continue

        params = (infile,outfile,force,' -c '.join(columns))
        cmd = "skim.py %s %s %s -c %s"%params
        if args.queue == 'local':
            if select: cmd += ' -s "%s"'%select
            submit = cmd
        else:
            if select: cmd += r' -s "\\\"%s\\\""'%select
            submit = "csub -o %s %s"%(logfile,cmd)
        if args.verbose:
            print submit
        subprocess.call(submit,shell=True)
        time.sleep(args.sleep)
        #break
