#!/usr/bin/env python
"""
Apply extinction correction to catalog files.
"""
import os
from os.path import basename
import yaml
import subprocess
import time
import glob

import fitsio

from ugali.utils.shell import mkdir

if __name__ == "__main__":
    from parser import Parser
    parser = Parser(description=__doc__)
    args = parser.parse_args()

    config = yaml.load(open(args.config))
    catdir = config['catdir']
    ebv    = config.get('ebv',None)
    bands  = config['bands']

    infiles = sorted(glob.glob(catdir+'/*.fits'))
    logdir = mkdir(os.path.join(catdir,'log'))

    for infile in infiles:
        if not args.force:
            fits = fitsio.FITS(infile)
            if 'EBV' in fits[1].get_colnames():
                print "Found column 'EBV'; skipping %s..."%(basename(infile))
                continue

        logbase = ('ebv_'+os.path.basename(infile)).replace('.fits','.log')     
        logfile = os.path.join(logdir,logbase)
        params = dict(infile=infile,force='-f' if args.force else '',
                      extinction = '-e %s'%ebv if ebv else '',
                      bands = '-b '+' -b '.join(bands))

        cmd = "extinction.py %(force)s %(extinction)s %(bands)s %(infile)s"%params

        submit = 'csub -o %s -n %s -q %s "%s"'%(logfile,args.njobs,args.queue,cmd)
        subprocess.call(submit,shell=True)
        time.sleep(args.sleep)
