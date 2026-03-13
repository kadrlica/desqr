#!/usr/bin/env python
"""
Update the tag value in the raw catalog skims.
"""
__author__ = "Alex Drlica-Wagner"
import os
import time
import subprocess
import yaml

import numpy as np
import pandas as pd
import fitsio

from desqr import utils
from desqr.utils import mkdir, multiproc, insert_columns
from desqr.logger import logger

def update_tag(row):
    filename = f"{dirname}/{row.band}/D{row.expnum:08d}_{row.band}_cat.fits"
    if not os.path.exists(filename):
        logger.warning(f"File not found: {filename}")
        
    cat = fitsio.read(filename, columns=['TAG'])
    cat['TAG'][:] = row.tag
    insert_columns(filename, cat, force=True)

if __name__ == "__main__":
    from desqr.parser import Parser
    parser = Parser(description=__doc__)
    parser.add_argument('--nproc', default=12, type=int)
    args = parser.parse_args()
    
    config = args.config
    dirname = config['rawdir']
    explist = config['explist']
    df = pd.read_csv(explist)

    if args.bands is None: args.bands = config['bands']
    df = df[np.in1d(df['band'], args.bands)].reset_index(drop=True)
    
    func = update_tag
    arglist = [(row,) for idx,row in df.iterrows()]
    results = utils.multiproc(func, arglist, processes=args.nproc)

    # If we want to do this on files containing multiple exposures (i.e., hpx)
    #for band in config['bands']:
    #    if args.bands and (band not in args.bands): continue
    #        
    #    logdir = mkdir(os.path.join(dirname, 'log'))
    #    logfile = os.path.join(logdir,'update_tag_%s.log'%band)
    # 
    #    # Argument list too long so use xargs...
    #    #cmd  = ' merge_column.py %s/%s/*.fits -m %s'%(dirname, band, explist)
    #    #cmd  = 'ls %s/%s/*.fits | xargs -n 1000 merge_column.py -m %s'%(dirname, band, explist)
    #    #cmd += ' -c %s -f --nproc %s'%(COLUMN, args.njobs)
    #    #cmd += ' -v' if args.verbose else ''
    # 
    #    #submit = 'csub -q %s -o %s %s'%(args.queue,logfile,cmd)
    #    submit = cmd
    #    if args.verbose: print(submit)
    #    if not args.dryrun: subprocess.call(submit, shell=True)
    #    time.sleep(args.sleep)
