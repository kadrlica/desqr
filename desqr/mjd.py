#!/usr/bin/env python
"""
Add MJD information to SE hpx files.
"""
__author__ = "Alex Drlica-Wagner"
import os
import shutil
import logging

import pandas as pd
import fitsio
import numpy as np

from utils import insert_columns

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infile',help='input file')
    parser.add_argument('outfile',nargs='?',help='output file')
    parser.add_argument('-m','--mjdfile',default='data/desdm_mjd.fits',
                        help='File containing MJD information')
    parser.add_argument('-f','--force',action='store_true',
                        help='overwrite output columns if they exist')
    parser.add_argument('-v','--verbose',action='store_true',
                        help='output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    outfile = args.outfile
    if args.outfile is None:
        outfile = args.infile

    if not os.path.exists(outfile):
        shutil.copy(args.infile,outfile)        

        
    logging.info("Loading %s..."%args.infile)
    expnum = pd.DataFrame(fitsio.read(args.infile,columns=['EXPNUM']).byteswap().newbyteorder())

    logging.info("Reading %s..."%args.mjdfile)
    mjd = pd.DataFrame(fitsio.read(args.mjdfile).byteswap().newbyteorder())

    merged = expnum.merge(mjd,left_on='EXPNUM',right_on='EXPNUM')
    assert (merged['EXPNUM'] == expnum['EXPNUM']).all()

    out = merged[['MJD_OBS','EXPTIME']].to_records(index=False)

    logging.info("Writing %s..."%outfile)
    insert_columns(outfile,out,ext=1,force=args.force)
