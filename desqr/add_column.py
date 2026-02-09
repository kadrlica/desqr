#!/usr/bin/env python
"""
Add a column to a fits file based on a formula.
"""
__author__ = "Alex Drlica-Wagner"

import numpy as np
import fitsio

from desqr import utils
from desqr.utils import add_column, multiproc
from desqr.logger import logger

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('filenames',nargs='+')
    parser.add_argument('-c','--column',default=None)
    parser.add_argument('--formula',default=None)
    parser.add_argument('-f','--force',action='store_true',
                        help='overwrite output columns if they exist')
    parser.add_argument('-n','--nproc', default=12, type=int,
                        help='number of processes to run')
    parser.add_argument('-v','--verbose',action='store_true',
                        help='output verbosity')
    args = parser.parse_args()

    if args.verbose: logger.setLevel(logger.DEBUG)

    print(args.column)
    print(utils.parse_formula(args.formula))

    func = utils.add_column
    arglist = [(f,args.column,args.formula,args.force) for f in args.filenames]
    results = utils.multiproc(func, arglist processes=args.nproc)

