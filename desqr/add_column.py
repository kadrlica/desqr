#!/usr/bin/env python
"""
Add a column to a fits file.
"""
__author__ = "Alex Drlica-Wagner"

import logging
from multiprocessing import Pool

import numpy as np
import fitsio

import utils
from utils import insert_columns, add_column


#def parse_formula(formula):
#    """ Figure out which columns to read. """
#    if not (('data[' in formula) and (']' in formula)):
#        msg = 'invalid formula'
#        raise Exception(msg)
# 
#    columns = []
#    for x in formula.split('data[')[1:]:
#        columns += x.split(']')[0].replace('"','').replace("'",'')
# 
#    return columns
#    
#def add_column(filename,column,formula,force=False):
#    columns = parse_formula(formula)
#    print("Running file: %s"%filename)
#    if args.verbose: print("  Reading columns: %s"%columns)
#    data = fitsio.read(filename,columns=columns)
# 
#    if args.verbose: print('  Evaluating formula: %s'%formula)
#    col = eval(formula)
# 
#    col = np.asarray(col,dtype=[(column,col.dtype)])
#    insert_columns(filename,col,force=force)
#    return True

def wrapper(args):
    return utils.add_column(*args)

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infiles',nargs='+')
    parser.add_argument('--column',default=None)
    parser.add_argument('--formula',default=None)
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-v','--verbose',action='store_true')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.getLogger().setLevel(level)

    print args.column
    print utils.parse_formula(args.formula)
    arguments = [(f,args.column,args.formula,args.force) for f in args.infiles]
    p = Pool(maxtasksperchild=1)
    out = p.map(wrapper,arguments)

