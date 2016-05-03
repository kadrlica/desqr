#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import multiprocessing

import fitsio
import numpy as np
import scipy.ndimage as nd
import pandas as pd

CCDS = np.arange(1,63)
COLUMNS = ['EXPNUM','CCDNUM','SPREAD_MODEL']

def bad_spread_model(filename,threshold=0.15):
    data = fitsio.read(filename,columns=COLUMNS)
    expnum = np.unique(data['EXPNUM'])
    if len(expnum) > 1:
        msg = 'Multiple exposures found in file: %s'%filename
        raise Exception(msg)
    expnum = np.asscalar(expnum)

    below = data['SPREAD_MODEL'] < -0.005
    total = np.ones(len(below))

    nbelow = nd.sum(below,data['CCDNUM'],CCDS)
    ntotal = nd.sum(total,data['CCDNUM'],CCDS)

    bad_ccds = CCDS[np.where(nbelow/ntotal > threshold)]

    dtype = [('expnum',int),('ccdnum',int)]

    out = np.recarray(len(bad_ccds),dtype=dtype)
    out['expnum'] = expnum
    out['ccdnum'] = bad_ccds
    return out

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infiles',nargs='+')
    parser.add_argument('-o','--outfile',default='ccd_blacklist.csv')
    parser.add_argument('--badfrac',default=0.15,type=float)
    args = parser.parse_args()
    
    data = None
    for filename in args.infiles:
        if data is None:
            data = bad_spread_model(filename,args.badfrac)
        else:
            data = np.append(data,bad_spread_model(filename,args.badfrac))
            
    pd.DataFrame(data).to_csv(args.outfile,index=False)
