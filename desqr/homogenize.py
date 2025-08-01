#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import os,sys
from collections import OrderedDict as odict

import pandas as pd
import numpy as np
import pylab as plt
import healpy as hp

BANDS = ['g','r','i','z']
NSIDE = 1024
DECAM = 1.1 #deg

def homogenize(data, percent=90, sum_thresh=10*0.3*90, num_thresh=10):
    """ Select a subset of exposures to yield uniform depth (in
    teff*texp and number of exposures) over the sky.

    Parameters
    ----------
    data    : array of exposure quantities
    percent : fraction overlap for calculating thresholds
    sum_thresh : threshold in the summed depth
    num_thresh : threshold in the number of exposures    

    Returns
    -------
    sum_exp,num_exp : exposures selected on depth and number
    """

    # First, sort the exposures by whatever quantity you want
    sortby = data['t_eff']*data['exptime']
    idx = np.argsort(sortby)[::-1]
    data = data[idx]

    # Break exposures into bands
    exposures = odict([(b,data[data['band']==b]) for b in BANDS])

    # Select exposure to give uniform depth
    sel_sum = odict([(b,np.zeros(len(e),dtype=bool)) for b,e in exposures.items()])
    # Select exposure to give uniform number of exposures
    sel_num = odict([(b,np.zeros(len(e),dtype=bool)) for b,e in exposures.items()])
    
    # Sum depth (tef*texp) in each hpx
    sum_skymaps = odict([(b,np.zeros(hp.nside2npix(NSIDE))) for b in BANDS])
    # Maximum depth (teff*texp) in each hpx
    max_skymaps = odict([(b,np.zeros(hp.nside2npix(NSIDE))) for b in BANDS])
    # Number of exposures overlapping each hpx
    num_skymaps = odict([(b,np.zeros(hp.nside2npix(NSIDE))) for b in BANDS])

    pixarea = hp.nside2pixarea(NSIDE,degrees=True)
    rad = np.radians(DECAM)

    for band in BANDS:
        print band
        _sum = sum_skymaps[band]
        _max = max_skymaps[band]
        _num = num_skymaps[band]

        exps = exposures[band]

        vec = hp.ang2vec(exps['ra'],exps['dec'],lonlat=True)
        # Loop through exposures
        for i,(v,d) in enumerate(zip(vec,exps)):
            print '\r%s/%s'%(i+1,len(vec)),
            sys.stdout.flush()
            pix = hp.query_disc(NSIDE,v,rad,inclusive=False,fact=4,nest=False)
            
            # Find effective exposure time that is achieved over
            # the specified fraction of the exposure e.g., where
            # more than `percent` of the pixels have a larger SUM(teff * exptime)
            depth = np.percentile(_sum[pix],100-percent)
            if depth < sum_thresh: 
                sel_sum[band][i] = True
                _sum[pix] += d['t_eff']*d['exptime']

            count = np.percentile(_num[pix],100-percent)
            if count < num_thresh:
                sel_num[band][i] = True
                _num[pix] += 1

            _max[pix]  = np.clip(_max[pix],d['t_eff']*d['exptime'],None)

        _sum[_sum==0] = np.nan
        _max[_max==0] = np.nan
        _num[_num==0] = np.nan
        print

    sum_exps = odict([(b,exposures[b][sel]) for b,sel in sel_sum.items()])
    num_exps = odict([(b,exposures[b][sel]) for b,sel in sel_num.items()])
    return sum_exps,num_exps

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('filename',default='delve_wide_r2_v2.csv')
    parser.add_argument('--num',default=10,type=int,
                        help="threshold on overlapping exposures")
    parser.add_argument('--sum',default=270,type=float,
                        help="threshold on total effective exposure time (s)")
    parser.add_argument('--percent',default=90,type=float,
                        help="exposure coverage percent")
    args = parser.parse_args()

    data = pd.read_csv(args.filename).to_records(index=False)
    names = [str(n) for n in data.dtype.names]
    for i,n in enumerate(names): 
        if n == 'telra': names[i] = 'ra'
        if n == 'teldec': names[i] = 'dec'
    data.dtype.names = names

    # Select a homogenious set of exposures
    sum_thresh = args.sum
    num_thresh = args.num
    percent    = args.percent
    print("sum_thresh: %g"%sum_thresh)
    print("num_thresh: %d"%num_thresh)
    sum_exps,num_exps = homogenize(data,percent=percent,
                                   sum_thresh=sum_thresh,num_thresh=num_thresh)

    sum_out = pd.concat(map(pd.DataFrame,sum_exps.values())).sort_values('expnum')
    num_out = pd.concat(map(pd.DataFrame,num_exps.values())).sort_values('expnum')

    outbase = os.path.splitext(os.path.basename(args.filename))[0]

    outfile = outbase + '_sum.csv'
    print("Writing %s..."%outfile)
    sum_out.to_csv(outfile,index=False)

    outfile = outbase + '_num.csv'
    print("Writing %s..."%outfile)
    num_out.to_csv(outfile,index=False)
