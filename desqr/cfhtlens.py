#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import os.path
from os.path import splitext
import glob
import warnings
import logging

import healpy as hp
import numpy as np
import pandas as pd
import fitsio

import astropy.wcs
import pylab as plt

EDGE = 8192
DATADIR = '/data/des40.b/data/cfhtlens/w4/mask'

def ang2msk(lon,lat):
    maskfiles = sorted(glob.glob(DATADIR+'/W4*_izrgu_finalmask_mosaic.fits'))

    values = np.zeros(len(lon),dtype=int)
    index = np.arange(len(lon),dtype=int)
    warnings.simplefilter('ignore',category=astropy.wcs.FITSFixedWarning)

    idxs = []

    for f in maskfiles:
        logging.debug(f)
        mask = fitsio.read(f)
        
        wcs = astropy.wcs.WCS(f)

        crval = wcs.wcs.crval
        radius = max(wcs.wcs.cdelt*np.array([wcs._naxis1,wcs._naxis2]))

        xpix,ypix = wcs.wcs_world2pix(lon,lat,0)

        xedge = np.where(mask[wcs._naxis1/2,:] < EDGE)[0]
        yedge = np.where(mask[:,wcs._naxis2/2] < EDGE)[0]
        XMIN,XMAX = xedge.min(),xedge.max()
        YMIN,YMAX = yedge.min(),yedge.max()

        sel =  (xpix > XMIN) & (xpix < XMAX) 
        sel &= (ypix > YMIN) & (ypix < YMAX)

        xpix = xpix[sel].astype(int)
        ypix = ypix[sel].astype(int)
        idx = index[sel]
        idxs.append(idx)

        val = mask[ypix,xpix]
        values[idx] |= val

    # Set objects outside the footprint to EDGE
    idxs = np.hstack(idxs)
    idx = index[~np.in1d(index,idxs)]
    values[idx] = EDGE

    return values
