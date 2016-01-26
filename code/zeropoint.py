#!/usr/bin/env python
"""
Utilities for calculating magnitudes from zeropoints.
"""
import os,shutil
from os.path import splitext
import glob
import numpy as np
np.seterr(divide='ignore')
import fitsio

from ugali.utils.logger import logger
import utils
from const import BADMAG,BADQSLR

DATACOLS = ['EXPNUM','CCDNUM','FLUX_PSF','FLUXERR_PSF','FLUX_AUTO','FLUXERR_AUTO']
ZPCOLS   = ['EXPNUM','CCDNUM','MAG_ZERO','QSLR_FLAG']
OUTCOLS  = ['MAG_PSF','MAGERR_PSF','MAG_AUTO','MAGERR_AUTO','MAG_ZERO','QSLR_FLAG']
 
def calc_mag(flux,zeropoint):
    mag = -2.5 * np.log10(flux) + zeropoint
    mag[~np.isfinite(mag)] = BADMAG
    return mag

def calc_magerr(flux,fluxerr):
    magerr = 2.5/np.log(10) * (fluxerr/flux)
    magerr[~np.isfinite(magerr)] = BADMAG
    return magerr

def read_zeropoint(zeropoint,blacklist=None):
    """ Read the zeropoints and apply the blacklist. """

    ext = os.path.splitext(zeropoint)[-1]
    if ext == '.csv':
        zp = np.recfromcsv(zeropoint)
    elif ext == '.fits':
        zp = fitsio.read(zeropoint,ext=1,columns=ZPCOLS)
    else:
        msg = "Unrecognized zeropoint extension: %s"%ext
        raise Exception(msg)
        
    # reject CCDs with MAG_ZERO < 0
    sel_zero = (zp['MAG_ZERO'] >= 0)
    if (~sel_zero).sum():
        msg = "Found %i CCDs with negative zeropoints."%(~sel_zero).sum()
        logger.warning(msg)
    # reject CCDs with QSLR_FLAGS > 0

    sel_flags = (zp['QSLR_FLAG'] == 0)
    if (~sel_flags).sum():
        msg = "Found %i CCDs with bad qSLR flags."%(~sel_flags).sum()
        logger.warning(msg)

    sel = sel_zero & sel_flags
    zp = zp[sel]
    
    if blacklist is None:
        return zp
    else:
        ext = os.path.splitext(blacklist)[-1]
        if ext == '.csv':
            bl = np.recfromcsv(blacklist)
        elif ext == '.fits':
            bl = fitsio.read(blacklist)
        else:
            msg = "Unrecognized blacklist extension: %s"%ext
            raise Exception(msg)

        # A hack to do the sorting
        zpval = zp['EXPNUM'] + zp['CCDNUM']/100.
        blval = bl['expnum'] + bl['ccdnum']/100.
        idx = np.in1d(zpval,blval)
        return zp[~idx]

def zeropoint(data,zp):
    out = np.recarray(len(data),dtype=[(n,float) for n in OUTCOLS])
    out.fill(BADMAG)
    out['QSLR_FLAG'].fill(BADQSLR)

    for expnum in np.unique(data['EXPNUM']):
        zidx = np.where(zp['EXPNUM'] == expnum)[0]
        if zidx.sum() == 0:
            logger.debug("No zeropoint for exposure: %i"%expnum)
            continue
        
        z = zp[zidx]
        z = z[np.argsort(z['CCDNUM'])]

        # There are some CCDs with no ZP
        didx = np.where(data['EXPNUM'] == expnum)[0]
        d = data[didx]
        inccd = np.in1d(d['CCDNUM'],np.unique(z['CCDNUM']))

        didx =didx[inccd]
        d = d[inccd]
        
        cidx = np.searchsorted(z['CCDNUM'],d['CCDNUM'])

        for x in ['PSF','AUTO']:
            out['MAG_'+x][didx] = calc_mag(d['FLUX_'+x],z['MAG_ZERO'][cidx])
            out['MAGERR_'+x][didx] = calc_magerr(d['FLUX_'+x],d['FLUXERR_'+x])
        out['QSLR_FLAG'][didx] = z['QSLR_FLAG'][cidx]
        out['MAG_ZERO'][didx] = z['MAG_ZERO'][cidx]
      
    return out

if __name__ == "__main__":
    import argparse
    description = "Apply zeropoints to an input file."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infile')
    parser.add_argument('zpfile')
    parser.add_argument('-b','--blacklist',default=None)
    parser.add_argument('-f','--force',action='store_true')
    opts = parser.parse_args()

    logger.info("Reading %s..."%opts.infile)
    data = fitsio.read(opts.infile,ext=1,columns=DATACOLS)
    zp = read_zeropoint(opts.zpfile,opts.blacklist)
    logger.info("Applying zeropoints %s..."%opts.zpfile)
    out = zeropoint(data,zp)

    # Writing...
    logger.info("Writing %s..."%opts.infile)
    utils.insert_columns(opts.infile,out,force=opts.force)
