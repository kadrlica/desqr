#!/usr/bin/env python
"""
Utilities for applying zeropoints.
"""
import os,shutil
from os.path import splitext
import glob
from collections import OrderedDict as odict
import logging

import numpy as np
import numpy.lib.recfunctions as rfn
np.seterr(divide='ignore')
import pandas as pd
import fitsio

form desqr import utils
from desqr.logger import logger
from desqr.const import BADMAG,BADZP

DATACOLS = ['EXPNUM','CCDNUM','FLUX_PSF','FLUXERR_PSF','FLUX_AUTO','FLUXERR_AUTO']
ZPCOLS   = ['EXPNUM','CCDNUM','MAG_ZERO','SIGMA_MAG_ZERO','FLAG']
OUTCOLS  = ['MAG_PSF','MAGERR_PSF','MAG_AUTO','MAGERR_AUTO','MAG_ZERO','SIGMA_MAG_ZERO','ZP_FLAG']

def calc_mag(flux,zeropoint):
    """ Calculate magnitude from flux and zeropoint using Pogson's law. """
    mag = -2.5 * np.log10(flux) + zeropoint
    mag[~np.isfinite(mag)] = BADMAG
    mag[(zeropoint>90)|(zeropoint<0)] = BADMAG
    return mag

def calc_magerr(flux,fluxerr):
    """ Calculate magnitude error from flux and flux error. """
    magerr = 2.5/np.log(10) * (fluxerr/flux)
    magerr[~np.isfinite(magerr)] = BADMAG
    return magerr

def read_blacklist(blacklist):
    """ Read the list of excluded exposures. 

    Parameters
    ----------
    blacklist: filename or list of filenames

    Returns
    -------
    bl : recarray of expnum, ccdnum to exclude
    """
    if not blacklist: return None

    if utils.isstring(blacklist):
        ext = os.path.splitext(blacklist)[-1]
        if blacklist.endswith(('.csv','.csv.gz')):
            bl = pd.read_csv(blacklist).to_records(index=False)
            bl.dtype.names = map(str.upper,bl.dtype.names)
        elif blacklist.endswith(('.fits','.fz','.fits.gz')):
            bl = fitsio.read(blacklist,trim_strings=True)
        else:
            msg = "Unrecognized blacklist extension: %s"%ext
            raise Exception(msg)
    elif hasattr(blacklist,'__iter__'):
        bl = rfn.stack_arrays([read_blacklist(f) for f in blacklist],
                              usemask=False,asrecarray=True)
    else:
        msg = 'Unrecognized blacklist: %s'%blacklist
        raise Exception(msg)

    bl.dtype.names = list(map(str.upper,bl.dtype.names))
    return bl

def read_zeropoint(zeropoint,blacklist=None):
    """ Read the zeropoints and apply the blacklist. """

    if zeropoint.endswith(('.csv','.csv.gz')):
        zp = pd.read_csv(zeropoint).to_records(index=False)
        zp.dtype.names = map(str.upper,zp.dtype.names)
    elif zeropoint.endswith(('.fits','.fz','.fits.gz')):
        zp = fitsio.read(zeropoint,ext=1,columns=ZPCOLS)
    else:
        msg = "Unrecognized zeropoint: %s"%ext
        raise Exception(msg)

    # Rename FLAG to ZP_FLAG
    if 'FLAG' in zp.dtype.names:
        msg = "Renaming FLAG to ZP_FLAG..."
        logger.info(msg)
        names = ['ZP_FLAG' if n=='FLAG' else n for n in zp.dtype.names]
        zp.dtype.names = names

    # Everything assumes zeropoints are positive
    # If necessary, change the sign of MAG_ZERO
    if (zp['MAG_ZERO'] < 0).sum()/float(len(zp)) > 0.90:
        msg = "Switching sign of MAG_ZERO."
        logger.info(msg)
        zp['MAG_ZERO'] *= -1
        
    # reject CCDs with MAG_ZERO < 0
    sel_zero = (zp['MAG_ZERO'] > 0)
    if (~sel_zero).sum():
        msg = "Found %i CCDs with negative zeropoints."%(~sel_zero).sum()
        logger.warning(msg)
        
    # reject CCDs with FLAGS != 0
    sel_flags = (zp['ZP_FLAG'] == 0)
    if (~sel_flags).sum():
        msg = "Found %i CCDs with bad ZP flags."%(~sel_flags).sum()
        logger.warning(msg)

    sel = sel_zero & sel_flags
    zp = zp[sel]
    
    if blacklist is None:
        return zp
    else:
        logger.info("Excluding CCDs from %s..."%blacklist)
        bl = read_blacklist(blacklist)
        # Unique value from EXPNUM, CCDNUM
        zpval = utils.uid(zp['EXPNUM'],zp['CCDNUM'])
        blval = utils.uid(bl['EXPNUM'],bl['CCDNUM'])
        idx = np.in1d(zpval,blval)
        return zp[~idx]

def zeropoint(data,zp):
    """
    Match zeropoints to objects based on expnum,ccdnum.

    ADW: This could be better done with pandas.DataFrame.merge

    Parameters:
    -----------
    data : the input catalog data
    zp   : the input zeropoints (expnum,ccdnum)

    Returns:
    --------
    out  : catalog of calibrated magnitude and zeropoints
    """
    # ADW: This could be made more efficient with pandas
    out = np.recarray(len(data),dtype=[(n,'f4') for n in OUTCOLS])
    out.fill(BADMAG)
    out['ZP_FLAG'].fill(BADZP)

    expnums = np.unique(data['EXPNUM'])    
    logger.info("Applying zeropoints to %i exposures..."%len(expnums))

    for expnum in expnums:
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
        
        # This is trouble if there are CCDs with multiple ZPs...
        cidx = np.searchsorted(z['CCDNUM'],d['CCDNUM'])

        for x in ['PSF','AUTO']:
            out['MAG_'+x][didx] = calc_mag(d['FLUX_'+x],z['MAG_ZERO'][cidx])
            out['MAGERR_'+x][didx] = calc_magerr(d['FLUX_'+x],d['FLUXERR_'+x])
        out['MAG_ZERO'][didx] = z['MAG_ZERO'][cidx]
        out['SIGMA_MAG_ZERO'][didx] = z['SIGMA_MAG_ZERO'][cidx]
        out['ZP_FLAG'][didx] = z['ZP_FLAG'][cidx]

    return out

if __name__ == "__main__":
    import argparse
    description = "Apply zeropoints to an input file."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infile')
    parser.add_argument('zpfile')
    parser.add_argument('-b','--blacklist',default=None,action='append')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-v','--verbose',action='store_true')
    args = parser.parse_args()

    if args.verbose: logger.setLevel(logger.DEBUG)

    logger.info("Reading %s..."%args.infile)
    data = fitsio.read(args.infile,ext=1,columns=DATACOLS)
    zp = read_zeropoint(args.zpfile,args.blacklist)
    logger.info("Applying zeropoints from %s..."%args.zpfile)
    out = zeropoint(data,zp)

    # Writing...
    logger.info("Writing %s..."%args.infile)
    utils.insert_columns(args.infile,out,force=args.force)
    logger.info("Done.")
