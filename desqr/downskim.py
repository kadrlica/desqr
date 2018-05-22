#!/usr/bin/env python
"""
Link the raw catalog files
"""
__author__ = "Alex Drlica-Wagner"

import os
import logging
import subprocess
import glob

import yaml
import pandas as pd
import numpy as np
import fitsio

from utils import mkdir
from const import BANDS, TAGS

#from archive.dexp import ObjectsTable

#bliss_temp = '/export/data/des50.b/data/BLISS/{dirname}/{expnum}/D00{expnum}_*_fullcat.fits'
bliss_des50 = '/data/des50.b/data/BLISS/{dirname}/{expnum}/D00{expnum}_*_fullcat.fits'
bliss_des60 = '/data/des60.b/data/BLISS/{dirname}/{expnum}/D00{expnum}_*_fullcat.fits'

outtemp = '%(unitname)s_%(band)s.fits'

DTYPES = [('FILENAME', 'S48'), ('PFW_ATTEMPT_ID', '>i8'), ('TAG', 'S13'), ('UNITNAME', 'S9'), ('REQNUM', '>i4'), ('ATTNUM', '>i2'), ('EXPNUM', '>i8'), ('CCDNUM', '>i2'), ('BAND', 'S1'), ('T_EFF', '>f4'), ('FWHM_WORLD', '>f4'), ('FLAGS', '>i2'), ('OBJECT_NUMBER', '>i4'), ('RA', '>f8'), ('DEC', '>f8'), ('FLUX_PSF', '>f4'), ('FLUXERR_PSF', '>f4'), ('FLUX_AUTO', '>f4'), ('FLUXERR_AUTO', '>f4'), ('A_IMAGE', '>f4'), ('B_IMAGE', '>f4'), ('THETA_IMAGE', '>f4'), ('CLASS_STAR', '>f4'), ('SPREAD_MODEL', '>f4'), ('SPREADERR_MODEL', '>f4'),('IMAFLAGS_ISO','>i2'),('MJD_OBS','>f8'),('EXPTIME','>f4')]

def get_dirname(expnum):
    return expnum//100 * 100

def bliss_selection(data):
    sel = data['FLUX_PSF'] > 0
    sel &= data['FLUX_AUTO'] > 0
    sel &= data['FLAGS'] < 4
    sel &= (data['IMAFLAGS_ISO'] & 2047) == 0 
    olderr = np.seterr(invalid='ignore',divide='ignore')
    sel &= 1.0857362*(data['FLUXERR_PSF']/data['FLUX_PSF']) < 0.5
    np.seterr(**olderr)
    return sel

maglites_selection = bliss_selection
    
def create_output(data,exp):
    out = np.recarray(len(data),dtype=DTYPES)
    for n in out.dtype.names:
        if n not in exp.dtype.names + data.dtype.names:
            msg = "Column not found: %s"%n
            raise Exception(msg)

        if n in exp.dtype.names:
            out[n] = exp[n]

        if n in data.dtype.names:
            out[n] = data[n]

    return out

def create_catalog(filename):
    fits = fitsio.FITS(filename)

    hdr = create_header(fits)
    data = fits['LDAC_OBJECTS'].read()

    cat = np.recarray(len(data),dtype=DTYPES)

    for name,dtype in DTYPES:
        if name not in data.dtype.names: continue
        cat[name] = data[name]
            
    cat['FILENAME'] = os.path.basename(filename)
    cat['PFW_ATTEMPT_ID'] = int(str(hdr['EXPNUM'])+'01')
    cat['TAG'] = TAG
    cat['UNITNAME'] = 'D%08d'%hdr['EXPNUM']
    cat['REQNUM'] = 1
    cat['ATTNUM'] = 1
    cat['EXPNUM'] = hdr['EXPNUM']
    cat['CCDNUM'] = hdr['CCDNUM']
    cat['BAND'] = hdr['BAND']
    cat['EXPTIME'] = hdr['EXPTIME']
    cat['MJD_OBS'] = hdr['MJD-OBS']
    cat['T_EFF'] = -1
    cat['RA'] = data['ALPHAWIN_J2000']
    cat['DEC'] = data['DELTAWIN_J2000']
    cat['OBJECT_NUMBER'] = data['NUMBER']
    cat['BAND'] = hdr['FILTER']
    return cat

def create_header(fits):
    from astropy.io.fits import Header

    if isinstance(fits,basestring):
        fits = fitsio.FITS(fits)

    hdrstr = '\n'.join(fits['LDAC_IMHEAD'].read()[0][0])
    return Header.fromstring(hdrstr,sep='\n')

def downskim(outfile,select,exp,force):
    """ Like download, but for skimming catalog files.
    """
    if os.path.exists(outfile) and not force:
        logging.warn('Found %s; skipping...'%outfile)
        return

    expnum = exp['EXPNUM']
    dirname = get_dirname(expnum)
    params = dict(expnum=expnum,dirname=dirname)
    
    filenames = glob.glob(bliss_des50.format(**params))
    if not filenames:
        filenames = glob.glob(bliss_des60.format(**params))
    if not filenames:
        msg = "File not found for: " + str(params)
        msg += '\n'+bliss_des50.format(**params)
        msg += '\n'+bliss_des60.format(**params)
        raise IOError(msg)

    data = []
    for f in filenames:
        data += [create_catalog(f)]
    data = np.hstack(data)

    sel = select(data)
    output = create_output(data[sel],exp)

    logging.info("Writing %s..."%outfile)
    fitsio.write(outfile,output)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('explist')
    parser.add_argument('outfile')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-e','--expnum',default=None,type=int)
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('-t','--tag',default='BLISS',choices=TAGS)
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.getLogger().setLevel(level)

    explist = pd.read_csv(args.explist).to_records()
    explist.dtype.names = map(str.upper,explist.dtype.names)
    exp = explist[explist['EXPNUM'] == args.expnum]

    if len(exp) == 0:
        raise ValueError("EXPNUM not found: %s"%args.expnum)
    elif len(exp) > 1:
        msg = "Multiple values for EXPNUM found: %s"%args.expnum
        for e in exp:
            msg += ("\n"+e)
        raise ValueError()
    exp = exp[0]

    TAG = args.tag
    if args.tag == 'BLISS':
        select = bliss_selection
    elif 'MAGLITES' in args.tag :
        select = maglites_selection
    else:
        raise Exception('Tag not found: %s'%args.tag)

    downskim(args.outfile,select,exp,args.force)
