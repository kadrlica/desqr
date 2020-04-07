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
bliss_des61 = '/data/des61.b/data/BLISS/{dirname}/{expnum}/D00{expnum}_*_fullcat.fits'

outtemp = '%(unitname)s_%(band)s.fits'

DTYPES = [('FILENAME', 'S48'), ('PFW_ATTEMPT_ID', '>i8'), ('TAG', 'S13'), ('UNITNAME', 'S9'), ('REQNUM', '>i4'), ('ATTNUM', '>i2'), ('EXPNUM', '>i8'), ('CCDNUM', '>i2'), ('BAND', 'S1'), ('T_EFF', '>f4'), ('FWHM_WORLD', '>f4'), ('FLAGS', '>i2'), ('OBJECT_NUMBER', '>i4'), ('RA', '>f8'), ('DEC', '>f8'), ('FLUX_PSF', '>f4'), ('FLUXERR_PSF', '>f4'), ('FLUX_AUTO', '>f4'), ('FLUXERR_AUTO', '>f4'), ('A_IMAGE', '>f4'), ('B_IMAGE', '>f4'), ('THETA_IMAGE', '>f4'), ('CLASS_STAR', '>f4'), ('SPREAD_MODEL', '>f4'), ('SPREADERR_MODEL', '>f4'),('IMAFLAGS_ISO','>i2'),('MJD_OBS','>f8'),('EXPTIME','>f4')]
TAG = ''

def get_dirname(expnum):
    """ Get base directory name """
    return expnum//100 * 100

def get_filenames(expnum):
    """ Get list of filenames """
    dirname = get_dirname(expnum)
    params = dict(expnum=expnum,dirname=dirname)
    
    for path in [bliss_des50,bliss_des60,bliss_des61]:
        filenames = glob.glob(path.format(**params))
        if filenames: break

    if not filenames:
        msg = "File not found for: " + str(params)
        msg += '\n'+bliss_des50.format(**params)
        msg += '\n'+bliss_des60.format(**params)
        msg += '\n'+bliss_des61.format(**params)
        raise IOError(msg)
    return filenames
    
def bliss_select(data):
    sel = data['FLUX_PSF'] > 0
    sel &= data['FLUX_AUTO'] > 0
    sel &= data['FLAGS'] < 4
    sel &= (data['IMAFLAGS_ISO'] & 2047) == 0 
    olderr = np.seterr(invalid='ignore',divide='ignore')
    sel &= 1.0857362*(data['FLUXERR_PSF']/data['FLUX_PSF']) < 0.5
    np.seterr(**olderr)
    return sel

maglites_select = bliss_select
delve_select = bliss_select
    
def create_output(data,exp,dtype=DTYPES):
    out = np.recarray(len(data),dtype=dtype)
    for n in out.dtype.names:
        if n not in exp.dtype.names + data.dtype.names:
            msg = "Column not found: %s"%n
            raise Exception(msg)

        # This has defaults set
        if n in data.dtype.names:
            out[n] = data[n]

        # Fill some columns from exp
        if n in ['PFW_ATTEMPT_ID','TAG','UNITNAME','REQNUM','ATTNUM','T_EFF','TAG']:
            out[n] = exp[n]

    return out

def create_catalog(filename,dtype=DTYPES,tag=TAG):
    fits = fitsio.FITS(filename)

    hdr = create_header(fits)
    data = fits['LDAC_OBJECTS'].read()

    cat = np.recarray(len(data),dtype=dtype)

    for name,dt in dtype:
        if '_APER_' in name:
            prefix,idx = name.rsplit('_',1)
            cat[name] = data[prefix][:,int(idx)]
            continue

        if name not in data.dtype.names: continue
        cat[name] = data[name]

            
    cat['FILENAME'] = os.path.basename(filename)
    cat['PFW_ATTEMPT_ID'] = int(str(hdr['EXPNUM'])+'01')
    cat['TAG'] = tag 
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


def create_downskim(filenames, select, exp, dtype=DTYPES, tag=TAG):
    """ Create a downskimmed catalog from filenames 

    Parameters
    ----------
    filenames : list of filenames
    select : selection function
    exp : exposure info

    Returns
    -------
    catalog
    """
    data = []
    for f in filenames:
        try: 
            data += [create_catalog(f,dtype,tag)]
        except IOError as e:
            # Sometimes the processing has failed for a specific file
            logging.warn(str(e))

    data = np.hstack(data)

    if select is None: sel = slice(None)
    else:              sel = select(data)
    output = create_output(data[sel],exp,dtype)

    return output
    
def downskim(outfile,select,exp,force):
    """ Like download, but for skimming catalog files.
    """
    if os.path.exists(outfile) and not force:
        logging.warn('Found %s; skipping...'%outfile)
        return

    expnum = exp['EXPNUM']
    filenames = get_filenames(expnum)
    output = create_downskim(filenames,select,exp,DTYPES,TAG)

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

    explist = pd.read_csv(args.explist,encoding='ascii').to_records()
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
    if TAG.startswith('BLISS'):
        select = bliss_select
    elif TAG.startswith('MAGLITES'):
        select = maglites_select
    elif TAG.startswith("DELVE"):
        select = delve_select
    else:
        raise Exception('Tag not found: %s'%args.tag)

    downskim(args.outfile,select,exp,args.force)
