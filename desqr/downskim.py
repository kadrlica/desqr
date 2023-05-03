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

try:
    from utils import mkdir, isstring
    from const import BANDS, TAGS
except ModuleNotFoundError:
    from .utils import mkdir, isstring
    from .const import BANDS, TAGS


# Order defined below
des70 = '/data/des70.c/data/BLISS/{dirname}/{expnum}'
des50 = '/data/des50.b/data/BLISS/{dirname}/{expnum}'
des60 = '/data/des60.b/data/BLISS/{dirname}/{expnum}'
des61 = '/data/des61.b/data/BLISS/{dirname}/{expnum}'
filebase = 'D{expnum:08d}_*_fullcat.fits'

PATHS = [os.path.join(dirname,filebase) for dirname in [des70,des50,des60,des61]]

#print("WARNING: NGC55 hack")
#PATHS = ['/home/s1/kadrlica/projects/delve/deep/v1/data/{dirname}/{expnum}/D00{expnum}_*-fullcat.fits']

# Old skim column names
DTYPES = [('FILENAME', 'S48'), ('PFW_ATTEMPT_ID', '>i8'), ('TAG', 'S13'), ('UNITNAME', 'S9'), ('REQNUM', '>i4'), ('ATTNUM', '>i2'), ('EXPNUM', '>i8'), ('CCDNUM', '>i2'), ('BAND', 'S1'), ('T_EFF', '>f4'), ('FWHM_WORLD', '>f4'), ('FLAGS', '>i2'), ('OBJECT_NUMBER', '>i4'), ('RA', '>f8'), ('DEC', '>f8'), ('FLUX_PSF', '>f4'), ('FLUXERR_PSF', '>f4'), ('FLUX_AUTO', '>f4'), ('FLUXERR_AUTO', '>f4'), ('A_IMAGE', '>f4'), ('B_IMAGE', '>f4'), ('THETA_IMAGE', '>f4'), ('CLASS_STAR', '>f4'), ('SPREAD_MODEL', '>f4'), ('SPREADERR_MODEL', '>f4'),('IMAFLAGS_ISO','>i2'),('MJD_OBS','>f8'),('EXPTIME','>f4')]
# New reduced column set...
#DTYPES = [('EXPNUM', '>i8'), ('CCDNUM', '>i2'), ('BAND', 'S1'), ('TAG', 'S13'), ('T_EFF', '>f4'), ('FWHM_WORLD', '>f4'), ('MJD_OBS','>f8'), ('EXPTIME','>f4'), ('OBJECT_NUMBER', '>i4'), ('RA', '>f8'), ('DEC', '>f8'), ('FLUX_PSF', '>f4'), ('FLUXERR_PSF', '>f4'), ('FLUX_AUTO', '>f4'), ('FLUXERR_AUTO', '>f4'), ('A_IMAGE', '>f4'), ('B_IMAGE', '>f4'), ('THETA_IMAGE', '>f4'), ('KRON_RADIUS', '>f4'), ('FLUX_RADIUS', '>f4'), ('CLASS_STAR', '>f4'), ('SPREAD_MODEL', '>f4'), ('SPREADERR_MODEL', '>f4'), ('FLAGS', '>i2'), ('IMAFLAGS_ISO','>i2')]

TAG = ''

def get_dirname(expnum):
    """ Get base directory name """
    return expnum//100 * 100

def get_filenames(expnum,path=None):
    """ Get list of filenames """
    dirname = get_dirname(expnum)
    params = dict(expnum=expnum,dirname=dirname)

    paths = [os.path.join(path,filebase)] if path is not None else PATHS

    for path in paths:
        filenames = glob.glob(path.format(**params))
        if filenames: 
            msg = "Found files in path: %s"%(os.path.dirname(filenames[0]))
            logging.debug(msg)
            break

    if not filenames:
        msg = "File not found for: " + str(params)
        for path in paths:
            msg += '\n'+path.format(**params)
        raise IOError(msg)
    return filenames
    
def bliss_select(data):
    """ Quick catalog object selection """
    # Require PSF and AUTO measurements
    sel = data['FLUX_PSF'] > 0
    sel &= data['FLUX_AUTO'] > 0
    # Not saturated, but blends and neighbors ok
    sel &= data['FLAGS'] < 4
    # No "bad" imaflags_iso
    # http://des-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=8814
    sel &= (data['IMAFLAGS_ISO'] & 2047) == 0 
    # MAGERR_PSF < 0.5 (PSF S/N > 2.17)
    olderr = np.seterr(invalid='ignore',divide='ignore')
    sel &= 1.0857362*(data['FLUXERR_PSF']/data['FLUX_PSF']) < 0.5
    np.seterr(**olderr)
    return sel

maglites_select = bliss_select
delve_select = bliss_select

def fgcm_select(data):
    """ Catalog selection for input to FGCM """
    sel = data['FLUX_PSF'] > 0
    sel &= data['FLAGS'] < 4
    sel &= (data['IMAFLAGS_ISO'] & 2047) == 0 

    ### ADW: We should have been selecting only stars
    #sel &= np.abs(data['SPREAD_MODEL']) < 0.01    

    olderr = np.seterr(invalid='ignore',divide='ignore')
    sel &= data['FLUX_PSF']/data['FLUXERR_PSF'] > 5
    np.seterr(**olderr)
    return sel
    
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
        if n in ['PFW_ATTEMPT_ID','TAG','UNITNAME','REQNUM','ATTNUM','T_EFF']:
            if n not in map(str.upper,exp.dtype.names): continue
            try:               out[n] = exp[n]
            except ValueError: out[n] = exp[n.lower()]

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

    # Mapping names...
    names = cat.dtype.names

    # Header names
    for name in ['EXPNUM','CCDNUM','BAND','EXPTIME']:
        cat[name] = hdr[name]
    if 'BAND' in names:     cat['BAND'] = hdr['FILTER']

    # Other names
    if 'FILENAME' in names: cat['FILENAME'] = os.path.basename(filename)
    if 'PFW_ATTEMPT_ID' in names: cat['PFW_ATTEMPT_ID'] = int(str(hdr['EXPNUM'])+'01')
    if 'TAG' in names:      cat['TAG'] = tag 
    if 'UNITNAME' in names: cat['UNITNAME'] = 'D%08d'%hdr['EXPNUM']
    if 'REQNUM' in names:   cat['REQNUM'] = 1
    if 'ATTNUM' in names:   cat['ATTNUM'] = 1
    if 'MJD_OBS' in names:  cat['MJD_OBS'] = hdr['MJD-OBS']
    if 'T_EFF' in names:    cat['T_EFF'] = -1
    if 'RA' in names:       cat['RA'] = data['ALPHAWIN_J2000']
    if 'DEC' in names:      cat['DEC'] = data['DELTAWIN_J2000']
    if 'OBJECT_NUMBER' in names: cat['OBJECT_NUMBER'] = data['NUMBER']

    return cat

def create_header(fits):
    from astropy.io.fits import Header

    if isstring(fits):
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
    path = None if 'PATH' not in exp.dtype.names else exp['PATH']
    filenames = get_filenames(expnum,path)
    output = create_downskim(filenames,select,exp,DTYPES,TAG)

    outdir = os.path.dirname(outfile)
    mkdir(outdir)

    logging.info("Writing %s..."%outfile)
    dirname = os.path.dirname(filenames[0])
    header = [dict(name='PATH',value=dirname,comment='original path')]
    fitsio.write(outfile,output,header=header,clobber=True)

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

    explist = pd.read_csv(args.explist,comment='#').to_records(index=False)
    explist.dtype.names = [str(n).upper() for n in explist.dtype.names]
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
