#!/usr/bin/env python
"""
Calculate extinction for an input catalog.
"""
__author__ = "Alex Drlica-Wagner"
__email__  = "kadrlica@fnal.gov"

import os
import shutil
from collections import OrderedDict as odict
import tempfile
import subprocess
import logging

import pandas as pd
import numpy as np
import fitsio
import healpy
import healpy as hp

# These are some packages that are hanging around
import ugali.utils.healpix as healpix
from ugali.utils.projector import cel2gal
import utils
from ugali.utils.logger import logger


BANDS = ['g','r','i','z','Y','y']

# Rv = 3.1
# Schlafly & Finkbeiner ApJ 737, 103 (2011)
SCHLAFLY11 = odict([
        ('g',3.237), 
        ('r',2.176), 
        ('i',1.595), 
        ('z',1.217), 
        ('Y',1.058), 
        ('y',1.058),
        ])

# Rv = 3.1
# A dictionary of A(Lambda)/E(B-V) values for a flat SED with Rv=3.10
# https://dessvn.cosmology.illinois.edu/websvn/desdm/devel/extinction
# https://dessvn.cosmology.illinois.edu/websvn/desdm/devel/extinction/trunk/python/extinction/extinction_tools.py
OLD_DESDM =  odict([
        # From filter scans on 02/2012
        ('g',3.704722),
        ('r',2.610357),
        ('i',1.947345),
        ('z',1.496843),
        ('Y',1.311188),
        ('y',1.311188),
])

DESDM = odict([
        # From filter scans on 03/2013
        ('u', 4.708272),
        ('g', 3.682995),
        ('r', 2.604808),
        ('i', 1.940133), 
        ('z', 1.450496),
        ('Y', 1.277421),
        ('y', 1.277421),
        ])

# R_lambda values for Y3A1 from:
# https://cdcvs.fnal.gov/redmine/projects/des-y3/wiki
Y3A1 = odict([
        ('g',3.185),
        ('r',2.140),
        ('i',1.571),
        ('z',1.198),
        ('Y',1.052),
        ('y',1.052),
        ])

OUTCOLS  = ['EXTINCTION','EBV']
FITSEXT = ('.fits','.fit','.fz')
CSVEXT = ('.csv')

def readfile(filename,columns=None):
    """Read data from CSV or FITS file.
    """
    base,ext =  os.path.splitext(filename)
    if ext in FITSEXT:
        return fitsio.read(args.infile,ext=1,columns=columns)
    elif ext in CSVEXT:
        return pd.read_csv(filename).to_records()
    else:
        msg = "Unrecognized file extension: %s"%ext
        raise Exception(msg)

def writefile(filename,data,force=False):
    """Write data to CSV or FITS file
    """
    base,ext =  os.path.splitext(filename)
    if ext in FITSEXT:
        utils.insert_columns(filename,data,force=args.force)        
    elif ext in CSVEXT:
        orig = pd.read_csv(filename).to_records()
        dtype = orig.dtype.descr + data.dtype.descr
        out = np.recarray(len(data),dtype=dtype)

        for col in orig.dtype.names:
            out[col] = orig[col]
        for col in data.dtype.names:
            out[col] = data[col]

        # Write to the CSV file using pandas
        pd.DataFrame(out).to_csv(filename,index=False)
    else:
        msg = "Unrecognized file extension: %s"%ext
        raise Exception(msg)

def ebv(ra,dec,ebvmap=None):
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    if not len(ra) == len(dec):
        msg = "Column lengths must match"
        raise Exception(msg)

    if ebvmap is None or ebvmap.lower() == 'sfd':
        # Download SFD map
        url = "http://lambda.gsfc.nasa.gov/data/foregrounds/SFD/lambda_sfd_ebv.fits"
        logger.info("Downloading %s..."%url)
        filename = tempfile.NamedTemporaryFile().name
        cmd = "wget %s -O %s"%(url,filename)
        subprocess.call(cmd,shell=True)
        ebvmap = healpy.read_map(filename)
        os.remove(filename)
    elif isinstance(ebvmap,basestring):
        logger.info("Loading %s..."%ebvmap)
        ebvmap = healpy.read_map(ebvmap)
    else:
        msg = "Unrecognized ebv: %s"%ebvmap
        raise Exception(msg)

    # The SFD map is in Galactic coordinates
    glon,glat = cel2gal(ra,dec)
    ebv = healpix.get_interp_val(ebvmap,glon,glat)
    return ebv

def extinction(ebv,band):
    """
    Calculate the extinction from the E(B-V) value and the band.
    
    ebv  : The dust value
    band : The DES band (string or array)
    """
    if isinstance(band,basestring):
        band = np.repeat(band,len(ebv))
        
    bands,inv = np.unique(band,return_inverse=True)
    values = np.array([DESDM[b] for b in bands])
    return values[inv] * ebv

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infile',help='input file')
    parser.add_argument('outfile',nargs='?',help='output file')
    parser.add_argument('-e','--ebv',default='SFD',
                        help='E(B-V) name or file')
    parser.add_argument('--ra',default='RA',
                        help='name of right ascension column')
    parser.add_argument('--dec',default='DEC',
                        help='name of declination column')
    parser.add_argument('-b','--bands',default=[],action='append',
                        help='bands to calculate extinction')
    parser.add_argument('--ext',default='EXTINCTION',
                        help='name of extinction output column')
    parser.add_argument('-f','--force',action='store_true',
                        help='overwrite output columns if they exist')
    parser.add_argument('-v','--verbose',action='store_true',
                        help='output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("Reading %s..."%args.infile)

    columns = [args.ra,args.dec]

    if len(args.bands) == 0:
        args.bands = BANDS
    elif len(args.bands) == 1 and args.bands[0] not in BANDS:
        columns += args.bands[0]
        bands = None
    else:
        bands = args.bands

    data = readfile(args.infile,columns)
    ra = data[args.ra]
    dec = data[args.dec]

    logger.info("Calculating E(B-V): %s..."%args.ebv)
    ebvval = ebv(ra,dec,args.ebv)

    values = [ebvval]
    dtypes = [('EBV','f4')]

    logger.info("Calculating extinction %s..."%args.ebv)
    for b in args.bands:
        if b in DESDM:
            band = np.repeat(b,len(ra))
            extname = args.ext + '_%s'%b.upper()
        else:
            band = data[b]
            extname = args.ext

        extval = extinction(ebvval,band)
        values.append(extval)
        dtypes.append((extname,'f4'))

    out = np.rec.fromarrays(values,dtype=dtypes)
    
    # Writing...
    if not args.outfile: 
        args.outfile = args.infile
    else: 
        logger.debug("Copying %s to %s..."%(args.infile,args.outfile))
        shutil.copy(args.infile,args.outfile)

    logger.info("Writing %s..."%args.outfile)
    writefile(args.outfile,out,force=args.force)
    logger.info("Done.")
