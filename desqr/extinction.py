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

import pandas as pd
import numpy as np
import fitsio
import healpy as hp

# These are some packages that are hanging around
import ugali.utils.healpix as healpix
from ugali.utils.projector import cel2gal

from desqr import utils
from desqr.logger import logger

BANDS = ['g','r','i','z','Y','y']

# DES delta(m)/E(B-V) at Rv = 3.1
# Schlafly & Finkbeiner ApJ 737, 103 (2011)
# https://iopscience.iop.org/article/10.1088/0004-637X/737/2/103/pdf
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

# R_lambda values for DES DR1 from:
# https://arxiv.org/abs/1801.03181
DR1 = odict([
        ('g',3.186),
        ('r',2.140),
        ('i',1.569),
        ('z',1.196),
        ('Y',1.048),
        ('y',1.048),
        ])
# DES DR2 is the same as DES DR1
# https://arxiv.org/abs/2101.05765
DR2 = DR1

# From Table 6 in Schlafly 2011 with Rv = 3.1
# http://iopscience.iop.org/article/10.1088/0004-637X/737/2/103/pdf
PS1 = odict([
        ('g', 3.172 ),
        ('r', 2.271 ),
        ('i', 1.682 ),
        ('z', 1.322 ),
        ('y', 1.087 ),
        ('w', 2.341 ),
        ])

OUTCOLS  = ['EXTINCTION','EBV']
FITSEXT = ('.fits','.fit','.fz')
CSVEXT = ('.csv')

def readfile(filename,columns=None,force=False):
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

@utils.ignore_warning(UserWarning)
def ebv(ra,dec,ebvmap=None):
    """Calculate E(B-V) value by interpolating a map."""
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    if not len(ra) == len(dec):
        msg = "Column lengths must match"
        raise Exception(msg)

    if ebvmap is None or (utils.isstring(ebvmap) and ebvmap.lower() == 'sfd'):
        # Download SFD map
        url = "http://lambda.gsfc.nasa.gov/data/foregrounds/SFD/lambda_sfd_ebv.fits"
        logger.info("Downloading %s..."%url)
        filename = tempfile.NamedTemporaryFile().name
        cmd = "wget %s -O %s"%(url,filename)
        subprocess.call(cmd,shell=True)
        ebvmap = hp.read_map(filename)
        os.remove(filename)
    elif utils.isstring(ebvmap):
        logger.info("Loading %s..."%ebvmap)
        ebvmap = hp.read_map(ebvmap)

    # The SFD map is in Galactic coordinates
    glon,glat = cel2gal(ra,dec)
    ebv = healpix.get_interp_val(ebvmap,glon,glat)
    return ebv

def extinction(ebv,band,coeff=None):
    """
    Calculate the extinction from the E(B-V) value and the band.

    Parameters
    -----------
    ebv  : The dust value E(B-V)
    band : The survey band (string or array)
    coeff: the reddening coefficients (R_b)

    Returns
    -------
    extinction : The A_b values [A_b = R_b * E(B-V)]
    """
    if utils.isstring(band):
        band = np.repeat(band,len(ebv))

    if coeff is None: coeff = DR1
        
    bands,inv = np.unique(band,return_inverse=True)
    values = np.array([coeff[b] for b in bands])
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
    parser.add_argument('--coeff',default='DR1',choices=['Y3A1','PS1','DESDM','DR1'],
                        help='reddening coefficients')
    parser.add_argument('-f','--force',action='store_true',
                        help='overwrite output columns if they exist')
    parser.add_argument('-v','--verbose',action='store_true',
                        help='output verbosity')
    args = parser.parse_args()

    if args.verbose: logger.setLevel(logger.DEBUG)

    if args.coeff == 'Y3A1':
        COEFF = Y3A1
    elif args.coeff == 'DESDM':
        COEFF = DESDM
    elif args.coeff == 'PS1':
        COEFF = PS1
    elif args.coeff == 'DR1':
        COEFF = DR1
    else:
        raise Exception('Unrecognized coefficient set: %s'%args.coeff)
    logger.info("Using %s coefficients."%args.coeff)

    columns = [args.ra,args.dec]
    if len(args.bands) == 0:
        args.bands = BANDS
    elif len(args.bands) == 1 and args.bands[0] not in BANDS:
        columns += args.bands[0]
        bands = None
    else:
        bands = args.bands

    logger.info("Reading %s..."%args.infile)
    data = readfile(args.infile,columns)
    ra = data[args.ra]
    dec = data[args.dec]

    logger.info("  Calculating E(B-V) from %s."%args.ebv)
    ebvval = ebv(ra,dec,args.ebv)

    values = [ebvval]
    dtypes = [('EBV','f4')]

    logger.info("  Calculating extinction %s..."%args.ebv)
    for b in args.bands:
        if b in DESDM:
            band = np.repeat(b,len(ra))
            extname = args.ext + '_%s'%b.upper()
        else:
            band = data[b]
            extname = args.ext

        extval = extinction(ebvval,band,coeff=COEFF)
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
