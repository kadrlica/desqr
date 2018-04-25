#!/usr/bin/env python
"""
Calculate extinction for an input catalog.
"""
__author__ = "Alex Drlica-Wagner"

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

# These are some packages that are hanging around
import ugali.utils.healpix as healpix
import utils

# Rv = 3.1
# Schlafly & Finkbeiner ApJ 737, 103 (2011)
SCHLAFLY11 = odict([
        ('g',3.237), 
        ('r',2.176), 
        ('i',1.595), 
        ('z',1.217), 
        ('Y',1.058), 
        ])

# Rv = 3.1
# https://dessvn.cosmology.illinois.edu/websvn/desdm/devel/extinction
DESEXT = odict([
        ('g',3.704722),
        ('r',2.610357),
        ('i',1.947345),
        ('z',1.496843),
        ('Y',1.311188),
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

def extinction(ra,dec,band,ebvmap=None):
    if isinstance(band,basestring):
        band = np.array(len(ra)*[band])

    if not len(ra) == len(dec) == len(band):
        msg = "Column lengths must match"
        raise Exception(msg)

    out = np.recarray(len(ra),dtype=[(n,'f2') for n in OUTCOLS])
    out.fill(0)

    if ebvmap is None or ebvmap.lower() == 'sfd':
        # Download SFD map
        url = "http://lambda.gsfc.nasa.gov/data/foregrounds/SFD/lambda_sfd_ebv.fits"
        logging.info("Downloading %s..."%url)
        filename = tempfile.NamedTemporaryFile().name
        cmd = "wget %s -O %s"%(url,filename)
        subprocess.call(cmd,shell=True)
        ebvmap = healpy.read_map(filename)
        os.remove(filename)
    elif isinstance(ebvmap,basestring) and os.path.exists(ebvmap):
        logging.info("Loading %s..."%ebvmap)
        ebvmap = healpy.read_map(ebvmap)

    ebv = healpix.get_interp_val(ebvmap,ra,dec)
    out['EBV'] = ebv
    
    bands,inv = np.unique(band,return_inverse=True)
    values = np.array([DESEXT[b] for b in bands])
    out['EXTINCTION'] = values[inv] * ebv

    return out

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infile',help='input file')
    parser.add_argument('outfile',nargs='?',help='output file [optional]')
    parser.add_argument('-e','--ebv',default='SFD',
                        help='E(B-V) name or file')
    parser.add_argument('--ra',default='RA',
                        help='name of right ascension column')
    parser.add_argument('--dec',default='DEC',
                        help='name of declination column')
    parser.add_argument('--band',default='BAND',
                        help='name of band column')
    parser.add_argument('-f','--force',action='store_true',
                        help='overwrite output columns if they exist')
    args = parser.parse_args()

    logging.getLogger().setLevel(logging.INFO)

    logging.info("Reading %s..."%args.infile)
    columns = [args.ra,args.dec,args.band]
    data = readfile(args.infile,columns)
    ra = data[args.ra]
    dec = data[args.dec]
    band = data[args.band]

    logging.info("Calculating extinction %s..."%args.ebv)
    out = extinction(ra,dec,band,args.ebv)
    
    # Writing...
    logging.info("Writing %s..."%args.infile)
    if not args.outfile: 
        args.outfile = args.infile
    else: 
        logging.debug("Copying %s to %s..."%(args.infile,args.outfile))
        shutil.copy(args.infile,args.outfile)
    writefile(args.outfile,out,force=args.force)
