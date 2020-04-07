#!/usr/bin/env python
import os,shutil
import glob
import numpy as np
import fitsio
import pandas as pd
import numpy.lib.recfunctions as recfuncs
import healpy

from ugali.utils.logger import logger
from ugali.utils.projector import cel2gal
#from ugali.utils.healpix import ang2pix
from ugali.utils.shell import mkdir

HPXBASE = 'hpx_%05d.fits'

# Alternative RA,DEC names
ALT_RADEC_COLUMNS = [
    ['RA','DEC'],
    ['RA','Dec'],
    ['ALPHAWIN_J2000','DELTAWIN_J2000'],
    ]

# Object to string conversion
OBJ = '|O'
STR = '|S30'

def ang2pix(nside, lon, lat, nest=False):
    """
    Input (lon, lat) in degrees instead of (theta, phi) in radians
    """
    theta = np.radians(90. - lat)
    phi = np.radians(lon)
    return healpy.ang2pix(nside, theta, phi, nest)

def read_csv_float32(filename):
    """
    Read a csv file at 32bit precision (excluding ra/dec columns

    Adapted from: https://stackoverflow.com/a/30495210/4075339
    """
    # Sample 100 rows of data to determine dtypes.
    df_test = pd.read_csv(filename, nrows=100, encoding='ascii')
     
    float_cols = np.array([c for c in df_test if df_test[c].dtype == "float64"])
    exclude = [name for pair in ALT_RADEC_COLUMNS for name in pair]
    float_cols = float_cols[~np.in1d(float_cols,exclude)]
    float32_cols = dict([(c, np.float32) for c in float_cols])

    return pd.read_csv(filename,encoding='ascii',engine='c',dtype=float32_cols).to_records(index=False)

def readfile(filename,float32=False):
    """
    Abstract file reading to deal with raw catalog files.
    """
    base,ext = os.path.splitext(filename)
    if np.char.endswith(filename,['.fits','.fits.gz','.fz']).sum():
        f = fitsio.FITS(filename,'r',upper=True)
        idx = 'LDAC_OBJECTS' if 'LDAC_OBJECTS' in f else 1
        nrows = f[idx].get_nrows()
        logger.info("%i objects found"%nrows)
        if not nrows: 
            f.close()
            return
         
        data = f[idx].read()
        f.close()
    elif np.char.endswith(filename, ['.csv','.csv.gz']).sum():
        # This is not perfect, but hey, it works...
        if float32:
            data = read_csv_float32(filename)
        else:
            data = pd.read_csv(filename,encoding='ascii').to_records(index=False)
        dtype = [(str(n).upper(),d if d!=OBJ else STR) for n,d in data.dtype.descr]
        data = data.astype(dtype) # not efficient...see float32 idea above
        nrows = len(data)
        logger.info("%i objects found"%nrows)
        if not nrows: 
            return

    names = list(data.dtype.names)
    if ('RA' not in names) and ('DEC' not in names):
        for ra,dec in np.char.upper(ALT_RADEC_COLUMNS):
            if (ra in names) and (dec in names):
                names[names.index(ra)] = 'RA'
                names[names.index(dec)] = 'DEC'
                data.dtype.names = names
                break
        else:
            msg = "No RA,DEC columns found."
            raise ValueError(msg)

    return data

def pixelize(infiles,outdir='hpx',outbase=HPXBASE,nside=16,gzip=False,force=False,
             float32=False):
    """
    Break catalog up into a set of healpix files.
    """

    mkdir(outdir)
    outfiles = glob.glob(outdir+'/*.fits')
    if len(outfiles) and not force:
        msg = "Found files: %s"%glob.glob(outdir+'/*.fits')
        raise Exception(msg)

    #if len(outfiles): 
    #    print("Removing existing files...")
    #    map(os.remove,outfiles)

    for ii,infile in enumerate(infiles):
        logger.info('(%i/%i) %s'%(ii+1, len(infiles), infile))
        data = readfile(infile,float32)
        if data is None: continue

        catalog_pix = ang2pix(nside,data['RA'],data['DEC'])
        #### Add object pixel (hack to get correct byte order)
        #object_pix = ang2pix(NSIDE_OBJ,data['RA'],data['DEC'],nest=True)
        #name = 'HPX%i'%NSIDE_OBJ; dtype = '>i4'
        #data = recfuncs.rec_append_fields(data,name,object_pix,dtype)

        for pix in np.unique(catalog_pix):
            logger.debug("Processing pixel %s"%pix)

            if 'BAND' in data.dtype.names:
                band = np.unique(data['BAND'])
                if len(band) != 1 and not force:
                    msg = "Found bands: %s"%band
                    raise Exception(msg)
                band = ','.join([b.strip() for b in band])
            else:
                band = 'None'

            outfile = os.path.join(outdir,outbase%pix)
            arr = data[catalog_pix == pix]

            if not os.path.exists(outfile):
                logger.debug("Creating %s"%outfile)
                out=fitsio.FITS(outfile,mode='rw')
                out.write(arr)
                out[1].write_key('COORDSYS','CEL',comment='Coordinate system')
                out[1].write_key('ORDERING','RING',comment='HEALPix ordering scheme')
                out[1].write_key('NSIDE',nside,comment='HEALPix nside')
                out[1].write_key('HPX',pix,comment='HEALPix pixel (RING)')
                out[1].write_key('BAND',band,comment='Photometric band')
            else:
                out=fitsio.FITS(outfile,mode='rw')
                out[1].append(arr)

            logger.debug("Writing %s"%outfile)
            out.close()

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('indir',help='input directory')
    parser.add_argument('outdir',help='output directory')
    parser.add_argument('--ra-dec',nargs=2,
                        help='names of input RA,DEC columns (case insensitive)')
    parser.add_argument('-o','--outbase',default=HPXBASE,
                        help='output file basename')
    parser.add_argument('-n','--nside',default=16,type=int,
                        help='output nside')
    parser.add_argument('--float32',action='store_true',
                        help='convert columns to float32 (excludes RA,DEC)')
    parser.add_argument('-f','--force',action='store_true',
                        help='overwrite existing files')
    parser.add_argument('-v','--verbose',action='store_true',
                        help='output verbosity')
    args = parser.parse_args()

    if args.verbose: logger.setLevel(logger.DEBUG)
    if args.ra_dec:
        ALT_RADEC_COLUMNS = [args.ra_dec] + ALT_RADEC_COLUMNS

    # Grab fits or csv files
    for ext in ['.fits','.fits.gz','.csv','.csv.gz']:
        infiles = sorted(glob.glob(args.indir+'/*'+ext))
        if infiles: break

    pixelize(infiles,args.outdir,args.outbase,nside=args.nside,
             force=args.force,float32=args.float32)
