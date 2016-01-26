#!/usr/bin/env python
import shutil
import os
from collections import OrderedDict as odict
import itertools

import numpy as np
import fitsio
import healpy

from ugali.utils.logger import logger

# Tools for working with the shell
def pwd():
    # Careful, won't work after a call to os.chdir...
    return os.environ['PWD']

def mkdir(dir):
    if not os.path.exists(dir):  os.makedirs(dir)
    return dir

def mkscratch():
    if os.path.exists('/scratch/'):    
        return(mkdir('/scratch/%s/'%os.environ['USER']))
    elif os.path.exists('/tmp/'):
        return(mkdir('/tmp/%s/'%os.environ['USER']))
    else:
        raise Exception('...')

def multiproc(func,args,kwargs):
    """
    Wrapper around `apply_async` to call `funcs`.

    func  : callable
    args  : list of tuples
    kwargs: dict
    """
    from multiprocessing import Pool

    nargs = len(args)
    results = nargs*[None]

    pool = Pool(maxtasksperchild=1)
    for i,arg in enumerate(args):
        def callback(result,i=i):
            results[i] = result
        res = pool.apply_async(func,arg,kwargs,callback=callback)
    pool.close()
    pool.join()
    return results

def found(filename):
    logger.warning("Found %s; skipping..."%filename)

def which(program):
    # http://stackoverflow.com/q/377017
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

# Various utility functions

def bfields(fields,bands):
    """ Create band-by-band fields"""
    if isinstance(fields,basestring): fields = [fields]
    if isinstance(bands,basestring): bands = [bands]
    bands = np.char.upper(np.char.strip(bands))
    out = [f+'_%s'%b for f,b in itertools.product(fields,bands)]
    if len(out) == 1: return out[0]
    else:             return out

def unique_pair(a,b):
    A = np.vstack([a,b]).T
    B = A[np.lexsort(A.T[::-1])]
    return B[np.concatenate(([True],np.any(B[1:]!=B[:-1],axis=1)))]

def mean_value(value,labels):
    index = np.unique(labels)
    mean_value = ndimage.mean(value,labels=labels,index=index)
    return mean_value

# Dealing with FITS files
def write(filename,data,header=None,force=False):
    if os.path.exists(filename) and not force:
        found(filename)
        return
    fitsio.write(filename,data,header=header,clobber=force)

def insert_columns(filename,data,ext=1,force=False):
    print filename
    fits = fitsio.FITS(filename,'rw')
    names = fits[ext].get_colnames()
    overlap = np.in1d(data.dtype.names,names)

    test = None
    if np.any(~overlap):
        idx = np.argmax(np.in1d(names,data.dtype.names))
        test = names[idx]
        orig = fits[ext].read(columns=[test])

    if np.any(overlap) and not force:
        logger.warning("Found overlapping columns; skipping...")
        return
    if len(data) != fits[ext].get_nrows():
        logger.warning("Number of rows does not match; skipping...")
        return
    for col in data.dtype.names:
        if col not in names:
            msg = "Inserting column: %s"%col
            logger.info(msg)
            fits[ext].insert_column(col,data[col])
        else:
            msg = "Found column %s"%col
            logger.warning(msg)
            fits[ext].write_column(col,data[col])

    fits.close()

    # It's already too late since the file has been written...
    if test is not None:
        new = fitsio.read(filename,ext=ext,columns=[test])
        if np.any(new != orig):
            msg = "Input and output do not match!"
            raise Exception(msg)


def insert_index(filename,data,ext=1,force=False):
    d = fitsio.read(filename,columns=['FILENAME','OBJECT_NUMBER'])    
    if np.any(d[['FILENAME','OBJECT_NUMBER']]!=data[['FILENAME','OBJECT_NUMBER']]):
        raise Exception("Mismatch in FILENAME or OBJECT_NUMBER")

    cols = list(data.dtype.names)
    cols.remove('FILENAME')
    cols.remove('OBJECT_NUMBER')
    insert_columns(filename,data[cols],ext,force)

def ccdnum2(infile,outfile=None,force=True):
    import pyfits

    f = pyfits.open(infile)
    names = np.char.upper(f[1].columns.names).tolist()

    if 'CCDNUM' in names:
        print 'CCDNUM column already found; skipping...'
        return

    if outfile is None: outfile = infile
    if os.path.exists(outfile) and not force:
        found(outfile)
        return

    if not 'FILENAME' in names:
        msg = 'FILENAME column not found.'
        raise Exception(msg)

    d = f[1].data
    ccdnum = np.array([a[2].strip('c') for a in np.char.split(d['filename'],'_')],
                      dtype=int)
    idx = names.index('EXPNUM');
    coldefs =  pyfits.ColDefs([
            pyfits.Column(name='CCDNUM',format='I',array=ccdnum)
            ])
    hdu = pyfits.BinTableHDU.from_columns(d.columns[:idx+1]+coldefs+d.columns[idx+1:])

    print "Writing %s..."%outfile
    hdu.writeto(outfile,clobber=force)


def ccdnum(infile,outfile=None,force=True):
    if outfile is None: outfile = infile
    if os.path.exists(outfile) and not force:
        found(outfile)
        return

    f = fitsio.FITS(outfile,'rw')
    names = np.char.upper(f[1].get_colnames()).tolist()

    if outfile != infile: shutil.copy(infile,outfile)

    if 'CCDNUM' in names:
        print 'CCDNUM column already found; skipping...'
        f.close()
        return

    if not 'FILENAME' in names:
        msg = 'FILENAME column not found.'
        f.close()
        raise Exception(msg)

    d = f[1].data
    d = f[1].read(columns=['FILENAME','EXPNUM'])
    ccdnum = np.array([a[2].strip('c') for a in np.char.split(d['FILENAME'],'_')],
                      dtype=int)
    idx = names.index('EXPNUM');

    print "Writing %s..."%outfile
    f[1].insert_column('CCDNUM',ccdnum,colnum=idx+1)
    f.close()

def load(args):
    infile,columns = args
    logger.debug("Loading %s..."%infile)
    return fitsio.read(infile,columns=columns)

def load_infiles(infiles,columns=None,multiproc=False):
    if isinstance(infiles,basestring):
        infiles = [infiles]

    logger.debug("Loading %s files..."%len(infiles))

    args = zip(infiles,len(infiles)*[columns])

    if multiproc:
        from multiprocessing import Pool
        p = Pool(maxtasksperchild=1)
        out = p.map(load,args)
    else:
        out = [load(arg) for arg in args]

    data = None
    for i,d in enumerate(out):
        if i == 0: 
            data = d
        elif d.dtype != data.dtype: 
            # ADW: Not really safe...
            data = np.append(data,d.astype(data.dtype))
        else:
            data = np.append(data,d)

    return data

def uid(expnum,ccdnum):
    return expnum + ccdnum/100.

def pix2ang(nside, pix, nest=False):
    """
    Return (lon, lat) in degrees instead of (theta, phi) in radians
    """
    theta, phi =  healpy.pix2ang(nside, pix, nest)
    lon = np.degrees(phi)
    lat = 90. - np.degrees(theta)                    
    return lon, lat

def ang2pix(nside, lon, lat, nest=False):
    """
    Input (lon, lat) in degrees instead of (theta, phi) in radians
    """
    theta = np.radians(90. - lat)
    phi = np.radians(lon)
    return healpy.ang2pix(nside, theta, phi, nest)


def get_vizier_catalog(ra,dec,radius=None,**kwargs):
    import warnings
    from astroquery.vizier import Vizier
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    kwargs.setdefault('row_limit',-1)
    coord = SkyCoord(ra*u.deg,dec*u.deg)
    radius = u.Quantity(radius,u.deg)
    vizier = Vizier(**kwargs)
    warnings.filterwarnings("ignore")
    tab = vizier.query_region(coord,radius)
    warnings.resetwarnings()
    return tab[0]

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    opts = parser.parse_args()


def print_problem(msg):
    import termcolor as color
    print color(msg,'red')
