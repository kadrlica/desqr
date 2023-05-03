#!/usr/bin/env python
import os
import time
import shutil
import itertools
import glob
from collections import OrderedDict as odict
from functools import wraps
import six
import warnings
import logging

import numpy as np
import numpy.lib.recfunctions as rfn
from numpy.lib import NumpyVersion
import fitsio
import healpy as hp

# Need to do better with the logging...
from ugali.utils.logger import logger
#import logging as logger

# Tools for working with the shell
def pwd():
    # Careful, won't work after a call to os.chdir...
    return os.environ['PWD']

def mkdir(path):
    # https://stackoverflow.com/a/600612/4075339
    import errno
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
    return path

def mkscratch():
    if os.path.exists('/scratch/'):    
        return(mkdir('/scratch/%s/'%os.environ['USER']))
    elif os.path.exists('/tmp/'):
        return(mkdir('/tmp/%s/'%os.environ['USER']))
    else:
        raise Exception('...')

def set_defaults(kwargs,defaults):
    for k,v in defaults.items():
        kwargs.setdefault(k,v)
    return kwargs
setdefaults=set_defaults

def mr_nice_guy(value=10):
    """ Play nice. """
    nice = os.nice(0)
    os.nice(value-nice)

def multiproc(func, arglist, kwargs={}, processes=20):
    """
    Wrapper around `apply_async` to call `funcs`.

    func      : callable
    arglist   : list of argument tuples
    kwargs    : dict
    processes : number of processes
    """

    nargs = len(arglist)
    results = nargs*[None]
    
    # Use all CPUs [os.cpu_count()]
    if processes <= 0: processes = None
    
    if processes == 1:
        results = [func(*args,**kwargs) for args in arglist]
    else:
        from multiprocessing import Pool
        pool = Pool(maxtasksperchild=1,processes=processes)
        for i,args in enumerate(arglist):
            def callback(result,i=i):
                results[i] = result
            res = pool.apply_async(func,args,kwargs,callback=callback)
        pool.close()
        pool.join()
        del pool

    return results

def found(filename):
    #logger = logging.getLogger()
    logger.warning("Found %s; skipping..."%filename)

def is_found(filename,force=False):
    """" Return boolean if file is found. """
    if os.path.exists(filename) and not force:
        found(filename)
        return True
    else:
        return False

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

def isstring(obj):
    """Python 2/3 compatible string check"""
    import six
    return isinstance(obj, six.string_types)

# Various utility functions
def bfield(field,band):
    """ Combine a single field and band. """
    return ('%s_%s'%(field,band.strip())).upper()

def bfields(fields,bands):
    """ Create band-by-band fields. """
    if np.isscalar(fields) and np.isscalar(bands):
        return bfield(fields,bands)
    fields = np.atleast_1d(fields)
    bands  = np.atleast_1d(bands)
    out = [bfield(f,b) for f,b in itertools.product(fields,bands)]
    return out

def unique_pair(a,b):
    A = np.vstack([a,b]).T
    B = A[np.lexsort(A.T[::-1])]
    return B[np.concatenate(([True],np.any(B[1:]!=B[:-1],axis=1)))]

def mean_value(value,labels):
    index = np.unique(labels)
    mean_value = ndimage.mean(value,labels=labels,index=index)
    return mean_value

# Dealing with FITS files
def write_fits(filename,data,header=None,force=False):
    if os.path.exists(filename) and not force:
        found(filename)
        return
    fitsio.write(filename,data,header=header,clobber=force)

write = write_fits

def insert_header_key(filename, name, value, comment=None):
    with fitsio.FITS(filename,'rw') as fits:
        fits[1].write_key(name, value, comment="my comment")
        

def insert_columns(filename,data,ext=1,force=False,check=True):
    """
    Insert columns from a recarray into a table in an existing FITS
    file. Input data must have the same length as existing table. By
    default, several corruption checks are run.
    
    Parameters:
    -----------
    filename : the FITS file to insert into
    data     : the recarray to insert
    ext      : the FITS extension (name or index)
    force    : overwrite columns if they already exist

    Returns:
    --------
    None
    """
    #logger = logging.getLogger()
    logger.info(filename)
    if not os.path.exists(filename):
        msg = "Requested file does not exist."
        raise IOError(msg)

    # ADW: Better to use the context manager

    with fitsio.FITS(filename,'rw') as fits:
        oldnames = fits[ext].get_colnames()
        newnames = data.dtype.names
        overlap = np.in1d(oldnames,newnames)
         
        if np.any(overlap) and not force:
            logger.warning("Found overlapping columns; skipping...")
            return
        if len(data) != fits[ext].get_nrows():
            logger.warning("Number of rows does not match; skipping...")
            return
         
        # Test that at least one other column not corrupted during insert
        test = None
        if np.any(~overlap) and check:
            idx = np.where(~overlap)[0][-1]
            test = oldnames[idx]
            orig = fits[ext].read(columns=[test])
         
        for col in newnames:
            if col not in oldnames:
                msg = "Inserting column %s"%col
                logger.info(msg)
                fits[ext].insert_column(col,data[col])
            else:
                msg = "Overwriting column %s"%col
                logger.warning(msg)
                fits[ext].write_column(col,data[col])

    # The file has already been written, but might as well know...
    if test and check:
        new = fitsio.read(filename,ext=ext,columns=[test])
        np.testing.assert_array_equal(new,orig,"Collateral damage: %s"%test)
        
        name = newnames[0]
        new = fitsio.read(filename,ext=ext,columns=[name])
        np.testing.assert_array_equal(new[name],data[name],
                                      "Input/output mismatch: %s"%name)


def insert_index(filename,data,ext=1,force=False):
    d = fitsio.read(filename,columns=['FILENAME','OBJECT_NUMBER'])    
    if np.any(d[['FILENAME','OBJECT_NUMBER']]!=data[['FILENAME','OBJECT_NUMBER']]):
        raise Exception("Mismatch in FILENAME or OBJECT_NUMBER")

    cols = list(data.dtype.names)
    cols.remove('FILENAME')
    cols.remove('OBJECT_NUMBER')
    insert_columns(filename,data[cols],ext,force)

def ccdnum(infile,outfile=None,force=True):
    """ Grab the CCDNUM from the filename.

    Parameters
    ----------
    infile  : input filename D{expnum:08d}_{band}_[c]{ccdnum}...
    outfile : output filename
    force   : overwrite existing CCDNUM column
    
    Returns
    -------
    None
    """
    
    if outfile is None: outfile = infile
    if os.path.exists(outfile) and not force:
        found(outfile)
        return

    f = fitsio.FITS(outfile,'rw')
    names = np.char.upper(f[1].get_colnames()).tolist()

    if outfile != infile: shutil.copy(infile,outfile)

    if 'CCDNUM' in names:
        print('CCDNUM column already found; skipping...')
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

    print("Writing %s..."%outfile)
    f[1].insert_column('CCDNUM',ccdnum,colnum=idx+1)
    f.close()

def load(args):
    infile,columns = args
    logger.debug("Loading %s..."%infile)
    return fitsio.read(infile,columns=columns)

def load_infiles(infiles,columns=None,multiproc=False):
    if isstring(infiles):
        infiles = [infiles]

    logger.debug("Loading %s files..."%len(infiles))

    args = zip(infiles,len(infiles)*[columns])

    if multiproc:
        from multiprocessing import Pool
        processes = multiproc if multiproc > 0 else None
        p = Pool(processes,maxtasksperchild=1)
        out = p.map(load,args)
    else:
        out = [load(arg) for arg in args]

    dtype = out[0].dtype
    for i,d in enumerate(out):
        if d.dtype != dtype: 
            # ADW: Not really safe...
            logger.warn("Casting input data to same type.")
            out[i] = d.astype(dtype)

    logger.debug('Concatenating arrays...')
    return np.concatenate(out)

def get_pixels(dirname,filename=None):
    """ Get array of exising pixel/filenames in a directory

    Parameters
    ----------
    dirname  : directory
    filename : filename base

    Returns
    -------
    array : healpix pixels and filenames
    """
    filenames = None

    if filename:
        filenames = glob.glob(os.path.join(dirname,filename))
        if not filenames: raise Exception("No files found.")

    if not filenames:
        filenames = glob.glob(os.path.join(dirname,'*.fits'))
    if not len(filenames):
        filenames = glob.glob(os.path.join(dirname,'*.fits.fz'))
    if not len(filenames):
        filenames = glob.glob(os.path.join(dirname,'*.csv'))

    filenames = sorted(filenames)
    basename = np.char.partition(filenames,'_')[:,0]
    pixels = np.char.partition(basename,'.')[:,0].astype(int)
    array = np.rec.fromarrays([pixels,filenames],names=['pixel','filename'])
    return array

def uid(expnum,ccdnum):
    """ Create unique float from EXPNUM, CCDNUM 
    
    Parameters
    ----------
    expnum : integer expnum 
    ccdnum : integer ccdnum (1-62)
    
    Returns
    -------
    uid   : float constructed from expnum + ccdnum/100.
    """
    return expnum + ccdnum/100.

def ruid(uid):
    """Extract EXPNUM, CCDNUM from unique float
    
    Parameters
    ----------
    uid   : float constructed from expnum + ccdnum/100.
    
    Returns
    -------
    expnum : integer expnum 
    ccdnum : integer ccdnum
    """
    uid = np.asarray(uid)
    expnum = np.floor(uid).astype(int)
    ccdnum = np.round(100 * (uid%1)).astype(int)
    return expnum,ccdnum

def file2pix(filenames):
    filenames = np.atleast_1d(filenames).astype(str)
    basenames = np.char.rpartition(filenames,'/')[:,-1]
    pixels = np.char.rpartition(np.char.partition(basenames,'.')[:,0],'_')[:,-1]
    return pixels.astype(int)

def rename_column(data, old, new):
    """
    Rename a column from old to new.
    """
    names = data.dtype.names

def check_formula(formula):
    """ Check that a formula is valid. """
    if not (('data[' in formula) and (']' in formula)):
        msg = 'Invalid formula:\n %s'%formula
        raise ValueError(msg)

def parse_formula(formula):
    """ Return the columns used in a formula. """
    check_formula(formula)

    columns = []
    for x in formula.split('data[')[1:]:
        col = x.split(']')[0].replace('"','').replace("'",'')
        columns += [col]

    return columns

def skim_data(infile,outfile,columns,select=None,force=False):
    """ Skim data from a FITS file. 

    ADW: Could this be replaced by a ftool?
    """
    logger.info("Running %s..."%infile)
    data,header=fitsio.read(infile,header=True,columns=columns)

    if select:
        check_formula(select)
        logger.debug("  Applying selection: %s"%select)
        data = data[eval(select)]

    if len(data) == 0:
        logger.debug("  No objects pass selection.")
    else:
        logger.debug("  Writing %s..."%outfile)
        write_fits(outfile,data,header=header,force=force)

def add_column(filename,column,formula,force=False):
    """ Add a column to a FITS file.

    ADW: Could this be replaced by a ftool?
    """
    columns = parse_formula(formula)
    logger.info("Running file: %s"%filename)
    logger.debug("  Reading columns: %s"%columns)
    data = fitsio.read(filename,columns=columns)

    logger.debug('  Evaluating formula: %s'%formula)
    col = eval(formula)

    col = np.asarray(col,dtype=[(column,col.dtype)])
    insert_columns(filename,col,force=force)
    return True

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

def get_local_catalog(ra,dec,radius,catalog,**kwargs):
    """ Get local catalog from disk 

    Parameters
    ----------
    ra : right ascension of catalog pixel (deg)
    dec: declination of catalog pixel (deg)
    catalog : catalog name
    
    Returns
    -------
    cat : object catalog
    """
    from ugali.utils.healpix import ang2disc
    name = catalog.lower()
    mapping = {}
    if name in ('gaia_dr1'):
        nside=32
        dirname = '/data/des40.b/data/gaia/dr1/healpix'
        basename = 'GaiaSource_%05d.fits'
        columns = ['SOURCE_ID','RA','DEC']
        #mapping.update({'RA':'_RAJ2000', 'DEC':'_DEJ2000'})
    elif name in ('gaia','gaia_dr2'): # Default...
        nside=32
        dirname = '/data/des40.b/data/gaia/dr2/healpix'
        basename = 'GaiaSource_%05d.fits'
        columns = ['SOURCE_ID','RA','DEC']
        #mapping.update({'RA':'_RAJ2000', 'DEC':'_DEJ2000'})
    elif name in ('gaia_edr3'):
        nside=32
        dirname = '/data/des40.b/data/gaia/edr3/healpix'
        basename = 'GaiaSource_%05d.fits'
        columns = ['SOURCE_ID','RA','DEC']
        #mapping.update({'RA':'_RAJ2000', 'DEC':'_DEJ2000'})
    elif 'atlas' in name:
        nside = 32
        dirname = '/data/des40.b/data/atlas-refcat2/healpix'
        basename = 'atlas-refcat2_%05d.fits'
        columns = ['OBJID','RA','DEC','G']
        mapping.update({})
    elif 'ps1' in name:
        nside = 32
        dirname = '/data/des40.b/data/ps1/dr1/healpix'
        basename = 'ps1_dr1_%05d.fits'
        columns = ['RA','DEC']
    elif 'decals' in name:
        nside = 32
        dirname = '/data/des40.b/data/decals/dr8/south_healpix'
        basename = 'decals-dr8-sweep_%05d.fits'
        columns = ['OBJID','RA','DEC']
    else:
        raise Exception('Unrecognized catalog: %s'%catalog)

    pixels = ang2disc(nside, ra, dec, radius, inclusive=True)
    filenames = [os.path.join(dirname,basename%p) for p in pixels]
    filenames = [f for f in filenames if os.path.exists(f)]
    cat = load_infiles(filenames,columns=columns)

    names = [mapping.get(n,n) for n in cat.dtype.names]
    cat.dtype.names = names
    return cat
        
def print_problem(msg):
    import termcolor as color
    print(color(msg,'red'))

def print_statistics(hpxmap):
    q = [5,16,50,84,95]
    p = np.nanpercentile(hpxmap[np.isfinite(hpxmap)],q)
    s68 = (p[3]-p[1])
    std = np.nanstd(hpxmap[np.isfinite(hpxmap)])
        
    print("  [%.0f,%.0f,%.0f]%%: %0.1f/%0.1f/%0.1f"%(q[0],q[2],q[4],p[0],p[2],p[4]))
    print("  68%% Interval: %0.1f"%(s68))
    print("  STDEV: %0.1f"%(std))

    return p,s68,std

def set_memory_limit(mlimit):
    """Set the (soft) memory limit for setrlimit.

    Parameters:
    -----------
    mlimit : soft memory limit (bytes)
                                                                                 
    Returns:
    --------
    soft, hard : memory limits (bytes)
    """
    import resource
    rsrc = resource.RLIMIT_AS
    resource.setrlimit(rsrc, (mlimit, mlimit))
    return resource.getrlimit(rsrc)

def verbose(func):
    """
    Decorator for timing functions

    Parameter
    ---------
    func : function to wrap
    
    Returns
    -------
    wrapper : wrapped function
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.debug("Running %r..."%(func.__name__))
        t0=time.time()
        ret = func(*args,**kwargs)
        logger.debug('%4.2fs'%(time.time()-t0))
        return ret

    return wrapper

def ignore_warning(warning=None):
    """
    Decorator to ignore a given warning occurring during method execution.
    https://stackoverflow.com/a/70292317

    Parameters
    ----------
    warning : Warning type or None
    
    Returns
    -------
    inner : wrapped function
    """

    def inner(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            with warnings.catch_warnings():
                if warning: warnings.filterwarnings("ignore", category=warning)
                else: warnings.filterwarnings("ignore")
                return func(*args, **kwargs)

        return wrapper

    return inner

# Some default healpix maps
def empty(nside,dtype=int): 
    return np.zeros(hp.nside2npix(nside),dtype=dtype)

def blank(nside,dtype=float): 
    return np.nan*np.ones(hp.nside2npix(nside),dtype=dtype)

def unseen(nside,dtype=float): 
    return hp.UNSEEN*np.ones(hp.nside2npix(nside),dtype=dtype)

def unstructure(arr):
    """Convert structured array to unstructured array.

    Poor stop-gap for recfunctions.structured_to_unstructured 
    in numpy >= 1.16.0

    Parameters
    ----------
    arr : Structured array to convert

    Returns
    -------
    unstructured : Unstructured array
    """
    if NumpyVersion(np.__version__) >= '1.16.0':
        return rfn.structured_to_unstructured(arr)
    else:
        out = rfn.repack_fields(arr.copy())
        dtype = out.dtype[0].str
        return out.view(dtype).reshape(len(out),-1)

def rec_append_fields(rec, names, arrs, dtypes=None):
    """
    Re-implement mlab.rec_append_fields for speed.

    Return a new record array with field names populated with data
    from arrays in *arrs*.  If appending a single field, then *names*,
    *arrs* and *dtypes* do not have to be lists. They can just be the
    values themselves.

    Parameters
    ----------
    rec    : recarray
    names  : names of new fields
    arrs   : values of new fields
    dtypes : dtypes of new fields

    Returns
    -------
    rec : recarray with fields appended
    """
    if (not isstring(names) and iterable(names) and len(names) and isstring(names[0])):
        if len(names) != len(arrs):
            raise ValueError("number of arrays do not match number of names")
    else:  # we have only 1 name and 1 array
        names = [names]
        arrs = [arrs]
    arrs = list(map(np.asarray, arrs))
    if dtypes is None:
        dtypes = [a.dtype for a in arrs]
    elif not iterable(dtypes):
        dtypes = [dtypes]
    if len(arrs) != len(dtypes):
        if len(dtypes) == 1:
            dtypes = dtypes * len(arrs)
        else:
            raise ValueError("dtypes must be None, a single dtype or a list")
    old_dtypes = rec.dtype.descr
    if six.PY2:
        old_dtypes = [(name.encode('utf-8'), dt) for name, dt in old_dtypes]
    newdtype = np.dtype(old_dtypes + list(zip(names, dtypes)))
    newrec = np.recarray(rec.shape, dtype=newdtype)
    for field in rec.dtype.fields:
        newrec[field] = rec[field]
    for name, arr in zip(names, arrs):
        newrec[name] = arr
    return newrec
