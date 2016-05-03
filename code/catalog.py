#!/usr/bin/env python
import os, sys
import time
from collections import OrderedDict as odict
from functools import wraps
import copy
import fitsio
import numpy as np
import numpy.lib.recfunctions as recfuncs
import healpy

from const import OBJECT_ID, UNIQUE_ID, BANDS, NSIDES, MINBANDS
from const import BADMAG
import utils
from utils import bfields, load_infiles
from pixelize import ang2pix

import scipy.ndimage as nd
from ugali.utils.logger import logger
import ugali.utils.projector

# Input columns
BAND = 'BAND'
IDX = [OBJECT_ID,BAND] + ['FILENAME','REQNUM','ATTNUM','OBJECT_NUMBER','TAG']
COORDS = ['RA','DEC']
HEALPIX = ['HPX%d'%nside for nside in NSIDES]
MAGPSF = ['MAG_PSF','MAGERR_PSF']
MAGAUTO = ['MAG_AUTO','MAGERR_AUTO']
MAGS = MAGPSF + MAGAUTO
SPREAD = ['SPREAD_MODEL','SPREADERR_MODEL']
CLASS = ['CLASS_STAR']
#FLAGS = ['FLAGS','QSLR_FLAG'] # Y2Q1
FLAGS = ['FLAGS']
NEPOCHS = ['NEPOCHS']
TEFF = ['T_EFF']
EXPNUM = ['EXPNUM','CCDNUM']
EXTINCTION = ['EXTINCTION']

BEST = MAGS + SPREAD + CLASS + FLAGS + EXPNUM + TEFF + EXTINCTION
INPUT_COLS = IDX + [BAND] + COORDS + BEST

# Output columns
WAVGMAGPSF  = ['WAVG_'+f for f in MAGPSF + ['MAGRMS_PSF']]
WAVGMAGAUTO = ['WAVG_'+f for f in MAGAUTO + ['MAGRMS_AUTO']]
WAVGMAGS   = WAVGMAGPSF + WAVGMAGAUTO
WAVGSPREAD = ['WAVG_'+f for f in SPREAD + ['SPREADRMS_MODEL']]
#              + ['MIN_SPREAD_MODEL','MIN_SPREADERR_MODEL']
#              + ['MAX_SPREAD_MODEL','MAX_SPREADERR_MODEL']]
### # Y2Q1
#WAVGCLASS =  ['WAVG_'+f for f in CLASS + ['CLASSRMS_STAR','CLASSMIN_STAR','CLASSMAX_STAR']] # Y2Q1
# Y2Q2
WAVGCLASS =  []
WAVGFLAGS = ['WAVG_'+f for f in FLAGS] 

WAVG = WAVGMAGS + WAVGSPREAD + WAVGFLAGS

OUTPUT_COLS = odict(
    [(i,('i8',-1)) for i in [OBJECT_ID]]
    + [(c,('f8',np.nan)) for c in COORDS]
    + [(i,('i8',-1)) for i in HEALPIX]
    + [(f,('i2',0)) for f in bfields(NEPOCHS,BANDS)]
    + [(f,('f4',BADMAG)) for f in bfields(MAGS,BANDS)]
    + [(f,('f4',-1)) for f in bfields(SPREAD[:1],BANDS)]
    + [(f,('f4',1)) for f in bfields(SPREAD[1:],BANDS)]
    + [(f,('f4',-1)) for f in bfields(CLASS,BANDS)]
    + [(f,('i2',99)) for f in bfields(FLAGS,BANDS)]
    + [(f,('f4',BADMAG)) for f in bfields(WAVGMAGS,BANDS)]
    + [(f,('f4',-1)) for f in bfields(WAVGSPREAD[:1],BANDS)]
    + [(f,('f4',1)) for f in bfields(WAVGSPREAD[1:],BANDS)]
    #+ [(f,('f4',-1)) for f in bfields(WAVGCLASS,BANDS)] # Y2Q1
    + [(f,('i2',99)) for f in bfields(WAVGFLAGS,BANDS)]
    + [(f,('i4',-1)) for f in bfields(EXPNUM,BANDS)]
    + [(f,('f4',-1)) for f in bfields(TEFF,BANDS)]
    + [(f,('f4',-1)) for f in bfields(EXTINCTION,BANDS)]
    )

def verbose(func):
    """
    Decorator for timing functions
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.debug("Running %r..."%(func.__name__))
        t0=time.time()
        ret = func(*args,**kwargs)
        logger.debug('%4.2fs'%(time.time()-t0))
        return ret

    return wrapper

def maximum_index(input,labels,index=None):
    if index is None: index = np.unique(labels)
    argmax = nd.maximum_position(input,labels,index)
    return np.asarray(argmax,dtype=int).reshape(len(argmax))

def minimum_index(input,labels,index=None):
    if index is None: index = np.unique(labels)
    argmin = nd.minimum_position(input,labels,index)
    return np.asarray(argmin,dtype=int).reshape(len(argmin))

def wavg(values,weights,labels=None,index=None,inverse=None,counts=None):
    # Weighted average:
    # http://en.wikipedia.org/wiki/Weighted_arithmetic_mean
    if labels is None: 
        labels = np.ones(len(values),dtype=int)

    if any(x is None for x in [index,inverse,counts]):
        index,inverse,counts = np.unique(labels,False,True,True)

    V1 = nd.sum(weights,labels=labels,index=index)
    V2 = nd.sum(weights**2,labels=labels,index=index)    

    wavg = nd.sum(weights*values,labels=labels,index=index)/V1

    wavg_err  = nd.maximum(np.sqrt(1/weights),labels=labels,index=index)

    wavg_rms  = nd.sum(weights*(values - wavg[inverse])**2,
                       labels=labels,index=index)

    err = np.seterr()
    np.seterr(all='ignore')
    wavg_rms *= (V1/(V1**2 - V2)) 
    wavg_rms  = np.sqrt(wavg_rms)
    np.seterr(**err)

    return wavg, np.where(counts > 1,wavg_rms,wavg_err)


@verbose
def coadd_coords(lon,lat,labels=None,index=None):
    """
    Calculated the median coordinates of each object.

    Returns
      out : Dictionary of average coordinates.
    """
    if labels is None: 
        labels = np.ones(len(lon),dtype=int)

    if index is None: index = np.unique(labels)

    x,y,z = healpy.rotator.dir2vec(lon,lat,lonlat=True)

    # Mean coordinates
    #x_out = nd.mean(x,labels=labels,index=index)
    #y_out = nd.mean(y,labels=labels,index=index)
    #z_out = nd.mean(z,labels=labels,index=index)

    # Median coordinates
    x_out = nd.median(x,labels=labels,index=index)
    y_out = nd.median(y,labels=labels,index=index)
    z_out = nd.median(z,labels=labels,index=index)

    # std coordinates
    #x_std = nd.standard_deviation(x,labels=labels,index=index)
    #y_std = nd.standard_deviation(y,labels=labels,index=index)
    #z_std = nd.standard_deviation(z,labels=labels,index=index)
    
    lon_out, lat_out = healpy.rotator.vec2dir(x_out,y_out,z_out,lonlat=True)

    # What about adding a dispersion?

    return lon_out % 360.,lat_out

@verbose
def coadd_healpix(lon,lat,nsides=NSIDES,nest=True):
    pix = [ang2pix(nside,lon,lat,nest=nest) for nside in NSIDES]
    return pix

@verbose
def coadd_mag(mag,magerr,labels,index=None,inverse=None,counts=None):
    if any(x is None for x in [index,inverse,counts]):
        index,inverse,counts = np.unique(labels,False,True,True)

    weights=1./magerr**2
    wavg_mag,wavg_magrms = wavg(mag,weights,labels,index,inverse,counts)

    # Magnitude errors in quadrature
    wavg_magerr = np.sqrt(nd.sum(magerr**2,labels=labels,index=index))/counts

    return wavg_mag,wavg_magerr,wavg_magrms

@verbose
def coadd_spread(spread,spreaderr,labels,index=None,inverse=None,counts=None):
    if any(x is None for x in [index,inverse,counts]):
        index,inverse,counts = np.unique(labels,False,True,True)

    weights=1./spreaderr**2

    wavg_spread,wavg_spreadrms = wavg(spread,weights,labels,index,inverse,counts)

    # Add spreaderr in quadrature
    wavg_spreaderr = np.sqrt(nd.sum(spreaderr**2,labels=labels,index=index))/counts
    ret = [wavg_spread,wavg_spreaderr,wavg_spreadrms]
    
    ### # extrema is only a little faster than positions individually
    ### extrema = nd.extrema(spreaderr,labels=labels,index=index)
    ### spread_min_idx = np.asarray(extrema[2],dtype=int)
    ### spread_max_idx = np.asarray(extrema[3],dtype=int)
    ###  
    ### min_spread = spread[spread_min_idx]
    ### min_spreaderr = extrema[0]
    ### ret += [min_spread,min_spreaderr]
    ### max_spread = spread[spread_max_idx]
    ### max_spreaderr = extrema[1]
    ### ret += [max_spread,max_spreaderr]
    return ret 

@verbose
def coadd_class(class_star,labels,index=None):
    if index is None: index = np.unique(labels)
    
    class_avg = nd.mean(class_star,labels=labels,index=index)
    class_std = nd.standard_deviation(class_star,labels=labels,index=index)
    class_min = nd.minimum(class_star,labels=labels,index=index)
    class_max = nd.maximum(class_star,labels=labels,index=index)
    return class_avg, class_std, class_min, class_max

@verbose
def coadd_flags(flags,labels,index=None):
    if index is None: index = np.unique(labels)

    # Add some checks if flags can't be cast to np.uint8
    # There is no numpy function for doing this...
    def array_or(flags):
        bin = flags.astype(dtype=np.uint8).reshape(len(flags),-1)
        bits = np.unpackbits(bin,axis=1)
        return np.packbits(np.any(bits,axis=0).astype(int))[0]
    
    flags_or = nd.labeled_comprehension(flags,labels,index,array_or,out_dtype=int,default=0)
    return flags_or

@verbose
def best_values(values,best,labels,index=None):
    if index is None: index = np.unique(labels)
    
    #argmax = nd.maximum_position(best,labels,index)
    argmax = maximum_index(best,labels,index)
    try:
        return [values[f][argmax] for f in values.dtype.names]
    except TypeError:
        return values[argmax]

def coadd_objects(data):
    """
    Create a unique object catalog.
    """
    unique_ids = np.unique(data[OBJECT_ID])
    nobjs = len(unique_ids)

    dtype = [(k,v[0]) for k,v in OUTPUT_COLS.items()]
    cat = np.recarray(nobjs,dtype=dtype)
    for k,v in OUTPUT_COLS.items():
        cat[k] = v[1]
    cat[OBJECT_ID] = unique_ids

    keys = copy.copy(data[IDX+EXPNUM])
    keys = recfuncs.rec_append_fields(keys,UNIQUE_ID,-np.ones(len(keys)),dtypes='>i8')
    # OBJECT_NUMBER has a different meaning (and type) in Y1A1 and Y2N.
    # Standardize it here (wouldn't be necessary if done on download).
    # ADW: Is this working properly?
    if not keys.dtype['OBJECT_NUMBER'] is not np.dtype('>i8'):
        keys['OBJECT_NUMBER'] = keys['OBJECT_NUMBER'].astype('>i8')

    x = coadd_coords(data['RA'],data['DEC'],data[OBJECT_ID],index=unique_ids)
    for i,f in enumerate(COORDS):
        cat[f] = x[i]

    x = coadd_healpix(cat[COORDS[0]],cat[COORDS[1]],NSIDES,nest=True)
    for i,f in enumerate(HEALPIX):
        cat[f] = x[i]

    for b in BANDS:
        logger.info("=== Creating %s-band catalog ==="%b)
        # This is an annoying feature of fitsio...
        band = '{0:<{1}}'.format(b,len(data[BAND][0]))

        sel = np.where(data[BAND]==band)[0]
        if len(sel) == 0:
            logger.warning("No objects found in %s-band."%b)
            continue

        d = data[sel]

        labels = d[OBJECT_ID]
        index,inverse,counts = np.unique(labels,False,True,True)
        
        # Match unique indices of this selection to the catalog indices
        idx = np.searchsorted(cat[OBJECT_ID],index)

        # Set the unique keys
        i = best_values(sel,d['T_EFF'],labels,index)
        keys[UNIQUE_ID][i] = keys[OBJECT_ID][i]
        
        for f in NEPOCHS:
            cat[bfields(f,b)][idx] = counts

        x = best_values(d[BEST],d['T_EFF'],labels,index)
        for i,f in enumerate(BEST):
            cat[bfields(f,b)][idx] = x[i]

        x = coadd_mag(d[MAGPSF[0]],d[MAGPSF[1]],labels,index,inverse,counts)
        for i,f in enumerate(WAVGMAGPSF):
            cat[bfields(f,b)][idx] = x[i]

        x = coadd_mag(d[MAGAUTO[0]],d[MAGAUTO[1]],labels,index,inverse,counts)
        for i,f in enumerate(WAVGMAGAUTO):
            cat[bfields(f,b)][idx] = x[i]

        x = coadd_spread(d[SPREAD[0]],d[SPREAD[1]],labels,index,inverse,counts)
        for i,f in enumerate(WAVGSPREAD):
            cat[bfields(f,b)][idx] = x[i]

        try:
            x = coadd_class(d['CLASS_STAR'],labels,index)
            for i,f in enumerate(WAVGCLASS):
                cat[bfields(f,b)][idx] = x[i]
        except ValueError,msg:
            logger.warning(msg)

        for i,(f,w) in enumerate(zip(FLAGS,WAVGFLAGS)):
            x = coadd_flags(d[f],labels,index)
            cat[bfields(w,b)][idx] = x
        
    return cat,keys

def good_objects(data):
    """
    Remove spurious single epoch detections.

    Cuts are made on:
    - 0 < MAG[ERR] < 99
    - -99 < SPREAD_MODEL < 99
    - 0 < SPREADERR_MODEL < 99
    - FLAGS < 4

    Parameters:
    -----------
    data : Input single epoch catalog

    Returns:
    --------
    data : Single epoch catalog with bad objecs removed.
    """
    nobjs = len(data)
    sel = np.ones(nobjs,dtype=bool)

    # Only calibrated objects with 0 < MAG < BADMAG where 
    # MAG = MAG[ERR]_[PSF/AUTO] 
    # (There are a small number of spurious measurements with MAG > BADMAG)
    mags = data[MAGS]
    dtype = mags.dtype[0].str
    sel &= np.all(mags.view(dtype).reshape((nobjs,-1))<BADMAG,axis=1)
    sel &= np.all(mags.view(dtype).reshape((nobjs,-1))>0,axis=1)


    # Objects with valid spread_model and spreaderr values
    # 99 is chosen as the max of NUMBER(7,5) data type
    sel &= (data['SPREAD_MODEL']>-99) & (data['SPREAD_MODEL']<99)
    sel &= (data['SPREADERR_MODEL']>0) & (data['SPREADERR_MODEL']<99)

    # Only objects without bad SExtractor flags
    sel &= (data['FLAGS'] < 4)
    #flags = data[['FLAGS']]
    #dtype = flags.dtype[0].str
    #sel &= np.all(flags.view(dtype).reshape((nobjs,-1)) < 4,axis=1)

    return data[sel]

def quality_cuts(cat,key=None):
    """
    Prune down the catalog by selecting on detection band, etc.
    """
    nobjs = len(cat)
    sel = np.zeros(nobjs,dtype=bool)

    # Objects with detections in the minimum number of bands.
    # Careful with this, the 'view' can be dangerous...
    columns = bfields(['MAG_PSF'],BANDS)
    dtype = cat.dtype[columns[0]]
    mags = cat[columns].view(dtype).reshape((cat.size,-1))
    sel |= (np.sum(mags < BADMAG,axis=1) >= MINBANDS)
    
    """
    # Objects with r,i,z
    mags = cat[bfields(['MAG_PSF','MAG_AUTO'],BANDS)].view(np.float).reshape((cat.size,-1))
    sel = np.any(mags!=BADMAG,axis=1)

    # Objects with g,r
    mags = cat[bfields(['MAG_PSF','MAG_AUTO'],['g','r'])].view(np.float).reshape((cat.size,-1))
    sel = np.all(mags!=BADMAG,axis=1)

    # Objects with any of r,i,z
    mags = cat[bfields(['MAG_PSF','MAG_AUTO'],['r','i','z'])].view(np.float).reshape((cat.size,-1))
    sel = np.any(mags!=BADMAG,axis=1)
    """

    if key is None:
        return cat[sel]
    else: 
        # Most objects accepted; key selection on the omitted IDs for speed
        ksel = ~np.in1d(key[OBJECT_ID],cat[OBJECT_ID][~sel])
        return cat[sel],key[ksel]

def check_keys(cat,key):
    """
    Check the number of non-unique objects against the number of
    measured magnitudes.
    """
    nkey = (keys['UNIQUE_ID'] >0).sum()

    columns = bfields(['MAG_PSF'],BANDS)
    dtype = cat.dtype[columns[0]]
    ncat = (cat[columns].view(dtype).reshape(len(cat),-1) < 90).sum()

    if ncat != nkey:
        msg = "Number of keys does not match catalog"
        raise ValueError(msg)

    return


if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infiles',nargs='+')
    parser.add_argument('-o','--outfile')
    parser.add_argument('-k','--keyfile')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('-b','--bands',default=None,action='append')
    parser.add_argument('--min-bands',default=None,type=int)
    opts = parser.parse_args()

    if vars(opts).get('verbose'): logger.setLevel(logger.DEBUG)
    if opts.bands: BANDS = opts.bands
    if opts.min_bands: MINBANDS = opts.min_bands

    if os.path.exists(opts.outfile) and not opts.force:
        logger.warning("Found %s; skipping..."%opts.outfile)
        sys.exit()

    logger.info("Loading files: %s"%opts.infiles)
    data = load_infiles(opts.infiles,INPUT_COLS)
    logger.info("All objects: %i"%len(data))

    good = good_objects(data)
    logger.info("Good objects: %i"%len(good))

    if len(good) == 0:
        logger.warning("No good objects found; exiting...")
        sys.exit()
        
    cat,key = coadd_objects(good)
    logger.info("Unique objects: %i"%len(cat))
 
    catalog,keys = quality_cuts(cat,key)
    check_keys(catalog,keys)
    logger.info("Quality objects: %i"%len(catalog))
    #print catalog[bfields(['MAG_PSF'],BANDS)]

    if opts.outfile and len(catalog):
        logger.info("Writing %s..."%opts.outfile)
        utils.write(opts.outfile,catalog,force=opts.force)

    if opts.keyfile and len(keys):
        logger.info("Writing %s..."%opts.keyfile)
        utils.write(opts.keyfile,keys,force=opts.force)
