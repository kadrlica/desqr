#!/usr/bin/env python
"""
Spatial matching algorithms.
"""
import os, sys
from collections import OrderedDict as odict
import gc

import fitsio
import numpy as np
import pandas as pd
import healpy as hp
import scipy
from scipy.spatial import cKDTree
import scipy.ndimage as nd

try:
    from matplotlib.mlab import rec_append_fields
except ImportError:
    from ugali.utils.mlab import rec_append_fields

from ugali.utils.projector import angsep

from desqr.utils import set_memory_limit, insert_columns
from desqr.const import ZEROSTR,OBJECT_ID
from desqr.split import split_qcat
from desqr.logger import logger

MATCHCOLS = ['RA','DEC','EXPNUM']

### def projector2(lon,lat):
###     from ugali.utils.projector import SphericalRotator
###     rotator = SphericalRotator(0,0)
###  
###     lon = np.asarray(lon)
###     lat = np.asarray(lat)
###  
###     x, y, z = rotator.cartesian(lon.ravel(),lat.ravel())
###     coords = np.empty((x.size, 3))
###     coords[:, 0] = x
###     coords[:, 1] = y
###     coords[:, 2] = z
###     return coords


def ang2vec(lon, lat):
    """Convert longitude and latitude into array of unit vectors.

    Parameters
    ----------
    lon : longitude (deg)
    lat : latitude (deg)

    Returns
    -------
    vec : 2D array of vectors
    """
    lonr = np.deg2rad(lon)
    latr = np.deg2rad(lat)
    coslat = np.cos(latr)
    return np.stack([
        np.atleast_1d(np.cos(lonr) * coslat),
        np.atleast_1d(np.sin(lonr) * coslat),
        np.atleast_1d(np.sin(latr)),
    ], axis=-1)

def projector(lon, lat):
    return ang2vec(lon,lat)

def centroid(lon,lat,stat='median',labels=None,index=None):
    if labels is None: 
        labels = np.ones(len(lon),dtype=int)

    if index is None: index = np.unique(labels)

    x,y,z = hp.rotator.dir2vec(lon,lat,lonlat=True)

    if stat == 'mean':
        x_out = nd.mean(x,labels=labels,index=index)
        y_out = nd.mean(y,labels=labels,index=index)
        z_out = nd.mean(z,labels=labels,index=index)
    elif stat == 'median':
        x_out = nd.median(x,labels=labels,index=index)
        y_out = nd.median(y,labels=labels,index=index)
        z_out = nd.median(z,labels=labels,index=index)
    else:
        msg = "Unrecognized stat: %s"%stat
        raise Exception(msg)

    lon_out, lat_out = hp.rotator.vec2dir(x_out,y_out,z_out,lonlat=True)

    return lon_out % 360.,lat_out


def match_to_self(lon, lat, radius, min_match=1, min_obs=1):
    """Match a list of lon/lat positions to itself.

    Parameters
    ----------
    lon : `np.ndarray`
        Array of longitude positions (degrees)
    lat : `np.ndarray`
        Array of latitude positions (degrees)
    radius : `float`
        Radius to match (degrees)
    min_obs : `int`, optional
        Minimum number of matches to count.

    Returns
    -------
    idx : `list` [`list` [`int`]]
        List of match indices for each object.
        Sorted from most matches to fewest, duplicates removed.
    """
    coords = ang2vec(lon, lat)
    tree = cKDTree(coords)
    angle = 2.0*np.sin(np.deg2rad(radius)/2.0)
    idx = tree.query_ball_tree(tree, angle, eps=0.0)
    # We start with the most matches, so we need this sorted
    len_arr = np.array([len(j) for j in idx])
    st = np.argsort(len_arr)[:: -1]
    match_id = np.full(len(idx), -1, dtype=np.int32)
    idx2 = []
    for j in st:
        if match_id[j] < 0 and len_arr[j] >= min_obs:
            match_id[idx[j]] = j
            idx2.append(idx[j])
    return idx2

def match_to_other(lon1, lat1, radius, lon2, lat2):
    """Match one list of lon/lat positions to another list of lon/lat positions.

    Parameters
    ----------
    lon1 : `np.ndarray`
        Array of longitude positions for list 1 (degrees)
    lat1 : `np.ndarray`
        Array of latitude positions for list 1 (degrees)
    radius : `float`
        Radius of matches (degrees)
    lon2 : `np.ndarray`
        Array of longitude positions for list 2 (degrees)
    lat2 : `np.ndarray`
        Array of latitude positions for list 2 (degrees)

    Returns
    -------
    idx : `list` [`list` [`int`]]
        Each row in idx corresponds to each position in lon1/lat1.
        The indices in the row correspond to the indices in lon2/lat2.
    """
    coords1 = ang2vec(lon1, lat1)
    coords2 = ang2vec(lon2, lat2)
    tree1 = cKDTree(coords1)
    # The second tree in the match does not need to be balanced, and
    # turning this off yields significantly faster runtime.
    tree2 = cKDTree(coords2, balanced_tree=False)
    angle = 2.0*np.sin(np.deg2rad(radius)/2.0)
    idx = tree1.query_ball_tree(tree2, angle, eps=0.0)
    return idx

def match_query(lon1,lat1,lon2,lat2,eps=0.01,n_jobs=1):
    """
    Perform a KDTree match after transforming to Cartesian coordinates.

    Parameters
    ----------
    lon1 : longitude from first array (deg)
    lat1 : latitude from first array (deg)
    lon2 : longitude from second array (deg)
    lat2 : latitude from second array (deg)
    eps  : precision (see `cKDTree.query`)
    n_jobs : Number of workers to use for processing (see `cKDTree.query`)

    Returns:
    --------
    idx1 : index into the first array
    idx2 : index into the second array
    ds   : angular seperation (deg)
    """
    coords1 = ang2vec(lon1, lat1)
    coords2 = ang2vec(lon2, lat2)
 
    tree = cKDTree(coords2)
    idx1 = np.arange(lon1.size) 

    if tuple(int(x) for x in scipy.__version__.split('.')[:2]) >= (1, 4):
        idx2 = tree.query(coords1, eps=eps, workers=n_jobs)[1]
    else:
        idx2 = tree.query(coords1, eps=eps, n_jobs=n_jobs)[1]

    ds = angsep(lon1, lat1, lon2[idx2], lat2[idx2])
    return np.atleast_1d(idx1), np.atleast_1d(idx2), np.atleast_1d(ds)

def match_ball_tree(lon,lat,radius=1.0):
    """ 
    Internal catalog match. Finds closest match to each object other
    than the object itself.

    Parameters
    ----------
    lon : longitude (deg)
    lat : latitude (deg)
    radius : match radius (arcsec)
    
    Returns
    -------
    match_id : 
    """
    if len(lon) != len(lat):
        msg = "Input lon and lat do not match"
        raise Exception(msg)
    if np.any(np.isnan(lon)):
        msg = "Invalid value found in lon: nan"
        raise ValueError(msg)

    nobjs = len(lon)
    logger.info("Matching %i objects."%nobjs)
     
    # First iteration...
    coords = ang2vec(lon, lat)
    #dist = (radius/3600.)*(np.pi/180.)
    dist = 2.0*np.sin(np.deg2rad(radius/3600.)/2.0)
    tree = cKDTree(coords)
    logger.info("Querying ball tree with radius %s arcsec..."%radius)
    idx = tree.query_ball_tree(tree,dist,eps=0.01)
    # Because memory...
    del tree, coords
    gc.collect()

    # This could probably be improved...
    logger.info("Filling matched IDs.")
    match_id = -1*np.ones(nobjs,dtype=int)
    # Ordering has a +/-1% effect on the number of matches
    # Starting with most match minimizes the number of unique objects
    # [831713 < 839700 < 847271]
    #np.arange(len(data)):
    #np.argsort([len(i) for i in idx]):
    for i in np.argsort([len(i) for i in idx])[::-1]:
        if match_id[i] >= 0: continue
        match_id[idx[i]] = i

    uid,inv = np.unique(match_id,return_inverse=True)
    match_id[:] = np.arange(len(uid))[inv]
    logger.info("Found %i unique objects."%len(uid))

    return match_id

def match_multi_stage(lon,lat,radius=1.0):
    """
    Multi-stage matching (not necessary for radius <= 0.5):
    1) Use match_ball_tree to self-match raw catalog
    2) Find the median location of the matched objects
    3) Use match_query to match between objects and medians
    4) Use match_ball_tree to self-match objects far from medians

    Parameters:
    -----------
    lon    : longitude coordinate (deg)
    lat    : latitude coordinate (deg)
    radius : matching radius (arcmin)
    """
    if len(lon) != len(lat):
        msg = "Input lon and lat do not match"
        raise Exception(msg)

    nobjs = len(lon)
    match_id = -1 * np.ones(nobjs,dtype=int)

    logger.info("Self-matching with ball tree...")
    ball1_id = match_ball_tree(lon,lat,radius)
    uid = np.unique(ball1_id)

    logger.info("Calculating match coordinates...")
    median_lon,median_lat = centroid(lon,lat,stat='median',labels=ball1_id,index=uid)

    logger.info("Matching against matches...")
    idx1,idx2,sep = match_query(lon,lat,median_lon,median_lat)
    matched = np.where(sep <= radius/3600.)[0]
    unmatched = np.where(sep > radius/3600.)[0]

    logger.info("Found %i matched matches."%len(matched))
    logger.info("Found %i unmatched matches."%len(unmatched))

    match_id[idx1[matched]] = uid[idx2[matched]]

    logger.info("Matching unmatched objects...")
    ball2_id = match_ball_tree(lon[unmatched],lat[unmatched],radius)
    ball2_id += (match_id.max()+1)
    match_id[idx1[unmatched]] = ball2_id

    if np.any(match_id == -1):
        msg = "Unmatched objects"
        raise Exception(msg)

    # Make sure the unique ids sequential
    uid,inv = np.unique(match_id,return_inverse=True)
    unique_id = np.arange(len(uid),dtype=int)[inv]
    return unique_id

def match_htm(data,radius=1.0):
    """ NOT IMPLEMENTED """
    logger.warning("'match_htm' not implemented")
    if True: return

    import esutil.htm
    htm = esutil.htm.HTM()
    kwargs = dict(radius=radius,maxmatch=-1)
    m = htm.match(data['RA'],data['DEC'],data['RA'],data['DEC'],**kwargs)
        
def match_exposures(data,radius=1.0):
    """ A matching algorithm that does not associate objects on the
    same exposure with each other. This algorithm ended up being too
    slow.

    Parameters:
    -----------
    data   : input data containing 'RA' and 'DEC'
    radius : matching radius (arcsec)

    Returns:
    --------
    match_id : the object matching id
    """

    expnums,counts = np.unique(data['EXPNUM'],return_counts=True)
    nobjs = len(data)
    nexp = len(expnums)
    logger.info("Found %i objects in %i exposures"%(nobjs,nexp))
    match_id = -1*np.ones(nobjs,dtype=int)

    for i,expnum in enumerate(expnums[counts.argsort()[::-1]]):
        exp_idx = np.where(data['EXPNUM'] == expnum)[0]

        if i == 0: 
            unique_idx = exp_idx
            match_id[unique_idx] = np.arange(len(unique_idx),dtype=int)

        m = match(data['RA'][exp_idx],data['DEC'][exp_idx],
                  data['RA'][unique_idx],data['DEC'][unique_idx])
        exp_match_idx,unique_match_idx,sep = m

        matched = np.where(sep <= radius/3600.)[0]
        unmatched = np.where(sep > radius/3600.)[0]
        logger.info("%i: EXPNUM=%i NOBJS=%i NMATCH=%i"%(i,expnum,len(exp_idx),len(matched)))
        # Assign the existing match_id to matched objects
        ii = unique_idx[unique_match_idx[matched]]
        jj = exp_idx[exp_match_idx[matched]]
        match_id[jj] = match_id[ii]

        # Create a new match_id for unmatched objects
        new_id = np.arange(len(unmatched),dtype=int) + np.max(match_id) + 1
        kk = exp_idx[exp_match_idx[unmatched]]
        match_id[kk] = new_id

        # Add the unmatched indices to list of unique indexes
        unique_idx = np.hstack([unique_idx,exp_idx[unmatched]])

    if np.any(match_id < 0):
        msg = "Unmatched objects found."
        raise Exception(msg)

    return match_id 

def split(data,match_id,objid=OBJECT_ID):
    """Split matched objects into pairs if consistently separated.
    
    Wraps around the `split_qcat` function from Eric Neilsen.
    
    Parameters:
    -----------
    data     : the input record array
    match_id : array of object match ID

    Returns:
    --------
    rec       : recarray of [OBJECT_ID,'PARENT_ID','SUB_ID','SPLIT_FLAG']
    """
    # Columns used by the split
    parid,subid,flag = ['PARENT_ID','SUB_ID','SPLIT_FLAG']
    dtype = [(objid,int),(parid,int),(subid,'i2'),(flag,'S15')]

    # Create a DataFrame and run the splitting
    df = pd.DataFrame(data.byteswap().newbyteorder())
    df[objid] = match_id
    df[parid] = match_id

    split = split_qcat(df)[[objid,parid,subid,flag]]

    if len(split) != len(df):
        msg = "Length of split objects differs from input."
        raise Exception(msg)
        
    # Increment the match_id for split objects
    objid_max = split[objid].max()
    subid_sel = (split[flag] == 'split') & (split[subid] > 0)
    # Assign new objid to each unique objid, subid pair
    group_id = split.loc[subid_sel].groupby([objid,subid]).ngroup()
    split.loc[group_id.index,objid] = objid_max + 1 + group_id.values

    # Convert back to recarray 
    rec = np.empty(len(split),dtype=dtype)
    for n in rec.dtype.names:
        rec[n] = split[n]

    # Number of split objects and new objects
    ndet = (split[flag] == 'split').sum()
    nobj = len(np.unique(group_id))
    logger.info("Split %i detections into %i new objects."%(ndet,nobj))

    return rec

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infiles',nargs='+')
    parser.add_argument('-r','--radius',default=1.0,type=float,
                        help='matching radius')
    parser.add_argument('-m','--mlimit',default=None,type=float,
                        help='memory limit (GB)')
    parser.add_argument('-o','--objid',default=OBJECT_ID,
                        help='name of output object identifier')
    parser.add_argument('-f','--force',action='store_true',
                        help='overwrite output columns if they exist')
    parser.add_argument('-s','--split',action='store_true',
                        help='split double objects')
    parser.add_argument('-v','--verbose',action='store_true',
                        help='output verbosity')
    args = parser.parse_args()
    
    if args.verbose: logger.setLevel(logger.DEBUG)
    
    if args.mlimit: 
        logger.info("Setting memory limit: %.1fGB"%(args.mlimit))
        soft,hard = set_memory_limit(args.mlimit*1024**3)
        logger.info("Memory limit: %.1fGB"%(soft/1024.**3))

    logger.info("Matching files: %s"%args.infiles)
    radius = args.radius
    fileidx = odict()

    data = []
    for i,f in enumerate(args.infiles):
        d,hdr = fitsio.read(f,header=True,columns=MATCHCOLS)
        if i == 0:
            imin = 0
            data = d
            pix = hdr.get('HPX')
        else:
            imin = len(data)
            data = np.append(data,d)
            if hdr.get('HPX') != pix:
                logger.warning('HEALPix pixels do not match')
        imax = len(data) 
        fileidx[f]=slice(imin,imax)

    zero_id = int(ZEROSTR%(pix,0))

    if len(data) == 1:
        match_id = np.array([0])
    else:
        match_id = match_ball_tree(data['RA'],data['DEC'],radius=args.radius)
        #match_id = match_multi_stage(data['RA'],data['DEC'],radius=args.radius)
    match_id += zero_id

    if args.split:
        # Run the split
        logger.info("Splitting objects...")
        out = split(data,match_id,objid=args.objid)
    else:
        # Just use the match_id
        out = np.rec.array(match_id,dtype=[(args.objid,match_id.dtype)],copy=False)

    # Write out columns
    for f,idx in fileidx.items():
        logger.info("Inserting column(s) into %s..."%f)
        insert_columns(f,out[idx],force=args.force)
     
