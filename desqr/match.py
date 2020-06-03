#!/usr/bin/env python
"""
Spatial matching algorithms.
"""
import os, sys
import fitsio
import numpy as np
import pandas as pd
from collections import OrderedDict as odict
import logging
import gc

import healpy
from scipy.spatial import cKDTree
import scipy.ndimage as nd
from matplotlib.mlab import rec_append_fields, rec_drop_fields

from ugali.utils.logger import logger
import ugali.utils.projector as proj

import utils
#from utils import logger
from const import ZEROSTR,OBJECT_ID
from split import split_qcat

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

def projector(lon,lat):
    return healpy.rotator.dir2vec(lon,lat,lonlat=True).T

def centroid(lon,lat,stat='median',labels=None,index=None):
    if labels is None: 
        labels = np.ones(len(lon),dtype=int)

    if index is None: index = np.unique(labels)

    x,y,z = healpy.rotator.dir2vec(lon,lat,lonlat=True)

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

    lon_out, lat_out = healpy.rotator.vec2dir(x_out,y_out,z_out,lonlat=True)

    return lon_out % 360.,lat_out

def match_query(lon1,lat1,lon2,lat2,eps=0.01):
    """
    Perform a KDTree match after transforming to Cartesian coordinates.

    Parameters
    ----------
    lon1 : longitude from first array (deg)
    lat1 : latitude from first array (deg)
    lon2 : longitude from second array (deg)
    lat2 : latitude from second array (deg)
    eps : precision (see `cKDTree.query`)

    Returns:
    --------
    idx1 : index into the first array
    idx2 : index into the second array
    ds   : angular seperation (deg)
    """
    coords1 = projector(lon1,lat1)
    coords2 = projector(lon2,lat2)
 
    tree = cKDTree(coords2)
    idx1 = np.arange(lon1.size) 
    idx2 = tree.query(coords1,eps=eps)[1]

    ds = proj.angsep(lon1, lat1, lon2[idx2], lat2[idx2])
    return idx1, idx2, ds

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
    coords = projector(lon,lat)
    dist = (radius/3600.)*(np.pi/180.)
    tree = cKDTree(coords)
    logger.info("Querying ball tree with radius %s arcsec..."%radius)
    idx = tree.query_ball_tree(tree,dist,eps=0.01)
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

def split(data,match_id):
    """Split matched objects into pairs if consistently separated.
    
    Parameters:
    -----------
    data     : the input record array
    match_id : the object match ID

    Returns:
    --------
    df       : a DataFrame containing 'QUICK_OBJECT_ID','SUB_ID', and 'SPLIT_FLAG'
    """
    logger.info("Splitting objects...")

    # Columns used by the split
    objid,subid,flag = ['QUICK_OBJECT_ID','SUB_ID','SPLIT_FLAG']
    data = data.byteswap().newbyteorder()

    # Create a DataFrame and run the splitting
    df=pd.DataFrame(rec_append_fields(data,objid,match_id))
    split = split_qcat(df)[[objid,subid,flag]]

    if len(split) != len(df):
        msg = "Length of split objects differs from input."
        raise Exception(msg)

    # Increment the match_id for split objects
    subsel = (split[subid] > 0)
    nsub = subsel.sum()
    split[objid][subsel] = split[objid].max() + np.arange(1,nsub+1)

    # Convert back to recarray 
    dtype = [(objid,int),(subid,'i2'),(flag,'S15')]
    rec = np.empty(len(split),dtype=dtype)
    for n in rec.dtype.names:
        rec[n] = split[n]

    logger.info("Split %i objects."%(nsub))

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
    parser.add_argument('-v','--verbose',action='store_true',
                        help='output verbosity')
    args = parser.parse_args()
    
    if args.verbose: logger.setLevel(logging.DEBUG)
    
    if args.mlimit: 
        logger.info("Setting memory limit: %.1fGB"%(args.mlimit))
        soft,hard = utils.set_memory_limit(args.mlimit*1024**3)
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

    # Run the split
    out = split(data,match_id)
    # Name the match_id column
    out.dtype.names = [args.objid] + list(out.dtype.names)[1:]

    ## Just use the match_id
    #out = np.rec.array(match_id,dtype=[(args.objid,match_id.dtype)],copy=False)

    # Write out columns
    for f,idx in fileidx.items():
        logger.info("Inserting column(s) into %s..."%f)
        utils.insert_columns(f,out[idx],force=args.force)


    ### uid,inv = np.unique(match_id,return_inverse=True)
    ### median_ra,median_dec = centroid(data['RA'],data['DEC'],stat='median',labels=match_id,index=uid)
    ### sep = proj.angsep(data['RA'],data['DEC'],median_ra[inv],median_dec[inv])
    ### import pylab as plt; plt.ion()
    ### plt.hist(sep*3600.,bins=100,log=True)
     
