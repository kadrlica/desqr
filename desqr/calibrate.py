#!/usr/bin/env python
"""
Derive photometric zeropoints.

Derived from Douglas Tucker's DELVE_tie_to_refcat2.py
"""
__author__ = "Alex Drlica-Wagner"

import os, sys
import logging
import subprocess
import glob
import datetime
from collections import OrderedDict as odict

import yaml
import pandas as pd
import numpy as np
import fitsio
from astropy.stats import sigma_clip

from utils import mkdir
from const import BANDS, TAGS
import downskim
import match
import utils

from ugali.utils.healpix import ang2disc

UID     = 'FILENAME'
FLUX    = 'FLUX_APER_08'
FLUXERR = 'FLUXERR_APER_08'
DTYPES = downskim.DTYPES + [(FLUX, '>f4'), (FLUXERR, '>f4')]
REFCAT_COLUMNS = ['OBJID','RA','DEC','G','DG','R','DR','I','DI','Z','DZ']

# Transformation coefficients
# mag_des = mag_ps + A[mag][0]*color_ps + A[mag][1]
A = odict([
        ('g', [+0.0994, -0.0076]), # [-0.2 < (g-r)_ps <= 1.2]
        ('r', [-0.1335, +0.0189]), # [-0.2 < (g-r)_ps <= 1.2]
        ('i', [-0.3407, +0.0026]), # [-0.2 < (i-z)_ps <= 0.3]
        ('z', [-0.2575, -0.0074]), # [-0.2 < (i-z)_ps <= 0.3]
        ('Y', [-0.6032, +0.0185]), # [-0.2 < (i-z)_ps <= 0.3]
        ])
# Updated based on fit to DES
A = odict([
        ('g', [+0.0994, -0.0076 - 0.0243]), # [-0.2 < (g-r)_ps <= 1.2]
        ('r', [-0.1335, +0.0189 + 0.0026]), # [-0.2 < (g-r)_ps <= 1.2]
        ('i', [-0.3407, +0.0026 - 0.0039]), # [-0.2 < (i-z)_ps <= 0.3]
        ('z', [-0.2575, -0.0074 - 0.0127]), # [-0.2 < (i-z)_ps <= 0.3]
        ('Y', [-0.6032, +0.0185]),          # [-0.2 < (i-z)_ps <= 0.3]
        ])


def derive_zeropoints(dataFrame, band):
    """ Modification of Douglas Tucker's `DELVE_tie_to_stds`.

    Parameters
    ----------
    dataFrame : pandas DataFrame of matched, merged objects containing
                observed measurements and reference catalog measurements
    band      : band to perform calibration on

    Returns
    -------
    dataFrame : output DataFrame containing zeropoint statistics
    """
    #NOTE: This function could be modified to use the dataFrame['BAND']
    #information to determine the reference magnitude for each
    #object. Not really necessary for single-exposure calibration, but
    #would remove the restriction on specifying `band` at input.

    logging.info("Starting to derive zeropoints...")
    logging.debug(datetime.datetime.now())

    aggFieldColName   = UID
    fluxObsColName    = FLUX
    fluxerrObsColName = FLUXERR

    # Transform ATLAS-REFCAT2 mags into DES mags for this filter band...
    magStdColName = '%s_des'%band
    magerrStdColName = '%serr_des'%band    
    if band is 'g':
        # g_des = g_ps + 0.0994*(g-r)_ps - 0.0076    [-0.2 < (g-r)_ps <= 1.2]
        dataFrame[magStdColName] = dataFrame['G']+\
            A[band][0]*(dataFrame['G']-dataFrame['R'])+A[band][1]
        dataFrame[magerrStdColName] = dataFrame['DG']  # temporary
        mask  = ( (dataFrame['G']-dataFrame['R']) > -0.2)
        mask &= ( (dataFrame['G']-dataFrame['R']) <= 1.2)
    elif band is 'r':
        # r_des = r_ps - 0.1335*(g-r)_ps + 0.0189    [-0.2 < (g-r)_ps <= 1.2]
        dataFrame[magStdColName] = dataFrame['R']+\
            A[band][0]*(dataFrame['G']-dataFrame['R'])+A[band][1]
        dataFrame[magerrStdColName] = dataFrame['DR']  # temporary
        mask  = ( (dataFrame['G']-dataFrame['R']) > -0.2)
        mask &= ( (dataFrame['G']-dataFrame['R']) <= 1.2)
    elif band is 'i':
        # i_des = i_ps - 0.3407*(i-z)_ps + 0.0026    [-0.2 < (i-z)_ps <= 0.3]
        dataFrame[magStdColName] = dataFrame['I']+\
            A[band][0]*(dataFrame['I']-dataFrame['Z'])+A[band][1]
        dataFrame[magerrStdColName] = dataFrame['DI']  # temporary
        mask  = ( (dataFrame['I']-dataFrame['Z']) > -0.2)
        mask &= ( (dataFrame['I']-dataFrame['Z']) <= 0.3)
    elif band is 'z':
        # z_des = z_ps - 0.2575*(i-z)_ps - 0.0074    [-0.2 < (i-z)_ps <= 0.3]
        dataFrame[magStdColName] = dataFrame['Z']+\
            A[band][0]*(dataFrame['I']-dataFrame['Z'])+A[band][1]
        dataFrame[magerrStdColName] = dataFrame['DZ']  # temporary
        mask  = ( (dataFrame['I']-dataFrame['Z']) > -0.2)
        mask &= ( (dataFrame['I']-dataFrame['Z']) <= 0.3)
    elif band is 'Y':
        # Y_des = z_ps - 0.6032*(i-z)_ps + 0.0185    [-0.2 < (i-z)_ps <= 0.3]
        dataFrame[magStdColName] = dataFrame['Z']+\
            A[band][0]*(dataFrame['I']-dataFrame['Z'])+A[band][1]
        dataFrame[magerrStdColName] = dataFrame['DZ']  # temporary
        mask  = ( (dataFrame['I']-dataFrame['Z']) > -0.2)
        mask &= ( (dataFrame['I']-dataFrame['Z']) <= 0.3)
    else:
        msg = "Unrecognized band: %s "%band
        raise ValueError(msg)

    dataFrame = dataFrame[mask].copy()

    # Add a 'MAG_OBS' column and a 'MAG_DIFF' column to the pandas DataFrame...
    dataFrame['MAG_OBS'] = -2.5*np.log10(dataFrame[fluxObsColName])
    dataFrame['MAG_DIFF'] = dataFrame[magStdColName]-dataFrame['MAG_OBS']

    ###############################################
    # Aggregate by aggFieldColName
    ###############################################

    # Make a copy of original dataFrame...
    df = dataFrame.copy()

    # Create an initial mask...
    mask1  = ( (df[magStdColName] >= 0.) & (df[magStdColName] <= 25.) )
    mask1 &= ( (df['FLUX_PSF'] > 10.) & (df['FLAGS'] < 2) & (np.abs(df['SPREAD_MODEL']) < 0.01) )
    if magerrStdColName != 'None':
        mask1 &= (df[magerrStdColName] < 0.1)

    magDiffGlobalMedian = df[mask1]['MAG_DIFF'].median()
    magDiffMin = magDiffGlobalMedian - 5.0
    magDiffMax = magDiffGlobalMedian + 5.0
    mask2 = ( (df['MAG_DIFF'] > magDiffMin) & (df['MAG_DIFF'] < magDiffMax) )
    mask = mask1 & mask2

    logging.info("Masking %s objects."%(~mask).sum())
    logging.info("%s objects remain."%mask.sum())

    if mask.sum() == 0:
        msg = "No unmasked objects!"
        raise RuntimeError(msg)

    # Iterate over the copy of dataFrame 3 times, removing outliers...
    #  We are using "Method 2/Group by item" from
    #  http://nbviewer.jupyter.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/07%20-%20Lesson.ipynb
    logging.info("Sigma-clipping...")
    nsig = 3.0
    for niter in range(3):
        logging.info("   iteration: %d...""" %niter)

        # make a copy of original df, and then delete the old one...
        newdf = df[mask].copy()
        del df

        # group by aggFieldColName...
        grpnewdf = newdf.groupby([aggFieldColName])

        # add/update new columns to newdf
        logging.debug("pre outlier: "+str(datetime.datetime.now()))
        newdf['Outlier']  = grpnewdf['MAG_DIFF'].transform( lambda x: abs(x-x.mean()) > nsig*x.std() )
        del grpnewdf
        logging.debug("post outlier: "+str(datetime.datetime.now()))

        nrows = newdf['MAG_DIFF'].size
        logging.info("  Number of rows remaining:  %d"%nrows)

        df = newdf
        mask = ( df['Outlier'] == False )

    # Perform pandas grouping/aggregating functions on sigma-clipped Data Frame...
    logging.info('Performing grouping/aggregating...')
    logging.debug("pre stats: "+str(datetime.datetime.now()))
    # These are unique columns (save space with float32)
    columns = [aggFieldColName,'EXPNUM','CCDNUM','BAND','EXPTIME','T_EFF']
    groupedDataFrame = df.groupby(columns)
    magZeroMedian = groupedDataFrame['MAG_DIFF'].median().astype(np.float32)
    magZeroMean   = groupedDataFrame['MAG_DIFF'].mean().astype(np.float32)
    magZeroStd    = groupedDataFrame['MAG_DIFF'].std().astype(np.float32)
    magZeroNum    = groupedDataFrame['MAG_DIFF'].count().astype(np.int32)
    magZeroErr    = (magZeroStd/np.sqrt(magZeroNum-1)).astype(np.float32)
    logging.debug("post stats: "+str(datetime.datetime.now()))

    # Rename these pandas series...
    magZeroMedian.name = 'MAG_ZERO_MEDIAN'
    magZeroMean.name   = 'MAG_ZERO_MEAN'
    magZeroStd.name    = 'MAG_ZERO_STD'
    magZeroNum.name    = 'MAG_ZERO_NUM'
    magZeroErr.name    = 'MAG_ZERO_MEAN_ERR'

    seriesList = []
    seriesList.extend([magZeroMedian, magZeroMean, magZeroStd, \
                       magZeroErr, magZeroNum])
    newDataFrame = pd.concat(seriesList, join='outer', axis=1 )

    logging.debug(datetime.datetime.now())
    return newDataFrame

def read_refcat(ra,dec,radius=1.5):
    """ Read the reference catalog in a localized region.

    Parameters
    ----------
    ra : ra (deg)
    dec : dec (deg)
    radius : regional radius (deg)

    Returns
    -------
    refcat : reference catalog columns
    """
    
    nside = 32
    pixels = ang2disc(nside, ra, dec, radius, inclusive=True)
    dirname = '/data/des40.b/data/atlas-refcat2/healpix'
    basename = 'atlas-refcat2_%05d.fits'
    filenames = [os.path.join(dirname,basename%p) for p in pixels]
    columns = REFCAT_COLUMNS
    refcat = utils.load_infiles(filenames,columns=columns)
    return refcat

def merge_refcat(catalog,refcat,angsep=1.0):
    """ Match and merge catalog and reference catalog.
    Both catalogs must contain the columns `RA` and `DEC`.

    Parameters
    ----------
    catalog : the observed catalog 
    refcat  : the reference catalog
    angsep  : the angular separation tolerence for matching (arcsec)

    Returns
    -------
    df : matched, merged pandas DataFrame
    """
    # Find index of closest match to each refcat object
    idx1,idx2,sep = match.match_query(refcat['RA'],refcat['DEC'],catalog['RA'],catalog['DEC'])

    # Matching tolerance
    angsep /= 3600. # (deg)
    sel = (sep < angsep)

    # Select sorted matched objects
    df1 = pd.DataFrame(catalog[idx2[sel]].byteswap().newbyteorder())
    df2 = pd.DataFrame(refcat[idx1[sel]].byteswap().newbyteorder())
    df2.rename(columns={'RA':'RA_REF','DEC':'DEC_REF'},inplace=True)

    # Return the merged data frame (note that columns will be duplicated)
    return pd.concat([df1,df2],axis=1)

def calibrate(outfile,select,exp,force):
    """ The calibration driver function.

    Parameters
    ----------
    outfile : output file name
    select  : selection function for skimming observed catalog
    exp     : exposure info
    force   : force overwrite of output file

    Returns
    -------
    zps     : zeropoints
    """
    if os.path.exists(outfile) and not force:
        logging.warn('Found %s; skipping...'%outfile)
        return

    expnum = exp['expnum']
    filenames = downskim.get_filenames(expnum)

    logging.info("Reading object catalog...")
    catalog = downskim.create_downskim(filenames,select,exp=exp,dtype=DTYPES)
    logging.info("  %s objects"%len(catalog))

    logging.info("Reading reference catalog...")
    refang = 1.2 # config param: reference catalog angular selection (deg)
    refcat = read_refcat(exp['telra'],exp['teldec'],refang)
    logging.info("  %s objects"%len(refcat))

    logging.info("Matching with reference catalog...")
    angsep = 1.0 # config param: angular separation for matching (arcsec)
    dataFrame = merge_refcat(catalog,refcat,angsep)
    logging.info("  %s objects"%len(dataFrame))

    # For debugging...
    if False:
        fname = os.path.join(mkdir('match'),os.path.basename(outfile))
        fname = fname.replace('.fits','.csv')
        dataFrame.to_csv(fname, index=False)

    logging.info("Deriving zeropoints...")
    output = derive_zeropoints(dataFrame,exp['band'])

    logging.info("Writing %s..."%outfile)
    ext = os.path.splitext(outfile)[-1]

    if ext in ('.fits','.fits.gz'):
        fitsio.write(outfile, output.to_records(), clobber=True)
    else:
        output.to_csv(outfile, float_format='%.4f')

    return output

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('explist', help='exposure list')
    parser.add_argument('outfile', help='output file')
    parser.add_argument('-f','--force',action='store_true',
                        help='force rerun')
    parser.add_argument('-e','--expnum',nargs='+',default=[],type=int,
                        action='append',
                        help='exposure number')
    parser.add_argument('-v','--verbose',action='store_true',
                        help='verbose output')
    parser.add_argument('-t','--tag',default=None,choices=TAGS,
                        help='data set tag')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.getLogger().setLevel(level)

    explist = pd.read_csv(args.explist).to_records(index=False)
    #explist.dtype.names = [str(n).upper() for n in explist.dtype.names]
    explist.dtype.names = [str(n).lower() for n in explist.dtype.names]

    #exp = explist[explist['EXPNUM'] == args.expnum]
    #if len(exp) == 0:
    #    raise ValueError("EXPNUM not found: %s"%args.expnum)
    #elif len(exp) > 1:
    #    msg = "Multiple values for EXPNUM found: %s"%args.expnum
    #    for e in exp: msg += ("\n"+e)
    #    raise ValueError(msg)
    #exp = exp[0]

    explist = explist[np.in1d(explist['expnum'], args.expnum)]
    if len(explist) == 0:
        raise ValueError("EXPNUM not found: %s"%args.expnum)

    status = 0
    for exp in explist:

        TAG = exp['tag'] if args.tag is None else args.tag

        if TAG.startswith('BLISS'):
            select = downskim.bliss_select
        elif TAG.startswith('MAGLITES'):
            select = downskim.maglites_select
        elif TAG.startswith("DELVE"):
            select = downskim.delve_select
        else:
            raise Exception('Tag not found: %s'%TAG)

        edict = dict(zip(exp.dtype.names,exp))
        outfile = args.outfile.format(**edict)

        try: 
            zps = calibrate(outfile,select,exp,args.force)
        except Exception as e:
            logging.error(str(e))
            status = 1

    sys.exit(status)
