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
from scipy import interpolate

from utils import mkdir
from const import BANDS, TAGS
import downskim
import match
import utils

from ugali.utils.healpix import ang2disc

SURVEYS = ['refcat2','desdr2']
UID     = 'FILENAME'
FLUX    = 'FLUX_PSF'
FLUXERR = 'FLUXERR_PSF'
DTYPES = downskim.DTYPES + [(FLUX, '>f4'), (FLUXERR, '>f4')]
REFCAT_COLUMNS = ['OBJID','RA','DEC','G','DG','R','DR','I','DI','Z','DZ']
DESDR2_COLUMNS = ['COADD_OBJECT_ID','RA','DEC',
                  'WAVG_MAG_PSF_G','WAVG_MAGERR_PSF_G',
                  'WAVG_MAG_PSF_R','WAVG_MAGERR_PSF_R',
                  'WAVG_MAG_PSF_I','WAVG_MAGERR_PSF_I',
                  'WAVG_MAG_PSF_Z','WAVG_MAGERR_PSF_Z']

INTERP_DIR = '/home/s1/kadrlica/projects/delve/calib/v4/interp'

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
REFCAT2 = odict([
        ('g', [+0.0994, -0.0076 - 0.0243]), # [-0.2 < (g-r)_ps <= 1.2]
        ('r', [-0.1335, +0.0189 + 0.0026]), # [-0.2 < (g-r)_ps <= 1.2]
        ('i', [-0.3407, +0.0026 - 0.0039]), # [-0.2 < (i-z)_ps <= 0.3]
        ('z', [-0.2575, -0.0074 - 0.0127]), # [-0.2 < (i-z)_ps <= 0.3]
        ('Y', [-0.6032, +0.0185]),          # [-0.2 < (i-z)_ps <= 0.3]
        ])

# Transformation to DES DR2 (NOOP)
DESDR2 = odict([
        ('g', [0.0, 0.0]), 
        ('r', [0.0, 0.0]), 
        ('i', [0.0, 0.0]), 
        ('z', [0.0, 0.0]), 
        ('Y', [0.0, 0.0]), 
        ])



# New from Chin Yi (Calibrated with des_ncsa delve_zps_v4.2.fits)
REFCAT2_INTERP = odict([
        ('g', ['transInterp.ref2_to_des.g_gr_ref2_North_v3.csv',  
               'transInterp.ref2_to_des.g_gr_ref2_South_v3.csv',
               +3.1-0.1, +2.4+1.2]),
        ('r', ['transInterp.ref2_to_des.r_gr_ref2_North_v3.csv',
               'transInterp.ref2_to_des.r_gr_ref2_South_v3.csv',
               -0.2, +1.8-0.3]),
        ('i', ['transInterp.ref2_to_des.i_iz_ref2_North_v3.csv',
               'transInterp.ref2_to_des.i_iz_ref2_South_v3.csv',
               -0.0-0.1, -0.5-0.4]),
        ('z', ['transInterp.ref2_to_des.z_iz_ref2_North_v3.csv',
               'transInterp.ref2_to_des.z_iz_ref2_South_v3.csv',
               -1.0+0.1, -3.8+1.2])])



def interp_to_des(dataFrame, band, interps, colnames):
    """ Transform to DES bandpass using interpolations.
    
    Paramters
    ---------
    dataFrame : input data frame containing magnitudes
    band      : output band

    Returns
    -------
    dataFrame : input DataFrame with new column added.
    mask      : mask for useable rows
    """
    logging.info("Performing interpolation transformation from reference catalog to DECam...")

    # Output column names
    magStdColName,magerrStdColName = colnames

    # Declination split
    split = -30 # Declination split (deg)

    # Offset to go from Refcat2 to DES
    offset = np.zeros(len(dataFrame), dtype=float)

    # Interpolation color
    if band in ['g','r']: 
        logging.info("Interpolating in (g-r) color...")
        color = dataFrame['G']-dataFrame['R']
    else:                 
        logging.info("Interpolating in (i-z) color...")
        color = dataFrame['I']-dataFrame['Z']

    # For interpolation
    kwargs = dict(bounds_error=False, fill_value=0., kind='linear')

    # Fill offset for objects in the north
    north_sel = dataFrame['DEC'] > split
    if north_sel.sum():
        filename=os.path.join('interp',interps[band][0])
        logging.info("North: %s"%filename)
        df_interp_north = pd.read_csv(filename)  
        interp_north = interpolate.interp1d(df_interp_north.bin_label.values.astype(float),
                                            df_interp_north.bin_median.values, **kwargs)
        offset[north_sel] = interp_north(color)[north_sel] - (interps[band][2]/1000)

    # Fill offset for objects in the south
    south_sel = dataFrame['DEC'] <= split
    if south_sel.sum():
        filename=os.path.join('interp',interps[band][1])
        logging.info("South: %s"%filename)
        df_interp_south = pd.read_csv(filename)  
        interp_south = interpolate.interp1d(df_interp_south.bin_label.values.astype(float),
                                            df_interp_south.bin_median.values, **kwargs)
        offset[south_sel] = interp_south(color)[south_sel] - (interps[band][3]/1000)

    # Calculate the magnitude
    if band is 'g':
        # Create linear interpolation of the median dmag vs. color bin calculated above...
        dataFrame[magStdColName] = dataFrame['G'] + offset
        dataFrame[magerrStdColName] = dataFrame['DG']  # temporary
        mask  = ( color > -0.2)
        mask &= ( color <= 1.0)
        mask &= ( dataFrame['DG'] <= 0.05)
        mask &= ( dataFrame['DR'] <= 0.05)
        mask &= ( dataFrame['DI'] <= 0.05)
        mask &= ( dataFrame['DZ'] <= 0.05)
    elif band is 'r':
        dataFrame[magStdColName] = dataFrame['R'] + offset
        dataFrame[magerrStdColName] = dataFrame['DR']  # temporary
        mask  = ( color > -0.2)
        mask &= ( color <= 1.0)
        mask &= ( dataFrame['DG'] <= 0.05)
        mask &= ( dataFrame['DR'] <= 0.05)
        mask &= ( dataFrame['DI'] <= 0.05)
        mask &= ( dataFrame['DZ'] <= 0.05)
    elif band is 'i':
        dataFrame[magStdColName] = dataFrame['I'] + offset
        dataFrame[magerrStdColName] = dataFrame['DI']  # temporary
        mask  = ( color > -0.2)
        mask &= ( color <= 0.3)
        mask &= ( dataFrame['DG'] <= 0.05)
        mask &= ( dataFrame['DR'] <= 0.05)
        mask &= ( dataFrame['DI'] <= 0.05)
        mask &= ( dataFrame['DZ'] <= 0.05)
    elif band is 'z':
        dataFrame[magStdColName] = dataFrame['Z'] + offset
        dataFrame[magerrStdColName] = dataFrame['DZ']  # temporary
        mask  = ( color > -0.2)
        mask &= ( color <= 0.3)
        mask &= ( dataFrame['DG'] <= 0.05)
        mask &= ( dataFrame['DR'] <= 0.05)
        mask &= ( dataFrame['DI'] <= 0.05)
        mask &= ( dataFrame['DZ'] <= 0.05)
    else:
        msg = "Unrecognized band: %s "%band
        raise ValueError(msg)

    return dataFrame, mask


def compare_methods(dataFrame,band):
    logging.info("Comparing transform methods...")
    logging.debug(datetime.datetime.now())

    magStdColName = '%s_des'%band
    magerrStdColName = '%serr_des'%band    

    dataFrame1, mask1 = linear_to_des(dataFrame.copy(), band, REFCAT2,
                                      colnames=[magStdColName,magerrStdColName])

    dataFrame2, mask2 = interp_to_des(dataFrame.copy(), band, REFCAT2_INTERP, 
                                      colnames=[magStdColName,magerrStdColName])

    return [dataFrame1, mask1], [dataFrame2, mask2]

def derive_zeropoints(dataFrame, band, survey='refcat2', transform='linear'):
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

    if   survey == 'refcat2': A = REFCAT2
    elif survey == 'desdr2':  A = DESDR2
    else: 
        raise Exception('Unrecognized survey: %s'%survey)        

    logging.info("Starting to derive zeropoints...")
    logging.debug(datetime.datetime.now())

    aggFieldColName   = UID
    fluxObsColName    = FLUX
    fluxerrObsColName = FLUXERR

    magStdColName = '%s_des'%band
    magerrStdColName = '%serr_des'%band    
    cols = [magStdColName,magerrStdColName]

    if transform=='linear':
        dataFrame, mask = linear_to_des(dataFrame.copy(), band, A, 
                                        colnames=cols)
    else:
        dataFrame, mask = interp_to_des(dataFrame.copy(), band, REFCAT2_INTERP, 
                                        colnames=cols)

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

    ## Selection of catalog magerr
    #if 'MAGERR_PSF' in df.columns:
    #    mask1 &= (df[FLUX]/df[FLUXERR]) > 20

    # Selection of reference catalog magerr
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
        grpnewdf = newdf.groupby(['EXPNUM', 'CCDNUM'])

        # add/update new columns to newdf
        newdf['Outlier']  = grpnewdf['MAG_DIFF'].transform( lambda x: abs(x-x.mean()) > nsig*x.std() )
        del grpnewdf

        nrows = newdf['MAG_DIFF'].size
        logging.info("  Number of rows remaining:  %d"%nrows)

        df = newdf
        mask = ( df['Outlier'] == False )

    # Perform pandas grouping/aggregating functions on sigma-clipped Data Frame...
    logging.info('Performing grouping/aggregating...')
    # These are unique columns (save space with float32)
    columns = ['EXPNUM','CCDNUM','BAND','EXPTIME']
    groupedDataFrame = df.groupby(columns)
    magZeroMedian = groupedDataFrame['MAG_DIFF'].median().astype(np.float32)
    magZeroMean   = groupedDataFrame['MAG_DIFF'].mean().astype(np.float32)
    magZeroStd    = groupedDataFrame['MAG_DIFF'].std().astype(np.float32)
    magZeroNum    = groupedDataFrame['MAG_DIFF'].count().astype(np.int32)
    magZeroErr    = (magZeroStd/np.sqrt(magZeroNum-1)).astype(np.float32)

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

def read_desdr2(ra,dec,radius=1.5):
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
    dirname = '/data/des40.b/data/y6a1/gold/1.1/healpix'
    basename = 'y6_gold_1_0_%05d.fits'
    filenames = [os.path.join(dirname,basename%p) for p in pixels]
    #filenames = [f for f in filenames if os.path.exists(f)]
    columns = DESDR2_COLUMNS
    refcat = utils.load_infiles(filenames,columns=columns)
    refcat = refcat[refcat['WAVG_MAG_PSF_R'] < 21]

    mapping = dict(zip(DESDR2_COLUMNS,REFCAT_COLUMNS))
    new_names = [mapping[n] for n in refcat.dtype.names]
    refcat.dtype.names = new_names

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

def calibrate(outfile,select,exp,survey,force,transform='linear'):
    """ The calibration driver function.

    Parameters
    ----------
    outfile : output file name
    select  : selection function for skimming observed catalog
    exp     : exposure info
    survey  : survey
    force   : force overwrite of output file
    transform : type of transformation [linear, interp]

    Returns
    -------
    zps     : zeropoints
    """
    if utils.is_found(outfile,force):
        return

    expnum = exp['expnum']
    path = None if 'path' not in exp.dtype.names else exp['path']
    filenames = expnum #filenames = downskim.get_filenames(expnum,path)

    logging.info("Reading object catalog...")
    logging.debug('  '+str(datetime.datetime.now()))
 #   catalog =  fitsio.read('/home/s1/chinyi/90adata/exposures/{}/D{:08d}_{}_zps.fits'.format(exp['band'],exp['expnum'],exp['band'],)) 
    catalog =  fitsio.read('/home/s1/chinyi/80bdata/for_others/Will_Cerny/AquaIII/raw_cat/{}/D{:08d}_{}_zps.fits'.format(exp['band'],exp['expnum'],exp['band'],))   
    #catalog = downskim.create_downskim(filenames,select,exp=exp,dtype=DTYPES)
    logging.debug('  '+str(datetime.datetime.now()))
    logging.info("  %s objects"%len(catalog))


    logging.info("Reading reference catalog...")
    refang = 1.2 # config param: reference catalog angular selection (deg)

    logging.debug('  '+str(datetime.datetime.now()))
    if survey == 'refcat2': 
        refcat = read_refcat(exp['telra'],exp['teldec'],refang)
    elif survey == 'desdr2':
        refcat = read_desdr2(exp['telra'],exp['teldec'],refang)
    else:
        raise Exception('Unrecognized survey: %s'%survey)
    logging.debug('  '+str(datetime.datetime.now()))
    logging.info("  %s objects"%len(refcat))

    logging.info("Matching with reference catalog...")
    #angsep = 1.0 # config param: angular separation for matching (arcsec)
    angsep = 0.5 # config param: angular separation for matching (arcsec)
    dataFrame = merge_refcat(catalog,refcat,angsep)
    logging.info("  %s objects"%len(dataFrame))

    # For debugging...
    if False:
        fname = os.path.join(mkdir('match'),os.path.basename(outfile))
        fname = fname.replace('.fits','.csv')
        dataFrame.to_csv(fname, index=False)

    logging.info("Deriving zeropoints...")
    output = derive_zeropoints(dataFrame,exp['band'],survey=survey, 
                               transform=transform)

    logging.info("Writing %s..."%outfile)
      
    header = [dict(name='SURVEY',value=survey,comment='reference catalog')]
    fitsio.write(outfile, output.to_records(), header=header, clobber=True)
   # ext = os.path.splitext(outfile)[-1]

    #if ext in ('.fits','.fits.gz'):
    #    dirname = os.path.dirname(filenames[0])
    #    header = [dict(name='PATH',value=dirname,comment='original path'),
    #              dict(name='SURVEY',value=survey,comment='reference catalog')]
    #    fitsio.write(outfile, output.to_records(), header=header, clobber=True)
    #else:
    #    output.to_csv(outfile, float_format='%.4f')

    return output

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('explist', help='exposure list')
    parser.add_argument('outfile', help='output file')
    parser.add_argument('-f','--force',action='store_true',
                        help='force rerun')
    parser.add_argument('-e','--expnum',nargs='+',default=[],type=int,
                        action='append',help='exposure number')
    parser.add_argument('-v','--verbose',action='store_true',
                        help='verbose output')
    parser.add_argument('-s','--survey',default='refcat2',choices=SURVEYS,
                        help='data set tag')
    parser.add_argument('-t','--tag',default=None,choices=TAGS,
                        help='data set tag')
    parser.add_argument('--transform',default='linear',
                        choices=['linear','interp'],
                        help='type of transformation')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.getLogger().setLevel(level)
    print("Pocca",args.explist)
    explist = pd.read_csv(args.explist).to_records(index=False)
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
    for i,exp in enumerate(explist):
        logging.info("(%s/%s): %s "%(i+1,len(explist),exp['expnum']))

        # Hack to allow for DELVE selects
        select = downskim.delve_select
        
        #TAG = exp['tag'] if args.tag is None else args.tag

        #if TAG.startswith('BLISS'):
        #    select = downskim.bliss_select
        #elif TAG.startswith('MAGLITES'):
        #    select = downskim.maglites_select
        #elif TAG.startswith("DELVE"):
        #    select = downskim.delve_select
        #else:
        #    raise Exception('Tag not found: %s'%TAG)

        edict = dict(zip(exp.dtype.names,exp))
        outfile = args.outfile.format(**edict)

        try: 
            zps = calibrate(outfile,select,exp,args.survey,args.force,transform=args.transform)
        except Exception as e:
            logging.error(str(e))
            status = 1

    sys.exit(status)
