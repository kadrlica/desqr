#!/usr/bin/env python
"""
Characterize photometric scatter.
"""
__author__ = "Alex Drlica-Wagner"
import os
from os.path import join
import glob
import yaml
from collections import OrderedDict as odict

import matplotlib
if os.getenv('TERM')=='screen' or not os.getenv('DISPLAY'):
    matplotlib.use('Agg')
import pylab as plt

import numpy as np
import scipy.ndimage as nd
import pandas as pd
import healpy as hp
import fitsio

from ugali.utils.healpix import ang2pix, pix2ang
import ugali.utils.projector
from ugali.utils.projector import angsep
from ugali.utils.shell import mkdir
from ugali.utils.logger import logger

import utils
from const import OBJECT_ID, UNIQUE_ID, BANDS, BADMAG
from utils import bfield, bfields, load_infiles
from utils import blank, rec_append_fields
from match import match_query

import plotting
from catalog import good_objects

TYPES = ['gaia','gaia_dr2','gaia_edr3','des','des_dr2','hpx_rms','wavg_rms']
COLUMNS = [OBJECT_ID,'RA','DEC']
GAIA_DIR = {
    'dr2':'/data/des40.b/data/gaia/dr2/healpix/',
    'edr3':'/data/des40.b/data/gaia/edr3/healpix/',
}

DES_BASE = {
    'dr2':'/data/des40.b/data/des/dr2/healpix/dr2_main_%05d.fits',
}

def plot_photometry(filename,outfile=None,survey='delve'):
    logger.info("Reading %s..."%filename)
    hpxmap = hp.read_map(filename,verbose=False)
    hpxmap *= 1000 # mmag

    cbar_kwargs = dict()
    hpxmap_kwargs = dict(xsize=2000)
    hist_kwargs = dict()

    if '_rms_' in filename:    
        label = 'MAG RMS (mmag)'
    elif '_gaia_' in filename: 
        label = r'$\Delta$(Gaia G) (mmag)'
        hpxmap_kwargs['vmin'] = vmin = -30
        hpxmap_kwargs['vmax'] = vmax = 30
        hist_kwargs['bins'] = np.linspace(vmin,vmax)
    elif '_des_' in filename:  
        label = r'$\Delta$(DES) (mmag)'
        hpxmap_kwargs['vmin'] = vmin = -0.05
        hpxmap_kwargs['vmax'] = vmax = 0.05
        hist_kwargs['bins'] = np.linspace(vmin,vmax)
        survey='des'
    else: 
        label = None

    cbar_kwargs['label'] = label

    fig,axes,smap = plotting.plot_hpxmap_hist(hpxmap,survey,cbar_kwargs,hpxmap_kwargs,
                                              hist_kwargs)
    #axes[0].annotate('%s band'%band, (0.05,0.93), xycoords='axes fraction')
    axes[1].set_xlabel(label)

    utils.print_statistics(hpxmap)

    if outfile is None: 
        outfile=os.path.basename(filename).split('.')[0]+'.png'
    logger.info("Writing %s..."%outfile)
    plt.savefig(outfile,bbox_inches='tight')

def get_des_catalog(hpx,columns=['RA','DEC','WAVG_MAG_PSF_G'],version='dr2'):
    """ Grab DES catalogs """
    basename = DES_BASE[version]

    pixels = [hpx]
    filenames = [basename%p for p in pixels]
    filenames = [f for f in filenames if os.path.exists(f)]
    if len(filenames):
        cat = load_infiles(filenames,columns=columns)
    else:
        cat = None
    return cat

def get_gaia_catalog(hpx,columns=['RA','DEC','PHOT_G_MEAN_FLUX'],version='dr2'):
    dirname = GAIA_DIR[version]
    basename = 'GaiaSource_%05d.fits'

    pixels = np.atleast_1d(hpx)
    filenames = [os.path.join(dirname,basename%p) for p in pixels]
    filenames = [f for f in filenames if os.path.exists(f)]
    cat = load_infiles(filenames,columns=columns)

    # Need to put the Gaia flux onto the AB system
    if version == 'edr3':
        #https://www.cosmos.esa.int/web/gaia/edr3-passbands
        gzp_g = 25.8010446445
    elif version == 'dr2':
        # For magnitudes published in Gaia DR2
        #https://www.cosmos.esa.int/web/gaia/iow_20180316
        gzp_g = 25.7933969562 

    mag_G = gzp_g - 2.5*np.log10(cat['PHOT_G_MEAN_FLUX'])

    cat = rec_append_fields(cat,'PHOT_G_MEAN_MAG_AB',mag_G)
    return cat

def gaia_transform(g, r, i, z, version='dr2'):
    """From Eli Rykoff for FGCM comparisons via Slack:
    https://darkenergysurvey.slack.com/archives/C016J5SDV9T/p1608182260343200

    This complex model behaves about as well as a random forest
    classifier for Gaia DR2.  The magnitude dependence of the
    transformation is huge because of background errors in Gaia DR2.

    The EDR3 transformations have a much smaller r-offset curvature
    due to background issues and should be able to go deeper.

    These transformations are valid for 0 < g - i < 1.5

    /nfs/slac/kipac/fs1/g/des/erykoff/des/y6a1/fgcm/run_v2.1.2/gedr3/gaia_superfitmodel3.py

    Parameters
    ----------
    g, r, i, z : DECam magnitudes in g,r,i,z  (AB system)
    version    : Version of the Gaia catalog to compare against

    Returns
    -------
    Gmag  : Predicted Gaia G-band magnitude based on DECam griz
    """

    magConst = 2.5 / np.log(10.0)
    lambdaStd = np.array([4790.28076172, 6403.26367188, 7802.49755859, 9158.77441406])
    i0Std = np.array([0.16008162, 0.18297842, 0.17169334, 0.1337308])
    i1Std = np.array([2.09808350e-05, -1.22070312e-04, -1.08942389e-04, -8.01086426e-05])
    i10Std = i1Std / i0Std
    fudgeFactors = np.array([0.25, 1.0, 1.0, 0.25])
    fudgeShift   = 0.0 # Additive shift in peak (mag)
    nBands = lambdaStd.size

    #g -= -46 * 1e-3
    #r -= 0   * 1e-3
    #i -= 38  * 1e-3
    #z -= 31  * 1e-3

    trainCat = np.rec.fromarrays([g,r,i,z],names=['g','r','i','z'])
    fluxg = 10**(g/-2.5)
    fluxr = 10**(r/-2.5)
    fluxi = 10**(i/-2.5)
    fluxz = 10**(z/-2.5)

    S = np.zeros((trainCat.size, nBands - 1), dtype='f8')

    S[:, 0] = (-1. / magConst) * (trainCat['r'] - trainCat['g']) / (lambdaStd[1] - lambdaStd[0])
    S[:, 1] = (-1. / magConst) * (trainCat['i'] - trainCat['r']) / (lambdaStd[2] - lambdaStd[1])
    S[:, 2] = (-1. / magConst) * (trainCat['z'] - trainCat['i']) / (lambdaStd[3] - lambdaStd[2])

    fnuPrime = np.zeros((trainCat.size, nBands))
    fnuPrime[:, 0] = S[:, 0] + fudgeFactors[0] * (S[:, 1] + S[:, 0])
    fnuPrime[:, 1] = fudgeFactors[1] * (S[:, 0] + S[:, 1]) / 2.0
    fnuPrime[:, 2] = fudgeFactors[2] * (S[:, 1] + S[:, 2]) / 2.0
    fnuPrime[:, 3] = S[:, 2] + fudgeFactors[3] * ((lambdaStd[3] - lambdaStd[2]) / (lambdaStd[3] - lambdaStd[1])) * (S[:, 2] - S[:, 1])

    if version=='dr2':
        # DR2 fit parameters
        pars = [ 1.43223290e+00,  1.50877061e+00, 8.43173013e-01, -5.99023967e-04, 
                 4.06188382e-01,  3.11181978e-01, 2.51002598e-01,  1.00000000e-05,  
                 4.94284725e-03,  1.80499806e-03]

    elif version=='edr3':
        # EDR3: Notice that the last two parameters (describing the
        # r-offset curvature due to background issues) are much smaller.
        pars = [ 2.61815727e+00,  2.69372875e+00,  1.45644592e+00, -5.99023051e-04,
                 3.97535324e-01,  3.15794343e-01,  2.55484718e-01,  1.00000000e-05,
                 8.30152817e-04, -3.57980758e-04]

        fudgeShift = 2.4e-3 # Additive shift in peak (mag)
    else:
        raise Exception("Unrecognized Gaia version: %s"%version)

    i10g,i10r,i10i,i10z = pars[0:4]
    kg,kr,ki,kz = pars[4:8]
    r1,r2 = pars[8:10]
    desFlux = (kg * fluxg * (1.0 + fnuPrime[:, 0] * i10g) +
               kr * fluxr * (1.0 + fnuPrime[:, 1] * i10r) +
               ki * fluxi * (1.0 + fnuPrime[:, 2] * i10i) +
               kz * fluxz * (1.0 + fnuPrime[:, 3] * i10z))
    mGDES = -2.5 * np.log10(desFlux)
    rMag = -2.5 * np.log10(fluxr)
    mGDES += r1 * (rMag - 17.0) + r2 * (rMag - 17.0)**2.
    mGDES += fudgeShift # Peak shift
    return mGDES

def gaia_match(filename,version='edr3',verbose=True):
    """ Match input catalog to Gaia catalog.

    Parameters
    ----------
    filename : catalog filename
    version  : Gaia catalog version ['dr2', 'edr3']

    Returns
    -------
    ra,dec,G_pred,G_gaia : object coordinates, predicted Gaia G, Gaia G
    """
    if not os.path.exists(filename): 
        msg = "Couldn't find %s"%filename
        raise Exception(msg)

    columns = ['RA','DEC']
    #spread,nepochs = ['SPREAD_MODEL_R','NEPOCHS_R']
    #mag_g,mag_r,mag_i,mag_z = bfields(['MAG_PSF'],['g','r','i','z'])
    spread,nepochs = ['WAVG_SPREAD_MODEL_R','NEPOCHS_R']
    mag_g,mag_r,mag_i,mag_z = bfields(['WAVG_MAG_PSF'],['g','r','i','z'])
    columns += [spread, nepochs, mag_g, mag_r, mag_i, mag_z]
                       
    # Hack to get pixel location
    hpx = int(filename.split('_')[-1].split('.')[0])
    if verbose:
        ra,dec = hp.pix2ang(32, hpx, lonlat=True)
        radius = np.degrees(hp.max_pixrad(32))
        msg = '%s (RA,DEC,RAD) = %.2f,%.2f,%.2f'%(os.path.basename(filename),ra,dec,radius)
        print(msg)

    #print("Getting coadd catalog...")
    cat = load_infiles([filename],columns)

    # Select stars with 16 < r < 20 and 0.5 < (g-i) < 1.5
    sel = (np.fabs(cat[spread])<0.002) & \
        (cat[mag_g]<90) & (cat[mag_r]<90) & (cat[mag_i]<90) & (cat[mag_z]<90) & \
        (cat[mag_r]>16) & (cat[mag_r]<20) & \
        ((cat[mag_g] - cat[mag_i]) > 0.5) & \
        ((cat[mag_g] - cat[mag_i]) < 1.5) & \
        (cat[nepochs] > 1)
    cat = cat[sel]

    if len(cat) == 0:
        msg = "WARNING: No objects passing selection in: %s"%filename
        print(msg)
        return np.array([]),np.array([]),np.array([]),np.array([])

    #msg = "Getting external catalog: %s"%catalog
    ext = get_gaia_catalog(hpx,version=version)
    ext = ext[ext['PHOT_G_MEAN_MAG_AB'] < 20]
    if len(ext) == 0:
        msg = "WARNING: No Gaia objects passing selection."
        print(msg)
        return np.array([],dtype=int), np.array([])

    m = match_query(cat['RA'],cat['DEC'],ext['RA'],ext['DEC'])

    # Use a fairly narrow matching radius
    cut = 0.5 # arcsec
    sel = m[-1]*3600. < cut

    cat_match = cat[m[0][sel]]
    ext_match = ext[m[1][sel]]

    cat_G = gaia_transform(cat_match[mag_g],cat_match[mag_r],cat_match[mag_i],cat_match[mag_z],
                           version=version)
    ext_G = ext_match['PHOT_G_MEAN_MAG_AB']

    return cat_match['RA'],cat_match['DEC'],cat_G,ext_G

def gaia_match_fgcm(filename,version='edr3',verbose=True):
    """ Match input catalog to Gaia catalog. Hack for FGCM

    Parameters
    ----------
    filename : catalog filename
    version  : Gaia catalog version ['dr2', 'edr3']

    Returns
    -------
    ra,dec,G_pred,G_gaia : object coordinates, predicted Gaia G, Gaia G
    """
    if not os.path.exists(filename): 
        msg = "Couldn't find %s"%filename
        raise Exception(msg)

    columns = ['RA','DEC']
    spread,nepochs = ['SPREAD_MODEL_R','NEPOCHS_R']
    mag_g,mag_r,mag_i,mag_z = bfields(['MAG_PSF'],['g','r','i','z'])
    #spread,nepochs = ['WAVG_SPREAD_MODEL_R','NEPOCHS_R']
    #mag_g,mag_r,mag_i,mag_z = bfields(['WAVG_MAG_PSF'],['g','r','i','z'])
    columns += [spread, nepochs, mag_g, mag_r, mag_i, mag_z]

    # Include ZPs
    columns += bfields('MAG_ZERO',['g','r','i','z'])
    columns += bfields('FGCM_ZPT',['g','r','i','z'])
                       
    # Hack to get pixel location
    hpx = int(filename.split('_')[-1].split('.')[0])
    if verbose:
        ra,dec = hp.pix2ang(32, hpx, lonlat=True)
        radius = np.degrees(hp.max_pixrad(32))
        msg = '%s (RA,DEC,RAD) = %.2f,%.2f,%.2f'%(os.path.basename(filename),ra,dec,radius)
        print(msg)

    #print("Getting coadd catalog...")
    cat = load_infiles([filename],columns)

    # Switch ZPs
    # mag = -2.5 * np.log10(flux) + zeropoint
    cat[mag_g] += cat['FGCM_ZPT_G'] - cat['MAG_ZERO_G']
    cat[mag_r] += cat['FGCM_ZPT_R'] - cat['MAG_ZERO_R']
    cat[mag_i] += cat['FGCM_ZPT_I'] - cat['MAG_ZERO_I']
    cat[mag_z] += cat['FGCM_ZPT_Z'] - cat['MAG_ZERO_Z']

    # Select stars with 16 < r < 20 and 0.5 < (g-i) < 1.5
    sel = (np.fabs(cat[spread])<0.002) & \
        (cat[mag_g]<90) & (cat[mag_r]<90) & \
        (cat[mag_i]<90) & (cat[mag_z]<90) & \
        (cat[mag_r]>16) & (cat[mag_r]<20) & \
        ((cat[mag_g] - cat[mag_i]) > 0.5) & \
        ((cat[mag_g] - cat[mag_i]) < 1.5) & \
        (cat[nepochs] > 1) 

    sel &= np.abs(cat['FGCM_ZPT_G'] - cat['MAG_ZERO_G']) < 0.1
    sel &= np.abs(cat['FGCM_ZPT_R'] - cat['MAG_ZERO_R']) < 0.1
    sel &= np.abs(cat['FGCM_ZPT_I'] - cat['MAG_ZERO_I']) < 0.1
    sel &= np.abs(cat['FGCM_ZPT_Z'] - cat['MAG_ZERO_Z']) < 0.1

    cat = cat[sel]

    if len(cat) == 0:
        msg = "WARNING: No objects passing selection in: %s"%filename
        print(msg)
        return np.array([]),np.array([]),np.array([]),np.array([])

    #msg = "Getting external catalog: %s"%catalog
    ext = get_gaia_catalog(hpx,version=version)
    ext = ext[ext['PHOT_G_MEAN_MAG_AB'] < 20]
    if len(ext) == 0:
        msg = "WARNING: No Gaia objects passing selection."
        print(msg)
        return np.array([],dtype=int), np.array([])

    m = match_query(cat['RA'],cat['DEC'],ext['RA'],ext['DEC'])

    # Use a fairly narrow matching radius
    cut = 0.5 # arcsec
    sel = m[-1]*3600. < cut

    cat_match = cat[m[0][sel]]
    ext_match = ext[m[1][sel]]

    cat_G = gaia_transform(cat_match[mag_g],cat_match[mag_r],cat_match[mag_i],cat_match[mag_z],
                           version=version)
    ext_G = ext_match['PHOT_G_MEAN_MAG_AB']

    return cat_match['RA'],cat_match['DEC'],cat_G,ext_G

def gaia_photometry(filename,nside=64,band=None,plot=False,version='edr3'):
    """ Calculate the median spread between DECam catalog and Gaia magnitudes.

    Parameters
    ----------
    filename : input filename (catalog file)
    nside    : output nside
    band     : band of interest [not used!]
    plot     : pause and plot
    version  : Gaia catalog version ['dr2', 'edr3']

    Returns
    -------
    upix, stat : unique pixel and value of statistic (median) in that pixel
    """
    if not os.path.exists(filename): 
        msg = "Couldn't find %s"%filename
        raise Exception(msg)

    ra,dec,cat_G,gaia_G = gaia_match(filename,version)
    #ra,dec,cat_G,gaia_G = gaia_match_fgcm(filename,version)

    if len(ra) == 0:
        msg = "WARNING: No objects passing selection."
        print(msg)
        return np.array([],dtype=int), np.array([])

    diff  = cat_G - gaia_G

    pix = hp.ang2pix(nside,ra,dec,lonlat=True)
    upix = np.unique(pix)
    stat = nd.median(diff,labels=pix,index=upix)

    if False:
        plt.figure()
        plt.hist(cat_G - ext_G)
        import pdb; pdb.set_trace()
        
    return upix,stat

def des_photometry(filename,nside=64,band=None,plot=False,version='dr2'):
    """ Calculate the median spread between catalog and DES DR2.

    Parameters
    ----------
    filename : input filename (catalog file)
    nside    : output nside
    band     : band of interest
    plot     : pause and plot

    Returns
    -------
    upix, stat : unique pixel and value of statistic (median) in that pixel
    """
    if not os.path.exists(filename): 
        msg = "Couldn't find %s"%filename
        raise Exception(msg)

    #columns = [OBJECT_ID,'RA','DEC']
    columns = ['RA','DEC']
    spread,nepochs = bfields(['WAVG_SPREAD_MODEL','NEPOCHS'],band)
    mag = bfield('WAVG_MAG_PSF',band)
    columns += [spread, nepochs, mag]

    # Hack to get pixel location
    hpx = int(filename.split('_')[-1].split('.')[0])
    ra,dec = hp.pix2ang(NSIDE, hpx, lonlat=True)
    radius = np.degrees(hp.max_pixrad(NSIDE))

    msg = '%s (RA,DEC,RAD) = %.2f,%.2f,%.2f'%(os.path.basename(filename),ra,dec,radius)
    print(msg)

    # Load the DES catalog
    ext = get_des_catalog(hpx,version=version,columns=columns)
    if ext is None:
        msg = "WARNING: No DES objects in pixel: %s"%hpx
        print(msg)
        return np.array([],dtype=int), np.array([])

    # Select stars with 16 < mag < 20
    sel = (np.fabs(ext[spread])<0.002) & \
          (ext[mag]>16) & (ext[mag]<20) & \
          (ext[nepochs] > 1)
    ext = ext[sel]

    if len(ext) == 0:
        msg = "WARNING: No DES objects passing selection in pixel: %s"%hpx
        print(msg)
        return np.array([],dtype=int), np.array([])

    # Load the test catalog
    cat = load_infiles([filename],columns)

    # Select stars with 16 < mag < 20
    sel = (np.fabs(cat[spread])<0.002) & \
          (cat[mag]>16) & (cat[mag]<20) & \
          (cat[nepochs] > 1)
    cat = cat[sel]

    if len(cat) == 0:
        msg = "WARNING: No objects passing selection in: %s"%filename
        print(msg)
        return np.array([],dtype=int), np.array([])

    m = match_query(cat['RA'],cat['DEC'],ext['RA'],ext['DEC'])
    cut = 0.5 # arcsec
    sel = m[-1]*3600. < cut

    cat_match = cat[m[0][sel]]
    ext_match = ext[m[1][sel]]

    # Difference between catalog and external DES catalog
    diff  = cat_match[mag] - ext_match[mag]

    pix = hp.ang2pix(nside,cat_match['RA'],cat_match['DEC'],lonlat=True)
    upix = np.unique(pix)
    stat = nd.median(diff,labels=pix,index=upix)

    if False:
        plt.figure()
        plt.hist(diff)
        import pdb; pdb.set_trace()
        
    return upix,stat


def wavg_rms_photometry(filename,nside=64,band=None,plot=False):
    """ Calculate the median WAVG_MAGRMS in pixels. 
    The WAVG RMS is the unbiased estimate of sample variance:
    https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Reliability_weights 

    Parameters
    ----------
    filename : input filename (catalog file)
    nside    : output nside
    band     : band of interest
    plot     : pause and plot

    Returns
    -------
    upix, stat : unique pixel and value of statistic (median) in that pixel
    """
    if not os.path.exists(filename): 
        msg = "Couldn't find %s"%filename
        raise Exception(msg)

    columns = ['RA','DEC']
    spread,nepochs = bfields(['WAVG_SPREAD_MODEL','NEPOCHS'],band)
    mag,magerr,magrms = bfields(['WAVG_MAG_PSF','WAVG_MAGERR_PSF','WAVG_MAGRMS_PSF'],band)
    columns += [spread, nepochs, mag, magerr, magrms]

    # Hack to get pixel location
    hpx = int(filename.split('_')[-1].split('.')[0])
    #hpx = ang2pix(NSIDE, cat['RA'], cat['DEC'])
    ra,dec = hp.pix2ang(NSIDE, hpx, lonlat=True)
    msg = '%s (RA,DEC) = %.2f,%.2f'%(os.path.basename(filename),ra,dec)
    print(msg)

    cat = load_infiles([filename],columns)

    # Select stars with 16 < r < 18
    #(cat[mag] > 16) & (cat[mag] < 18) &\
    sel = (np.fabs(cat[spread]) < 0.002) & \
          (cat[mag] > 16) & (cat[mag] < 17) &\
          (cat[magrms] < 90) &\
          (cat[nepochs] > 1)
    cat = cat[sel]

    if len(cat) == 0:
        msg = "WARNING: No objects passing selection in: %s"%filename
        print(msg)
        return np.array([],dtype=int), np.array([])

    pix = hp.ang2pix(nside,cat['RA'],cat['DEC'],lonlat=True)
    upix = np.unique(pix)
    stat = nd.median(cat[magrms],labels=pix,index=upix)

    if False:
        plt.figure()
        plt.hist(cat[magrms],bins=50)
        import pdb; pdb.set_trace()
        
    return upix,stat

def hpx_rms_photometry(filename,nside=64,band=None,plot=False):
    """ Calculate the unweighted MAG RMS in pixels. 
    The RMS is the unweighted standard deviation of the measured magnitudes.

    Parameters
    ----------
    filename : input filename (hpx file)
    nside    : output nside
    band     : band of interest
    plot     : pause and plot

    Returns
    -------
    upix, stat : unique pixel and value of statistic (median) in that pixel
    """
    if not os.path.exists(filename): 
        msg = "Couldn't find %s"%filename
        raise Exception(msg)

    columns = ['RA','DEC']
    objid,spread,mag = ['QUICK_OBJECT_ID','SPREAD_MODEL','MAG_PSF']
    columns += [objid, spread, mag]

    # Hack to get pixel location
    hpx = int(filename.split('_')[-1].split('.')[0])
    #hpx = ang2pix(NSIDE, cat['RA'], cat['DEC'])
    ra,dec = hp.pix2ang(NSIDE, hpx, lonlat=True)
    msg = '%s (RA,DEC) = %.2f,%.2f'%(os.path.basename(filename),ra,dec)
    print(msg)

    cat = load_infiles([filename],columns)

    # Select stars with 16 < r < 18
    sel = (np.fabs(cat[spread]) < 0.002) & \
          (cat[mag] > 16) & (cat[mag] < 18)
    cat = cat[sel]

    if len(cat) == 0:
        msg = "WARNING: No objects passing selection in: %s"%filename
        print(msg)
        return np.array([],dtype=int), np.array([])

    # Calculate "per object" quantities
    df = pd.DataFrame(cat.byteswap().newbyteorder())
    grp = df.groupby(objid)
    ra = grp['RA'].median()
    dec = grp['DEC'].median()
    cts = grp[mag].count()
    rms = grp[mag].std()

    pix = hp.ang2pix(nside,ra,dec,lonlat=True)
    upix = np.unique(pix)
    stat = nd.median(rms,labels=pix,index=upix)

    if False:
        plt.figure()
        plt.hist(rms,bins=50,range=(np.nanmin(rms),np.nanmax(rms)))
        import pdb; pdb.set_trace()
        
    return upix,stat

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config',help='configuration file')
    parser.add_argument('-b','--band',default='r',choices=BANDS+['griz'])
    parser.add_argument('-n','--nside',default=128,type=int)
    parser.add_argument('-o','--outbase',default='photo')
    parser.add_argument('-p','--pix',default=None,type=int,action='append')
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('--type',choices=TYPES,default='gaia_edr3')
    parser.add_argument('--nproc',default=4,type=int)
    args = parser.parse_args()

    print('Calculating photometric offsets...')
    band = args.band

    config = yaml.safe_load(open(args.config))
    OBJECT_ID = config.get('objid',OBJECT_ID)
    NSIDE = config['nside']
    catdir  = config['catdir']
    catbase = config['catbase']
    filebase = os.path.join(catdir,catbase)

    # Output nside
    nside = args.nside
    #outdir = mkdir('release/photometry')
    outdir  = './'
    outbase = args.outbase

    kwargs = dict(nside=args.nside,band=band)
    if args.type in ('gaia','gaia_edr3'):
        func = gaia_photometry
        kwargs['version'] = 'edr3'
        outbase += '_gaia_%(version)s'%kwargs
        print("Calculating offsets to Gaia %(version)s..."%kwargs)
    elif args.type == 'gaia_dr2':
        func = gaia_photometry
        kwargs['version'] = 'dr2'
        outbase += '_gaia_%(version)s'%kwargs
        print("Calculating offsets to Gaia %(version)s..."%kwargs)
    elif args.type == 'des_dr2':
        func = des_photometry
        kwargs['version'] = 'dr2'
        outbase += '_des_%(version)s_%(band)s'%kwargs
        print("Calculating offsets to DES %(version)s..."%kwargs)
    elif args.type == 'wavg_rms':
        func = wavg_rms_photometry
        outbase += '_%s_%s'%(args.type,band)
        print("Calculating %s in %s band..."%(args.type,band))
    elif args.type == 'hpx_rms':
        filebase = os.path.join(config['hpxdir'],'{band:s}',config['hpxbase'])
        filebase = filebase.format(band=band)
        func = hpx_rms_photometry
        outbase += '_%s_%s'%(args.type,band)
        print("Calculating %s in %s band..."%(args.type,band))
    else:
        msg = "Unrecognized type: %s"%args.type
        raise Exception(msg)

    if args.pix is not None:
        pixels = args.pix
    else:
        pixels = np.arange(hp.nside2npix(NSIDE))

    # Remove this...
    if False:
        RA = [30,35]
        DEC = [-25,-20]
        ra, dec = hp.pix2ang(NSIDE,pixels,lonlat=True)
        sel = (ra>RA[0]) & (ra<RA[1]) & (dec>DEC[0]) & (dec<DEC[1])
        pixels = pixels[sel]
     
    if len(pixels) == 0:
        msg = "Invalid pixel: %s"%args.pix
        raise Exception(msg)
     
    filenames = [filebase%p for p in pixels]
    filenames = [f for f in filenames if os.path.exists(f)]

    if not len(filenames):
        msg = "No valid files found."
        raise Exception(msg)
    
    # Args must be tuple
    print("Processing %i files..."%len(filenames))

    arglist = zip(filenames)

    results = utils.multiproc(func,arglist,kwargs,processes=args.nproc)
     
    hpxmap = blank(nside)
    if None in results:
        print("WARNING: %i processes failed..."%results.count(None))

    for pix,stat in [r for r in results if r is not None]:
        hpxmap[pix] = stat
     
    hpxmap = np.ma.MaskedArray(hpxmap,np.isnan(hpxmap),fill_value=np.nan)

    outfile = os.path.join(outdir,outbase+'_n%i.fits.gz'%nside)
    print("Writing %s..."%outfile)
    hp.write_map(outfile,hpxmap,overwrite=True)

    q = [5,50,95]
    p = np.percentile(hpxmap.compressed(),q)
    print("Global Photometric Percentiles:")
    print('%s (%s%%)'%(p,q))

    print("Plotting %s..."%outfile)
    pngfile = outfile.replace('.fits.gz','.png')
    plot_photometry(outfile,pngfile)
    plt.ion()
