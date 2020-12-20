#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import glob
import os
from os.path import join
import yaml
import matplotlib
if os.getenv('TERM')=='screen' or not os.getenv('DISPLAY'):
    matplotlib.use('Agg')
from collections import OrderedDict as odict
from multiprocessing import Pool

import fitsio
import numpy as np
import scipy.ndimage as nd
import pylab as plt
import healpy as hp
import matplotlib.colors as colors

from ugali.utils.healpix import ang2pix, pix2ang
import ugali.utils.projector
from ugali.utils.projector import angsep
from ugali.utils.shell import mkdir

import utils
from const import OBJECT_ID, UNIQUE_ID, BANDS, BADMAG, NSIDES
from utils import bfields, load_infiles, get_vizier_catalog, get_local_catalog
from match import match_query

import plotting
from footprint import blank
from catalog import good_objects

COLUMNS = [OBJECT_ID,'RA','DEC']

def get_gaia_catalog(hpx,columns=['RA','DEC','PHOT_G_MEAN_FLUX']):
    dirname = '/data/des40.b/data/gaia/dr2/healpix/'
    basename = 'GaiaSource_%05d.fits'

    pixels = [hpx]
    filenames = [os.path.join(dirname,basename%p) for p in pixels]
    filenames = [f for f in filenames if os.path.exists(f)]
    cat = load_infiles(filenames,columns=columns)
    return cat

#def gaia_transform(g,r,i):
#    """ DES DR1 transformation to Gaia G. """
#    G = r - 0.10020396 + 0.14969636 * (g-i) - 0.01253971 * (g-i)**2 - 0.03451075 * (g-i)**3
#    return G

def gaia_transform(g,r,i,z):
    """ From Eli... This complex model behaves about as well as a random forest classifier for Gaia DR2.
    The magnitude dependence of the transformation is huge because of background errors in Gaia DR2.

    This is valid for 0 < g - i < 1.5
    """

    magConst = 2.5 / np.log(10.0)
    lambdaStd = np.array([4790.28076172, 6403.26367188, 7802.49755859, 9158.77441406])
    i0Std = np.array([0.16008162, 0.18297842, 0.17169334, 0.1337308])
    i1Std = np.array([2.09808350e-05, -1.22070312e-04, -1.08942389e-04, -8.01086426e-05])
    i10Std = i1Std / i0Std
    fudgeFactors = np.array([0.25, 1.0, 1.0, 0.25])
    nBands = lambdaStd.size

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

    # DR2 fit parameters
    pars = [ 1.43223290e+00,  1.50877061e+00, 8.43173013e-01, -5.99023967e-04, 
             4.06188382e-01,  3.11181978e-01, 2.51002598e-01,  1.00000000e-05,  
             4.94284725e-03,  1.80499806e-03]

    # EDR3: Notice that the last two parameters (describing the
    # r-offset curvature due to background issues) are much smaller.
    #pars = [ 2.61815727e+00,  2.69372875e+00,  1.45644592e+00, -5.99023051e-04,
    #         3.97535324e-01,  3.15794343e-01,  2.55484718e-01,  1.00000000e-05,
    #         8.30152817e-04, -3.57980758e-04]

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
    return mGDES

def gaia_photometry(catfile,nside=64,band=None,plot=False):
    if not os.path.exists(catfile): 
        msg = "Couldn't find %s"%catfile
        raise Exception(msg)

    #columns = [OBJECT_ID,'RA','DEC']
    columns = ['RA','DEC']
    spread,nepochs = ['WAVG_SPREAD_MODEL_R','NEPOCHS_R']
    mag_g,mag_r,mag_i,mag_z = bfields(['MAG_PSF'],['g','r','i','z'])
    #mag_g,mag_r,mag_i,mag_z = bfields(['WAVG_MAG_PSF'],['g','r','i','z'])
    columns += [spread, nepochs, mag_g, mag_r, mag_i, mag_z]

    # Hack to get pixel location
    hpx = int(catfile.split('_')[-1].split('.')[0])
    #hpx = ang2pix(NSIDE, cat['RA'], cat['DEC'])
    ra,dec = pix2ang(NSIDE, hpx)
    radius = np.degrees(hp.max_pixrad(NSIDE))

    msg = '%s (RA,DEC,RAD) = %.2f,%.2f,%.2f'%(os.path.basename(catfile),ra,dec,radius)
    print(msg)

    #print "Getting coadd catalog: DES"
    cat = load_infiles([catfile],columns)
    # Select stars with 16 < r < 20 and 0.0 < (g-i) < 1.5
    sel = (np.fabs(cat[spread])<0.002) & \
        (cat[mag_g]<90) & (cat[mag_r]<90) & (cat[mag_i]<90) & (cat[mag_z]<90) & \
        (cat[mag_r]>16) & (cat[mag_r]<20) & \
        ((cat[mag_g] - cat[mag_i]) > 0.0) & \
        ((cat[mag_g] - cat[mag_i]) < 1.5) & \
        (cat[nepochs] > 1)
    cat = cat[sel]

    if len(cat) == 0:
        msg = "WARNING: No objects passing selection in: %s"%catfile
        print(msg)
        return np.array([],dtype=int), np.array([])

    #msg = "Getting external catalog: %s"%catalog
    ext = get_gaia_catalog(hpx)

    m = match_query(cat['RA'],cat['DEC'],ext['RA'],ext['DEC'])

    # Use a fairly wide matching radius (2 arcsec)
    cut = 1.0
    sel = m[-1]*3600. < cut

    cat_match = cat[m[0][sel]]
    ext_match = ext[m[1][sel]]

    cat_G = gaia_transform(cat_match[mag_g],cat_match[mag_r],cat_match[mag_i],cat_match[mag_z])
    # Need to put Gaia on AB system
    ext_G = -2.5 * np.log10(ext_match['PHOT_G_MEAN_FLUX']) + 25.7934
    diff  = cat_G - ext_G

    pix = ang2pix(nside,cat_match['RA'],cat_match['DEC'])
    upix = np.unique(pix)
    stat = nd.median(diff,labels=pix,index=upix)

    if False:
        plt.figure()
        plt.hist(cat_G - ext_G)
        import pdb; pdb.set_trace()
        
    return upix,stat


if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config',help='Configuration file.')
    parser.add_argument('-o','--outbase',default='photo')
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('-n','--nside',default=128,type=int)
    parser.add_argument('-p','--pix',default=None,type=int,action='append')
    parser.add_argument('-b','--band',default=['r'],action='append',choices=BANDS)
    opts = parser.parse_args()

    config = yaml.load(open(opts.config))
    OBJECT_ID = config.get('objid',OBJECT_ID)
    NSIDE = config['nside']

    nside = opts.nside
    npix = hp.nside2npix(config['nside'])
    catdir = config['catdir']
    #catdir = '/data/des40.b/data/des/y3a2/gold/v2.2/healpix'

    outdir = mkdir('release/photometry')
    outbase = opts.outbase
    outbase += '_gaia'

    hpxmaps = odict()

    for band in opts.band:
        print('Calculating photometry for %s-band...'%band)

        if opts.pix is not None:
            pixels = sorted([p for p in opts.pix if len(glob.glob(catdir+'/*%05d.fits'%p))])
        else:
            pixels = sorted([p for p in range(npix) if len(glob.glob(catdir+'/*%05d.fits'%p))])
         
        if len(pixels) == 0:
            msg = "Invalid pixel: %s"%opts.pix
            raise Exception(msg)
         
        catfiles = [glob.glob(catdir+'/*%05d.fits'%p)[0] for p in pixels]
        #catfiles = ['cat/y1a1_gold_1_0_2_01376.fits']
         
        # Args must be tuple
        print("Processing %i files..."%len(catfiles))

        print('Calculating photometric offsets...')
        args = zip(catfiles)
        kwargs = dict(nside=opts.nside,band=band)
        func = gaia_photometry

        results = utils.multiproc(func,args,kwargs)
        #results = [func(*a,**kwargs) for a in args]
         
        hpxmap = blank(nside)
         
        if None in results:
            print("WARNING: %i processes failed..."%results.count(None))
        for pix,stat in [r for r in results if r is not None]:
            hpxmap[pix] = stat
         
        hpxmap = np.ma.MaskedArray(hpxmap,np.isnan(hpxmap),fill_value=np.nan)
        hpxmaps[band] = hpxmap

        outfile = join(outdir,outbase+'_%s_n%i.fits'%(band,nside))
        print("Writing %s..."%outfile)
        hp.write_map(outfile,hpxmap,overwrite=True)

        q = [5,50,95]
        p = np.percentile(hpxmap.compressed(),q)
        print("Global Median Photometry:")
        print('%s (%s%%)'%(p,q))
    plt.ion()
