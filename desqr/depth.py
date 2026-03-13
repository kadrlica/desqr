#!/usr/bin/env python
"""
Characterize photometric depth over the survey.
"""
import os,sys
from os.path import join
import glob
import time
import subprocess

import matplotlib
try:             os.environ['DISPLAY']
except KeyError: matplotlib.use('Agg')
matplotlib.use('Agg')
import pylab as plt

import numpy as np
import scipy.ndimage as nd
import scipy.stats
import yaml
import fitsio
import healpy as hp

import ugali.utils.projector
from ugali.utils.projector import angsep

from const import OBJECT_ID, UNIQUE_ID, BANDS, BADMAG, NSIDES
from footprint import blank, empty
import plotting
from plotting import draw_footprint, draw_peak
import utils
from utils import mkdir, logger
from utils import bfields, load_infiles

NSIDE = 1024

def draw_maglim_hist(hpxmap,**kwargs):
    DeprecationWarning()
    maglims = hpxmap[hpxmap > 0]
    kwargs.setdefault('bins',np.linspace(maglims.min(),maglims.max(),100))
    kwargs.setdefault('histtype','step')
    kwargs.setdefault('normed',True)
    kwargs.setdefault('lw',1.5)

    ax = plt.gca()
    n,b,p = ax.hist(maglims,**kwargs)
    ax.set_xlabel('Magnitude Limit')
    ax.set_ylabel('Normalized Pixel Counts')

    ax2 = ax.twinx()
    plt.hist(maglims,cumulative=-1,color='r',**kwargs)
    ax2.set_ylabel('Cumulative', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    quantiles = [5,50,95]
    percentiles = np.percentile(maglims,quantiles)
    for q,p in zip(quantiles,percentiles):
        draw_peak(p,color='r',label='%.2f (%g%%)'%(p,100-q))

    ax.set_xlim(kwargs['bins'].min(),kwargs['bins'].max())
    ax2.set_ylim(0,1)
    return quantiles,percentiles

def maglim_range(hpxmap):
    DeprecationWarning()
    pix = np.where(hpxmap > 0)
    sort = np.sort(hpxmap[pix])
    vmin = sort[int(len(sort) * 0.005)]
    vmax = sort[int(len(sort) * 0.995)]
    return vmin, vmax

def draw_maglim_pixel(hpxmap,**kwargs):
    DeprecationWarning()
    nside = hp.npix2nside(len(hpxmap))
    pix = np.where(hpxmap > 0)
    if len(pix[0]) == 0:
        logger.warning("No maglims found")
        return

    ra,dec = hp.pix2ang(nside,pix,lonlat=True)
    ra_center,dec_center = np.median(ra),np.median(dec)

    vmin,vmax = maglim_range(hpxmap)

    kwargs.setdefault('rot',(ra_center, dec_center, 0.))
    kwargs.setdefault('min',vmin)
    kwargs.setdefault('max',vmax)
    hp.gnomview(hpxmap,**kwargs)

def draw_magerr(mag,magerr,**kwargs):
    DeprecationWarning()
    kwargs.setdefault('fmt','o')
    bins = kwargs.pop('bins',np.linspace(16,25,100))
                      
    ax = plt.gca()
    labels = np.digitize(mag,bins)
    index = np.unique(labels)
    centers = ((bins[1:]+bins[:-1])/2.)[index]
    median = nd.median(magerr,labels=labels,index=index)
    mean = nd.mean(magerr,labels=labels,index=index)
    std = nd.standard_deviation(magerr,labels=labels,index=index)
     
    ax.errorbar(centers,mean,yerr=std,**kwargs)
    return median,mean,std

def plot_depth(filename,outfile=None,survey='delve'):
    logger.info("Reading %s..."%filename)
    hpxmap = hp.read_map(filename)

    label = r'Magnitude Limit (mag)'
    cbar_kwargs = dict(label=label)
    hpxmap_kwargs = dict(xsize=2000)
    fig,axes,smap = plotting.plot_hpxmap_hist(hpxmap,survey,cbar_kwargs,hpxmap_kwargs)
    #axes[0].annotate('%s band'%band, (0.05,0.93), xycoords='axes fraction')
    axes[1].set_xlabel(label)

    if outfile is None: 
        outfile=os.path.basename(filename).split('.')[0]+'.png'
    logger.info("Writing %s..."%outfile)
    plt.savefig(outfile,bbox_inches='tight')


def calculate_depth(filename, nside=NSIDE, snr=10.0,
                    mag='MAG_PSF', magerr='MAGERR_PSF', spread='WAVG_SPREAD_MODEL',
                    stargal=None):
    """ Calculate the photometric depth with arbitrary magnitude and object selection. 

    Parameters
    ----------
    filename : input filename
    nside : nside of output map
    snr : signal-to-noise of maglim
    mag : magnitude variable prefix
    magerr : magnitude error variable prefix
    spread : star-galaxy classifier
    stargal : object selection (None = none; 0 = stars; 1 = galaxies)
    """
    MAGS = bfields(mag,BANDS)
    MAGERRS = bfields(magerr,BANDS)
    SPREADS = bfields(spread,BANDS)

    # From propagation of errors:
    # mag = -2.5 * log10(flux)
    # magerr = -2.5/ln(10) * fluxerr/flux
    mag_snr = (2.5/np.log(10)) / (snr)

    logger.info(filename)
    ret = dict()
    for band,mag,magerr,spread in zip(BANDS,MAGS,MAGERRS,SPREADS):
        columns = ['RA','DEC'] + [mag, magerr, spread]
        data = fitsio.read(filename,columns=columns)

        h, edges = np.histogram(data[mag], bins=np.arange(17, 30, 0.1))
        mag_bright = edges[np.argmax(h)] - 3.
        mag_faint = edges.max()

        cut = (data[mag] > mag_bright) & (data[mag] < mag_faint)
        if stargal == 0:
            cut &= (np.fabs(data[spread]) < 0.002)
        elif stargal == 1:
            cut &= (data[spread] > 0.003)

        d = data[cut]
        if len(d) < 2:
            logger.warning("Insufficent objects in %s-band"%band)
            ret[band] = [np.array([],dtype=int),np.array([])]            
            continue

        pix = hp.ang2pix(nside,d['RA'],d['DEC'],lonlat=True)

        match = ugali.utils.projector.match(d['RA'],d['DEC'],d['RA'],d['DEC'],nnearest=2)

        delta_mag = d[mag][match[1]] - d[mag][match[0]]
        delta_log_magerr = np.log10(d[magerr][match[1]]) - np.log10(d[magerr][match[0]])

        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = delta_log_magerr / delta_mag
            cut_nan_inf = np.isfinite(ratio) & (delta_mag > 0.5)

        if cut_nan_inf.sum() < 2:
            logger.warning("Insufficent objects in %s-band"%band)
            ret[band] = [np.array([],dtype=int),np.array([])]            
            continue

        kde = scipy.stats.gaussian_kde(ratio[cut_nan_inf])

        values = np.linspace(0., 1., 1000)
        kde_values = kde.evaluate(values)
        slope = values[np.argmax(kde_values)]
        
        maglims = d[mag] - ((np.log10(d[magerr]) - np.log10(mag_snr)) / slope)

        hpx = np.unique(pix)
        maglim = nd.median(maglims,labels=pix,index=hpx)
        ret[band] = [hpx,maglim]
        
    return ret

def depth_psf(filename, nside=NSIDE, snr=10.0):
    return calculate_depth(filename, nside, snr,
                           mag='MAG_PSF', magerr='MAGERR_PSF', spread='WAVG_SPREAD_MODEL',
                           stargal=0)

def depth_auto(filename, nside=NSIDE, snr=10.0):
    return calculate_depth(filename, nside, snr,
                           mag='MAG_AUTO', magerr='MAGERR_AUTO', spread='WAVG_SPREAD_MODEL',
                           stargal=None)

def depth_psf_mag(filename, nside=NSIDE, snr=10.0):
    return calculate_depth(filename, nside, snr,
                           mag='PSF_MAG_APER_8', magerr='PSF_MAG_ERR_APER_8', spread='SPREAD_MODEL',
                           stargal=0)

def depth_bdf_mag(filename, nside=NSIDE, snr=10.0):
    return calculate_depth(filename, nside, snr,
                           mag='BDF_MAG', magerr='BDF_MAG', spread='SPREAD_MODEL',
                           stargal=None)


def teff(infile,nside=NSIDE,mode='median'):
    TEFFS = bfields('TEFF',BANDS)

    logger.info(infile)
    ret = dict()
    for band,teff in zip(BANDS,TEFFS):
        data = fitsio.read(infile,columns=['RA','DEC',teff])
        pix = hp.ang2pix(nside,data['RA'],data['DEC'],lonlat=True)
        hpx = np.unique(pix)        

        if mode.lower() == 'median':
            teff_value = nd.median(data[teff],labels=pix,index=hpx)
        elif mode.lower() == 'mean':
            teff_value = nd.mean(data[teff],labels=pix,index=hpx)
        else:
            msg = 'Unrecognized mode: %s'%mode
            raise Exception(msg)

        ret[band] = [hpx,maglim]
    return ret


if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('configfile',help='config file')
    parser.add_argument('-n','--nside',default=NSIDE,type=int,
                        help='output nside')
    parser.add_argument('-v','--verbose',action='store_true',
                        help='output verbosity')
    parser.add_argument('--type',default='psf',
                        choices=['psf','auto','psf_mag','bdf_mag'],
                        help='which mag to use')
    parser.add_argument('--snr',default=10,type=float,
                        help='which snr to calculate')
    parser.add_argument('--nproc',default=4,type=int)
    args = parser.parse_args()

    if args.verbose: logger.setLevel(logger.INFO)
    
    config = yaml.safe_load(open(args.configfile))
    BANDS = config['bands']
    NSIDE = args.nside
    survey = config.get('survey')
    outdir = mkdir('release/depth')

    filenames = sorted(glob.glob(config['catdir'] + '/*_*.fits'))
    #filenames = filenames[:100]
    arglist = list(zip(filenames,))
    kwargs = dict(nside=args.nside,snr=args.snr)

    if args.type == 'psf':
        func = depth_psf
    elif args.type == 'auto': 
        func = depth_auto
    elif args.type == 'psf_mag':
        func = depth_psf_mag
    elif args.type == 'bdf_mag': 
        func = depth_bdf_mag

    else: raise Exception('Unrecognized magnitude: %s'%args.type)

    out = utils.multiproc(func,arglist,kwargs,processes=args.nproc)

    basename = 'maglim_%s_%isig'%(args.type,args.snr)
    hpxmaps = dict()
    for b in BANDS:
        logger.info("Filling %s-band..."%b)
        hpxmaps[b] = hpxmap = blank(NSIDE)
        for i,maglims in enumerate(out):
            logger.info(str(filenames[i]))
            hpxmap[maglims[b][0]] = maglims[b][1]

        outfile = join(outdir,basename+'_%s_n%i.fits.gz'%(b,NSIDE))
        logger.info("Writing %s..."%outfile)
        hp.write_map(outfile,hpxmap,overwrite=True)

    for b in BANDS:
        outfile = join(outdir,basename+'_%s_n%i.fits.gz'%(b,NSIDE))
        logger.info("Plotting %s..."%outfile)
        pngfile = outfile.replace('.fits.gz','.png')
        plot_depth(outfile,pngfile,survey=survey)



###def depth_psf(infile,nside=NSIDE,snr=10.):
###    """ MAG_PSF depth with star/galaxy selection """
###    #MAGS = bfields('MAG_PSF',BANDS)
###    #MAGERRS = bfields('MAGERR_PSF',BANDS)
###    #SPREADS = bfields('WAVG_SPREAD_MODEL',BANDS)
###    MAGS = bfields('WAVG_MAG_PSF',BANDS)
###    MAGERRS = bfields('WAVG_MAGERR_PSF',BANDS)
###    SPREADS = bfields('WAVG_SPREAD_MODEL',BANDS)
### 
###    # From propagation of errors:
###    # mag = -2.5 * log10(flux)
###    # magerr = -2.5/ln(10) * fluxerr/flux
###    mag_snr = (2.5/np.log(10)) / (snr)
### 
###    logger.info(infile)
###    ret = dict()
###    for band,mag,magerr,spread in zip(BANDS,MAGS,MAGERRS,SPREADS):
###        data = fitsio.read(infile,columns=['RA','DEC',mag,magerr,spread])
### 
###        h, edges = np.histogram(data[mag], bins=np.arange(17, 30, 0.1))
###        mag_bright_end = edges[np.argmax(h)] - 3.
### 
###        cut = (np.fabs(data[spread]) < 0.002) & (data[mag] > mag_bright_end) & (data[mag] < 30.)
### 
###        d = data[cut]
###        if len(d) < 2:
###            logger.warning("Insufficent objects in %s-band"%band)
###            ret[band] = [np.array([],dtype=int),np.array([])]            
###            continue
### 
###        pix = hp.ang2pix(nside,d['RA'],d['DEC'],lonlat=True)
### 
###        match = ugali.utils.projector.match(d['RA'],d['DEC'],d['RA'],d['DEC'],nnearest=2)
### 
###        delta_mag = d[mag][match[1]] - d[mag][match[0]]
###        delta_log_magerr = np.log10(d[magerr][match[1]]) - np.log10(d[magerr][match[0]])
### 
###        old = np.seterr(divide='ignore',invalid='ignore')
###        ratio = delta_log_magerr / delta_mag
###        cut_nan_inf = np.isfinite(ratio) & (delta_mag > 0.5)
###        np.seterr(**old)
### 
###        if cut_nan_inf.sum() < 2:
###            logger.warning("Insufficent objects in %s-band"%band)
###            ret[band] = [np.array([],dtype=int),np.array([])]            
###            continue
### 
###        kde = scipy.stats.gaussian_kde(ratio[cut_nan_inf])
### 
###        values = np.linspace(0., 1., 1000)
###        kde_values = kde.evaluate(values)
###        slope = values[np.argmax(kde_values)]
###        
###        maglims = d[mag] - ((np.log10(d[magerr]) - np.log10(mag_snr)) / slope)
### 
###        hpx = np.unique(pix)
###        maglim = nd.median(maglims,labels=pix,index=hpx)
###        ret[band] = [hpx,maglim]
###    return ret
### 
###def depth_auto(infile,nside=NSIDE,snr=10.):
###    """ MAG_AUTO depth without star/galaxy selection """
### 
###    MAGS = bfields('MAG_AUTO',BANDS)
###    MAGERRS = bfields('MAGERR_AUTO',BANDS)
### 
###    # From propagation of errors:
###    # mag = -2.5 * log10(flux)
###    # magerr = -2.5/ln(10) * fluxerr/flux
###    mag_snr = (2.5/np.log(10)) / (snr)
### 
###    logger.info(infile)
###    ret = dict()
###    for band,mag,magerr in zip(BANDS,MAGS,MAGERRS):
###        data = fitsio.read(infile,columns=['RA','DEC',mag,magerr])
### 
###        h, edges = np.histogram(data[mag], bins=np.arange(17, 30, 0.1))
###        mag_bright = edges[np.argmax(h)] - 3.
###        mag_faint = edges.max()
###        
###        cut = (data[mag] > mag_bright) & (data[mag] < mag_faint)
### 
###        d = data[cut]
###        if len(d) < 2:
###            logger.warning("Insufficent objects in %s-band"%band)
###            ret[band] = [np.array([],dtype=int),np.array([])]            
###            continue
### 
###        pix = hp.ang2pix(nside,d['RA'],d['DEC'],lonlat=True)
### 
###        match = ugali.utils.projector.match(d['RA'],d['DEC'],d['RA'],d['DEC'],nnearest=2)
### 
###        delta_mag = d[mag][match[1]] - d[mag][match[0]]
###        delta_log_magerr = np.log10(d[magerr][match[1]]) - np.log10(d[magerr][match[0]])
### 
###        old = np.seterr(divide='ignore',invalid='ignore')
###        ratio = delta_log_magerr / delta_mag
###        cut_nan_inf = np.isfinite(ratio) & (delta_mag > 0.5)
###        np.seterr(**old)
### 
###        if cut_nan_inf.sum() < 2:
###            logger.warning("Insufficent objects in %s-band"%band)
###            ret[band] = [np.array([],dtype=int),np.array([])]            
###            continue
### 
###        kde = scipy.stats.gaussian_kde(ratio[cut_nan_inf])
### 
###        values = np.linspace(0., 1., 1000)
###        kde_values = kde.evaluate(values)
###        slope = values[np.argmax(kde_values)]
###        
###        maglims = d[mag] - ((np.log10(d[magerr]) - np.log10(mag_snr)) / slope)
### 
###        hpx = np.unique(pix)
###        maglim = nd.median(maglims,labels=pix,index=hpx)
###        ret[band] = [hpx,maglim]
###    return ret

