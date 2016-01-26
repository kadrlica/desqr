#!/usr/bin/env python
import glob
import os
import time
from os.path import join
import subprocess
import matplotlib
try:             os.environ['DISPLAY']
except KeyError: matplotlib.use('Agg')

import fitsio
import numpy as np
import scipy.ndimage as nd
import scipy.stats
import pylab as plt
import matplotlib.colors as colors
import healpy

from multiprocessing import Pool

from ugali.utils.healpix import ang2pix,pix2ang
from ugali.utils.projector import angsep
from ugali.utils.shell import mkdir
import ugali.utils.projector

from const import OBJECT_ID, UNIQUE_ID, BANDS, BADMAG, NSIDES
from utils import bfields, load_infiles
from footprint import blank,unseen
import footprint
from plotting import draw_footprint, draw_peak

NSIDE = 1024

def draw_maglim_hist(skymap,**kwargs):
    maglims = skymap[skymap > 0]
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

def maglim_range(skymap):
    pix = np.where(skymap > 0)
    sort = np.sort(skymap[pix])
    vmin = sort[int(len(sort) * 0.005)]
    vmax = sort[int(len(sort) * 0.995)]
    return vmin, vmax

def draw_maglim_pixel(skymap,**kwargs):
    nside = healpy.npix2nside(len(skymap))
    pix = np.where(skymap > 0)
    if len(pix[0]) == 0:
        print "No maglims found"
        return

    ra,dec = pix2ang(nside,pix)
    ra_center,dec_center = np.median(ra),np.median(dec)

    vmin,vmax = maglim_range(skymap)

    kwargs.setdefault('rot',(ra_center, dec_center, 0.))
    kwargs.setdefault('min',vmin)
    kwargs.setdefault('max',vmax)
    healpy.gnomview(skymap,**kwargs)

def draw_maglim_footprint(skymap,**kwargs):
    pass

def draw_magerr(mag,magerr,**kwargs):
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


def depth(infile,nside=NSIDE,signal_to_noise=10.):
    MAGS = bfields('MAG_PSF',BANDS)
    MAGERRS = bfields('MAGERR_PSF',BANDS)
    SPREADS = bfields('WAVG_SPREAD_MODEL',BANDS)

    # From propagation of errors:
    # mag = -2.5 * log10(flux)
    # magerr = -2.5/ln(10) * fluxerr/flux
    mag_snr = (2.5/np.log(10)) / (signal_to_noise)

    print infile
    ret = dict()
    for band,mag,magerr,spread in zip(BANDS,MAGS,MAGERRS,SPREADS):
        data = fitsio.read(infile,columns=['RA','DEC',mag,magerr,spread])

        h, edges = np.histogram(data[mag], bins=np.arange(17, 30, 0.1))
        mag_bright_end = edges[np.argmax(h)] - 3.

        cut = (np.fabs(data[spread]) < 0.002) & (data[mag] > mag_bright_end) & (data[mag] < 30.)

        d = data[cut]
        if len(d) < 2:
            print "WARNING: Insufficent objects in %s-band"%band
            ret[band] = [np.array([],dtype=int),np.array([])]            
            continue

        pix = ang2pix(nside,d['RA'],d['DEC'])

        match = ugali.utils.projector.match(d['RA'],d['DEC'],d['RA'],d['DEC'],nnearest=2)

        delta_mag = d[mag][match[1]] - d[mag][match[0]]
        delta_log_magerr = np.log10(d[magerr][match[1]]) - np.log10(d[magerr][match[0]])

        old = np.seterr(divide='ignore',invalid='ignore')
        ratio = delta_log_magerr / delta_mag
        cut_nan_inf = np.isfinite(ratio) & (delta_mag > 0.5)
        np.seterr(**old)

        if cut_nan_inf.sum() < 2:
            print "WARNING: Insufficent objects in %s-band"%band
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

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    opts = parser.parse_args()

    outdir = mkdir('release/depth')

    infiles = sorted(glob.glob('cat/cat_hpx_*.fits'))
    p = Pool(maxtasksperchild=1)
    out = p.map(depth,infiles)

    skymaps = dict()
    for b in BANDS:
        skymap = blank(NSIDE)
        for maglims in out:
            skymap[maglims[b][0]] = maglims[b][1]
        skymaps[b] = np.ma.MaskedArray(skymap,np.isnan(skymap),fill_value=np.nan)
        outfile = join(outdir,'y2q1_maglim_%s_n%i_ring.fits'%(b,NSIDE))
        healpy.write_map(outfile,skymaps[b].data)
        subprocess.call('gzip -f %s'%outfile,shell=True)
        
    out = dict()
    outstr = '|_. Band |_. Footprint |_. Distribution |_. Magnitude Limit |\n'
    template = '|_. %(band)s |{{thumbnail(%(map)s, size=300)}}|{{thumbnail(%(hist)s, size=300)}}|_. %(maglim)s |\n'
     
    for b in BANDS:
        skymap = skymaps[b]
        out['band'] = b

        plt.figure()
        vmin,vmax = maglim_range(skymap)
        im = footprint.draw_footprint(skymap,vmin=vmin,vmax=vmax)
        plt.colorbar(im,label=r'Magnitude Limit (mag)')
        plt.title(r'10$\sigma$ Limiting Magnigude (%s-band)'%b)
        outfile = join(outdir,'y2q1_maglim_%s_n%i_car.png'%(b,NSIDE))
        plt.savefig(outfile,bbox_inches='tight')
        out['map'] = os.path.basename(outfile)

        plt.figure()
        q,p = draw_maglim_hist(skymap)
        plt.title('Magnitude Limits (%s-band)'%b)
        plt.legend(loc='upper left')
        outfile = join(outdir,'y2q1_maglim_%s_hist.png'%b)
        plt.savefig(outfile,bbox_inches='tight')
        out['hist'] = os.path.basename(outfile)
        out['maglim'] = '/'.join(['%.2f'%_p for _p in p])
        outstr += template%out
    print outstr
