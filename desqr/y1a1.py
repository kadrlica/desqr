#!/usr/bin/env python
import os
from os.path import join
import glob
import fitsio
from collections import OrderedDict as odict
import matplotlib
try:             os.environ['DISPLAY']
except KeyError: matplotlib.use('Agg')

from multiprocessing import Pool

import pylab as plt
import numpy as np
import scipy.stats
from scipy.stats import norm
import healpy

from utils import bfields
from const import BANDS,COLORS
from footprint import blank
import footprint
import plotting

import ugali.utils.projector
from ugali.utils.shell import mkdir


COLUMNS = ['RA','DEC']+bfields(['WAVG_MAG_PSF','WAVG_MAGERR_PSF'],BANDS)
Y1A1_COLUMNS = COLUMNS
Y2Q1_COLUMNS = COLUMNS + ['WAVG_SPREAD_MODEL_R']
NSIDE = 16

def residuals(pixel,nside=NSIDE,plot=False,y1a1dir='y1a1/v1/hpx/',y2q1dir='cat/'):
    y1a1file = os.path.join(y1a1dir,'cat_hpx_%05d.fits'%pixel)
    y2q1file = os.path.join(y2q1dir,'cat_hpx_%05d.fits'%pixel)
    print y1a1file,y2q1file

    y1a1 = fitsio.read(y1a1file,columns=Y1A1_COLUMNS)
    y2q1 = fitsio.read(y2q1file,columns=Y2Q1_COLUMNS)
    y2q1 = y2q1[(y2q1['WAVG_SPREAD_MODEL_R'] < 0.002)]

    #hpx = ang2pix(nside,y1a1['RA'],y1a1['DEC'])    
    #y1a1 = recfuncs.rec_append_fields('HPX',hpx,dtypes=int)
    # 
    #hpx = ang2pix(nside,y2q1['RA'],y2q1['DEC'])    
    #y2q1 = recfuncs.rec_append_fields('HPX',hpx,dtypes=int)

    #if plot:
    #    
    #    fig,axes = plt.subplots(1,2,figsize=(12,5))

    kwargs = dict(histtype='step',lw=1.5,normed=True)
    ret = odict()
    for band in BANDS:
        mag,magerr = bfields(['WAVG_MAG_PSF','WAVG_MAGERR_PSF'],band)
        color = COLORS[band]

        y1 = y1a1[(y1a1[mag] < 22)&(y1a1[mag] > 17)]
        y2 = y2q1[(y2q1[mag] < 22)&(y2q1[mag] > 17) ]       

        match = ugali.utils.projector.match(y1['RA'],y1['DEC'],y2['RA'],y2['DEC'])
        sepsec = 3600.*match[-1]
        sel = (sepsec < 1.0)
        idx1 = match[0][sel]
        idx2 = match[1][sel]

        y1_match = y1[idx1]
        y2_match = y2[idx2]
        
        res = (y2[mag][idx2] - y1[mag][idx1])
        res_clip,lo,hi = scipy.stats.sigmaclip(res,5,5)
        
        mu,sigma = norm.fit(res_clip)
        median = np.median(res_clip)
        ret[band] = (median,mu,sigma)

        if plot:
            bins = np.linspace(-0.1,0.1,100); centers = (bins[1:]+bins[:-1])/2.
            kwargs['bins'] = bins
            axes[0].hist(res,color=color,**kwargs)
            mu,sigma = norm.fit(res[(res>bins.min())&(res<bins.max())])
            label = r'$%s\ (\mu=%.2f,\sigma=%.2f)$'%(band,mu,sigma)
            axes[0].plot(centers,norm.pdf(centers,mu,sigma),color=color,label=label)

    if plot:
        plt.sca(axes[0])
        plt.legend(fontsize=8)
        plt.sca(axes[1])
        plt.legend(fontsize=8)

    return pixel,ret

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('-p','--pix',default=None,type=int,action='append')
    parser.add_argument('-n','--nside',default=16,type=int)
    opts = parser.parse_args()

    nside = opts.nside
    npix = healpy.nside2npix(16)
    outdir = mkdir('release/compare_gold')
    y1a1dir = 'y1a1/v1/hpx'
    y2q1dir = 'cat/'
    if opts.pix is not None:
        pixels = sorted([p for p in opts.pix if len(glob.glob(y1a1dir+'/*%05d.fits'%p))])
    else:
        pixels = sorted([p for p in range(npix) if len(glob.glob(y1a1dir+'/*%05d.fits'%p)) and len(glob.glob(y2q1dir+'/*%05d.fits'%p)) ])

    if len(pixels) == 0:
        msg = "Invalid pixel: %s"%opts.pix
        raise Exception(msg)
    
    args = [pix for pix in pixels]
    p = Pool(maxtasksperchild=1)
    out = p.map(residuals,args)

median_skymaps = odict()
mean_skymaps = odict()
std_skymaps = odict()
for band in BANDS:    
    median_skymap = blank(nside)
    mean_skymap = blank(nside)
    std_skymap = blank(nside)
    for pix,val in out:
        median_skymap[pix] = val[band][0]
        mean_skymap[pix] = val[band][1]
        std_skymap[pix] = val[band][2]
    median_skymaps[band] = median_skymap
    mean_skymaps[band] = mean_skymap
    std_skymaps[band] = std_skymap



for band in BANDS:    
    plt.figure()
    im = plotting.draw_footprint(median_skymaps[band])
    plt.colorbar(im,label=r'Median Offset (mag)')
    plt.title(r'Median Magnitude Offset (%s-band)'%band)
    outfile = join(outdir,'y2q1_y1a1_gold_median_%s_n%i_car.png'%(band,nside))
    plt.savefig(outfile,bbox_inches='tight')
    plt.figure()
    im = plotting.draw_footprint(mean_skymaps[band])
    plt.colorbar(im,label=r'Mean Offset (mag)')
    plt.title(r'Mean Magnitude Offset (%s-band)'%band)
    outfile = join(outdir,'y2q1_y1a1_gold_mean_%s_n%i_car.png'%(band,nside))
    plt.savefig(outfile,bbox_inches='tight')
    plt.figure()
    im = plotting.draw_footprint(std_skymaps[band])
    plt.colorbar(im,label=r'Standard Deviation (mag)')
    plt.title(r'Standard Deviation Magnitude Offset (%s-band)'%band)
    outfile = join(outdir,'y2q1_y1a1_gold_std_%s_n%i_car.png'%(band,nside))
    plt.savefig(outfile,bbox_inches='tight')
