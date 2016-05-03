#!/usr/bin/env python
import glob
import os
from os.path import join
import matplotlib
try:             os.environ['DISPLAY']
except KeyError: matplotlib.use('Agg')

import fitsio
import numpy as np
import scipy.ndimage as nd
import pylab as plt
import matplotlib.colors as colors
import healpy

from ugali.utils.healpix import ang2pix,pix2ang
from ugali.utils.shell import mkdir

from const import OBJECT_ID, UNIQUE_ID, BANDS, BADMAG, NSIDES
from utils import bfields, load_infiles

def draw_peak(peak,**kwargs):
    kwargs.setdefault('ls','--')
    kwargs.setdefault('label','%.1f '%(peak))
    ax = plt.gca()
    ax.axvline(peak,**kwargs)

def draw_hist(skymap,**kwargs):
    if isinstance(skymap,np.ma.MaskedArray):
        pix = np.where(~skymap.mask)
    else:
        pix = np.where((np.isfinite(skymap)) & (skymap !=healpy.UNSEEN))

    data = skymap[pix]
    kwargs.setdefault('bins',np.linspace(data.min(),data.max(),100))
    kwargs.setdefault('histtype','step')
    kwargs.setdefault('normed',True)
    kwargs.setdefault('lw',1.5)

    ax = plt.gca()
    n,b,p = ax.hist(data,**kwargs)

    ax2 = ax.twinx()
    plt.hist(data,cumulative=-1,color='r',**kwargs)
    ax2.set_ylabel('Cumulative', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    ax2.set_ylim(0,1)

    plt.sca(ax)
    quantiles = [5,50,95]
    percentiles = np.percentile(data,quantiles)
    for q,p in zip(quantiles,percentiles):
        draw_peak(p,color='r',label='%.1f (%g%%)'%(p,100-q))

    ax.set_xlim(kwargs['bins'].min(),kwargs['bins'].max())

    return quantiles,percentiles

def draw_footprint(skymap,proj='car',**kwargs):

    if not isinstance(skymap,np.ma.MaskedArray):
        mask = ~np.isfinite(skymap) | (skymap==healpy.UNSEEN)
        skymap = np.ma.MaskedArray(skymap,mask=mask)
    pix = np.where(~skymap.mask)

    vmin,vmax = np.percentile(skymap[pix],[0.5,99.5])

    kwargs.setdefault('vmin',vmin)
    kwargs.setdefault('vmax',vmax)
    kwargs.setdefault('rasterized',True)
    kwargs.setdefault('cmap','jet')

    nside = healpy.npix2nside(len(skymap))

    xmin,xmax = -70,110
    ymin,ymax = -70,10
    # desgw
    #xmin,xmax = 30,180
    #ymin,ymax = -80,-50
    # maglites
    xmin,xmax = 80,280
    ymin,ymax = -90,-50

    delta = 0.1
    xx,yy = np.meshgrid(np.arange(xmin,xmax,delta),np.arange(ymin,ymax,delta))
    pp = ang2pix(nside,xx,yy)

    ax = plt.gca()
    im = ax.pcolormesh(xx[::-1],yy,skymap[pp],**kwargs)
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('DEC (deg)')
    #ax.set_xticks([90,60,30,0,-30,-60])
    #ax.set_yticks([-70,-50,-30,-10,10])
    #ax.set_xticks(np.arange(180,-180-30,-30))
    #ax.set_yticks(np.arange(-90,90+20,20))
    ax.set_xlim(xmax,xmin)
    ax.set_ylim(ymin,ymax)
    ax.grid()
    return im

def draw_pixel(skymap,**kwargs):
    if isinstance(skymap,np.ma.MaskedArray):
        pix = np.where(~skymap.mask)
    else:
        pix = np.where((np.isfinite(skymap)) & (skymap !=healpy.UNSEEN))

    vmin,vmax = np.percentile(skymap[pix],[0.5,99.5])
    kwargs.setdefault('vmin',vmin)
    kwargs.setdefault('vmax',vmax)
    kwargs.setdefault('rasterized',True)
    kwargs.setdefault('cmap','jet')

    nside = healpy.npix2nside(len(skymap))
    pixrad = np.degrees(healpy.max_pixrad(nside))
    ra,dec = pix2ang(nside,pix)

    xmin,xmax = ra.min()-pixrad,ra.max()+pixrad
    ymin,ymax = dec.min()-pixrad,dec.max()+pixrad

    delta = 0.01
    xx,yy = np.meshgrid(np.arange(xmin,xmax,delta),np.arange(ymin,ymax,delta))
    pp = ang2pix(nside,xx,yy)

    ax = plt.gca()
    im = ax.pcolormesh(xx[::-1],yy,skymap[pp],**kwargs)
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('DEC (deg)')
    ax.set_xlim(xmax,xmin)
    ax.set_ylim(ymin,ymax)
    ax.grid()
    return im


if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    opts = parser.parse_args()
