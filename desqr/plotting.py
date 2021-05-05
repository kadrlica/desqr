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
import healpy as hp

from ugali.utils.shell import mkdir

from const import OBJECT_ID, UNIQUE_ID, BANDS, BADMAG, NSIDES
from utils import bfields, load_infiles

def draw_peak(peak,**kwargs):
    kwargs.setdefault('ls','--')
    kwargs.setdefault('label','%.1f '%(peak))
    ax = plt.gca()
    ax.axvline(peak,**kwargs)

def draw_hist(hpxmap,**kwargs):
    if isinstance(hpxmap,np.ma.MaskedArray):
        pix = np.where(~hpxmap.mask)
    else:
        pix = np.where((np.isfinite(hpxmap)) & (hpxmap !=hp.UNSEEN))

    data = hpxmap[pix]
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

def draw_footprint(hpxmap,proj='car',**kwargs):
    """
    Draw plot of footprint.
    """
    if not isinstance(hpxmap,np.ma.MaskedArray):
        mask = ~np.isfinite(hpxmap) | (hpxmap==hp.UNSEEN)
        hpxmap = np.ma.MaskedArray(hpxmap,mask=mask)
    pix = np.where(~hpxmap.mask)

    vmin,vmax = np.percentile(hpxmap[pix],[0.5,99.5])
    #vmin,vmax = np.percentile(hpxmap[pix],[0.5,99])

    kwargs.setdefault('vmin',vmin)
    kwargs.setdefault('vmax',vmax)
    kwargs.setdefault('rasterized',True)
    kwargs.setdefault('cmap','viridis')
    extent = kwargs.pop('extent',[180,-180,-90,20])
    xmin,xmax,ymin,ymax = extent

    nside = hp.npix2nside(len(hpxmap))

    steps = 500
    xx,yy = np.meshgrid(np.linspace(xmin,xmax,steps),
                        np.linspace(ymin,ymax,steps))
    pp = hp.ang2pix(nside,xx,yy,lonlat=True)

    ax = plt.gca()
    im = ax.pcolormesh(xx[::-1],yy,hpxmap[pp],**kwargs)
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('DEC (deg)')
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    ax.grid(ls=':',color='black',lw=0.5)
    return im, xx, yy, hpxmap[pp]

def draw_skymap(hpxmap,**kwargs):
    from skymap.survey import SurveySkymap
    smap = SurveySkymap()
    return smap.draw_hpxmap(hpxmap,**kwargs)

def draw_des(hpxmap,**kwargs):
    """ 
    Draw DES footprint:
       110 > RA > -70, -70 < DEC < 10
    """
    from skymap.survey import DESSkymap
    smap = DESSkymap()
    smap.draw_des()
    return smap.draw_hpxmap(hpxmap,**kwargs)

def draw_maglites(hpxmap,**kwargs):
    """ 
    Draw MagLiteS footprint:
       280 > RA > 80, -90 < DEC < -50
    """
    from skymap.survey import MaglitesSkymap
    smap = MaglitesSkymap()
    smap.draw_maglites()
    return smap.draw_hpxmap(**kwargs)

def draw_bliss(hpxmap,**kwargs):
    """ 
    Draw BLISS footprint:
       360 > RA > 120, -60 < DEC < -10
    """
    from skymap.survey import BlissSkymap
    smap = BlissSkymap()
    smap.draw_bliss()
    return smap.draw_hpxmap(hpxmap,**kwargs)

def draw_delve(hpxmap,**kwargs):
    """ 
    Draw BLISS footprint:
       360 > RA > 120, -60 < DEC < -10
    """
    from skymap.survey import DelveSkymap
    smap = DelveSkymap()
    smap.draw_milky_way()
    return smap.draw_hpxmap(hpxmap,**kwargs)

def draw_desgw(hpxmap,proj='car',**kwargs):
    kwargs.setdefault('extent',[180,30,-80,-50])
    return draw_footprint(hpxmap,proj,**kwargs)

def draw_survey(hpxmap,survey=None,**kwargs):
    kwargs.setdefault('cmap','viridis')
    if survey == 'des': return draw_des(hpxmap,**kwargs)
    elif survey == 'maglites': return draw_maglites(hpxmap,**kwargs)
    elif survey == 'bliss': return draw_bliss(hpxmap,**kwargs)
    elif survey == 'delve': return draw_bliss(hpxmap,**kwargs)
    else: return draw_skymap(hpxmap,**kwargs)

def draw_pixel(hpxmap,**kwargs):
    if isinstance(hpxmap,np.ma.MaskedArray):
        pix = np.where(~hpxmap.mask)
    else:
        pix = np.where((np.isfinite(hpxmap)) & (hpxmap !=hp.UNSEEN))

    vmin,vmax = np.percentile(hpxmap[pix],[0.5,99.5])
    kwargs.setdefault('vmin',vmin)
    kwargs.setdefault('vmax',vmax)
    kwargs.setdefault('rasterized',True)
    kwargs.setdefault('cmap','jet')

    nside = hp.npix2nside(len(hpxmap))
    pixrad = np.degrees(hp.max_pixrad(nside))
    ra,dec = hp.pix2ang(nside,pix,lonlat=True)

    xmin,xmax = ra.min()-pixrad,ra.max()+pixrad
    ymin,ymax = dec.min()-pixrad,dec.max()+pixrad

    delta = 0.01
    xx,yy = np.meshgrid(np.arange(xmin,xmax,delta),np.arange(ymin,ymax,delta))
    pp = hp.ang2pix(nside,xx,yy,lonlat=True)

    ax = plt.gca()
    im = ax.pcolormesh(xx[::-1],yy,hpxmap[pp],**kwargs)
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('DEC (deg)')
    ax.set_xlim(xmax,xmin)
    ax.set_ylim(ymin,ymax)
    ax.grid(ls=':',color='black',lw=0.5)
    return im


if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    opts = parser.parse_args()
