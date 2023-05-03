#!/usr/bin/env python
import glob
import os
from os.path import join
import matplotlib
if os.getenv('TERM')=='screen' or not os.getenv('DISPLAY'):
    matplotlib.use('Agg')

import fitsio
import numpy as np
import scipy.ndimage as nd
import pylab as plt
import matplotlib.colors as colors
import healpy as hp

from mpl_toolkits.axisartist import Subplot
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from ugali.utils.shell import mkdir

try:
    from desqr.const import OBJECT_ID, UNIQUE_ID, BANDS, BADMAG, NSIDES
    from desqr.utils import bfields, load_infiles, setdefaults, isstring
except ModuleNotFoundError:
    from .const import OBJECT_ID, UNIQUE_ID, BANDS, BADMAG, NSIDES
    from .utils import bfields, load_infiles, setdefaults, isstring

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

def draw_peak_hist(hpxmap,**kwargs):
    ax = plt.gca()

    if isinstance(hpxmap,np.ma.MaskedArray):
        pix = np.where(~hpxmap.mask)
    else:
        pix = np.where((np.isfinite(hpxmap)) & (hpxmap != hp.UNSEEN))

    data = hpxmap[pix]

    vmin = kwargs.pop('vmin',np.percentile(data,q=0.1))
    vmax = kwargs.pop('vmax',np.percentile(data,q=99.9))
    nbins = kwargs.pop('nbins',100)
    defaults = dict(bins=np.linspace(vmin,vmax,nbins),
                    histtype='step',density=True,lw=1.5,
                    peak=False,quantiles=False,color='k')
    setdefaults(kwargs,defaults)

    do_peak = kwargs.pop('peak')
    do_quantiles = kwargs.pop('quantiles')
    do_overflow = kwargs.pop('overflow',False)
    # Deal with bug: https://github.com/matplotlib/matplotlib/issues/6448/
    if do_overflow:
        data = np.clip(data,kwargs['bins'].min(),kwargs['bins'].max())
    else:
        data = data[(data > kwargs['bins'].min())&(data < kwargs['bins'].max())]

    n,b,p = ax.hist(data,**kwargs)

    ret = dict()
    peak = ((b[1:]+b[:-1])/2.)[np.argmax(n)]
    ret['peak'] = peak
    if do_peak:
        draw_peak(peak,color='k',label='%.1f'%(peak))

    ret['mean'] = np.mean(data)
    ret['std']  = np.std(data)

    quantiles = [5,16,50,84,95]
    percentiles = np.percentile(data,quantiles)
    ret['quantiles']   = quantiles
    ret['percentiles'] = percentiles
    for p,q in zip(percentiles,quantiles):
        ret['q%02d'%q] = p

    if do_quantiles:
        for q,p in zip(quantiles,percentiles):
            draw_peak(p,color='r',label='%.1f (%g%%)'%(p,100-q))

    ax.set_xlim(kwargs['bins'].min(),kwargs['bins'].max())
    return ret

def create_hpxmap_hist_figure():
    #fig = plt.figure(figsize=(10.5,3.8))
    fig = plt.figure(figsize=(12.0,3.8))
    gridspec=plt.GridSpec(1, 3)
    gridspec.update(left=0.07,right=0.91,bottom=0.15,top=0.95,wspace=0.08)
    return fig, gridspec

def plot_hpxmap_hist(hpxmap,survey=None,
                     cbar_kwargs=dict(),hpxmap_kwargs=dict(),hist_kwargs=dict()):
    """Plot two-panel figure with skymap and histogram.
    
    Parameters
    ----------
    hpxmap : healpix map
    survey : survey configuration
    cbar_kwargs : kwargs passed to draw_inset_colorbar
    hpxmap_kwargs : kwargs passed to draw_hpxmap
    hist_kwargs : kwargs passed to draw_peak_hist
    
    Returns
    -------
    fig,[ax1,ax2],smap
    """
    from skymap.survey import SurveyMcBryde

    hist_defaults = dict(peak=True)
    setdefaults(hist_kwargs,hist_defaults)

    cbar_defaults = dict()
    if survey != 'des':
        cbar_defaults['loc'] = 'upper center'
    setdefaults(cbar_kwargs,cbar_defaults)

    if isstring(hpxmap):
        hpxmap = hp.read_map(f)

    fig,gridspec = create_hpxmap_hist_figure()
    ax1 = Subplot(fig,gridspec[0:2])
    fig.add_subplot(ax1)
    plt.sca(ax1)

    smap,im = draw_survey(hpxmap,survey,**hpxmap_kwargs)
    smap.draw_inset_colorbar(**cbar_kwargs)
    smap.draw_milky_way()

    #ax1.axis['right'].major_ticklabels.set_visible(False)
    #ax1.axis['top'].major_ticklabels.set_visible(False)
    ax1.axis[:].set_visible(False)

    #ax1.set_axis_off()
    
    ax2 = Subplot(fig,gridspec[2])
    fig.add_subplot(ax2)
    plt.sca(ax2)
    ret = draw_peak_hist(hpxmap,**hist_kwargs)
    ax2.yaxis.set_major_locator(MaxNLocator(6,prune='both'))
    ax2.xaxis.set_major_locator(MaxNLocator(5))
    ax2.axis['left'].major_ticklabels.set_visible(False)
    ax2.axis['right'].major_ticklabels.set_visible(True)
    ax2.axis['right'].label.set_visible(True)
    ax2.axis['right'].label.set_text(r'Normalized Area')
    ax2.axis['bottom'].label.set_visible(True)

    #plt.subplots_adjust(bottom=0.15,top=0.95,wspace=0.1,right=0.90)

    return fig,[ax1,ax2],smap


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
    """ Draw a skymap hpxmap """
    from skymap.survey import SurveySkymap
    smap = SurveySkymap()
    return smap,smap.draw_hpxmap(hpxmap,**kwargs)

def draw_des(hpxmap,**kwargs):
    """ Draw DES footprint:
       110 > RA > -70, -70 < DEC < 10
    """
    from skymap.survey import DESSkymap
    smap = DESSkymap()
    smap.draw_des()
    return smap,smap.draw_hpxmap(hpxmap,**kwargs)

def draw_maglites(hpxmap,**kwargs):
    """ Draw MagLiteS footprint:
       280 > RA > 80, -90 < DEC < -50
    """
    from skymap.survey import MaglitesSkymap
    smap = MaglitesSkymap()
    smap.draw_maglites()
    return smap,smap.draw_hpxmap(**kwargs)

def draw_bliss(hpxmap,**kwargs):
    """ Draw BLISS footprint:
       360 > RA > 120, -60 < DEC < -10
    """
    from skymap.survey import BlissSkymap
    smap = BlissSkymap()
    smap.draw_bliss()
    return smap,smap.draw_hpxmap(hpxmap,**kwargs)

def draw_delve(hpxmap,**kwargs):
    """ Draw DELVE footprint
    """
    from skymap.survey import SurveyMcBryde
    nside = hp.npix2nside(len(hpxmap))
    vec = hp.ang2vec(180, -30, lonlat=True)
    pix = hp.query_disc(nside, vec, np.radians(1.0))
    val = hpxmap[pix]
    kw = dict(meridians=False, parallels=False)

    if np.isnan(val).all() or (val == hp.UNSEEN).all():
        kw['lon_0'] = kwargs.get('lon_0',0)
        smap = SurveyMcBryde(**kw)
    else:
        kw['lon_0'] = kwargs.get('lon_0',180)
        smap = SurveyMcBryde(**kw)

    smap.draw_meridians(fontsize=10)
    smap.draw_parallels(fontsize=10)

    return smap,smap.draw_hpxmap(hpxmap,**kwargs)

def draw_desgw(hpxmap,proj='car',**kwargs):
    kwargs.setdefault('extent',[180,30,-80,-50])
    return draw_footprint(hpxmap,proj,**kwargs)

def draw_survey(hpxmap,survey=None,**kwargs):
    """ Draw a survey skymap """
    kwargs.setdefault('cmap','viridis')
    if survey == 'des': return draw_des(hpxmap,**kwargs)
    elif survey == 'maglites': return draw_maglites(hpxmap,**kwargs)
    elif survey == 'bliss': return draw_bliss(hpxmap,**kwargs)
    elif survey == 'delve': return draw_delve(hpxmap,**kwargs)
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
