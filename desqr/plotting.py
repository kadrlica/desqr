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

from desqr import utils
from desqr.const import OBJECT_ID, UNIQUE_ID, BANDS, BADMAG, NSIDES
from desqr.utils import bfields, load_infiles, setdefaults, isstring
from desqr.utils import mkdir, calc_statistics

viridis_w = plt.cm.viridis.copy()
viridis_w.set_under('none')
gray_w = plt.cm.gray.copy()
gray_w.set_under('none')
    
def draw_peak(peak,**kwargs):
    kwargs.setdefault('ls','--')
    kwargs.setdefault('label',f'{peak:.1f}')
    ax = plt.gca()
    ax.axvline(peak,**kwargs)

def draw_peak_hist(hpxmap,**kwargs):
    """ Draw the histogram and annotate peak. """
    ax = plt.gca()

    data = utils.masked_array(hpxmap)
    vmin,vmax = utils.calc_statistics(hpxmap, q=[0.1, 99.9])['p']
    vmin = kwargs.pop('vmin', vmin)
    vmax = kwargs.pop('vmax', vmax)
    nbins = kwargs.pop('nbins',100)
    defaults = dict(bins=np.linspace(vmin,vmax,nbins),
                    histtype='step',density=True,lw=1.5,
                    peak=False,quantiles=False,color='k')
    setdefaults(kwargs,defaults)

    do_stats = kwargs.pop('stats')
    do_peak = kwargs.pop('peak')
    do_quantiles = kwargs.pop('quantiles')
    do_overflow = kwargs.pop('overflow',True)
    # Deal with bug: https://github.com/matplotlib/matplotlib/issues/6448/
    if do_overflow:
        data = np.clip(data, kwargs['bins'].min(), kwargs['bins'].max())
    else:
        data = data[(data > kwargs['bins'].min()) & (data < kwargs['bins'].max())]

    # Plot the histogram
    ax.hist(data.compressed(), **kwargs)

    # Calculate statistics on the histogram
    stats = utils.calc_statistics(data, q=[5,16,50,84,95], bins=kwargs['bins'])

    if do_peak:
        draw_peak(stats['peak'], color='k', label=f"{stats['peak']:.1f}")
        
    if do_quantiles:
        for q, p in zip(stats['q'], stats['p']):
            draw_peak(p, color='gray', ls=':', label=f'{p:.1f} ({100-q:g}%)')

    if do_stats:
        text = ""
        for s in ['peak', 'mean', 'median', 'std', 's68']:
            text += f"{s}: {stats[s]:.1f}\n"
        ax.annotate(text.rstrip(), (0.7,0.95), xycoords='axes fraction',
                    fontsize=8, ha='left', va='top')
            
    if not np.isnan(kwargs['bins']).all():
        ax.set_xlim(np.nanmin(kwargs['bins']),np.nanmax(kwargs['bins']))
        
    return stats

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

    smap,im = draw_survey(hpxmap, survey, **hpxmap_kwargs)
    smap.draw_inset_colorbar(**cbar_kwargs)
    smap.draw_milky_way()

    #ax1.axis['right'].major_ticklabels.set_visible(False)
    #ax1.axis['top'].major_ticklabels.set_visible(False)
    ax1.axis[:].set_visible(False)

    #ax1.set_axis_off()
    
    ax2 = Subplot(fig,gridspec[2])
    fig.add_subplot(ax2)
    plt.sca(ax2)
    stats = draw_peak_hist(hpxmap,**hist_kwargs)
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
    return smap,smap.draw_hpxmap(hpxmap,**kwargs)

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
    # Need to handle masked and unmasked arrays; so convert to ma
    val = np.ma.masked_invalid(np.ma.masked_values(val, hp.UNSEEN))

    kw = dict(meridians=False, parallels=False)
    if val.mask.all():
        kw['lon_0'] = kwargs.get('lon_0',0)
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
