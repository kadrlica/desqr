#!/usr/bin/env python
import glob
import os
import matplotlib
try:             os.environ['DISPLAY']
except KeyError: matplotlib.use('Agg')

import fitsio
import numpy as np
import scipy.ndimage as nd
import pylab as plt
import matplotlib.colors as colors
import healpy

from ugali.utils.healpix import ang2pix
from ugali.utils.projector import angsep
from ugali.utils.shell import mkdir
import catalog
from utils import bfields
from const import BANDS, OBJECT_ID

def plot_peak(peak,text,**kwargs):
    ax = plt.gca()
    ax.axvline(peak,**kwargs)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x = peak +  0.05*(xlim[1]-xlim[0])
    y = ylim[1] - 0.075*(ylim[1] - ylim[0])
    ax.annotate(text,(x,y),xycoords='data',va='center')

def plot_skymap(infiles,nside=1024):
    counts = np.zeros(healpy.nside2npix(nside),dtype=int)
    skymap = healpy.UNSEEN * np.ones(healpy.nside2npix(nside))

    for i,f in enumerate(infiles):
        if i%10 == 0: print "(%i/%i)"%(i,len(infiles))
        d = fitsio.read(f,columns=['RA','DEC','MAG_PSF_G'])
        sel = (d['MAG_PSF_G'] < 23)
        d = d[sel]
        pix,cts = np.unique(ang2pix(nside,d['RA'],d['DEC']),return_counts=True)
        counts[pix] += cts

    idx = np.where(counts>0)
    pixarea = healpy.nside2pixarea(nside,degrees=True)
    skymap[idx] = np.log10(counts[idx]/pixarea)
    im = healpy.mollview(skymap,title='Object Density (g < 23)',return_projected_map = True)
    healpy.graticule()
    
    return counts,im

def check_nan(data):
    print "Checking for NaNs..."
    for name in data.dtype.names:
        nans = np.isnan(coadd[name]).sum()
        if nans: print name, nans

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-p','--pix',default=2585,type=int)
    parser.add_argument('-s','--skymap',action='store_true')
    opts = parser.parse_args()

    pltdir = mkdir('plots/')
    
    if opts.skymap:
        print "Creating skymap..."
        infiles = glob.glob('cat/*/cat_hpx*.fits')
        cts,im = plot_skymap(infiles)
        plt.savefig(pltdir+'skymap.png',bbox_inches='tight')

    infiles = glob.glob('hpx/*/hpx*%i.fits'%opts.pix)
    print "Loading SE files: %s"%infiles
    data = catalog.load_infiles(infiles,catalog.INPUT_COLS+['EXPNUM'])
    good = catalog.good_objects(data)
    
    infile = glob.glob('cat/cat_hpx_*%i.fits'%opts.pix)[0]
    print "Loading COADD file: %s"%infile
    coadd = fitsio.read(infile)

    check_nan(coadd)

    se = good[np.in1d(good[OBJECT_ID],coadd[OBJECT_ID])]

    uid,inv,cts = np.unique(se[OBJECT_ID],False,True,True)
    
    if not (uid == coadd[OBJECT_ID]).all():
        raise Exception("Object IDs do not match")

    kwargs = dict(bins=100,histtype='step',lw=1.5)

    # Spatial distribution
    maglim = 30

    plt.figure()
    sel = (se['MAG_PSF'] < maglim) & (se['BAND'] == 'g')
    plt.hist2d(se['RA'][sel],se['DEC'][sel],bins=250,norm=colors.LogNorm())
    plt.colorbar(label='log(Counts)')
    plt.title('HPX %05d (g < 30)'%opts.pix)
    plt.xlabel('RA (deg)')
    plt.ylabel('DEC (deg)')
    plt.savefig(pltdir+'spatial_se_%05d.png'%opts.pix)
     
    plt.figure()
    sel = coadd['MAG_PSF_G'] < maglim
    plt.hist2d(coadd['RA'][sel],coadd['DEC'][sel],bins=250,norm=colors.LogNorm())
    plt.colorbar(label='log(Counts)')
    plt.title('HPX %05d (g < %s)'%(opts.pix,maglim))
    plt.xlabel('RA (deg)')
    plt.ylabel('DEC (deg)')
    plt.savefig(pltdir+'spatial_coadd_%05d.png'%opts.pix)
     
    # Matching and separation
    sep = angsep(coadd['RA'][inv],coadd['DEC'][inv],se['RA'],se['DEC'])
     
    plt.figure()
    n,b,p = plt.hist(sep*3600.,**kwargs)
    c = (b[:-1]+b[1:])/2.
    peak = c[np.argmax(n)]; text = '%.1f (mas)'%(peak*1e3)
    plot_peak(peak,text,lw=2,ls='--',c='r')
    plt.xlabel("Separation (arcsec)")
    plt.ylabel("Counts")
    plt.savefig(pltdir+'angsep_%05d.png'%opts.pix,bbox_inches='tight')
     
    plt.figure()
    plt.hist(sep*3600.,log=True,**kwargs)
    c = (b[:-1]+b[1:])/2.
    peak = c[np.argmax(n)]; text = '%.1f (mas)'%(peak*1e3)
    plot_peak(peak,text,lw=2,ls='--',c='r')
    plt.xlabel("Separation (arcsec)")
    plt.ylabel("Counts")
    plt.savefig(pltdir+'log_angsep_%05d.png'%opts.pix,bbox_inches='tight')
     
    # Spread Model
    kwargs.update(dict(bins=np.linspace(-0.01,0.04,100)))
     
    plt.figure()
    for x in bfields(['SPREAD_MODEL','WAVG_SPREAD_MODEL'],['R','I']):
        plt.hist(coadd[x],label=x,**kwargs)
     
    plt.xlabel('SPREAD')
    plt.ylabel('Counts')
    plt.legend()
    plt.savefig(pltdir+'spread_%05d.png'%opts.pix,bbox_inches='tight')
     
    plt.figure()
    plt.hist(coadd['SPREAD_MODEL_R']+3*coadd['SPREADERR_MODEL_R'],
             label='SPREAD_R+3*ERR',**kwargs)
    plt.hist(coadd['WAVG_SPREAD_MODEL_R']+3*coadd['WAVG_SPREADRMS_MODEL_R'],
             label='WAVG_SPREAD_R+3*RMS',**kwargs)
    plt.hist(coadd['WAVG_SPREAD_MODEL_R']+3*coadd['WAVG_SPREADERR_MODEL_R'],
             label='WAVG_SPREAD_R+3*ERR',**kwargs)
    plt.xlabel('SPREAD_R + 3*ERR')
    plt.ylabel('Counts')
    plt.legend()
    plt.savefig(pltdir+'spread_3err_r_%05d.png'%opts.pix,bbox_inches='tight')
     
    for mmin,mmax in [[16,20],[20,21],[21,22],[22,23],[23,26]]:
        sel = (coadd['WAVG_MAG_PSF_R'] > mmin) & (coadd['WAVG_MAG_PSF_R'] <= mmax)
        d = coadd[sel]
        plt.figure()
        plt.plot(np.nan,np.nan,'ow',label='%g < WAVG_MAG_PSF < %g'%(mmin,mmax))
        plt.hist(d['SPREAD_MODEL_R']+3*d['SPREADERR_MODEL_R'],
                 label='SPREAD_R+3*ERR',**kwargs)
        plt.hist(d['WAVG_SPREAD_MODEL_R']+3*d['WAVG_SPREADRMS_MODEL_R'],
                 label='WAVG_SPREAD_R+3*RMS',**kwargs)
        plt.hist(d['WAVG_SPREAD_MODEL_R']+3*d['WAVG_SPREADERR_MODEL_R'],
                 label='WAVG_SPREAD_R+3*ERR',**kwargs)
        plt.xlabel('SPREAD_R + 3*ERR')
        plt.ylabel('Counts')
        plt.legend()
        plt.savefig(pltdir+'spread_3err_r_%g_%g_%05d.png'%(mmin,mmax,opts.pix),
                    bbox_inches='tight')
     
    kwargs.update(dict(bins=np.linspace(-0.001,0.01,100)))
    plt.figure()
    for x in bfields(['SPREADERR_MODEL','WAVG_SPREADERR_MODEL','WAVG_SPREADRMS_MODEL'],['R','I']):
        plt.hist(coadd[x],label=x,**kwargs)
     
    plt.xlabel('SPREADERR')
    plt.ylabel('Counts')
    plt.legend()
    plt.savefig(pltdir+'spreaderr_%05d.png'%opts.pix,bbox_inches='tight')
     
     
    # exposure-to-exposure variation in spread_model
    mmin,mmax = 16,23
    plt.figure()
    sel = (se['BAND'] == 'r')&(se['MAG_PSF'] > mmin)&(se['MAG_PSF'] <= mmax)
    s = se[sel]
    expnum,expinv,expcts = np.unique(s['EXPNUM'],False,True,True)
    argsort = np.argsort(expcts)
    kwargs.update(bins=np.linspace(-0.01,0.04,100),normed=True,lw=2.5)
    plt.hist(s['SPREAD_MODEL'],label='ALL',zorder=10,**kwargs)
    kwargs.update(lw=1.5)
    for e in expnum[argsort][-5:]:
        plt.hist(s['SPREAD_MODEL'][s['EXPNUM']==e],label='EXPNUM=%i'%e,**kwargs)
    plt.legend()
    plt.title('Spread Model per Exposure (%g < r < %s)'%(mmin,mmax))
    plt.ylabel('Counts')
    plt.xlabel('SPREAD_MODEL_R')
    plt.savefig(pltdir+'spread_model_exp_%g_%g_%05d.png'%(mmin,mmax,opts.pix),bbox_inches='tight')
     
    VARS = [['SPREAD_MODEL_R','SPREADERR_MODEL_R'],
            ['WAVG_SPREAD_MODEL_R','SPREADERR_MODEL_R'],
            ['WAVG_SPREAD_MODEL_R','WAVG_SPREADERR_MODEL_R'],
            ['WAVG_SPREAD_MODEL_R','WAVG_SPREADRMS_MODEL_R'],
            ]
    bins = [np.linspace(-0.01,0.01,100),np.linspace(0,0.005,100)]
    for spread,err in VARS:
        plt.figure()
        plt.hist2d(coadd[spread],coadd[err],bins=bins,norm=colors.LogNorm())
        plt.xlabel(spread)
        plt.ylabel(err)
     
    # Magnitude
    kwargs.update(dict(bins=np.linspace(16,26,100)))
    plt.figure()
    for i,v in enumerate(['WAVG_MAG_PSF','MAG_PSF']):
        for c,b in zip(['g','r','gold','indigo','gray'],BANDS):
            var = v+'_%s'%b.upper()
            alpha = 0.5 if i > 0 else 1.0
            plt.hist(coadd[var],alpha=alpha,label=var,color=c,**kwargs)
        plt.legend(loc='upper left',fontsize=10)
     
    plt.savefig(pltdir+'mag_%05d.png'%opts.pix,bbox_inches='tight')
     
    for c,b in zip(['g','r','gold','indigo','gray'],BANDS):
        plt.figure()
        for i,v in enumerate(['WAVG_MAG_PSF','MAG_PSF']):
            var = v+'_%s'%b.upper()
            err = var.replace('MAG_','MAGERR_')
            bins = kwargs['bins']
            labels = np.digitize(coadd[var],bins)
            index = np.unique(labels)
            medians = nd.median(coadd[err],labels=labels,index=index)
            alpha = 0.5 if i > 0 else 1.0
            mag = bins[index[:-1]]; magerr = medians[:-1]
            plt.plot(mag,magerr,'-o',c=c,alpha=alpha,label=var)
            argmax = np.argmax(magerr > 0.1)
            maglim = (mag[argmax] + mag[argmax-1])/2.
            if i==1: plot_peak(maglim,'%.1f'%maglim,color=c,lw=1.5,ls='--')
            plt.axhline(0.1,color='gray',lw=1.5,ls='--')
        plt.xlabel('MAG_PSF')
        plt.xlabel('MAGERR_PSF')
        plt.legend(loc='upper left',fontsize=10)
        plt.savefig(pltdir+'magerr_%s_%05d.png'%(b,opts.pix),bbox_inches='tight')
