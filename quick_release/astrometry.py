#!/usr/bin/env python
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
import matplotlib.colors as colors
import healpy

from ugali.utils.healpix import ang2pix, pix2ang
import ugali.utils.projector
from ugali.utils.projector import angsep
from ugali.utils.shell import mkdir

import utils
from const import OBJECT_ID, UNIQUE_ID, BANDS, BADMAG, NSIDES
from utils import bfields, load_infiles, get_vizier_catalog
from match import match_query

import footprint
import plotting
from footprint import blank
from catalog import good_objects

COLUMNS = [OBJECT_ID,'RA','DEC']

CATALOGS = odict([
        ('2MASS',dict(catalog='II/246/out',
                      columns=['2MASS','_RAJ2000','_DEJ2000','Jmag'])
                      #column_filters={})
         ),
        ('UCAC4',dict(catalog='I/322A/out',
                      columns=['UCAC4','_RAJ2000','_DEJ2000','Jmag'])
                      #,'pmRA','e_pmRA','pmDE','e_pmDE'],
                      #column_filters={})
         ),
        ('GAIA',dict(catalog='I/337/gaia',
                     columns=['Source','_RAJ2000','_DEJ2000','<Gmag>'])
                     #,'pmRA','e_pmRA','pmDE','e_pmDE'],
                     #column_filters={})
         ),

        ])

def hist_peak(bins,num):
    centers = (bins[:-1]+bins[1:])/2.
    return centers[np.argmax(num)]

def draw_angsep(sep,**kwargs):
    """ Angular separation in mas
    """
    kwargs.setdefault('bins',np.linspace(0,250,101))
    kwargs.setdefault('histtype','step')
    kwargs.setdefault('lw',1.5)
    kwargs.setdefault('normed',True)

    ax = plt.gca()
    n,b,p = ax.hist(sep,**kwargs)
    ax.set_xlabel("Angular Separation (mas)")
    if kwargs.get('normed'):
        ax.set_ylabel("Normalized Counts")
    else:
        ax.set_ylabel("Counts")

    peak = hist_peak(b,n)
    plotting.draw_peak(peak,color='b',label='Peak = %.0f mas'%peak)
    median = np.median(sep)
    plotting.draw_peak(median,color='r',label='Median = %.0f mas'%median)
    rms = np.sqrt(np.sum(sep**2)/len(sep))
    plotting.draw_peak(rms,color='g',label='RMS = %.0f mas'%rms)

    plt.legend()

    return n,b

def draw_astrometry_pixel(skymap,**kwargs):
    """ skymap in units of mas """
    im = plotting.draw_pixel(skymap,**kwargs)
    ax = plt.gca()
    plt.colorbar(im,label=r'Median Angular Separation (mas)')

def draw_astrometry_footprint(skymap,**kwargs):
    """ skymap in units of mas """
    im = plotting.draw_des(skymap,**kwargs)
    ax = plt.gca()
    plt.colorbar(im,label=r'Median Angular Separation (mas)')

def draw_astrometry_hist(skymap,**kwargs):
    """ skymap in units of mas """
    kwargs.setdefault('bins',np.linspace(0,200,101))
    q,p = plotting.draw_hist(skymap,**kwargs)
    plt.legend(loc='upper right')
    ax = plt.gca()
    ax.set_xlabel('Median Separation (mas)')

def external_astrometry(catfile,catalog='2MASS',nside=64,band='r',plot=False):
    nice = os.nice(0)
    os.nice(10-nice)

    if not os.path.exists(catfile): 
        msg = "Couldn't find %s"%catfile
        raise Exception(msg)

    columns = [OBJECT_ID,'RA','DEC']
    spread,mag,nepochs = bfields(['WAVG_SPREAD_MODEL','MAG_PSF','NEPOCHS'],band)
    columns += [spread,mag,nepochs]

    # Hack to get pixel location
    hpx = int(catfile.split('_')[-1].split('.')[0])
    #hpx = ang2pix(NSIDE, cat['RA'], cat['DEC'])
    ra,dec = pix2ang(NSIDE, hpx)
    radius = np.degrees(healpy.max_pixrad(NSIDE))
    print os.path.basename(catfile), 
    print '(RA,DEC,RAD) = %.2f,%.2f,%.2f'%(ra,dec,radius)

    #print "Getting coadd catalog: DES"
    cat = load_infiles([catfile],columns)
    # Select stars with 17 < mag < 21
    sel = (np.fabs(cat[spread])<0.002) & \
        (cat[mag]>15) & \
        (cat[mag]<21) & \
        (cat[nepochs] > 1)
    cat = cat[sel]

    #print "Getting external catalog: %s"%catalog
    ext = get_vizier_catalog(ra,dec,radius,**CATALOGS[catalog])

    m = match_query(cat['RA'],cat['DEC'],ext['_RAJ2000'],ext['_DEJ2000'])

    # Use a fairly wide matching radius (2 arcsec)
    cut = 2.0
    sel = m[-1]*3600. < cut
    sepdeg = m[-1][sel]
    sepsec = m[-1][sel] * 3600.
    sepmas = sepsec * 1000.

    sep = sepmas

    pix = ang2pix(nside,cat['RA'][sel],cat['DEC'][sel])
    upix = np.unique(pix)
    peak = nd.median(sep,labels=pix,index=upix)

    if plot:
        plt.figure()
        draw_angsep(sep,bins=np.linspace(0,cut*1000.,101))
        if isinstance(plot,basestring):
            outfile = plot
            plt.savefig(outfile,bbox_inches='tight')

    #return cat,ext,m
    return upix,peak

def internal_astrometry(catfile,hpxfile,nside=128,band='r',plot=False):
    nice = os.nice(0)
    os.nice(10-nice)

    #print catfile,hpxfile,nside

    #catfile = glob.glob('cat/*_%05d.fits'%pix)[0]
    if not os.path.exists(catfile): 
        msg = "Couldn't find %s"%catfile
        raise Exception(msg)

    columns = [OBJECT_ID,'RA','DEC']

    spread,mag,nepochs = bfields(['WAVG_SPREAD_MODEL','MAG_PSF','NEPOCHS'],band)
    columns += [spread,mag,nepochs]

    cat = load_infiles([catfile],columns)
    # Select stars with 17 < mag < 21
    sel = (np.fabs(cat[spread])<0.002) & \
        (cat[mag]>17) & \
        (cat[mag]<21) & \
        (cat[nepochs] > 1)
    cat = cat[sel]

    if len(cat) == 0:
        print "WARNING: No objects passing selection in: %s"%catfile
        return np.array([],dtype=int), np.array([])

    #hpxfiles = glob.glob('hpx/%s/*_%05d.fits'%(band,pix))
    hpx = load_infiles(hpxfile, [OBJECT_ID, 'RA', 'DEC'])
    hpx = hpx[np.in1d(hpx[OBJECT_ID],cat[OBJECT_ID])]

    if len(hpx) == 0:
        print "WARNING: No matched objects in: %s"%hpxfile
        return np.array([],dtype=int), np.array([])
        
    #keyfile = 'key/key_hpx_%05d.fits'%pix
    #key = load_infiles([keyfile],[OBJECT_ID,'FILENAME','OBJECT_NUMBER'])
    #key = key[np.in1d(key[OBJECT_ID],cat[OBJECT_ID])]
    # 
    #key_id = np.char.add(key['FILENAME'],key['OBJECT_NUMBER'].astype(str))
    #hpx_id = np.char.add(hpx['FILENAME'],hpx['OBJECT_NUMBER'].astype(str))
    # 
    #hpx = hpx[np.in1d(hpx_id,key_id)]

    uid,inv,cts = np.unique(hpx[OBJECT_ID],False,True,True)

    # Make sure that the order matches between coadd and se.
    if not np.all(uid == cat[OBJECT_ID]):
        cat = cat[np.in1d(cat[OBJECT_ID],hpx[OBJECT_ID])]
    if not np.all(uid == cat[OBJECT_ID]):
        cat = cat[np.argsort(cat[OBJECT_ID])]
    assert np.all(uid == cat[OBJECT_ID])
    
    ra,dec = cat['RA'][inv],cat['DEC'][inv]

    sepdeg = angsep(ra,dec,hpx['RA'],hpx['DEC'])
    sepsec = sepdeg * 3600.
    sepmas = sepsec * 1000.
    sel = [sepsec > 1e-5]
    sep = sepmas[sel]

    pix = ang2pix(nside,ra[sel],dec[sel])
    upix = np.unique(pix)
    peak = nd.median(sep,labels=pix,index=upix)

    if plot:
        plt.figure()
        draw_angsep(sep)
        if isinstance(plot,basestring):
            outfile = plot
            plt.savefig(outfile,bbox_inches='tight')

    return upix,peak

def distance(args,plot=False):
    nice = os.nice(0)
    os.nice(10-nice)

    pix,nside = args
    catfile = glob.glob('cat/*_%05d.fits'%pix)[0]
    if not os.path.exists(catfile): 
        msg = "Couldn't find %s"%catfile
        raise Exception(msg)

    print catfile

    columns = [OBJECT_ID,'RA','DEC']

    spread,mag,nepochs = bfields(['WAVG_SPREAD_MODEL','MAG_PSF','NEPOCHS'],band)
    columns += [spread,mag,nepochs]

    cat = load_infiles([catfile],columns)
    sel = (np.fabs(cat[spread])<0.002)&(cat[mag]>16)&(cat[mag]<22)&(cat[nepochs] > 1)
    cat = cat[sel]

    if len(cat) == 0:
        print "WARNING: No catalog objects passing selection"
        return np.array([],dtype=int), np.array([])

    ra,dec = cat['RA'],cat['DEC']

    m = ugali.utils.projector.match(ra,dec,ra,dec,nnearest=2)
    sep = m[-1]

    hpx = ang2pix(nside,ra,dec)
    peak = nd.median(sep,labels=hpx,index=np.unique(hpx))

    if plot:
        plt.figure()
        draw_angsep(sep)
        if isinstance(plot,basestring):
            outfile = plot
            plt.savefig(outfile,bbox_inches='tight')

    return hpx,peak

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config',help='Configuration file.')
    parser.add_argument('-o','--outbase',default='y3q2_astrom')
    parser.add_argument('-t','--title',default='Astrometric Precision')
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('-n','--nside',default=128,type=int)
    parser.add_argument('-p','--pix',default=None,type=int,action='append')
    parser.add_argument('-b','--band',action='append',choices=BANDS)
    parser.add_argument('--type',choices=['internal','external'],default='external')
    parser.add_argument('--catalog',choices=CATALOGS.keys(),default='UCAC4')
    opts = parser.parse_args()

    if not opts.band: opts.band = BANDS

    config = yaml.load(open(opts.config))
    OBJECT_ID = config.get('objid',OBJECT_ID)
    NSIDE = config['nside']

    nside = opts.nside
    npix = healpy.nside2npix(config['nside'])
    catdir = config['catdir']

    outdir = mkdir('release/astrometry')
    outbase = opts.outbase
    outbase += '_%s'%(opts.type[:3])
    if opts.type == 'external':
        outbase += '_%s'%(opts.catalog.lower())

    skymaps = odict()

    for band in opts.band:
        print 'Calculating astrometry for %s-band...'%band
        hpxdir = join(config['hpxdir'],band)
         
        if opts.pix is not None:
            pixels = sorted([p for p in opts.pix if len(glob.glob(catdir+'/*%05d.fits'%p))])
        else:
            pixels = sorted([p for p in range(npix) if (len(glob.glob(catdir+'/*%05d.fits'%p)) and len(glob.glob(hpxdir+'/*%05d.fits'%p))) ])
         
        if len(pixels) == 0:
            msg = "Invalid pixel: %s"%opts.pix
            raise Exception(msg)
         
        catfiles = [glob.glob(catdir+'/*%05d.fits'%p)[0] for p in pixels]
        hpxfiles = [glob.glob(hpxdir+'/*%05d.fits'%p)[0] for p in pixels]
         
        #catfiles = ['cat/y1a1_gold_1_0_2_01376.fits']
        #hpxfiles = ['hpx/i/hpx_i_01376.fits']
         
        # Args must be tuple
        print "Processing %i files..."%len(catfiles)
         
        if opts.type == 'internal':
            args = zip(catfiles,hpxfiles)
            kwargs = dict(nside=opts.nside,band=band)
            astrometry = internal_astrometry
        elif opts.type == 'external':
            args = zip(catfiles,)
            kwargs = dict(nside=opts.nside,band=band,catalog=opts.catalog.upper())
            astrometry = external_astrometry
        else:
            raise Exception()
         
        results = utils.multiproc(astrometry,args,kwargs)
        #results = [astrometry(*args[-1],**kwargs)]
         
        #args = [(pix,nside) for pix in pixels][:10]
        #p = Pool(maxtasksperchild=1,processes=20)
        #out = p.map(astrometry,args)
         
        skymap = blank(nside)
         
        if None in results:
            print "WARNING: %i processes failed..."%results.count(None)
        for pix,peak in [r for r in results if r is not None]:
            skymap[pix] = peak
         
        skymap = np.ma.MaskedArray(skymap,np.isnan(skymap),fill_value=np.nan)
         
        skymaps[band] = skymap

        q = [5,50,95]
        p = np.percentile(skymap.compressed(),q)
        print "Global Median Angular Separation:"
        print '%s (%s%%)'%(p,q)
         
    print '|_. Band |_. Median Separation Map |_. Median Separation Distribution |_. Astrometric Precision |'
    template = '|_. %(band)s |{{thumbnail(%(map)s, size=300)}}|{{thumbnail(%(hist)s, size=300)}}|_. %(value)s mas |'
            
    for band,skymap in skymaps.items():
        out = dict(band=band)
         
        if opts.pix is not None:
            pix = opts.pix[0]
         
            plt.figure()
            draw_astrometry_pixel(skymap)
            plt.title(opts.title+" (%i; %s-band)"%(pix,band))
            outfile = join(outdir,outbase+'_%s_hpx%05d_n%i_car.png'%(band,pix,nside))
            plt.savefig(outfile,bbox_inches='tight')
            out['map'] = os.path.basename(outfile)
         
            kwargs['plot']=True
            utils.multiproc(astrometry,args,kwargs)
            plt.title(opts.title+" (%i; %s-band)"%(pix,band))
            outfile = join(outdir,outbase+'_%s_hpx%05d_n%i_hist.png'%(band,pix,nside))
            plt.savefig(outfile,bbox_inches='tight')
            out['hist'] = os.path.basename(outfile)
         
            outfile = join(outdir,outbase+'_%s_hpx%05d_n%i.fits'%(band,pix,nside))
            healpy.write_map(outfile,skymap)
         
        else:
            plt.figure()
            draw_astrometry_footprint(skymap)
            plt.title(opts.title+" (%s-band)"%band)
            outfile = join(outdir,outbase+'_%s_n%i_car.png'%(band,nside))
            plt.savefig(outfile,bbox_inches='tight')
            out['map'] = os.path.basename(outfile)
         
            plt.figure()
            draw_astrometry_hist(skymap,bins=np.linspace(0,500,101))
            plt.title(opts.title+" (%s-band)"%band)
            outfile = join(outdir,outbase+'_%s_n%i_hist.png'%(band,nside))
            plt.savefig(outfile,bbox_inches='tight')
            out['hist'] = os.path.basename(outfile)
         
            outfile = join(outdir,outbase+'_%s_n%i.fits'%(band,nside))
            healpy.write_map(outfile,skymap)


        q = [5,50,95]
        p = np.percentile(skymap.compressed(),q)
        out['value'] = '/'.join(['%.0f'%_p for _p in p])
        out['percentile'] = '/'.join(['%.0f'%_q for _q in q])
        outstr = template%out
        print outstr

    plt.ion()
