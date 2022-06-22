#!/usr/bin/env python
"""
Characterize astrometric scatter.
"""
import glob
import os
from os.path import join
import logging
from collections import OrderedDict as odict

import yaml
import matplotlib
if os.getenv('TERM')=='screen' or not os.getenv('DISPLAY'):
    matplotlib.use('Agg')
import pylab as plt

import numpy as np
import scipy.ndimage as nd
import healpy as hp
import fitsio

import ugali.utils.projector
from ugali.utils.projector import angsep

import utils
import plotting
from const import OBJECT_ID, UNIQUE_ID, BANDS, BADMAG, NSIDES
from utils import bfields, load_infiles, get_vizier_catalog, get_local_catalog
from utils import mkdir, blank, mr_nice_guy
from match import match_query

COLUMNS = [OBJECT_ID,'RA','DEC']
TYPES = ['internal','hpxrms','gaia_dr1','gaia_dr2','gaia_edr3']

# Catalogs to be accessed from Vizier
VIZIER = odict([
        ('2MASS',dict(catalog='II/246/out',
                      columns=['2MASS','_RAJ2000','_DEJ2000','Jmag'])
                      #column_filters={})
         ),
        ('UCAC4',dict(catalog='I/322A/out',
                      columns=['UCAC4','_RAJ2000','_DEJ2000','Jmag'])
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
    kwargs.setdefault('density',True)

    ax = plt.gca()
    n,b,p = ax.hist(sep,**kwargs)
    ax.set_xlabel("Angular Separation (mas)")
    if kwargs.get('density'):
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

def draw_astrometry_pixel(hpxmap,**kwargs):
    """ hpxmap in units of mas """
    im = plotting.draw_pixel(hpxmap,**kwargs)
    ax = plt.gca()
    label=r'Median Angular Separation (mas)'
    plt.colorbar(im,label=label,
                 orientation='horizontal',fraction=0.05,shrink=0.75)

def draw_astrometry_footprint(hpxmap,survey=None,**kwargs):
    """ Draw skymap of astrometric offsets.

    Parameters
    ----------
    hpxmap : healpix map of astrometric offsets (mas)
    survey : survey map
    kwargs : passed to plotting.draw_survey

    Returns
    -------
    None
    """
    im = plotting.draw_survey(hpxmap,survey,**kwargs)[0]
    ax = plt.gca()
    label=r'Median Angular Separation (mas)'
    plt.colorbar(im,label=label,
                 orientation='horizontal',fraction=0.05,shrink=0.75)

def draw_astrometry_hist(hpxmap,**kwargs):
    """ hpxmap in units of mas """
    kwargs.setdefault('bins',np.linspace(0,200,101))
    q,p = plotting.draw_hist(hpxmap,**kwargs)
    plt.legend(loc='upper right')
    ax = plt.gca()
    ax.set_xlabel('Median Separation (mas)')

def plot_astrometry(filename,outfile=None,survey='delve'):
    """ Plot the astrometry over the footprint """
    print("Reading %s..."%filename)
    hpxmap = hp.read_map(filename,verbose=False)

    label = r'Median Angular Separation (mas)'
    cbar_kwargs = dict(label=label)
    hpxmap_kwargs = dict(xsize=1000)
    hist_kwargs = dict(bins=np.linspace(1,100,51))
    fig,axes,smap = plotting.plot_hpxmap_hist(hpxmap,survey,cbar_kwargs,hpxmap_kwargs,hist_kwargs)
    axes[1].set_xlabel(label)

    if outfile is None: 
        outfile=os.path.basename(filename).split('.')[0]+'.png'
    print("Writing %s..."%outfile)
    plt.savefig(outfile,bbox_inches='tight')

def external_astrometry(catfile,nside=128,band='r',catalog='gaia_dr2',plot=False):
    """ Calculate the median astrometric spread between DECam catalog and an external catalog.

    Parameters
    ----------
    filename : input filename (catalog file)
    nside    : output nside
    band     : band of interest ['g','r','i','z']
    plot     : pause and plot

    Returns
    -------
    pix, stat : healpix pixels and statistics (medians)
    """
    #mr_nice_guy()

    if not os.path.exists(catfile): 
        msg = "Couldn't find %s"%catfile
        raise Exception(msg)

    # Get pixel from filename (fragile)
    pixel = int(catfile.split('_')[-1].split('.')[0])
    ra,dec = hp.pix2ang(NSIDE, pixel, lonlat=True)
    radius = np.degrees(hp.max_pixrad(NSIDE))
    msg = '%s: (RA,DEC) = %.2f,%.2f'%(os.path.basename(catfile),ra,dec)
    print(msg)

    # Load catalog
    columns = [OBJECT_ID,'RA','DEC']
    spread,mag,nepochs = bfields(['WAVG_SPREAD_MODEL','MAG_PSF','NEPOCHS'],band)
    columns += [spread,mag,nepochs]

    cat = load_infiles([catfile],columns)
    # Select stars with 16 < mag < 19
    sel = (np.fabs(cat[spread])<0.002) & \
          (cat[mag]>16) & \
          (cat[mag]<19) & \
          (cat[nepochs] > 1)
    cat = cat[sel]

    if len(cat) == 0:
        msg = "WARNING: No objects passing selection in: %s"%catfile
        print(msg)
        return np.array([],dtype=int), np.array([])

    #print "Getting external catalog: %s"%catalog
    if catalog in list(VIZIER.keys()):
        ext = get_vizier_catalog(ra,dec,radius,**CATALOGS[catalog])
    else:
        ext = get_local_catalog(ra,dec,radius,catalog)

    m = match_query(cat['RA'],cat['DEC'],ext['RA'],ext['DEC'])

    # Use a fairly wide matching radius (2 arcsec)
    cut = 2.0
    sel = m[-1]*3600. < cut
    sepdeg = m[-1][sel]
    sepsec = m[-1][sel] * 3600.
    sepmas = sepsec * 1000.

    sep = sepmas

    pix = hp.ang2pix(nside,cat['RA'][sel],cat['DEC'][sel],lonlat=True)
    upix = np.unique(pix)
    try:
        peak = nd.median(sep,labels=pix,index=upix)
    except ValueError as e:
        msg = "WARNING: Failed to calculate peak for: %s"%catfile
        print(msg)
        peak = np.nan*np.ones(len(upix))
        
    if plot:
        plt.figure()
        draw_angsep(sep,bins=np.linspace(0,cut*1000.,101))
        if isinstance(plot,basestring):
            outfile = plot
            plt.savefig(outfile,bbox_inches='tight')

    #return cat,ext,m
    return upix,peak

def internal_astrometry(catfile,nside=128,band='r',plot=False):
    """Calculate internal relative astrometry.

    Parameters
    ----------
    catfile : merged catalog file
    nside   : nside for calculation
    band    : band to use ['g','r','i','z','griz']
    plot    : plot output

    Returns
    -------
    pix,stat : healpix pixels and output statistics
    """
    #mr_nice_guy()

    if not os.path.exists(catfile): 
        msg = "Couldn't find %s"%catfile
        raise Exception(msg)

    # Get pixel from filename (fragile)
    pixel = int(catfile.split('_')[-1].split('.')[0])
    ra,dec = hp.pix2ang(NSIDE, pixel, lonlat=True)
    msg = '%s: (RA,DEC) = %.2f,%.2f'%(os.path.basename(catfile),ra,dec)
    print(msg)

    # Build hpx file names (fragile)
    b = '*' if band=='griz' else band
    hpxfile = os.path.join(hpxdir, b, hpxbase.format(band=b))
    hpxfiles = glob.glob(hpxfile%pixel)
    
    if not len(hpxfiles): 
        print("WARNING: No matched hpx files: %s"%hpxfiles)
        return np.array([],dtype=int), np.array([])

    columns = [OBJECT_ID,'RA','DEC']
    # Only load one band for catalog
    b = 'r' if band=='griz' else band
    spread,mag,nepochs = bfields(['WAVG_SPREAD_MODEL','MAG_PSF','NEPOCHS'],b)
    columns += [spread,mag,nepochs]

    cat = load_infiles([catfile],columns)
    # Select stars with 17 < mag < 21
    sel = (np.fabs(cat[spread])<0.002) & \
        (cat[mag]>17) & \
        (cat[mag]<21) & \
        (cat[nepochs] > 1)
    cat = cat[sel]

    if len(cat) == 0:
        print("WARNING: No objects passing selection in: %s"%catfile)
        return np.array([],dtype=int), np.array([])

    # Load hpx files
    hpx = load_infiles(hpxfiles, [OBJECT_ID, 'RA', 'DEC'])
    hpx = hpx[np.in1d(hpx[OBJECT_ID],cat[OBJECT_ID])]

    if len(hpx) == 0:
        print("WARNING: No matched objects in: %s"%hpxfile)
        return np.array([],dtype=int), np.array([])

    # Make sure that the order matches between coadd and single epoch.
    uid,inv,cts = np.unique(hpx[OBJECT_ID],False,True,True)

    if not np.all(uid == cat[OBJECT_ID]):
        cat = cat[np.in1d(cat[OBJECT_ID],hpx[OBJECT_ID])]
    if not np.all(uid == cat[OBJECT_ID]):
        cat = cat[np.argsort(cat[OBJECT_ID])]
    assert np.all(uid == cat[OBJECT_ID])
    
    ra,dec = cat['RA'][inv],cat['DEC'][inv]

    sepdeg = angsep(ra,dec,hpx['RA'],hpx['DEC'])
    sepsec = sepdeg * 3600.
    sepmas = sepsec * 1000.
    sel = [sepsec > 1e-5] # remove same objects
    sep = sepmas[sel]

    pix = hp.ang2pix(nside,ra[sel],dec[sel],lonlat=True)
    upix = np.unique(pix)
    peak = nd.median(sep,labels=pix,index=upix)

    if plot:
        plt.figure()
        draw_angsep(sep)
        if isinstance(plot,basestring):
            outfile = plot
            plt.savefig(outfile,bbox_inches='tight')

    return upix,peak

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config',help='configuration file')
    parser.add_argument('-b','--band',default='r',choices=BANDS+['griz'])
    parser.add_argument('-n','--nside',default=128,type=int)
    parser.add_argument('-o','--outbase',default='astrom')
    parser.add_argument('-p','--pix',default=None,type=int,action='append')
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('--type',choices=TYPES,default='gaia_dr2')
    parser.add_argument('--nproc',default=4,type=int)
    args = parser.parse_args()

    if args.verbose: logging.getLogger().setLevel(logging.DEBUG)
    print("Calculating astrometric offsets ...")
    band = args.band

    config = yaml.safe_load(open(args.config))
    OBJECT_ID = config.get('objid',OBJECT_ID)
    NSIDE = config['nside']
    NPIX = hp.nside2npix(NSIDE)
    survey=config.get('survey')
    catdir = config['catdir']
    catbase = config['catbase']
    hpxdir = config['hpxdir']
    hpxbase = config['hpxbase']

    # Output nside
    nside = args.nside
    outdir = mkdir('release/astrometry')
    outbase = args.outbase

    # Find healpix pixels
    if args.pix is not None:
        pixels = args.pix
    else:
        pixels = np.arange(hp.nside2npix(NSIDE))

    if len(pixels) == 0:
        msg = "Invalid pixel: %s"%args.pix
        raise Exception(msg)
     
    filenames = [os.path.join(catdir,catbase%p) for p in pixels]
    filenames = [f for f in filenames if os.path.exists(f)]
    #filenames = filenames[1000:1250]

    if not len(filenames):
        msg = "No valid files found."
        raise Exception(msg)

    # Args must be tuple
    print("Processing %i files..."%len(filenames))

    kwargs = dict(nside=args.nside,band=band)
    if args.type in ('hpx_rms'):
        func    = internal_astrometry
        arglist = list(zip(filenames,))
        outbase += '_%s_%s'%(args.type,band)
        print('Calculating internal astrometry in %(band)s...'%kwargs)
    else:
        func    = external_astrometry
        arglist = list(zip(filenames,))
        kwargs['catalog'] = args.type
        outbase += '_%s_%s'%(args.type,band)
        print('Calculating external astrometry vs. %(catalog)s...'%kwargs)

    # Dispatch the jobs
    results = utils.multiproc(func,arglist,kwargs,processes=args.nproc)

    # Assemble the output map
    hpxmap = blank(nside)
    if None in results:
        print("WARNING: %i processes failed..."%results.count(None))

    for pix,peak in [r for r in results if r is not None]:
        hpxmap[pix] = peak
     
    hpxmap = np.ma.MaskedArray(hpxmap,np.isnan(hpxmap),fill_value=np.nan)

    outfile = join(outdir,outbase+'_n%i.fits.gz'%(nside))
    print("Writing %s..."%outfile)
    hp.write_map(outfile,hpxmap,overwrite=True)

    q = [5,50,95]
    p = np.percentile(hpxmap.compressed(),q)
    print("Global Astrometric Percentiles:")
    print('%s (%s%%)'%(p,q))

    print("Plotting %s..."%outfile)
    pngfile = outfile.replace('.fits.gz','.png')
    plot_astrometry(outfile,pngfile,survey=survey)

    plt.ion()
        
