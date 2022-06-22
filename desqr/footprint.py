#!/usr/bin/env python
"""
Count number of objects passing selections over the footprint.
"""
import os
from os.path import join
import glob
from collections import OrderedDict as odict
from multiprocessing import Pool
import functools

import matplotlib
if os.getenv('TERM')=='screen' or not os.getenv('DISPLAY'):
    matplotlib.use('Agg')
from matplotlib.colors import LogNorm
import pylab as plt
import numpy as np
import scipy.ndimage as nd

import yaml
import healpy as hp
import fitsio

from const import OBJECT_ID, UNIQUE_ID, BANDS, BADMAG, NSIDES
from utils import bfield, bfields, load_infiles, mkdir
from utils import empty, blank, unseen
import utils
import plotting

MAGVAR = 'WAVG_MAG_PSF'
EXTVAR = 'EXTENDED_CLASS'
REDVAR = 'EXTINCTION'

# band-dependent magnitude limit
MAGLIM = odict([
    ['g',23.0],
    ['r',23.0],
    ['i',22.5], # was 22.0
    ['z',22.0],
    ['Y',21.0],
])

# star-galaxy classification
EXTCLASS = odict([
    ['star', [0,1]],
    ['gal' , [2,3]],
    [True  , range(4)],
    [None  , None],
])

# color cuts (if color not in this list, selects all)
COLOR = odict([
    # [color, [[cmin, cmax], [mmin, mmax]] ]
    ['gr', [[0.1,0.3], [17,20]] ],
    ['gi', [[0.2,0.4], [17,20]] ],
])

# Coverage Footprint
def create_columns(bands=BANDS,extclass=False,deredden=False):
    """ Create the list of column names to load. """
    COLUMNS  = ['RA','DEC']
    COLUMNS += bfields(MAGVAR,bands)
    if extclass: 
        COLUMNS += bfields(EXTVAR,bands)
    if deredden:
        COLUMNS += bfields(REDVAR,bands)
    return COLUMNS
    
def select_mag(data, var, mag_range=None):
    """ Select objects based on maglim cut.

    Parameters
    ----------
    data      : input catalog
    var       : magnitude  variable name
    mag_range : allowed magnitude range [min, max]

    Returns
    -------
    selection
    """
    # NOOP; select everything
    if mag_range is None:
        return np.ones(len(data), dtype=bool)
    elif np.isscalar(mag_range):
        return data[var] < mag_range
    elif len(mag_range) == 2:
        return (data[var] > mag_range[0]) & (data[var] < mag_range[1])
    else:
        raise Exception()

def select_color(data, var1, var2, color_range=None):
    """ Select objects based on a color cut. 

    color = data[var1] - data[var2]

    Parameters
    ----------
    data  : input catalog
    var1  : first magnitude variable name
    var2  : second magnitude variable name
    color_range : color range [min, max]

    Returns
    -------
    selection
    """
    # NOOP; select everything
    if color_range is None:
        return np.ones(len(data),dtype=bool)
    
    color = data[var1] - data[var2]
    return (color > color_range[0]) & (color < color_range[1])

def select_class(data, var, values):
    """ Select objects based on discrete classes (star/galaxy).

    Parameters
    ----------
    data   : input catalog
    var    : classification variable
    values : list of allowed values

    Returns
    -------
    selection
    """
    if values is None:
        return np.ones(len(data), dtype=bool)
    else:
        return np.in1d(data[var], values)


def selection_function(data, bands, mag=True, color=False, extclass=None):
    """ Build a function to do a specific selection on mag and class.

    Parameters
    ----------
    bands       : bands to analyze
    mag         : perform magnitude selection
    color       : perform color selection
    extclass    : extended class cut ['star','gal']
    
    
    Returns
    -------
    func : selection function
    """
    sel = np.ones(len(data),dtype=bool)

    # NOOP; select everything
    if (bands == 'any') or (bands is None): 
        return sel

    # Perform magnitude selection
    if mag:
        for b in bands:
            mag_range = MAGLIM[b]
            sel &= select_mag(data, bfield(MAGVAR,b), mag_range)

    # Perform color selection
    if color:
        color_range, mag_range = COLOR[bands]

        if color_range is not None:
            magvar1 = bfield(MAGVAR,bands[0])
            magvar2 = bfield(MAGVAR,bands[1])
            sel &= select_color(data, magvar1, magvar2, color_range)

        if mag_range is not None:
            sel &= select_mag(data, magvar1, mag_range)

    # Perform extended class selection
    if extclass:
        for b in bands:
            extvalues = EXTCLASS[extclass]
            sel &= select_class(data, bfield(EXTVAR,b), extvalues)
            
    return sel

def count(filename, deredden=False):
    """Count the number of objects passing a selection.

    Parameters
    ----------
    filename : file containing the data
    deredden : flag to deredden magnitudes

    Returns
    -------
    ret : dict with [pix,cts] for each selection
    """
    print(filename)
    data = fitsio.read(filename,columns=COLUMNS)

    # Deredden the magnitudes
    if deredden:
        for c in COLUMNS:
            if not c.startswith(MAGVAR): continue
            data[c] = data[c] - data[bfield(EXTVAR,c[-1])]

    pixels = hp.ang2pix(nside,data['RA'],data['DEC'],lonlat=True)
    ret = dict()
    for name,(selfn,counts) in skymaps.items():
        pix,cts = np.unique(pixels[selfn(data)],return_counts=True)
        ret[name] = [pix,cts]
    return ret

def plot_footprint(filename,outfile=None,survey='delve',log=False):
    """ Plot the object density over the footprint """
    print("Reading %s..."%filename)

    counts = hp.read_map(filename,verbose=False)
    nside = hp.get_nside(counts)
    pixarea = hp.nside2pixarea(nside,degrees=True)*3600 #arcmin2
    hpxmap = counts/pixarea
    hpxmap[counts <= 0] = np.nan

    label = r'Object Density (arcmin$^{-2}$)'
    cbar_kwargs = dict(label=label)
    hpxmap_kwargs = dict(xsize=1000)
    if log: hpxmap_kwargs.update(norm=LogNorm())
    hist_kwargs = dict(bins=np.linspace(1/pixarea,100/pixarea,51))
    fig,axes,smap = plotting.plot_hpxmap_hist(hpxmap,survey,cbar_kwargs,hpxmap_kwargs,hist_kwargs)
    axes[1].set_xlabel(label)

    if outfile is None: 
        outfile=os.path.basename(filename).split('.')[0]+'.png'
    print("Writing %s..."%outfile)
    plt.savefig(outfile,bbox_inches='tight')
    return fig,axes,smap

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config',help='configuration file')
    parser.add_argument('-b','--band',dest='bands',action='append',default=None)
    parser.add_argument('-c','--color',action='store_true')
    parser.add_argument('-d','--deredden',action='store_true')
    parser.add_argument('-e','--extclass',default=None,choices=['star','gal'])
    parser.add_argument('-m','--mag',action='store_true')
    parser.add_argument('-n','--nside',default=1024,type=int)
    parser.add_argument('-p','--plot',action='store_true')
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('--nproc',default=4,type=int)

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))
    nside = args.nside
    survey = config['survey']
    catdir = config['catdir']

    pixarea = hp.nside2pixarea(nside,degrees=True)
    if not args.extclass:
        outdir = mkdir('release/footprint')
        outbase = 'footprint_%s_n%i.fits.gz'
    else:
        outdir = mkdir('release/stargal')
        outbase = 'footprint_{}_%s_n%i.fits.gz'.format(args.extclass)

    bands = config['bands']
    if not args.bands: 
        args.bands = bands + [''.join(bands)] + ['any']

    print("Preparing to count objects...")
    print("  bands: %s"%args.bands)
    print("  color: %s"%args.color)
    print("  extclass: %s"%args.extclass)
    print("  deredden: %s"%args.deredden)

    COLUMNS = create_columns(bands,extclass=args.extclass,deredden=args.deredden)
    print("  columns: %s"%COLUMNS)

    skymaps = odict()
    for i,k in enumerate(args.bands):
        print("Adding %s-band selection..."%k)
        #skymaps[k] = [SELECT[k], empty(nside)]
        #fn = selection_function(k, extclass=EXTCLASS[args.extclass])
        kwargs = dict(bands=k, mag=args.mag, color=args.color, extclass=args.extclass)
        fn = functools.partial(selection_function, **kwargs)
        skymaps[k] = [fn, empty(nside)]
            
    filenames = sorted(glob.glob(catdir+'/cat_hpx_*.fits'))
    #filenames = filenames[:2000]
    print("Processing %i files..."%len(filenames))

    # Launch the jobs
    func    = count
    arglist = list(zip(filenames))
    kwargs  = dict(deredden=args.deredden)
    results = utils.multiproc(func,arglist,kwargs,processes=args.nproc)

    if None in results:
        print("WARNING: %i processes failed..."%results.count(None))

    print("Filling counts...")
    for name, (selfn,counts) in skymaps.items():
        for res in results: 
            if res is None: continue
            pix,cts = res[name]
            counts[pix] = cts

    print("Writing maps...")
    filenames = []
    for name,(selfn,counts) in skymaps.items():
        outfile = join(outdir,outbase%(name,nside))

        print("Writing %s..."%outfile)
        if os.path.exists(outfile): os.remove(outfile)
        hp.write_map(outfile,counts,dtype=int)
        filenames.append(outfile)

    if args.plot:
        print("Plotting maps...")
        for name in skymaps.keys():
            outfile = join(outdir,outbase%(name,nside))
            pngfile = outfile.replace('.fits.gz','.png')
            plot_footprint(outfile,pngfile,survey=survey)
