#!/usr/bin/env python
"""
Map the footprint.
"""
import glob
import os
from os.path import join
import matplotlib
if not os.getenv('DISPLAY'):    matplotlib.use('Agg')
if os.getenv('TERM')=='screen': matplotlib.use('Agg')
from collections import OrderedDict as odict
from multiprocessing import Pool

import fitsio
import numpy as np
import scipy.ndimage as nd
import pylab as plt
import matplotlib.colors as colors
import healpy as hp
import yaml

from const import OBJECT_ID, UNIQUE_ID, BANDS, BADMAG, NSIDES
from utils import bfields, load_infiles, mkdir
from utils import ang2pix,pix2ang
import plotting

# Coverage Footprint
def create_columns(bands=BANDS):
    COLUMNS = ['RA','DEC']+bfields('WAVG_MAG_PSF',bands)
    #COLUMNS = ['RA','DEC']+bfields('WAVG_MAG_PSF',BANDS)+bfields('NEPOCHS',BANDS)
    return COLUMNS
    
def empty(nside): 
    return np.zeros(hp.nside2npix(nside),dtype=int)

def blank(nside): 
    return np.nan*np.ones(hp.nside2npix(nside),dtype=int)

def unseen(nside): 
    return hp.UNSEEN*np.ones(hp.nside2npix(nside),dtype=int)

def select(data,mag,maglim=90):
    return data[mag] < maglim

def count(filename):
    """Count the number of objects passing selections"""
    print(filename)
    data = fitsio.read(filename,columns=COLUMNS)
    pixels = ang2pix(nside,data['RA'],data['DEC'])
    ret = dict()
    for name,(sel,counts) in skymaps.items():
        pix,cts = np.unique(pixels[sel(data)],return_counts=True)
        ret[name] = [pix,cts]
    return ret

def plot(filename,survey=None,title=None):
    """ Plot output file """
    counts = hp.read_map(filename)
    pixarea = hp.get_nside(counts)
    hpxmap = np.ma.MaskedArray(np.log10(counts/pixarea), counts <= 0, fill_value=np.nan)

    plt.figure(figsize=(10, 6.18))
    if survey == 'des':
        from skymap.survey import DESSkymap
        smap = DESSkymap()
        smap.draw_des()
        #im = plotting.draw_des(skymap,cmap='viridis')
    elif survey == 'maglites':
        from skymap.survey import MaglitesSkymap
        smap = MaglitesSkymap()
        smap.draw_maglites()
        #smap.draw_lmc(); smap.draw_smc()
        #im = plotting.draw_maglites(skymap,cmap='viridis')
    elif survey == 'bliss':
        from skymap.survey import BlissSkymap
        smap = BlissSkymap()
        smap.draw_bliss()
        #im = plotting.draw_bliss(skymap,cmap='viridis')
    elif survey == 'delve':
        from skymap.survey import SurveyMcBryde
        smap = SurveyMcBryde(lon_0 = 180)
    else:
        from skymap.survey import SurveySkymap
        smap = SurveySkymap()
        #im = plotting.draw_footprint(skymap,cmap='viridis')
        
    im,_lon,_lat,_val = smap.draw_hpxmap(hpxmap,cmap='viridis')
    plt.colorbar(im,label=r'$\log_{10}({\rm Density}\ [\deg^{-2}])$')
    plt.suptitle("Coverage (%s)"%title)

    outfile = filename.replace('.fits.gz','.png')
    print("Writing %s..."%outfile)
    plt.savefig(outfile,bbox_inches='tight')
    

sel_g     = lambda x: select(x,'WAVG_MAG_PSF_G')
sel_r     = lambda x: select(x,'WAVG_MAG_PSF_R')
sel_i     = lambda x: select(x,'WAVG_MAG_PSF_I')
sel_z     = lambda x: select(x,'WAVG_MAG_PSF_Z')
sel_Y     = lambda x: select(x,'WAVG_MAG_PSF_Y')
sel_g     = lambda x: select(x,'WAVG_MAG_PSF_G',23.0)
sel_r     = lambda x: select(x,'WAVG_MAG_PSF_R',23.0) 
sel_i     = lambda x: select(x,'WAVG_MAG_PSF_I',22.0) 
sel_z     = lambda x: select(x,'WAVG_MAG_PSF_Z',22.0) 
sel_Y     = lambda x: select(x,'WAVG_MAG_PSF_Y',21.0) 

sel_any   = lambda x: np.ones(len(x),dtype=bool)
sel_gr    = lambda x: sel_g(x) & sel_r(x)
sel_gi    = lambda x: sel_g(x) & sel_i(x)
sel_iz    = lambda x: sel_i(x) & sel_z(x)
sel_griz  = lambda x: sel_g(x) & sel_r(x) & sel_i(x) & sel_z(x)
sel_grizY = lambda x: sel_g(x) & sel_r(x) & sel_i(x) & sel_z(x) & sel_Y(x)
sel_all   = sel_grizY

SELECT = odict([
    ['any' ,sel_any  ],
    ['g'   ,sel_g    ],
    ['r'   ,sel_r    ],
    ['i'   ,sel_i    ],
    ['z'   ,sel_z    ],
    ['Y'   ,sel_Y   ],
    ['gr'  ,sel_gr  ],
    ['gi'  ,sel_gi  ],
    ['iz'  ,sel_iz  ],
    ['griz' ,sel_griz ],
    ['grizY' ,sel_grizY ],
    ])

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config')
    parser.add_argument('-n','--nside',default=1024,type=int)
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('-s','--survey',default='des')
    parser.add_argument('-b','--band',action='append',default=None)
    parser.add_argument('-p','--plot',action='store_true')

    args = parser.parse_args()

    config = yaml.load(open(args.config))
    nside = args.nside
    pixarea = hp.nside2pixarea(nside,degrees=True)
    outdir = mkdir('release/footprint')
    bands = config['bands']
    if not args.band: 
        args.band = ['any'] + bands + [''.join(bands)]
    catdir = config['catdir']

    COLUMNS = create_columns(bands)
    skymaps = odict()
    for i,k in enumerate(args.band):
        print("Adding %s-band selection..."%k)
        skymaps[k] = [SELECT[k],empty(nside)]

    print("Selecting data...")
    infiles = sorted(glob.glob(catdir+'/cat_hpx_*.fits'))
    pool = Pool(maxtasksperchild=1,processes=30)
    results = pool.map(count,infiles)

    print("Filling counts...")
    for name, (sel,counts) in skymaps.items():
        for res in results:
            pix,cts = res[name]
            counts[pix] = cts

    outstr = '|_. Band |_. Footprint|_. Area (deg^2) |\n'
    template = '|_. %(band)s |{{thumbnail(%(map)s, size=300)}}|_. %(area).2f |\n'

    out = dict()
    for name,(sel,counts) in skymaps.items():
        if counts.sum() == 0:
            print("No counts for selection: %s"%name)
            continue

        out['band'] = name
        out['area'] = pixarea * (counts > 0).sum() 
        out['map'] = ''

        outfile = join(outdir,'footprint_%s_n%i_equ.fits.gz'%(name,nside))
        print("Writing %s..."%outfile)
        if os.path.exists(outfile): os.remove(outfile)
        hp.write_map(outfile,counts,dtype=int)

        if args.plot:
            pngfile = plot(outfile,survey=args.survey,title=name)
            out['map'] = pngfile

        outstr += template%out

    print(outstr)

    if args.plot:
        kwargs = dict(bins=np.linspace(1,100,100),histtype='step',lw=1.5)
        plt.figure()
        for name,(sel,counts) in skymaps.items():
            plt.hist(counts,label=name,**kwargs)
         
        plt.legend(fontsize=10)    
        plt.xlabel('Counts per Pixel')
        plt.ylabel('Number of Pixels')
        outfile = join(outdir,'footprint_hist_n%i.png'%nside)
        plt.savefig(outfile,bbox_inches='tight')
        plt.ion()
