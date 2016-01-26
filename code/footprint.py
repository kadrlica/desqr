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

from const import OBJECT_ID, UNIQUE_ID, BANDS, BADMAG, NSIDES
from utils import bfields, load_infiles, mkdir
from utils import ang2pix,pix2ang
import plotting

# Coverage Footprint

COLUMNS = ['RA','DEC']+bfields('WAVG_MAG_PSF',BANDS)

def empty(nside): 
    return np.zeros(healpy.nside2npix(nside),dtype=int)

def blank(nside): 
    return np.nan*np.ones(healpy.nside2npix(nside),dtype=int)

def unseen(nside): 
    return healpy.UNSEEN*np.ones(healpy.nside2npix(nside),dtype=int)

def select(data,mag,maglim=90):
    return data[mag] < maglim

sel_all   = lambda x: np.ones(len(x),dtype=bool)
sel_g     = lambda x: select(x,'WAVG_MAG_PSF_G')
sel_r     = lambda x: select(x,'WAVG_MAG_PSF_R') 
sel_i     = lambda x: select(x,'WAVG_MAG_PSF_I') 
sel_z     = lambda x: select(x,'WAVG_MAG_PSF_Z') 
sel_Y     = lambda x: select(x,'WAVG_MAG_PSF_Y') 
sel_gr    = lambda x: sel_g(x) & sel_r(x)
sel_iz    = lambda x: sel_i(x) & sel_z(x)
sel_grizY = lambda x: sel_g(x) & sel_r(x) & sel_i(x) & sel_z(x) & sel_Y(x)

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-n','--nside',default=1024,type=int)
    parser.add_argument('-v','--verbose',action='store_true')
    opts = parser.parse_args()

    nside = opts.nside
    pixarea = healpy.nside2pixarea(nside,degrees=True)
    plotdir = mkdir('release/footprint')

    skymaps = [
        ['any'  ,sel_all  ],
        #['g'    ,sel_g    ],
        #['r'    ,sel_r    ],
        ['i'    ,sel_i    ],
        ['z'    ,sel_z    ],
        #['Y'    ,sel_Y    ],
        #['gr'   ,sel_gr   ],
        ['iz'   ,sel_iz   ],
        #['all',sel_grizY],
        ]

    for i,(n,s) in enumerate(skymaps):
        skymaps[i] += [empty(nside)]

    infiles = sorted(glob.glob('cat/cat_hpx_*.fits'))
    for f in infiles:
        if opts.verbose: print f
        data = fitsio.read(f,columns=COLUMNS)
        pixels = ang2pix(nside,data['RA'],data['DEC'])
        for name,sel,counts in skymaps:
            s = sel(data)
            pix,cts = np.unique(pixels[s],return_counts=True)
            counts[pix] += cts

    outstr = '|_. Band |_. Footprint|_. Area (deg^2) |\n'
    template = '|_. %(band)s |{{thumbnail(%(map)s, size=300}}|_. %(area).2f |\n'
    out = dict()
    for name,sel,counts in skymaps:
        out['band'] = name
        density = counts/pixarea
        skymap = np.ma.MaskedArray(np.log10(density),counts == 0,fill_value=np.nan)
        #skymap = np.ma.MaskedArray(counts.view(float),counts == 0,fill_value=np.nan)
        area = pixarea * (counts > 0).sum() 
        out['area'] = area
        plt.figure()
        im = plotting.draw_footprint(skymap,cmap='jet')
        plt.colorbar(im,label=r'$\log_{10}({\rm Density}\ [\deg^{-2}])$')
        #plt.text(-15,-15,r'${\rm Area} = %.0f\ \deg^2$'%area,fontsize=16)
        plt.title("Y2Q1 Coverage (%s)"%name)
        outfile = join(plotdir,'y2q1_footprint_%s_n%i_car.png'%(name,nside))
        print "Writing %s..."%outfile
        plt.savefig(outfile,bbox_inches='tight')
        out['map'] = os.path.basename(outfile)
        outfile = join(plotdir,'y2q1_footprint_%s_n%i_equ.fits'%(name,nside))
        print "Writing %s..."%outfile
        #healpy.write_map(outfile,skymap.data,dtype=int)
        outstr += template%out
     
    print outstr

    kwargs = dict(bins=np.linspace(1,100,100),histtype='step',lw=1.5)
    plt.figure()
    for name,sel,counts in skymaps:
        plt.hist(counts,label=name,**kwargs)

    plt.legend(fontsize=10)    
    plt.xlabel('Number of Pixels')
    plt.ylabel('Counts per Pixel')
    outfile = join(plotdir,'y2q1_footprint_hist.png')
    plt.savefig(outfile,bbox_inches='tight')
