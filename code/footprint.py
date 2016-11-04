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
#COLUMNS = ['RA','DEC']+bfields('WAVG_MAG_PSF',BANDS)+bfields('NEPOCHS',BANDS)

def empty(nside): 
    return np.zeros(healpy.nside2npix(nside),dtype=int)

def blank(nside): 
    return np.nan*np.ones(healpy.nside2npix(nside),dtype=int)

def unseen(nside): 
    return healpy.UNSEEN*np.ones(healpy.nside2npix(nside),dtype=int)

def select(data,mag,maglim=90):
    return data[mag] < maglim

def count(filename):
    print filename
    data = fitsio.read(filename,columns=COLUMNS)
    pixels = ang2pix(nside,data['RA'],data['DEC'])
    ret = dict()
    for name,sel,counts in skymaps:
        s = sel(data)
        pix,cts = np.unique(pixels[s],return_counts=True)
        ret[name] = [pix,cts]
    return ret

sel_g     = lambda x: select(x,'WAVG_MAG_PSF_G')
sel_r     = lambda x: select(x,'WAVG_MAG_PSF_R') 
sel_i     = lambda x: select(x,'WAVG_MAG_PSF_I') 
sel_z     = lambda x: select(x,'WAVG_MAG_PSF_Z') 
sel_Y     = lambda x: select(x,'WAVG_MAG_PSF_Y') 
sel_any   = lambda x: np.ones(len(x),dtype=bool)
sel_gr    = lambda x: sel_g(x) & sel_r(x)
sel_iz    = lambda x: sel_i(x) & sel_z(x)
sel_griz  = lambda x: sel_g(x) & sel_r(x) & sel_i(x) & sel_z(x)
sel_grizY = lambda x: sel_g(x) & sel_r(x) & sel_i(x) & sel_z(x) & sel_Y(x)
sel_all   = sel_grizY
#sel_g     = lambda x: select(x,'WAVG_MAG_PSF_G',21.0)
#sel_r     = lambda x: select(x,'WAVG_MAG_PSF_R',21.0) 
#sel_i     = lambda x: select(x,'WAVG_MAG_PSF_I',21.0) 
#sel_z     = lambda x: select(x,'WAVG_MAG_PSF_Z',21.0) 
#sel_Y     = lambda x: select(x,'WAVG_MAG_PSF_Y',21.0) 

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config',nargs='?')
    parser.add_argument('-n','--nside',default=1024,type=int)
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('-s','--survey',default='des')
    args = parser.parse_args()

    nside = args.nside
    pixarea = healpy.nside2pixarea(nside,degrees=True)
    plotdir = mkdir('release/footprint')

    skymaps = [
        ['any' ,sel_any ],
        ['all' ,sel_all ],
        ['g'   ,sel_g   ],
        ['r'   ,sel_r   ],
        ['i'   ,sel_i   ],
        ['z'   ,sel_z   ],
        ['Y'   ,sel_Y   ],
        #['gr'  ,sel_gr  ],
        #['iz'  ,sel_iz  ],
        ]

    for i,(n,s) in enumerate(skymaps):
        skymaps[i] += [empty(nside)]

    infiles = sorted(glob.glob('cat/cat_hpx_*.fits'))
    for f in infiles:
        if args.verbose: print f
        data = fitsio.read(f,columns=COLUMNS)
        #data = data[(data['NEPOCHS_G'] == 1)]
        pixels = ang2pix(nside,data['RA'],data['DEC'])
        for name,sel,counts in skymaps:
            s = sel(data)
            pix,cts = np.unique(pixels[s],return_counts=True)
            counts[pix] += cts

    #p = Pool(maxtasksperchild=1)
    #out = p.map(footprint,infiles)
    #for name,sel,counts in skymaps:
    #    out

    #data = load_infiles(infiles,columns=COLUMNS)
    #pixels = ang2pix(nside,data['RA'],data['DEC'])

    outstr = '|_. Band |_. Footprint|_. Area (deg^2) |\n'
    template = '|_. %(band)s |{{thumbnail(%(map)s, size=300)}}|_. %(area).2f |\n'

    out = dict()
    for name,sel,counts in skymaps:
        out['band'] = name
        density = counts/pixarea
        skymap = np.ma.MaskedArray(np.log10(density),counts == 0,fill_value=np.nan)
        area = pixarea * (counts > 0).sum() 
        out['area'] = area
        plt.figure()
        if args.survey == 'des':
            im = plotting.draw_des(skymap,cmap='viridis')
        elif args.survey == 'maglites':
            im = plotting.draw_maglites(skymap,cmap='viridis')
        else:
            im = plotting.draw_footprint(skymap,cmap='viridis')
        plt.colorbar(im,label=r'$\log_{10}({\rm Density}\ [\deg^{-2}])$')
        plt.title("Coverage (%s)"%name)
        outfile = join(plotdir,'footprint_%s_n%i_car.png'%(name,nside))
        print "Writing %s..."%outfile
        plt.savefig(outfile,bbox_inches='tight')
        out['map'] = os.path.basename(outfile)
        outfile = join(plotdir,'footprint_%s_n%i_equ.fits'%(name,nside))
        print "Writing %s..."%outfile
        healpy.write_map(outfile,counts,dtype=int)
        outstr += template%out
     

    print outstr

    kwargs = dict(bins=np.linspace(1,100,100),histtype='step',lw=1.5)
    plt.figure()
    for name,sel,counts in skymaps:
        plt.hist(counts,label=name,**kwargs)

    plt.legend(fontsize=10)    
    plt.xlabel('Counts per Pixel')
    plt.ylabel('Number of Pixels')
    outfile = join(plotdir,'y2q1_footprint_hist_n%i.png'%nside)
    plt.savefig(outfile,bbox_inches='tight')
