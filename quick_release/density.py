#!/usr/bin/env python
from os.path import exists, join
from multiprocessing import Pool
import glob

import fitsio
import numpy as np
import pylab as plt
import healpy as hp

from utils import pix2ang, ang2pix
from utils import load_infiles
from plotting import draw_footprint

NSIDE = 128
MAG_AUTO = 'WAVG_MAG_AUTO_I'
#MAG_AUTO = 'MAG_AUTO_I'
MAG_PSF = 'WAVG_MAG_PSF_I'
MAG = MAG_AUTO
SPREAD = 'WAVG_SPREAD_MODEL_I'
COLUMNS = [
    'HPX2048',
    MAG,
    SPREAD,
    #'WAVG_CLASS_STAR_I',
    ]
OBJTYPE='Galaxies'

def define_cut(data):
    """ 
    This is a function for defining cuts to apply to the data
    before calculating the density.
    """
    ## This is a non-op for the time being
    #cut = slice(None)

    # Magnitude cut
    #cut = (data[MAG]>21)&(data[MAG] < 22)
    cut = (data[MAG] > 19.5)&(data[MAG] < 20.5)
    #cut = (data[MAG] > 20)&(data[MAG] < 21)
    #cut = (data[MAG] > 22)

    # Cut for stars
    #cut &= (np.abs(data[SPREAD]) < 0.003)

    # Cut for galaxies
    cut &= (data[SPREAD] > 0.005)
    #cut &= (np.abs(data['WAVG_CLASS_STAR_I']) < 0.2)

    return cut

def density(infile, nside=NSIDE):
    """
    This is the function called by the multiprocessing pool to
    calculate the density for each infile.
    """
    area = hp.nside2pixarea(nside,degrees=True)
    scale = int(np.log(2048/NSIDE)/np.log(2))

    data = fitsio.read(infile,columns=COLUMNS)
    cut = define_cut(data)
    
    hpx = data['HPX2048'][cut]//4**scale
    idx,cts = np.unique(hpx,return_counts=True)

    return [idx,cts/area]

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    opts = parser.parse_args()

    dirname = './cat'
    infiles = sorted(glob.glob(join(dirname,'cat_hpx_*.fits')))

    print "Running %i files..."%len(infiles)
    p = Pool(maxtasksperchild=1)
    out = p.map(density,infiles)

    idx,density = np.hstack(out)    
    idx = idx.astype(int)
    skymap = hp.UNSEEN*np.ones(hp.nside2npix(NSIDE))
    skymap[idx] = density
    skymap = hp.reorder(skymap,'NEST','RING')

fig,ax = plt.subplots(1,2,figsize=(16,6))
plt.sca(ax[0])
im = draw_footprint(skymap)
plt.colorbar(im,label=r'${\rm Density\ (deg^{-2})}$')
plt.sca(ax[1])
im = draw_footprint(np.log10(skymap))
plt.colorbar(im,label=r'${\rm \log_{10}(Density\ (deg^{-2}))}$')
fig.suptitle('%s Density'%OBJTYPE)
plt.savefig('density_%s.png'%OBJTYPE.lower(),bbox_inches='tight')

plt.ion()
plt.show()
