#!/usr/bin/env python
"""
Checking and comparing photometric calibration zeropoints.
"""
import glob
import os
from os.path import join
import subprocess
import matplotlib
try:             os.environ['DISPLAY']
except KeyError: matplotlib.use('Agg')

import fitsio
import numpy as np
import scipy.ndimage as nd
import pylab as plt
import matplotlib.colors as colors
import healpy
import numpy.lib.recfunctions as recfn

from multiprocessing import Pool

from ugali.utils.healpix import ang2pix
from ugali.utils.projector import angsep
from ugali.utils.shell import mkdir

from const import OBJECT_ID, UNIQUE_ID, BANDS, BADMAG, NSIDES
from utils import bfields, load_infiles, uid
from footprint import blank
import footprint
import plotting
import download

COLORS = [
    ('g-r',['g','r']),
    ('r-i',['r','i']),
    ('i-z',['i','z']),
    ('z-Y',['z','Y']),
]
    
if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    opts = parser.parse_args()

    nside = 128
    HPX = 'HPX%i'%nside
    outdir = mkdir('release/calibration')

    gcmfile = 'y1a1_gcm_v0.fits'
    if not os.path.exists(gcmfile):
        query = download.gcm_query()
        download.download(gcmfile,query,section='dessci')

    print "Loading %s..."%gcmfile
    gcm = fitsio.read(gcmfile)
    names = ['MAG_ZERO',HPX]
    values = [gcm['ZEROPOINT']+25,ang2pix(nside,gcm['RA'],gcm['DEC'])]
    gcm = recfn.append_fields(gcm,names,values,usemask=False,asrecarray=True)
    gcm['BAND'] = np.char.strip(gcm['BAND'])
    
    
    qslrfile = 'y2n_y1a1_qslr_v6.fits'
    if not os.path.exists(qslrfile):
        query = download.qslr_query()
        download.download(qslrfile,query,section='dessci')
    print "Loading %s..."%qslrfile
    qslr = fitsio.read(qslrfile)
    names = [HPX]
    values = [ang2pix(nside,qslr['RA_MEAN'],qslr['DEC_MEAN'])]
    qslr = recfn.append_fields(qslr,names,values,usemask=False,asrecarray=True)
    qslr = qslr[qslr['QSLR_FLAG']==0]
        
    gcm_uid = uid(gcm['EXPNUM'],gcm['CCDNUM'])
    qslr_uid = uid(qslr['EXPNUM'],qslr['CCDNUM'])
     
    gcm = gcm[np.in1d(gcm_uid,qslr_uid)]
    qslr = qslr[np.in1d(qslr_uid,gcm_uid)]

    gcm = gcm[np.lexsort((gcm['CCDNUM'],gcm['EXPNUM']))]
    qslr = qslr[np.lexsort((qslr['CCDNUM'],qslr['EXPNUM']))]
     
    if (gcm['EXPNUM']!=qslr['EXPNUM']).any() \
            or (gcm['CCDNUM']!=qslr['CCDNUM']).any():
        msg = "GCM and qSLR do not match"
        raise Exception(msg)
     
gcm = gcm[np.lexsort((gcm['CCDNUM'],gcm['EXPNUM']))]
qslr = qslr[np.lexsort((qslr['CCDNUM'],qslr['EXPNUM']))]

if (gcm['EXPNUM']!=qslr['EXPNUM']).any() \
        or (gcm['CCDNUM']!=qslr['CCDNUM']).any():
    msg = "GCM and qSLR do not match"
    raise Exception(msg)

dzero = gcm['MAG_ZERO']-qslr['MAG_ZERO']
sep = angsep(gcm['RA'],gcm['DEC'],qslr['RA_MEAN'],qslr['DEC_MEAN'])

plt.figure()
plotting.draw_hist(dzero,normed=False)
plt.legend(loc='upper right')
plt.xlabel(r'${\rm ZP_{GCM} - ZP_{qSLR}}')
plt.ylabel("Number of CCDs")



gcm_bands = dict()
qslr_bands = dict()
diff_bands = dict()
for b in BANDS:
    zp_g = gcm['MAG_ZERO'][gcm['BAND'] == b]
    zp_q = qslr['MAG_ZERO'][qslr['BAND'] == b]
    delta = zp_g - zp_q
    labels = qslr[HPX][qslr['BAND'] == b]
    index = np.unique(labels)
    for data,out in [(zp_g,gcm_bands),(zp_q,qslr_bands),(delta,diff_bands)]:
        skymap = blank(nside)
        skymap[index] = nd.median(data,labels=labels,index=index)   
        out[b] = skymap

diff_colors = dict()
for name,(b1,b2) in COLORS:
    gcm_color = gcm_bands[b1] - gcm_bands[b2]
    qslr_color = qslr_bands[b1] - qslr_bands[b2]
    skymap = gcm_color - qslr_color
    diff_colors[name] = skymap


for b in BANDS:
    skymap = diff_bands[b]
    plt.figure()
    im = plotting.draw_footprint(skymap)
    plt.colorbar(im,label=r'${\rm ZP_{GCM} - ZP_{qSLR}}$')
    plt.title('Median Zeropoint Offset (%s-band)'%(b))
    outfile = join(outdir,'y2q1_zeropoint_%s_n%i_car.png'%(b,nside))
    plt.savefig(outfile,bbox_inches='tight')

for name,(b1,b2) in COLORS:
    skymap = diff_colors[name]
    plt.figure()
    im = plotting.draw_footprint(skymap)
    plt.colorbar(im,label=r'${\rm (%s)_{GCM} - (%s)_{qSLR}}$'%(name,name))
    plt.title('Color Offset (%s)'%(name))
    outfile = join(outdir,'y2q1_colors_%s_n%i_car.png'%(name,nside))
    plt.savefig(outfile,bbox_inches='tight')


    plt.figure()
    im = plotting.draw_hist(skymap[~np.isnan(skymap)])
    plt.title('Color Offset (%s)'%(name))
    plt.xlabel(r'${\rm (%s)_{GCM} - (%s)_{qSLR}}$'%(name,name))
    plt.ylabel('Normalized Number of CCDs')
    plt.legend(loc='upper left')
    outfile = join(outdir,'y2q1_colors_%s_n%i_hist.png'%(name,nside))
    plt.savefig(outfile,bbox_inches='tight')
