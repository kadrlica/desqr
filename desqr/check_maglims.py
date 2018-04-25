import numpy
import healpy
import pylab

import ugali.utils.skymap
import ugali.utils.healpix

pylab.ion()

infile = 'maglim_hpx_02662.fits'
pix_16 = 2662
ra_center, dec_center = ugali.utils.healpix.pixToAng(16, pix_16)

bands = ['G', 'R', 'I', 'Z', 'Y']
for band in bands:
    m_maglim = ugali.utils.skymap.readSparseHealpixMap(infile, 'MAGLIM_%s'%(band))
    m_maglim[m_maglim < -1.] = healpy.UNSEEN

    cut = (m_maglim != healpy.UNSEEN)
    maglim_sort = numpy.sort(m_maglim[cut])
    maglim_vmin = maglim_sort[int(len(maglim_sort) * 0.005)]
    maglim_vmax = maglim_sort[int(len(maglim_sort) * 0.995)]
    
    title = '%05i_%s'%(pix_16, band)

    reso = 0.25
    xsize = int(8. * 60. / reso)
    #pylab.figure(figsize=(6, 6))
    healpy.gnomview(m_maglim, rot=(ra_center, dec_center, 0.), xsize=xsize, reso=reso, unit='MAGLIM_%s'%(band), title=title, min=maglim_vmin, max=maglim_vmax)
    
