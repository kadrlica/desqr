import numpy
import healpy
import pyfits
import scipy.stats
import pylab

import ugali.utils.healpix
import ugali.utils.projector
import ugali.utils.skymap

pylab.ion()

############################################################

def maglim(infile, outfile=None, nside=2**10, signal_to_noise=10., save=False, plot=False):
    reader = pyfits.open(infile)
    data = reader[1].data
    reader.close()

    pix_16 = ugali.utils.healpix.angToPix(16, numpy.median(data['RA']), numpy.median(data['DEC']))
    ra_center, dec_center = ugali.utils.healpix.pixToAng(16, pix_16)
    proj = ugali.utils.projector.Projector(ra_center, dec_center)

    pix_nside = ugali.utils.healpix.subpixel(pix_16, 16, nside)
    dict_maglim = {}

    bands  = ['G', 'R', 'I', 'Z', 'Y']
    for band in bands:
        # Select an appropriate bright magnitude for the analysis
        h, edges = numpy.histogram(data['WAVG_MAG_PSF_%s'%band], bins=numpy.arange(15, 30, 0.1))
        mag_bright_end = edges[numpy.argmax(h)] - 3.
        
        cut = (numpy.fabs(data['WAVG_SPREAD_MODEL_%s'%(band)]) < 0.002) & (data['WAVG_MAG_PSF_%s'%band] > mag_bright_end) & (data['WAVG_MAG_PSF_%s'%(band)] < 30.)
        x, y = proj.sphereToImage(data['RA'][cut], data['DEC'][cut])
        results = ugali.utils.projector.match(data['RA'][cut], data['DEC'][cut], data['RA'][cut], data['DEC'][cut], nnearest=2)

        delta_mag = data['MAG_PSF_%s'%(band)][cut][results[1]] - data['MAG_PSF_%s'%(band)][cut][results[0]]
        delta_log_magerr = numpy.log10(data['MAGERR_PSF_%s'%(band)][cut][results[1]]) - numpy.log10(data['MAGERR_PSF_%s'%(band)][cut][results[0]])

        cut_nan_inf = ~numpy.isnan(delta_log_magerr / delta_mag) & ~numpy.isinf(delta_log_magerr / delta_mag) & (delta_mag > 0.5)
    
        kde = scipy.stats.gaussian_kde(delta_log_magerr[cut_nan_inf] / delta_mag[cut_nan_inf])

        values = numpy.linspace(0., 1., 1000)
        kde_values = kde.evaluate(values)
        slope = values[numpy.argmax(kde_values)]
        
        maglim = data['MAG_PSF_%s'%(band)][cut] - ((numpy.log10(data['MAGERR_PSF_%s'%(band)][cut]) - numpy.log10(1. / signal_to_noise)) / slope)

        # HEALPix results 

        npix = healpy.nside2npix(nside)
        pix = ugali.utils.healpix.angToPix(nside, data['RA'][cut], data['DEC'][cut])
        m = numpy.histogram(pix, bins=numpy.arange(0, npix + 1))[0].astype(float)
        m[m == 0.] = healpy.UNSEEN

        for pix_select in numpy.nonzero(m > 0)[0]:
            cut_pix_select = (pix == pix_select)
            #if numpy.sum(cut_pix_select) >= 10:
            if numpy.sum(cut_pix_select) > 0:
                m[pix_select] = numpy.median(maglim[cut_pix_select])
            else:
                m[pix_select] = healpy.UNSEEN

        dict_maglim['MAGLIM_%s'%(band)] = m[pix_nside]

    if outfile is None:
        return pix_nside,dict_maglim
    else:
        ugali.utils.skymap.writeSparseHealpixMap(pix_nside, dict_maglim, nside, outfile)

############################################################

# Testing
#infile = '/project/kicp/bechtol/des/mw_substructure/y2n/data/catalog/v6/hpx/cat_hpx_02662.fits'
#outfile = 'maglim_hpx_02662.fits'
# 
#maglim(infile, outfile)

############################################################

