#!/usr/bin/env python
import fitsio
import numpy as np
import pylab as plt

from const import BANDS, COLORS
from zeropoint import calc_mag,calc_magerr
from depth import draw_magerr

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    opts = parser.parse_args()

    survey = np.recfromcsv('y2n_y1a1_survey_v3.csv')
    zp = fitsio.read('y2n_y1a1_qslr_v6.fits')

    teff = 2
    
    for band in BANDS:
        kwargs = dict(histtype='step',lw=1.5,bins=np.linspace(16,25,100),color=COLORS[band])    
        y1a1sel = (survey['tag'] == 'Y1A1_FINALCUT')&(survey['band'] == band)
        y2nsel = (survey['tag'] == 'Y2N_FIRSTCUT')&(survey['band'] == band)

        sel = np.argmin(np.abs(survey[y1a1sel]['t_eff'] - teff))
        y1a1exp = survey[y1a1sel][sel]['expnum']
        y1a1file = 'raw/%s/D%08d_%s_cat.fits'%(band,y1a1exp,band)
        print 'Y1A1:',y1a1file
        y1a1 = fitsio.read(y1a1file)
        y1a1zp = np.median(zp['MAG_ZERO'][zp['EXPNUM']==y1a1exp])
        y1a1mag = calc_mag(y1a1['FLUX_PSF'],y1a1zp)
        y1a1magerr = calc_magerr(y1a1['FLUX_PSF'],y1a1['FLUXERR_PSF'])

        sel = np.argmin(np.abs(survey[y2nsel]['t_eff'] - teff))
        y2nexp = survey[y2nsel][sel]['expnum']
        y2nfile = 'raw/%s/D%08d_%s_cat.fits'%(band,y2nexp,band)
        print 'Y2N:',y2nfile
        y2n = fitsio.read(y2nfile)
        y2nzp = np.median(zp['MAG_ZERO'][zp['EXPNUM']==y2nexp])
        y2nmag = calc_mag(y2n['FLUX_PSF'],y2nzp)
        y2nmagerr = calc_magerr(y2n['FLUX_PSF'],y2n['FLUXERR_PSF'])

        plt.hist(y1a1mag,label='Y1A1',**kwargs)
        plt.hist(y2nmag,alpha=0.5,label='Y2N',**kwargs)
        plt.legend(fontsize=10,loc='upper left')

        plt.figure()
        draw_magerr(y1a1mag,y1a1magerr,color=kwargs['color'],label='Y1A1')
        draw_magerr(y2nmag,y2nmagerr,color=kwargs['color'],alpha=0.5,label='Y2N')
        plt.legend(fontsize=10,loc='upper left')

        break
plt.ion()
