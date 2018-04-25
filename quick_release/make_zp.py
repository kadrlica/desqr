#!/usr/bin/env python
"""
Generic python script.
"""
import numpy as np
import fitsio

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infile',default='desgw_survey_v1.csv')
    parser.add_argument('outfile',default='desgw_zp_v1.fits')
    args = parser.parse_args()

    data = np.recfromcsv(args.infile)
    
    dtype = [('FILENAME','S60'),('EXPNUM','i4'),('CCDNUM','i4'),('BAND','S1'),
             ('MAG_ZERO','f4'),('SIGMA_MAG_ZERO','f4'),
             ('RA_MEAN','f8'),('DEC_MEAN','f8'),
             ('NCALSTARS','i4'),('QSLR_FLAG','i2')]

    out = np.ones(len(data)*62,dtype=dtype)
    out['FILENAME'][:] = ''
    out['EXPNUM'][:] = data['expnum'].repeat(62)
    out['CCDNUM'][:] = np.tile(np.arange(1,63),len(data))
    out['BAND'][:] = data['band'].repeat(62)
    out['MAG_ZERO'][out['BAND']=='i'] = 30.0
    out['MAG_ZERO'][out['BAND']=='z'] = 29.67
    out['SIGMA_MAG_ZERO'][:] = np.nan
    out['RA_MEAN'][:] = data['telra'].repeat(62)
    out['DEC_MEAN'][:] = data['teldec'].repeat(62)
    out['NCALSTARS'][:] = 0
    out['QSLR_FLAG'][:] = 0

    outfile = args.outfile
    fitsio.write(outfile,out,clobber=True)
