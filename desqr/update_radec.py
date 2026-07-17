#!/usr/bin/env python
"""
Add/update RA, Dec columns to catalog using Phase3 solutions, no color terms.

Run with:
cmd  = 'ls %s/%s/*.fits | xargs -n 100 update_radec.py -f'%(dirname, band)
"""


### Notes to the user community (= Alex):
### Updated for the new WCS solutions and pixmappy, 13 Jun 2026

### This script will do whatever list of expnums you place on its command line.  It will
### print a message to stdout if it can't find an input catalog, and also print out some progress
### information when it succeeds.  It produces an 'radec' array that is inserted into the RA, DEC
### column of the input (overwriting if force = True).

import os, sys
import argparse

import numpy as np
from astropy.table import Table
import fitsio

import pixmappy as pm
from desqr import utils 
from desqr.logger import logger

DEGREE = np.pi / 180.   # Degree, in radians
DEFAULT_COLOR = 0.61    # Color used if no g-i available

# Read in all the Delve map components
# Calibration file path: export CAL_PATH=/home/s1/garyb/PM
#exposures_file='/home/s1/garyb/PM/delveTest.hdf5'
#exposures_file='/data/delve01.b/data/garyb/Results/delveExposures.hdf5'
#exposures_file='/data/delve01.b/data/garyb/Results/delveExposures6.0.hdf5'
exposures_file='/data/delve01.b/data/garyb/Results/delveExposures6.1.hdf5'
dmc = pm.DelveMaps(exposures_file=exposures_file)
#dmc = None

def get_expnum(filename):
    return int(os.path.basename(filename).split('_')[0].strip('D'))

def get_pixmappy_radec(filename):
    """Map objects in the input table from their pixel coordiates to ICRS RA,Dec
    using Phase 3 WCS models, ignoring color corrections.  Table is updated."""

    expnum = get_expnum(filename)

    if dmc is None:
        logger.warning(f'Pixmappy map not found.')
        return

    if not os.path.exists(filename):
        logger.warning(f'Catalog not found for exposure: {expnum:08d}')
        return

    if expnum not in dmc.exptab['expnum']:
        logger.warning(f"Solution not found for exposure: {expnum}")
        return
    
    # Use fitsio to close file
    columns = ['CCDNUM', 'XWIN_IMAGE', 'YWIN_IMAGE']
    cat = fitsio.read(filename, columns=columns)
    
    # Empty array to hold RA, Dec computed:
    rd = np.zeros((len(cat),2), dtype=float)
    for ccdnum in np.unique(cat['CCDNUM']):
        detpos = pm.ccdnum2detpos[ccdnum]
        use = np.where(cat['CCDNUM']==ccdnum)[0]

        try:
            mm = dmc.getDelveWCS(expnum, detpos)
        except:
            ### Currently returning None if solution does not exist for this exposure
            return

        tmp = mm.toSky(cat['XWIN_IMAGE'][use], cat['YWIN_IMAGE'][use], np.ones(len(use), dtype=float)*DEFAULT_COLOR)
        rd[use] = np.array(tmp).T
        # End of chip loop.

    radec = np.rec.fromarrays(rd.T, names=['RA','DEC'])
    return radec

def get_original_radec(filename):
    """ Get the original RA,Dec from the SCAMP solution."""
    radec = fitsio.read(filename, columns=['ALPHAWIN_J2000', 'DELTAWIN_J2000'])
    radec.dtype.names = ['RA','DEC']
    return radec

def update_radec(filename, force=False):
    """Calculate the new RA, Dec and insert into existing file.

    Parameters
    ----------
    filename : the catalog fits file
    force : overwrite columns if they already exist

    Returns
    -------
    radec : the output columns
    """
    radec = get_pixmappy_radec(filename)

    if radec is None:
        logger.warning("No RA,DEC returned by pixmappy; using original from SCAMP...")
        radec = get_original_radec(filename)

    radec['RA'] = radec['RA'] % 360.0
    utils.insert_columns(filename, radec, force=force)

    return radec

if __name__=='__main__':
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('filenames', nargs='+', help='input filenames')
    parser.add_argument('-c','--cache-clear',default=100,type=int,
                        help='number of exposures to run before clearing pixmappy cache')
    parser.add_argument('-f','--force',action='store_true',
                        help='overwrite output columns if they exist')
    parser.add_argument('-v','--verbose',action='store_true',
                        help='output verbosity')
    args = parser.parse_args()

    if args.verbose: logger.setLevel(logger.DEBUG)

    tab = Table.read(exposures_file)
    
    for i, filename in enumerate(args.filenames):
        
        radec = update_radec(filename, force=args.force)

        # Clear pixmappy caches once in a while
        if (i % args.cache_clear == 0) and dmc:
            logger.info("Clearing pixmappy cache...")
            dmc.clearCache()
