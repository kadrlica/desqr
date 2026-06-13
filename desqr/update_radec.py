#!/usr/bin/env python
"""
Add/update RA, Dec columns to catalog using Phase3 solutions, no color terms.

Run with:
cmd  = 'ls %s/%s/*.fits | xargs -n 100 update_radec.py -f'%(dirname, band)
"""


### Notes to the user community (= Alex):
### To run this code you'll need to have the pixmappy package installed.  It includes some data files
### that are accessed to define tree rings and such.  A copy currently lives in the ~garyb/pixmappy
### directory on the FNAL machines.
###  
### Newer WCS information is stored in files
### expoData3.hdf5
### epochShifts.fits
###  
### that are in ~garyb/PM directory.  You can alter the code below to point to someplace for them.
### Also you can see here templates for where the input catalogs should be found (`catfile`)
### and for the name that you want the output files to be placed.
###  
### This script will do whatever list of expnums you place on its command line.  It will
### print a message to stdout if it can't find an input catalog, and also print out some progress
### information when it succeeds.  The outputs will be new copies of the catalog files, with an
### 'radec' column added giving each object's estimate RA, Dec in degrees.
###  
### Any existing file of the same name as desired output will be overwritten.
###  
### # Format with band,expnum:
### #catfile = '/data/delve01.b/data/chinyi/dr3qr_skims/cat/delve_ncsa/{0:s}/D{1:08d}_{0:s}_cat.fits'
###  
### # Format with expnum:
### outfile = 'test{:08d}_cat.fits'

import os, sys
import argparse

import numpy as np
from scipy.interpolate import interp1d
from astropy.table import Table
from astropy.time import Time
import fitsio

import pixmappy as pm
from desqr import utils 
from desqr.logger import logger

DEGREE = np.pi / 180.   # Degree, in radians
DEFAULT_COLOR = 0.61    # Color used if no g-i available

# Read in table of exposure data
#exposureDataPath = '~garyb/PM/expoData3.hdf5'
#exposureDataPath = '~garyb/PM/expoData_v3.hdf5'
#exposureDataPath = '~garyb/PM/exposureData3.hdf5'
exposureDataPath = '/data/delve01.b/data/garyb/Results/delveExposures.hdf5'
shiftDataPath = '~garyb/PM/epochShifts.fits'

exposureData = Table.read(exposureDataPath)

# Read in all the DECam map components
despmc = pm.DESMaps()

# Use the Y1 polynomials until this day:
startY6 = 20180619    # Begin using Y6 polynomials (old epoch 20180829) at new epoch 20180619

class ShiftFinder:
    def __init__(self, path=shiftDataPath):
        '''Class that will generate Linear PixelMap appropriate to
        a given DECam CCD on a given night of observing. Input
        path to the shifts file created from DES+DELVE data on
        construction.  Then call this object with the (MJD or NITE)
        of the observation and the CCDNUM.

        The returned Linear instance will have a name like
        `20121101/S14` with the indicated nite being
        the night of observations.
        
        Also returns the epoch value and the polynomial set to use.'''
        
        self.tab = Table.read(path)
        self.DECAM_MJD0 = self.tab.meta['MJD0']
        self.starts = self.tab['startDay']
    def __call__(self, mjd, ccdnum):
        if mjd > 20000000:
            # This is a NITE, not an MJD
            nite = int(mjd)
            # Convert NITE notation to the DECam day count
            day =  np.floor(Time('{:04d}-{:02d}-{:02d}'.format(nite//10000, (nite//100)%100, nite%100)).mjd - self.DECAM_MJD0)
        elif mjd<10000:
            # This is a day count
            day = mjd
            ymd = Time(day+self.DECAM_MJD0, format='mjd').ymdhms
            nite =  int(ymd[0]*10000 + ymd[1]*100 + ymd[2])            
        else:
            day = np.floor(mjd - self.DECAM_MJD0 - 0.7)  # Obs of a given NITE are ~0.9-1.5 days past 00:00 UT of the NITE
            ymd = Time(day+self.DECAM_MJD0, format='mjd').ymdhms
            nite =  int(ymd[0]*10000 + ymd[1]*100 + ymd[2])
        index = np.searchsorted(self.starts, day, side='right') - 1
        if index<0:
            raise(ValueError('Requested MJD '+str(mjd) + ' is before epochs begin'))
        startDay = self.starts[index]
        epoch = self.tab['epoch'][index]
        detpos = pm.ccdnum2detpos[ccdnum]
        name = '{:08d}/{:s}'.format(nite,detpos)

        dt = day - self.tab['d0'][index]
        affine = self.tab['affine'][index,ccdnum,:,:] @ np.array([1,dt])
            
        # PixMappy class wants flattened 2x3
        mm = {'Type':'Linear', 'name':name, 'Coefficients':affine.T.flatten()}
        poly = 'Y1'
        if nite >=  20180619:  # start of Y6
            poly = 'Y6'
        
        return mm, epoch, poly
    
sf = ShiftFinder()

class LCg(pm.PixelMap):

    @staticmethod
    def type():
        return 'LCg'

    def __init__(self, name, **kwargs):

        '''PixelMap that makes the recomputed lateral color correction for
        g band (9 Feb 2026)'''
        # These are tabulations of the corrections.
        # Locations of reference points for new shifts
        super(LCg,self).__init__(name)

        self.includeOld = True   # Combine new tweaks with old polynomial
        rmid = np.array([0.14265489, 0.24708551, 0.31898602, 0.37742935, 0.42796466,
           0.47313273, 0.51434951, 0.5525    , 0.58818117, 0.62181823,
           0.65372682, 0.6841488 , 0.71327443, 0.74125653, 0.76822007,
           0.79426879, 0.81948993, 0.84395769, 0.8677358 , 0.88515029,
           0.89657206, 0.90785014, 0.91898982, 0.92999608, 0.9408736 ,
           0.95162679, 0.96225982, 0.97277663, 0.98318096, 0.99347632,
           1.00366609, 1.01375344, 1.0237414 , 1.03363285, 1.04343054,
           1.05313708, 1.06275497, 1.07228659, 1.08173424, 1.09110008,
           1.1003862 ])
        # Values of new radial shifts (converted to degrees)
        dr_new = np.array([0.00029848, 0.00045575, 0.00099359, 0.00064308, 0.0008012 ,
           0.00122591, 0.00104458, 0.00095653, 0.00086855, 0.00112881,
           0.00113355, 0.00101216, 0.0009162 , 0.0010048 , 0.00117232,
           0.00154478, 0.00178469, 0.00159779, 0.00146685, 0.00177162,
           0.00207413, 0.00256249, 0.00298204, 0.00349814, 0.00386177,
           0.00444939, 0.00509988, 0.00570929, 0.00680259, 0.00760426,
           0.00807254, 0.00848433, 0.00938779, 0.01006415, 0.01018218,
           0.01041629, 0.01055002, 0.00978171, 0.00851019, 0.00698859,
           0.00297378]) / 3600.
        # Coefficients of old correction polynomial, (r, r^3, r^5)
        self.old_poly = np.array([-1.83206613e-07,  1.54312254e-10, -7.02153565e-14])
    
        # Build LUT that uses r^2 as x and dr/r as y
        x = np.concatenate([[0,], rmid*rmid])
        y = np.concatenate([[0,], dr_new / rmid])
        self.rfunc = interp1d(x,y, kind='linear', bounds_error=False,
                            fill_value=(y[0], y[-1]))

        # Now color factors
        cmid = np.arange(0.1,4,0.2)
        y = np.array([-1.78058871e-01, -7.88957863e-02, -2.67864600e-02, -6.70515556e-03,
       -8.94105499e-04,  1.22985517e-04, -7.61278786e-04,  1.03024975e-03,
        3.97318282e-03,  1.11527501e-02,  2.63817935e-02,  4.91845151e-02,
        1.06382047e-01,  2.08982573e-01,  3.65628275e-01,  5.62582647e-01,
        7.90394250e-01,  1.01409756e+00,  1.20539754e+00,  1.33343773e+00])

        self.cfunc = interp1d(cmid, y, kind='linear', bounds_error=False,
                            fill_value=(y[0], y[-1]))
        self.oldRef = 0.61

    def __call__(self, u, v, c):
        '''Apply tweaks to DECam pixel positions
        :param x,y: uv plane coordinate arrays about optic axis, degrees
        :param c:   g-i color of source(s).
        :returns: u, v tweaked positions.
        '''
        if np.array(c).ndim>0:
            uv = np.stack([u,v],axis=-1)
            rsq = np.sum(uv*uv, axis=-1)
            dr = self.cfunc(c) * self.rfunc(rsq)
            if self.includeOld:
                dr += (c-self.oldRef)*(self.old_poly[0] + self.old_poly[1]*rsq + self.old_poly[2]*rsq*rsq)
            uv = uv* (1 + dr)[:,np.newaxis]
            return uv[:,0], uv[:,1]
        else:
            # work properly for single scalar input
            rsq = u*u+v*v
            dr = self.cfunc(c) * self.rfunc(rsq)
            if self.includeOld:
                dr += (c-self.oldRef)*(self.old_poly[0] + self.old_poly[1]*rsq + self.old_poly[2]*rsq*rsq)
            return u*(1 + dr),v*(1+dr)



class DCR(pm.PixelMap):
    @staticmethod
    def type():
        return 'DCR'

    def __init__(self, name, **kwargs):
        '''Tranformation of uv coordinates with new DCRs.
        Arguments:
        `band`: 'g', 'r', 'i', or 'z'
        `airmass`:   sec(z)
        `parallactic`:  parallactic angle, radians'''

        super(DCR,self).__init__(name)

        if 'band' not in kwargs or 'airmass' not in kwargs or 'parallactic' not in kwargs:
            raise ValueError('Missing arguments for DCR PixelMap')
        
        airmass = kwargs['airmass']
        parallactic = kwargs['parallactic']
        band = kwargs['band']
        
        # These are tabulations of the corrections.
        # Locations of reference points for new shifts
        self.uv0 = np.sqrt(airmass*airmass-1) * np.array([np.sin(parallactic), np.cos(parallactic)])
        
        self.includeOld = True  # False would omit the older corrections
        
        # DCR corrections derived in notebook
        dcr_c = {'g': np.array([ 0.01406037,  0.01076252,  0.00730764,  0.00383608,  0.00110561,
        -0.00072597, -0.00309405, -0.0052118 , -0.00696111, -0.00911936,
        -0.01182188, -0.01441902, -0.01823026, -0.02295002, -0.02804703,
        -0.03301929, -0.03774145, -0.04195041, -0.04537747, -0.04815899]),
             'r': np.array([ 0.00135103,  0.00053695, -0.00041605, -0.0014421 , -0.00207438,
        -0.00233447, -0.0025306 , -0.00246591, -0.00224002, -0.00229327,
        -0.0022802 , -0.00200004, -0.00191176, -0.00185283, -0.00186362,
        -0.00173601, -0.00136048, -0.00096188, -0.0007111 , -0.00054123]),
             'i': np.array([-2.22035193e-05, -4.58244048e-04, -8.31121853e-04, -1.32070648e-03,
        -1.63553360e-03, -1.75022417e-03, -1.85652269e-03, -1.73563109e-03,
        -1.53945681e-03, -1.41544631e-03, -1.22796424e-03, -8.43919448e-04,
        -5.94254775e-04, -3.73200138e-04, -1.93397461e-04,  1.17475431e-04,
         4.56881245e-04,  7.55475017e-04,  9.94815993e-04,  1.13631860e-03]),
             'z': np.array([ 6.97339148e-04,  5.02633413e-04,  2.36852346e-04,  3.78280505e-05,
        -1.06232380e-04, -1.86430189e-04, -2.33292462e-04, -2.52334759e-04,
        -2.70898005e-04, -3.26641379e-04, -3.82809476e-04, -4.17721174e-04,
        -4.78793049e-04, -5.36618191e-04, -5.99623980e-04, -6.42545220e-04,
        -6.60334158e-04, -6.35937805e-04, -6.16684307e-04, -5.87527546e-04])}

        # Now color factors
        cmid = np.arange(0.1,4,0.2)
        # Shift cfunc to be zero at g-i=1.1
        self.cfunc = interp1d(cmid, (dcr_c[band]-dcr_c[band][5]) / 3600.,
                              kind='linear', bounds_error=False, fill_value='extrapolate')
        # Previous solutions
        dcrConstant = {'g':45.0, 'r':8.4, 'i':3.2, 'z':1.4, 'Y':1.1}  
        ## coeffs = {'g':39.2, 'r':8.4, 'i':3.2, 'z':1.4}
        self.oldRef = 0.61
        # Convert mas/mag to degree/mag
        self.oldFactor = dcrConstant[band] / 3600. / 1000.  

    def __call__(self, u, v, c):
        '''Apply DCR correction to (u,v) coordinates (in degrees)'''
        
        if np.array(c).ndim>0:
            duv = self.cfunc(c)[:,np.newaxis] * self.uv0
            if self.includeOld:
                duv += ((c-self.oldRef)*self.oldFactor)[:,np.newaxis] * self.uv0
            return u+duv[:,0], v+duv[:,1]
        else:
            # Work properly for scalar inputs
            duv = self.cfunc(c) * self.uv0
            if self.includeOld:
                duv += self.oldFactor * (c-self.oldRef) * self.uv0
            return u+duv[0], v+duv[1]
 
# Add these to PixelMapCollection atoms
pm.PixelMapCollection.addAtom(LCg)
pm.PixelMapCollection.addAtom(DCR)


def cubic(uv, coeffs):
    '''Evaluate 2d cubic polynomial transformation'''
    x = uv[:,0]
    y = uv[:,1]
    A = np.stack( [np.ones_like(x),
                   x,y,
                   x**2, x*y, y**2,
                   x**3, x**2*y, x*y**2, y**3],
                  axis=-1)
    return A @ coeffs



## Douglas's color transformations
def rz2gi(rz):
    return np.where(rz<0.8, 0.182 + 2.125*rz, 1.182+0.863*rz)  # -0.7 < rz < 2.9

def ri2gi(ri):
    return np.where(ri<0.5, 0.127 + 3.458*ri, 1.196+1.226*ri)

def gr2gi(gr):
    return np.where(gr<1.2, -0.046 + 1.399*gr, -3.631+4.228*gr) # gr>1.2 not great.

def bprp2gi(bprp):
    # My fit to his g,i vs G,bprp
    # at https://des.ncsa.illinois.edu/releases/dr2/dr2-docs/dr2-interpolations
    return np.where(bprp<2., -0.447+1.309*bprp, 0.237+0.967*bprp) 


def run_filename(filename):
    '''Map objects in the input table from their pixel coordiates to ICRS RA,Dec
    using Phase 3 WCS models, ignoring color corrections.  Table is updated.'''

    expnum = int(os.path.basename(filename).split('_')[0].strip('D'))

    # Find the exposure in data table
    iexp = np.where(exposureData['expnum']==expnum)[0]
    if len(iexp)!=1:
        print('Did not find exposure {:d} in exposureData'.format(expnum))
        return 
    iexp = iexp[0]
    
    print('Running',expnum)
    band = exposureData['band'][iexp]
    nite = exposureData['nite'][iexp]
    pole = exposureData['pole'][iexp]
    # Create the projection
    expoProject = pm.Gnomonic(*pole)

    mjdmid = exposureData['mjdmid'][iexp]
        
    if not os.path.exists(filename):
        print('Did not find catalog for exposure {:08d}'.format(expnum))
        return

    # Use fitsio to close file
    columns = ['CCDNUM', 'XWIN_IMAGE', 'YWIN_IMAGE']
    cat = fitsio.read(filename, columns=columns)

    if False:
        # Skip DCR for now
        airmass = exposureData['airmass'][iexp]
        parallactic = exposureData['parallactic'][iexp]

        # Define the dcr PixelMap for this exposure
        dcrname = 'D{:07d}/dcr'.format(expnum)
        dcrmap = {'Type':'DCR','band':band,'airmass':airmass, 'parallactic':parallactic}
        despmc.update({dcrname:dcrmap})  

    # Set up the cubic exposure transformation
    coeffs = exposureData['coeffs'][iexp]

    # Now build PixelMap from pixel xy to exposure uv system,
    ccdnums = np.unique(cat['CCDNUM'])

    # Empty array to hold RA, Dec computed:
    rd = np.zeros((len(cat),2), dtype=float)
    for ccdnum in ccdnums:
        detpos = pm.ccdnum2detpos[ccdnum]
        use = np.where(cat['CCDNUM']==ccdnum)[0]

        # Aquire ccdshift and build exposure solution,
        # including new shift 
        shift, epoch, poly = sf(mjdmid, ccdnum)
        if band=='g':
            elements = ['{:s}/{:s}/rings'.format(band,detpos),
                        '{:s}/{:s}/lowedge'.format(band,detpos),
                        '{:s}/{:s}/highedge'.format(band,detpos),
                        '{:s}/{:s}/{:s}/poly'.format(band,poly,detpos),
                        shift['name'],
                        'g/color2']  ###,
                        ###    dcrname]
            if not despmc.hasMap('g/color2'):
                despmc.update({'g/color2':{'Type':'LCg'}})
        else:
            elements = ['{:s}/{:s}/rings'.format(band,detpos),
                        '{:s}/{:s}/lowedge'.format(band,detpos),
                        '{:s}/{:s}/highedge'.format(band,detpos),
                        '{:s}/{:s}/{:s}/poly'.format(band,poly,detpos),
                        shift['name'],
                        '{:s}/color'.format(band)] ###,
                        ###dcrname]
        if not despmc.hasMap(shift['name']):
            # Add the shift map to the collection if it's new
            name = shift.pop('name')
            despmc.update({name:shift})

        # Assign exposure to a DECam solution epoch
        mapname = 'D{:07d}/{:s}'.format(expnum,detpos)
        despmc.update({mapname:{'Type':'Composite','Elements':elements}})
            
        # Instantiate the map
        mm = despmc.getMap(mapname)
        uv = mm(cat['XWIN_IMAGE'][use], cat['YWIN_IMAGE'][use], np.ones(len(use),dtype=float)*DEFAULT_COLOR)
          ###gi[use])
        uv = np.stack(uv, axis=-1)
        # Apply cubic transformation
        uv += cubic(uv, coeffs)

        # Project to sky
        tmp = np.array(expoProject.toSky(*uv.T))
        rd[use] = tmp.T

        # End of chip loop.
    return rd

def update_radec(filename):
    """Calculate the new RA, Dec and insert into existing file."""
    radec = run_filename(filename)

    if radec is None:
        logger.warning("No RA,DEC returned; skipping")
        return
    
    radec = np.rec.fromarrays(radec.T, names=['RA','DEC'])
    utils.insert_columns(filename, radec, force=args.force)

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

    for i, filename in enumerate(args.filenames):
        radec = update_radec(filename)

        # Clear pixmappy caches once in a while
        if i % args.cache_clear == 0:
            logger.info("Clearing pixmappy cache...")
            despmc.clearCache()

    # Run with:
    #cmd  = 'ls %s/%s/*.fits | xargs -n 100 update_radec.py -f'%(dirname, band)
