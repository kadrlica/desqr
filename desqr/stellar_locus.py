#!/usr/bin/env python
""" Run stellar locus regression validation. """
import os, glob
import datetime
import yaml
import warnings
warnings.filterwarnings("ignore")

import matplotlib
if os.getenv('TERM').startswith('screen') or not os.getenv('DISPLAY'):
    matplotlib.use('Agg')
import pylab as plt

import numpy as np
import scipy.stats as scipyStats
from scipy.stats import median_abs_deviation
import fitsio
import healpy as hp
import pandas as pd

from desqr import utils
from desqr.utils import load_infiles, mkdir
from desqr.utils import blank, bfield, bfields
from desqr.classify import bdf_extended_class, wavg_extended_class
from desqr.extinction import ebv, extinction
from desqr.logger import logger
from desqr.const import BANDS
from desqr import plotting

def load_slr_hpxmap(filename, param='wperp'):
    """ Load healpix map from slr file. """
    fits = fitsio.FITS(filename)
    nside = fits[1].read_header()['NSIDE']
    data = fits[1].read()

    hpxmap = blank(nside)
    hpxmap[data['pix']] = data[param]
    hpxmap = utils.masked_array(hpxmap)
    return hpxmap

def plot_slr(filename,outfile=None,survey='delve',param='wperp'):
    logger.info("Reading %s..."%filename)
    
    try:
        hpxmap = hp.read_map(filename)
    except ValueError:
        hpxmap = load_slr_hpxmap(filename, param)
    hpxmap = utils.masked_array(hpxmap)
    hpxmap *= 1000 # mmag

    cbar_kwargs = dict()
    hpxmap_kwargs = dict(xsize=2000)
    hist_kwargs = dict()

    if param == 'wperp' or '_wperp_' in filename:    
        label = 'slr wperp (mmag)'
    elif param == 'slope' or '_slope_' in filename:
        label = 'slr slope'
    elif param == 'intercept' or '_intercept_' in filename:
        label = 'slr intercept'
    else: 
        label = None

    cbar_kwargs['label'] = label

    fig,axes,smap = plotting.plot_hpxmap_hist(hpxmap,survey,cbar_kwargs,hpxmap_kwargs,
                                              hist_kwargs)
    #axes[0].annotate('%s band'%band, (0.05,0.93), xycoords='axes fraction')
    axes[1].set_xlabel(label)

    utils.print_statistics(hpxmap, unit='mmag')

    if outfile is None: 
        outfile=os.path.basename(filename).split('.')[0]+'.png'
    logger.info("Writing %s..."%outfile)
    plt.savefig(outfile,bbox_inches='tight')

# Everything below this is copied directly from pipe_analysis/utils.py.
# Should we move all those functions here once pipe_analysis is rewritten?

def calcP1P2(mags, coeffs):
    # P1 =A′u+B′g+C′r+D′i+E′z+F′
    # P2′=Au+Bg+Cr+Di+Ez+F
    p1p2 = (
        float(coeffs[0]) * mags[0]
        + float(coeffs[1]) * mags[1]
        + float(coeffs[2]) * mags[2]
        + float(coeffs[3]) * mags[3]
        + float(coeffs[4]) * mags[4]
        + float(coeffs[5])
    )
    return p1p2

def p1CoeffsFromP2x0y0(p2Coeffs, x0, y0):
    """Compute Ivezic P1 coefficients using the P2 coeffs and origin (x0, y0)
    Reference: Ivezic et al. 2004 (2004AN....325..583I)
    theta = arctan(mP1), where mP1 is the slope of the equivalent straight
                         line (the P1 line) from the P2 coeffs in the (x, y)
                         coordinate system and x = c1 - c2, y = c2 - c3
    P1 = cos(theta)*c1 + ((sin(theta) - cos(theta))*c2 - sin(theta)*c3 + deltaP1
    P1 = 0 at x0, y0 ==> deltaP1 = -cos(theta)*x0 - sin(theta)*y0
    Parameters
    ----------
    p2Coeffs : `list` of `float`
       List of the four P2 coefficients from which, along with the origin point
       (``x0``, ``y0``), to compute/derive the associated P1 coefficients.
    x0, y0 : `float`
       Coordinates at which to set P1 = 0 (i.e. the P1/P2 axis origin).
    Returns
    -------
    p1Coeffs: `list` of `float`
       The four P1 coefficients.
    """
    mP1 = p2Coeffs[0] / p2Coeffs[2]
    cosTheta = np.cos(np.arctan(mP1))
    sinTheta = np.sin(np.arctan(mP1))
    deltaP1 = -cosTheta * x0 - sinTheta * y0
    p1Coeffs = [cosTheta, sinTheta - cosTheta, -sinTheta, deltaP1]

    return p1Coeffs


def p2p1CoeffsFromLinearFit(m, b, x0, y0):
    """Derive the Ivezic et al. 2004 P2 and P1 equations based on linear fit
    Where the linear fit is to the given region in color-color space.
    Reference: Ivezic et al. 2004 (2004AN....325..583I)
    For y = m*x + b fit, where x = c1 - c2 and y = c2 - c3,
    P2 = (-m*c1 + (m + 1)*c2 - c3 - b)/sqrt(m**2 + 1)
    P2norm = P2/sqrt[(m**2 + (m + 1)**2 + 1**2)/(m**2 + 1)]
    P1 = cos(theta)*x + sin(theta)*y + deltaP1, theta = arctan(m)
    P1 = cos(theta)*(c1 - c2) + sin(theta)*(c2 - c3) + deltaP1
    P1 = cos(theta)*c1 + ((sin(theta) - cos(theta))*c2 - sin(theta)*c3 + deltaP1
    P1 = 0 at x0, y0 ==> deltaP1 = -cos(theta)*x0 - sin(theta)*y0
    Parameters
    ----------
    m : `float`
       Slope of line to convert.
    b : `float`
       Intercept of line to convert.
    x0, y0 : `float`
       Coordinates at which to set P1 = 0.
    Returns
    -------
    result : `lsst.pipe.base.Struct`
       Result struct with components:
       - ``p2Coeffs`` : four P2 equation coefficents (`list` of `float`).
       - ``p1Coeffs`` : four P1 equation coefficents (`list` of `float`).
    """
    # Compute Ivezic P2 coefficients using the linear fit slope and intercept
    scaleFact = np.sqrt(m ** 2 + 1.0)
    p2Coeffs = [-m / scaleFact, (m + 1.0) / scaleFact, -1.0 / scaleFact, -b / scaleFact]
    p2Norm = 0.0
    for coeff in p2Coeffs[:-1]:  # Omit the constant normalization term
        p2Norm += coeff ** 2
    p2Norm = np.sqrt(p2Norm)
    p2Coeffs /= p2Norm

    # Compute Ivezic P1 coefficients equation using the linear fit slope and
    # point (x0, y0) as the origin
    p1Coeffs = p1CoeffsFromP2x0y0(p2Coeffs, x0, y0)

    return {"p2Coeffs":p2Coeffs, "p1Coeffs":p1Coeffs,}


def calcQuartileClippedStats(dataArray, nSigmaToClip=3.0):
    """Calculate the quartile-based clipped statistics of a data array.
    The difference between quartiles[2] and quartiles[0] is the interquartile
    distance.  0.74*interquartileDistance is an estimate of standard deviation
    so, in the case that ``dataArray`` has an approximately Gaussian
    distribution, this is equivalent to nSigma clipping.
    Parameters
    ----------
    dataArray : `list` or `numpy.ndarray` of `float`
        List or array containing the values for which the quartile-based
        clipped statistics are to be calculated.
    nSigmaToClip : `float`, optional
        Number of \"sigma\" outside of which to clip data when computing the
        statistics.
    Returns
    -------
    result : `lsst.pipe.base.Struct`
        The quartile-based clipped statistics with ``nSigmaToClip`` clipping.
        Atributes are:
        ``median``
            The median of the full ``dataArray`` (`float`).
        ``mean``
            The quartile-based clipped mean (`float`).
        ``stdDev``
            The quartile-based clipped standard deviation (`float`).
        ``rms``
            The quartile-based clipped root-mean-squared (`float`).
        ``clipValue``
            The value outside of which to clip the data before computing the
            statistics (`float`).
        ``goodArray``
            A boolean array indicating which data points in ``dataArray`` were
            used in the calculation of the statistics, where `False` indicates
            a clipped datapoint (`numpy.ndarray` of `bool`).
    """
    quartiles = np.percentile(dataArray, [25, 50, 75])
    assert len(quartiles) == 3
    median = quartiles[1]
    interQuartileDistance = quartiles[2] - quartiles[0]
    clipValue = nSigmaToClip * 0.74 * interQuartileDistance
    good = np.logical_not(np.abs(dataArray - median) > clipValue)
    quartileClippedMean = dataArray[good].mean()
    quartileClippedStdDev = dataArray[good].std()
    quartileClippedRms = np.sqrt(np.mean(dataArray[good] ** 2))

    return {
        "median":median,
        "mean":quartileClippedMean,
        "stdDev":quartileClippedStdDev,
        "rms":quartileClippedRms,
        "clipValue":clipValue,
        "goodArray":good,}

def stellarLocusResid(gmags,rmags,imags,gr,ri):    
    okfitcolors = (
        (gr < 1.1)
        & (gr > 0.3)
        & (np.abs(ri) < 1.0)
        & np.abs(gmags < 20)
        & np.abs(rmags < 20)
        & np.abs(imags < 20)
    )
    if okfitcolors.sum() == 0:
        return None
    
    # Eventually switch to using orthogonal regression instead of linear (as in pipe-analysis)?

    slope, intercept, r_value, p_value, std_err = scipyStats.linregress(
        gr[okfitcolors], ri[okfitcolors]
    )
    
    p2p1coeffs = p2p1CoeffsFromLinearFit(slope, intercept, 0.3, slope * 0.3 + intercept)
    p1coeffs = p2p1coeffs["p1Coeffs"].copy()
    # hack to put the zeros in for u, z coeffs
    p1coeffs.insert(0, 0.0)
    p1coeffs.insert(4, 0.0)
    p2coeffs = list(p2p1coeffs["p2Coeffs"].copy())
    p2coeffs.insert(0, 0.0)
    p2coeffs.insert(4, 0.0)
    umags = np.zeros(len(gmags))
    zmags = np.zeros(len(gmags))
    p1_fit = calcP1P2([umags, gmags, rmags, imags, zmags], p1coeffs)
    p2_fit = calcP1P2([umags, gmags, rmags, imags, zmags], p2coeffs)
    okp1_fit = (p1_fit < 0.6) & (p1_fit > -0.2)
    
    if okp1_fit.sum() < 10:
        return None
    
    # Do a second iteration, removing large (>3 sigma) outliers in p2:
    clippedStats = calcQuartileClippedStats(p2_fit[okp1_fit], 3.0)
    keep = np.abs(p2_fit) < clippedStats["clipValue"]

    if (okfitcolors & keep).sum() < 10: 
        return None 
    
    slope, intercept, r_value, p_value, std_err = scipyStats.linregress(
        gr[okfitcolors & keep], ri[okfitcolors & keep]
    )
    
    p2p1coeffs = p2p1CoeffsFromLinearFit(slope, intercept, 0.3, slope * 0.3 + intercept)
    p1coeffs = p2p1coeffs["p1Coeffs"].copy()
    # hack to put the zeros in for u, z coeffs
    p1coeffs.insert(0, 0.0)
    p1coeffs.insert(4, 0.0)
    p2coeffs = list(p2p1coeffs["p2Coeffs"].copy())
    p2coeffs.insert(0, 0.0)
    p2coeffs.insert(4, 0.0)
    p1_fit = calcP1P2([umags, gmags, rmags, imags, zmags], p1coeffs)
    p2_fit = calcP1P2([umags, gmags, rmags, imags, zmags], p2coeffs)
    okp1_fit = (p1_fit < 0.6) & (p1_fit > -0.2)

    return p1_fit[okp1_fit], p2_fit[okp1_fit], p1coeffs, p2coeffs, slope, intercept


def calcWperp(gmags,rmags,imags,gr,ri):
    # given a set of gri mags return the rms width of the stellar locus in mag
    try:
        p1,p2,p1c,p2c=stellarLocusResid(gmags,rmags,imags,gr,ri)
    except:
           return np.nan
    if len(p2) >= 10:
        #wperp=median_absolute_deviation(p2)
        wperp=median_abs_deviation(p2, scale='normal')
        return wperp
    else:
        #print("here2")
        return np.nan

def filter_catalog(catalog, mag_base="MAG_PSF"):
    # maglims set from https://cdcvs.fnal.gov/redmine/projects/des-y3/wiki/Y3A2_Empirical_Stellar_Locus
    # could try a brighter faint cut also want to use "WAVG_MAG_PSF_%s"
    
    magDict={}
    colorDict={}
    bands=["g","r","i"]
    colors = [['g', 'r'],
              ['r', 'i']]

    try:
        extclass = bdf_extended_class(catalog)
    except:
        extclass = wavg_extended_class(catalog['WAVG_SPREAD_MODEL_G'], catalog['WAVG_SPREADERR_MODEL_G'])

    initialSel =  extclass == 0
    initialSel &= (catalog[bfield(mag_base,"G")] > 16) & (catalog[bfield(mag_base,"G")] < 22) 
    initialSel &= (catalog[bfield(mag_base,"R")] > 16) & (catalog[bfield(mag_base,"R")] < 21) 
    initialSel &= (catalog[bfield(mag_base,"I")] > 16) & (catalog[bfield(mag_base,"I")] < 21)
    
    catalog = catalog[initialSel].copy()
    sfd_path = config['ebv']
    ebvval = ebv(catalog["RA"], catalog["DEC"], sfd_path)

    for band in bands:
        extval=extinction(ebvval,band.lower())
        magDict[band]=catalog[bfield(mag_base,band)] - extval
        
    for color in colors:
        colorDict[color[0]+color[1]]=magDict[color[0]]-magDict[color[1]]
    
    return catalog,magDict,colorDict,initialSel

def calcWAndLine(gmags,rmags,imags,gr,ri):
    # given a set of gri mags return the rms width of the stellar locus in mag
    try:
        p1,p2,p1c,p2c, slope, intercept = stellarLocusResid(gmags,rmags,imags,gr,ri)
    except:
           return np.nan,np.nan,np.nan, np.nan
    if len(p2) >= 10:
        wperp=median_abs_deviation(p2, scale='normal')
        return wperp, slope, intercept, len(p2)
    else:
        return np.nan,np.nan,np.nan, len(p2)
    
    
def apply_calcWAndLine(args):
    infile,bands,colors,nside = args
    logger.info(f"Loading {infile}...")
    data = load_infiles(infile, columns=COLUMNS)
    
    catpix=hp.ang2pix(nside,data["RA"],data["DEC"], lonlat=True)
    
    cat,magDict,colorDict,mask=filter_catalog(data, mag_base=mag_base)
    
    pix_unique = np.unique(catpix)
    ncols=5

    if len(cat) < 1:
        out_arr=np.ones((ncols,len(pix_unique)))
        out_arr[0,:]=pix_unique
        out_arr[1,:]=np.unique(catpix, return_counts=True)[1]
        out_arr[2:,:]=np.nan
        return out_arr
    
    wPerpVals = []
    slopeVals = []
    interceptVals = []
    catlenVals = []
    
    catpix=catpix[mask]
    
    for pix in pix_unique:
        sel=(catpix==pix)
        if sel.sum() > 10:
            wperp,slope,intercept,catlen=calcWAndLine(*[magDict[b][sel] for b in bands],
                                       *[colorDict[c[0]+c[1]][sel] for c in colors])
            wPerpVals.append(wperp)
            slopeVals.append(slope)
            interceptVals.append(intercept)
            catlenVals.append(catlen)
            
        else:
            wPerpVals.append(np.nan)
            slopeVals.append(np.nan)
            interceptVals.append(np.nan)
            catlenVals.append(sel.sum())
            
    out_arr=np.ones((ncols,len(pix_unique)))
    out_arr[0,:]=pix_unique
    out_arr[1,:]=catlenVals
    out_arr[2,:]=wPerpVals
    out_arr[3,:]=slopeVals
    out_arr[4,:]=interceptVals
    return out_arr

def calcWAndLine_files(infiles,bands,colors,nside,processes=1):
    """ Load multiple input files.
    
    Parameters:
    -----------
    infiles   : list of input fits files
    distmods : array of distance moduli to step through
    iso_func : interpolated isochrone
    processes : number of cores to execute on
    
    Returns:
    --------
    data : return numpy array 
    """
    infiles = np.atleast_1d(infiles)
    logger.debug("Loading {len(infiles)} files...")

    arglist = [ ((f, bands, colors, nside), ) for f in infiles]
    out = utils.multiproc(apply_calcWAndLine, arglist, processes=processes)
    
    out=np.hstack(out)
    #check for duplicates
    s = np.sort(out[0], axis=None)
    
    if (s[1:] == s[:-1]).sum() > 0:
        print((s[1:] == s[:-1]).sum())
        msg = "Duplicated healpix in index."
        raise Exception(msg)

    out = pd.DataFrame({
        'pix':       out[0].astype(int),
        'nobjects':  out[1].astype(int),
        'wperp':     out[2],
        'slope':     out[3],
        'intercept': out[4],
    })
    return out


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('configfile',help='configuration file')
    parser.add_argument('-t','--type', dest='mag_base',
                        choices=["best","wavg","mag_psf", "mag_bdf"], default="best",
                        help='type of mag',)
    parser.add_argument('-n','--nside', help='nside', default=128, type=int)
    parser.add_argument('--nproc', default=20, type=int)
    parser.add_argument('-o','--outbase',default='slr')
    parser.add_argument('-p','--pix',default=None,type=int,action='append')
    parser.add_argument('-v','--verbose',action='store_true')
    args = parser.parse_args()

    if args.verbose: logger.setLevel(logger.DEBUG)

    config = yaml.safe_load(open(args.configfile))
    NSIDE = config['nside']
    survey = config.get('survey')
    catbase = config['catbase']
    catdir = config['catdir']
    filebase = os.path.join(catdir,catbase)

    outdir = mkdir('release/slr')

    nside=args.nside
    logger.info(f"nside: {nside}")

    if args.mag_base=="best":
        mag_base="MAG_PSF"
    elif args.mag_base=="wavg":
        mag_base="WAVG_MAG_PSF"
    elif args.mag_base=="psf_mag":
        mag_base="PSF_MAG"
    elif args.mag_base=="bdf_mag":
        mag_base="BDF_MAG"
    logger.info(f"mag: {mag_base}")

    bands = ["g","r","i"]
    colors = [['g', 'r'],
              ['r', 'i']]

    COLUMNS = ["RA","DEC"]
    COLUMNS += bfields(mag_base, bands)
    COLUMNS += ['WAVG_SPREAD_MODEL_G', 'WAVG_SPREADERR_MODEL_G']

    if args.pix is not None:
        pixels = args.pix
    else:
        pixels = np.arange(hp.nside2npix(NSIDE))

    if len(pixels) == 0:
        msg = "Invalid pixel: %s"%args.pix
        raise Exception(msg)
     
    #filenames = sorted(glob.glob(f'{catdir}/*.fits'))
    filenames = [filebase%p for p in pixels]
    filenames = [f for f in filenames if os.path.exists(f)]
    # This LMC hpx failed
    filenames = [f for f in filenames if not f.endswith('_11759.fits')]
    logger.info(f"Running {len(filenames)} files...")
          
    out=calcWAndLine_files(filenames, bands, colors, nside, processes=args.nproc)

    outbase = args.outbase + f"_all_{args.mag_base}_n{nside}.fits.gz"
    filename = os.path.join(outdir, outbase)
    logger.info(f"Writing {filename}...")
    fitsio.write(filename, out.to_records(index=False),
                 header={'NSIDE':nside}, clobber=True)

    for param in ['wperp']:
        hpxmap = load_slr_hpxmap(filename, param)
        outbase = args.outbase + f"_{param}_{args.mag_base}_n{nside}.fits.gz"
        outfile = os.path.join(outdir, outbase)

        logger.info(f"Writing {outfile}...")
        hp.write_map(outfile,hpxmap,overwrite=True)

        print(f"Global SLR {param}:")
        utils.print_statistics(hpxmap*1000, q=[5,50,95])

        print("Plotting %s..."%outfile)
        pngfile = outfile.replace('.fits.gz','.png')
        plot_slr(outfile,pngfile,survey=survey)
        plt.ion()

