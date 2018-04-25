#!/usr/bin/env python
import sys
import os
from os.path import join, basename
import glob
from termcolor import colored as color
import logging

import yaml
import fitsio
import numpy as np

from const import BANDS,OBJECT_ID,UNIQUE_ID
from utils import bfields

SECTIONS = ['download','pixelize','zeropoint','match','catalog']
OK  = color('OK','green')
BAD = color('BAD','red')

def print_running(idx,total,indent=0,step=10):
    """
    'idx' is indexed starting from 1.
    """
    i = idx
    idt = indent*' '
    if i == 1 or i % step == 0 or i == total: 
        sys.stdout.write("\r%s(%i/%i): "%(idt,i,total))
        sys.stdout.flush()

def dir_exists(dirname):
    if not os.path.exists(dirname):
        msg = "No directory:"
        print color(msg,'red')
        print dirname
        return False
    return True
    

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config')
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('-b','--band',dest='bands',default=None,
                        action='append',choices=BANDS)
    parser.add_argument('-s','--section',dest='sections',default=None,
                        action='append',choices=SECTIONS)
    opts = parser.parse_args()

    config = yaml.load(open(opts.config))
    sections = opts.sections if opts.sections else SECTIONS

    bands  = opts.bands if opts.bands else config['bands']
    rawdir = config['rawdir']
    hpxdir = config['hpxdir']

    # First, check the download

    RAWCOUNT = None
    HPXCOUNT = None
    CATCOUNT = None

    for band in bands:
        print '\n'+80*'-'
        print "Checking band: %s"%(band)
        print 2*' '+str(sections)

        for section in sections:
            print "\nChecking '%s'..."%section

            ##############################
            if section == 'download':

                dirname = join(config['rawdir'],band)
                if not dir_exists(dirname): continue

                files = sorted(glob.glob(dirname+'/*.fits'))
                nfiles = len(files)

                print "  Filenames: "
                exp = np.recfromcsv(config['explist'])
                unit = np.array([basename(f).split('_')[0] for f in files])
                missing = ~np.in1d(exp['unitname'][exp['band']==band],unit)
                extra = ~np.in1d(unit,exp['unitname'][exp['band']==band])

                print_running(nfiles,nfiles,indent=4)
                if missing.sum():
                    msg = "Missing %s file(s):"%missing.sum()
                    print color(msg,'red')
                    print 4*' '+'%s'%exp['unitname'][missing]
                elif extra.sum():
                    msg = "  Found %s extra file(s):"%extra.sum()
                    print color(msg,'red')
                    print 4*' '+'%s'%unit[extra]
                else: print OK
                
                RAWCOUNT = 0
                print "  Objects:"
                for i,f in enumerate(files):
                    print_running(i+1,nfiles,indent=4,step=1)
                    fits = fitsio.FITS(f,'r')
                    num = fits[1].get_nrows()
                    ra = fits[1].read(columns='RA')
                    RAWCOUNT += num
                    bad = False
                    if num < 1e4:
                        msg = "Found %i objects in %s"%(num,f)
                        print color(msg,'yellow')
                        bad = True
                    elif np.any( (ra<0) | (ra>360) ):
                        msg = "Bad RA value in %s"%f
                        print color(msg,'red')
                        bad = True
                    fits.close()
                if not bad: print OK
                print 4*' '+"Number of Objects: %i"%RAWCOUNT

                print 2*' '+"Fluxes:"
                for i,f in enumerate(files):
                    print_running(i+1,nfiles,indent=4,step=1)
                    fits = fitsio.FITS(f,'r')

                    auto = fits[1].read(columns='FLUX_AUTO')
                    sel = ((auto < 0) | (auto > 1e9))
                    if np.any(sel):
                        msg = "Bad FLUX_AUTO value in %s"%f
                        print color(msg,'red')
                        print 4*' '+str(auto[sel])
                        bad = True
                if not bad: print OK

            ##############################
            if section == 'pixelize':
                dirname = join(config['hpxdir'],band)
                if not dir_exists(dirname): continue

                files = sorted(glob.glob(dirname+'/*.fits'))
                nfiles = len(files)

                HPXCOUNT = 0
                print 2*' '+"Objects:"
                for i,f in enumerate(files):
                    print_running(i+1,nfiles,indent=4,step=1)
                    fits = fitsio.FITS(f,'r')

                    num = fits[1].get_nrows()
                    HPXCOUNT += num

                    ra = fits[1].read(columns='RA')
                    sel = np.isnan(ra) | (ra<0) | (ra>360)

                    bad = False
                    if np.any( sel ):
                        msg = "Bad RA value in %s"%f
                        print color(msg,'red')
                        print 4*' '+str(ra[sel])
                        bad = True
                    fits.close()
                if not bad: print OK
                print 4*' '+"Number of Objects: %i"%HPXCOUNT

                print 4*' '+"raw: %i | hpx: %i"%(RAWCOUNT,HPXCOUNT),
                if RAWCOUNT is not None and RAWCOUNT != HPXCOUNT:
                    print BAD
                else:
                    print OK

            ##############################
            if section == 'zeropoint':
                dirname = join(config['hpxdir'],band)
                if not dir_exists(dirname): continue

                files = sorted(glob.glob(dirname+'/*.fits'))
                nfiles = len(files)
                print "  Zeropoints:"
                for i,f in enumerate(files):
                    print_running(i+1,nfiles,indent=4,step=1)
                    fits = fitsio.FITS(f,'r')
                    zp = fits[1].read(columns='MAG_ZERO')
                    sel = np.isnan(zp) | (zp > 35) | (zp < 25)

                    bad = False
                    if np.any( sel ):
                        msg = "Bad zeropoint in %s"%f
                        print 4*' '+str(zp[sel])
                        print color(msg,'red')
                        bad = True
                    fits.close()
                if not bad: 
                    msg = 'OK'
                    print color(msg,'green')

                print "  Magnitudes:"
                for i,f in enumerate(files):
                    print_running(i+1,nfiles,indent=4,step=1)
                    fits = fitsio.FITS(f,'r')
                    auto = fits[1].read(columns='MAG_AUTO')
                    sel = np.isnan(auto) | (auto < 5) | ((auto > 30) & (auto != 99))
                    bad = False
                    if np.any( sel ):
                        msg = "Unusual MAG_AUTO in %s"%f
                        print color(msg,'yellow')
                        print 4*' '+str(auto[sel])
                        bad = True
                    fits.close()
                if not bad: print OK
            
            ##############################
            if section == 'match':
                dirname = join(config['hpxdir'],band)
                if not dir_exists(dirname): continue
                
                files = sorted(glob.glob(dirname+'/*.fits'))
                nfiles = len(files)

                print "  Object IDs:"
                bad = False
                for i,f in enumerate(files):
                    print_running(i,nfiles,indent=4,step=1)
                    objid = fitsio.read(f,columns=OBJECT_ID)
                    if np.any(np.isnan(objid) | (objid < 0)):
                        msg = "Bad %s in %s"%(OBJECT_ID,f)
                        print color(msg,'red')
                        bad = True
                if not bad: print OK
                    
            ##############################
            if section == 'catalog':
                catdir = config['catdir']
                if not dir_exists(catdir): continue
                keydir = config['keydir']
                if not dir_exists(keydir): continue
                
                files = sorted(glob.glob(catdir+'/*.fits'))
                nfiles = len(files)

                print "  Matching:"
                CATCOUNT = 0
                bad = False
                for i,f in enumerate(files):
                    print_running(i+1,nfiles,indent=4,step=1)

                    epochs = bfields('NEPOCHS',bands)
                    epoch  = bfields('NEPOCHS',band)
                    coords = ['RA','DEC']
                    #spread = bfields('SPREAD_MODEL',bands)
                    columns = epochs + coords
                    data = fitsio.read(f,columns=columns)

                    CATCOUNT += len(data)

                    sel = np.isnan(data['RA'])|np.isnan(data['DEC'])
                    if np.any(sel):
                        msg = "RA,DEC NaN in %s"%f
                        print color(msg,'red')
                        print 4*' '+str(data[['RA','DEC']][sel])
                        bad = True

                    nobs = float((data[epoch] > 0).sum())
                    for b,e in zip(bands,epochs):
                        if e == epoch: continue
                        frac = ((data[e] > 0) & (data[epoch]>0)).sum()/nobs
                        if frac < 0.5 and nobs > 1e4:
                            msg = "Poor %s-%s match in %s"%(band,b,f)
                            print color(msg,'yellow')
                            print 4*' '+'Match fraction = %.2f'%frac
                if not bad: print OK
                print 4*' '+"Number of Objects: %i"%CATCOUNT
