#!/usr/bin/env python
import sys
import os
from os.path import join, basename
import glob
from termcolor import colored as color
import logging
from multiprocessing import Pool
from multiprocessing import Process, Value, Lock

import yaml
import fitsio
import numpy as np

from const import BANDS,OBJECT_ID,UNIQUE_ID,BADMAG
from utils import bfields

SECTIONS = ['download','pixelize','zeropoint','match','catalog','nan']
OK  = color('OK','green')
FAIL = color('FAIL','red')
BAD = color('BAD','red')
   
def bad_values_str(values):
    opts = np.get_printoptions()
    np.set_printoptions(precision=4,threshold=25,edgeitems=4)
    ret = repr(values)
    np.set_printoptions(**opts)
    return ret

def print_running(idx,total,indent=0,step=10):
    """
    'idx' is indexed starting from 1.
    """
    i = idx
    idt = indent*' '
    if i == 1 or i % step == 0 or i == total: 
        sys.stdout.write("\r%s(%i/%i): "%(idt,i,total))
        sys.stdout.flush()

def print_wheel(idx,total,indent=0,step=10):
    i = idx
    idt = indent*' '
    if 0 <= i % 10 <= 1: sys.stdout.write("\r%s(*): "%(idt))
    if 2 <= i % 10 <= 3: sys.stdout.write("\r%s(-): "%(idt))
    if 4 <= i % 10 <= 5: sys.stdout.write("\r%s(\): "%(idt))
    if 6 <= i % 10 <= 7: sys.stdout.write("\r%s(|): "%(idt))
    if 8 <= i % 10 <= 9: sys.stdout.write("\r%s(/): "%(idt))
    sys.stdout.flush()

def print_points(idx,total,indent=0,step=10):
    i = idx
    idt = indent*' '
    n = i % 10
    sys.stdout.write("\r%s(%-9s): "%(idt,n*'.'))
    sys.stdout.flush()

def dir_exists(dirname):
    if not os.path.exists(dirname):
        msg = "No directory:"
        print color(msg,'red')
        print dirname
        return False
    return True

# These pool functions should be moved to utils
global counter
counter = Value('i',0)

def init_counter(args):
    """ Initialize the counter for later use """
    global counter
    counter = args

def run_pool(func, args, **kwargs):
    """ Initialize the pool with a shared counter """
    global counter
    counter = Value('i',0)
    pool = Pool(processes=30,initializer=init_counter, initargs=(counter,),**kwargs)
    return pool.map(func,args)

def check_files(explist,files):
    exp = explist
    nfiles = len(files)
    unit = np.array([basename(f).split('_')[0] for f in files])
    missing = ~np.in1d(exp['unitname'],unit)
    extra = ~np.in1d(unit,exp['unitname'])
    
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

    print_running(nfiles,nfiles,indent=4)
    if missing.sum():
        msg = "Missing %s file(s):"%missing.sum()
        print color(msg,'red')
        print 4*' '+'%s'%exp['unitname'][missing]
    elif extra.sum():
        msg = "Extra %s file(s):"%extra.sum()
        print color(msg,'red')
        print 4*' '+'%s'%unit[extra]
    else: print OK

def count_objects(args):
    global counter
    with counter.get_lock(): counter.value += 1

    f,nfiles,band = args
    print_running(counter.value,nfiles,indent=4,step=1)

    fits = fitsio.FITS(f,'r')
    num = fits[1].get_nrows()
    hpx = (fits[1].read_header().get('HPX') is not None)
    fits.close()

    # Ignore low numbers of objects in pixelized files
    if (not hpx) and (num < 5e3 if band=='Y' else num < 1e4):
        msg = "Only %i objects in %s"%(num,f)
        print color(msg,'yellow')

    return num

def check_columns(args,columns=None,select=None,msg=None):
    """ Abstract base function for checking a column"""
    global counter
    f,nfiles,band = args

    with counter.get_lock(): counter.value += 1
    print_running(counter.value,nfiles,indent=4,step=1)

    try: 
        data = fitsio.read(f,columns=columns)
    except ValueError as e:
        msg = "Couldn't read %(columns)s from %(filename)s"
        msg = msg%dict(columns=columns,filename=f)
        print color(msg,'red')
        return True
    
    sel = select(data)
    if np.any(sel):
        #if not msg: msg = "Bad %(columns)s value in %(filename)s"
        #msg = msg%dict(columns=columns,filename=f)
        if msg: 
            print msg%dict(columns=columns,filename=f)
        if not isinstance(sel,(bool,type(np.bool))):
            msg = "Bad %(columns)s value in %(filename)s"
            print msg%dict(columns=columns,filename=f)
            print 4*' '+bad_values_str(data[sel])
        return True

    return False

def check_ra(args):
    kwargs = dict(columns='RA', select=lambda x: (x<0) | (x>360))
    return check_columns(args,**kwargs)

def check_flux(args):
    kwargs = dict(columns='FLUX_AUTO', select=lambda x: (x<0) | (x>1e9))
    return check_columns(args,**kwargs)

def check_ccdnum(args):
    kwargs = dict(columns='CCDNUM', select=lambda x: (x<1) | (x>62))
    return check_columns(args,**kwargs)

def check_zeropoint(args):
    def select(x):
        olderr = np.seterr(invalid='ignore')
        sel = (~np.isfinite(x)) | (x < 5)
        sel |= (x>35) & ~((x==BADMAG) | (x==999))
        np.seterr(**olderr)
        return sel
    kwargs = dict(columns='MAG_ZERO',select=select)
    return check_columns(args,**kwargs)

def check_magnitude(args):
    def select(x):
        olderr = np.seterr(invalid='ignore')
        sel = (~np.isfinite(x)) | (x < 5) 
        sel |= ((x > 35) & (x != BADMAG))
        np.seterr(**olderr)
        return sel
    kwargs = dict(columns='MAG_AUTO',select=select)
    return check_columns(args,**kwargs)

def check_objid(args):
    kwargs = dict(columns=OBJECT_ID,select=lambda x: (~np.isfinite(x)) | (x < 0))
    return check_columns(args,**kwargs)

def check_match(args):
    def select(x):
        nobjs = float(len(x))
        frac = (x > 0).sum()/nobjs
        bad = (frac < 0.1) and (nobjs > 1e4)
        if bad:
            msg = 'Match fraction = %.2f;'%frac
            print color(msg,'yellow'),
        return bad
    msg = color("Poor match in %(filename)s",'yellow')
    kwargs = dict(columns=bfields('NEPOCHS',band),select=select,msg=msg)
    return check_columns(args,**kwargs)

def check_nan(args):
    kwargs = dict(select = lambda x: bool(np.any(np.isnan(x.view('>f4')))) )

    kwargs.update(columns=bfields(['WAVG_MAGRMS_AUTO','WAVG_MAGRMS_PSF'],BANDS),
                  msg="Bad MAGRMS value in %(filename)s")
    ret = check_columns(args,**kwargs)

    kwargs.update(columns=bfields(['WAVG_SPREADRMS_MODEL'],BANDS),
                  msg="Bad SPREADRMS value in %(filename)s")
    ret |= check_columns(args,**kwargs)

    return ret

def check_header_keys(args)
    global counter
    with counter.get_lock(): counter.value += 1
    keys = np.atleast_1d(['RA','DEC'])
    f,nfiles,band = args
    print_running(counter.value,nfiles,indent=4,step=1)

    fits = fitsio.FITS(f,'r')
    names = fits[1].get_colnames()
    fits.close()
    return keys[~np.in1d(keys,names)]

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

                RAWCOUNT = 0

                files = sorted(glob.glob(dirname+'/*.fits'))
                nfiles = len(files)
                args = [(f,nfiles,band) for f in files]

                print 2*' '+"Files:"
                exp = np.recfromcsv(config['explist'])
                exp = exp[exp['band']==band]
                check_files(exp,files)

                print 2*' '+"Objects:"
                out = run_pool(count_objects,args)
                RAWCOUNT = np.sum(out)
                print RAWCOUNT, OK
                
                print 2*' '+"CCDNUM:"
                out = run_pool(check_ccdnum,args)
                print FAIL if np.any(out) else OK

                print 2*' '+"RA:"
                out = run_pool(check_ra,args)
                print FAIL if np.any(out) else OK
                
                print 2*' '+"Fluxes:"
                out = run_pool(check_flux,args)
                print FAIL if np.any(out) else OK

                # Check that the columns and dtypes are the same in
                # all files...

            ##############################
            if section == 'pixelize':
                dirname = join(config['hpxdir'],band)
                if not dir_exists(dirname): continue

                HPXCOUNT = 0

                files = sorted(glob.glob(dirname+'/*.fits'))
                nfiles = len(files)
                args = [(f,nfiles,band) for f in files]

                print 2*' '+"Objects:"
                out = run_pool(count_objects,args)
                HPXCOUNT = np.sum(out)
                print HPXCOUNT

                if RAWCOUNT is not None:
                    print 4*' '+"raw: %i | hpx: %i"%(RAWCOUNT,HPXCOUNT),
                    print FAIL if RAWCOUNT != HPXCOUNT else OK

                print 2*' '+"RA:"
                out = run_pool(check_ra,args)
                print FAIL if np.any(out) else OK

            ##############################
            if section == 'zeropoint':
                dirname = join(config['hpxdir'],band)
                if not dir_exists(dirname): continue

                files = sorted(glob.glob(dirname+'/*.fits'))
                nfiles = len(files)
                args = [(f,nfiles,band) for f in files]

                print 2*' '+"Zeropoints:"
                out = run_pool(check_zeropoint,args)
                print FAIL if np.any(out) else OK

                print 2*' '+"Magnitudes:"
                out = run_pool(check_magnitude,args)
                print FAIL if np.any(out) else OK
            
            ##############################
            if section == 'match':
                dirname = join(config['hpxdir'],band)
                if not dir_exists(dirname): continue
                
                files = sorted(glob.glob(dirname+'/*.fits'))
                nfiles = len(files)
                args = [(f,nfiles,band) for f in files]

                print "  Object IDs:"
                out = run_pool(check_objid,args)
                print FAIL if np.any(out) else OK
                    
            ##############################
            if section == 'catalog':
                catdir = config['catdir']
                if not dir_exists(catdir): continue
                keydir = config['keydir']
                if not dir_exists(keydir): continue
                
                files = sorted(glob.glob(catdir+'/*.fits'))
                nfiles = len(files)
                args = [(f,nfiles,band) for f in files]

                CATCOUNT = 0

                print 2*' '+"Objects:"
                out = run_pool(count_objects,args)
                CATCOUNT = np.sum(out)
                print CATCOUNT, OK
                 
                print 2*' '+"RA:"
                out = run_pool(check_ra,args)
                print FAIL if np.any(out) else OK

                print 2*' '+"Match Fraction:"
                out = run_pool(check_match,args)
                print FAIL if np.any(out) else OK

    if 'nan' in sections:
        catdir = config['catdir']
        if not dir_exists(catdir): pass
        keydir = config['keydir']
        if not dir_exists(keydir): pass
        
        files = sorted(glob.glob(catdir+'/*.fits'))
        nfiles = len(files)
        args = [(f,nfiles,band) for f in files]

        print 2*' '+"NaN:"
        out = run_pool(check_nan,args)
        print FAIL if np.any(out) else OK

    print "\nDone."
                
