#!/usr/bin/env python
import os
import yaml
import subprocess
import time
import glob

import numpy as np
import healpy
from ugali.utils.shell import mkdir
from utils import found

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-s','--sleep',default=2,type=float)
    parser.add_argument('-q','--queue',default='condor')
    opts = parser.parse_args()
    
    force = '-f' if opts.force else ''

    config = yaml.load(open(opts.config))
    hpxdir = config['hpxdir']
    catdir = mkdir(config['catdir'])
    keydir = mkdir(config['keydir'])
    logdir = mkdir(os.path.join(catdir,'log'))

    for pix in np.arange(healpy.nside2npix(config['nside'])):
        infiles = glob.glob(hpxdir+'/*/*%05d*.fits'%pix)
        basename = 'hpx_%05d.fits'%pix
        outfile = os.path.join(catdir,'cat_'+basename)
        keyfile = os.path.join(keydir,'key_'+basename)
        logfile = os.path.join(logdir,basename.replace('.fits','.log'))

        if len(infiles) == 0: continue

        if os.path.exists(outfile) and not opts.force:
            found(outfile)
            continue

        minbands = config.get('minbands')
        minbands = '--min-bands %s'%minbands if minbands else ''
        bands = ' '.join(['-b %s'%b for b in config.get('bands',[])])
        params=(' '.join(infiles),outfile,keyfile,bands,minbands,force)
        cmd = 'catalog.py -v %s -o %s -k %s %s %s %s'%params

        if opts.queue == 'local':
            print cmd
            submit = cmd
        else:
            submit = 'csub -o %s %s'%(logfile,cmd)
        subprocess.call(submit,shell=True)
        if opts.queue != 'local': time.sleep(opts.sleep)
