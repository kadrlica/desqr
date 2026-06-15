#!/usr/bin/env python
"""
Assemble multi-band unique catalog.
"""
import os
import yaml
import subprocess
import time
import glob

import numpy as np
import healpy as hp

from desqr.utils import is_found, mkdir

if __name__ == "__main__":
    from desqr.parser import Parser
    parser = Parser()
    parser.set_defaults(njobs=24)
    parser.set_defaults(mlimit=50) # GB
    args = parser.parse_args()

    print("Running catalog creation...")
    
    config = args.config
    hpxdir = config['hpxdir']
    catdir = mkdir(config['catdir'])
    keydir = mkdir(config['keydir'])
    logdir = mkdir(os.path.join(catdir,'log'))

    if args.pix: pixels = args.pix
    else: pixels = np.arange(hp.nside2npix(config['nside']))

    for i,pix in enumerate(pixels):
              
        infiles = glob.glob(hpxdir+'/*/*%05d.fits'%pix)
        catfile = os.path.join(catdir,config['catbase']%pix)
        keyfile = os.path.join(keydir,config['keybase']%pix)
        logfile = os.path.join(logdir,os.path.basename(catfile).replace('.fits','.log'))

        if len(infiles) == 0: continue
        if is_found(catfile,args.force): continue
        ra,dec = hp.pix2ang(config['nside'], pix, lonlat=True)
        print(f"({i+1}/{len(pixels)}): RA, Dec, Hpx = {ra:.2f}, {dec:.2f}, {i}")

        minbands = config.get('minbands')
        minbands = '--min-bands %s'%minbands if minbands else ''
        minepochs = config.get('minepochs')
        minepochs = '--min-epochs %s'%minepochs if minepochs else ''
        ebv = config.get('ebv',None)
        ebv = '--ebv %s'%ebv if ebv else ''

        force = '-f' if args.force else ''
        bands = ' '.join(['-b %s'%b for b in config.get('bands',[])])
        params=(' '.join(infiles),catfile,keyfile,bands,minbands,minepochs,ebv,force)
        cmd = 'catalog.py -v %s -o %s -k %s %s %s %s %s %s'%params

        #if args.queue == 'local':
        #    print(cmd)
        #    submit = cmd
        #else:
        #    submit = 'csub -o %s -n %s %s'%(logfile,args.njobs,cmd)
        submit = f"csub -q {args.queue} -o {logfile} -n {args.njobs} {cmd}"

        subprocess.call(submit,shell=True)
        if args.queue != 'local': time.sleep(args.sleep)
