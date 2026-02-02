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

from utils import is_found, mkdir

if __name__ == "__main__":
    from parser import Parser
    parser = Parser()
    args = parser.parse_args()
    
    config = args.config
    hpxdir = config['hpxdir']
    catdir = mkdir(config['catdir'])
    keydir = mkdir(config['keydir'])
    logdir = mkdir(os.path.join(catdir,'log'))

    if args.pix: pixels = args.pix
    else: pixels = np.arange(hp.nside2npix(config['nside']))

    for i,pix in enumerate(pixels):
        print("(%s/%s): %s"%(i+1,len(pixels), pix))
              
        infiles = glob.glob(hpxdir+'/*/*%05d.fits'%pix)
        catfile = os.path.join(catdir,config['catbase']%pix)
        keyfile = os.path.join(keydir,config['keybase']%pix)
        logfile = os.path.join(logdir,os.path.basename(catfile).replace('.fits','.log'))

        if len(infiles) == 0: continue
        if is_found(catfile,args.force): continue

        minbands = config.get('minbands')
        minbands = '--min-bands %s'%minbands if minbands else ''

        ebv = config.get('ebv',None)
        ebv = '--ebv %s'%ebv if ebv else ''

        force = '-f' if args.force else ''
        bands = ' '.join(['-b %s'%b for b in config.get('bands',[])])
        params=(' '.join(infiles),catfile,keyfile,bands,minbands,ebv,force)
        cmd = 'catalog.py -v %s -o %s -k %s %s %s %s %s'%params

        if args.queue == 'local':
            print(cmd)
            submit = cmd
        else:
            #submit = 'csub -o %s %s'%(logfile,cmd)
            submit = 'csub -o %s -n %s %s'%(logfile,args.njobs,cmd)

        subprocess.call(submit,shell=True)
        if args.queue != 'local': time.sleep(args.sleep)
