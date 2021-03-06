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
import healpy
from ugali.utils.shell import mkdir
from utils import found

if __name__ == "__main__":
    from parser import Parser
    parser = Parser()
    parser.add_argument('-p','--pix',default=None,action='append',type=int,
                        help='pixels to submit')
    args = parser.parse_args()
    
    force = '-f' if args.force else ''

    config = yaml.load(open(args.config))
    hpxdir = config['hpxdir']
    catdir = mkdir(config['catdir'])
    keydir = mkdir(config['keydir'])
    logdir = mkdir(os.path.join(catdir,'log'))

    pixels = args.pix
    if pixels is None:
        pixels = np.arange(healpy.nside2npix(config['nside']))

    for pix in pixels:
        infiles = glob.glob(hpxdir+'/*/*%05d*.fits'%pix)
        basename = 'hpx_%05d.fits'%pix
        outfile = os.path.join(catdir,'cat_'+basename)
        keyfile = os.path.join(keydir,'key_'+basename)
        logfile = os.path.join(logdir,basename.replace('.fits','.log'))

        if len(infiles) == 0: continue

        if os.path.exists(outfile) and not args.force:
            found(outfile)
            continue

        minbands = config.get('minbands')
        minbands = '--min-bands %s'%minbands if minbands else ''
        bands = ' '.join(['-b %s'%b for b in config.get('bands',[])])
        params=(' '.join(infiles),outfile,keyfile,bands,minbands,force)
        cmd = 'catalog.py -v %s -o %s -k %s %s %s %s'%params

        if args.queue == 'local':
            print cmd
            submit = cmd
        else:
            #submit = 'csub -o %s %s'%(logfile,cmd)
            submit = 'csub -o %s -n %s %s'%(logfile,args.njobs,cmd)

        subprocess.call(submit,shell=True)
        if args.queue != 'local': time.sleep(args.sleep)
