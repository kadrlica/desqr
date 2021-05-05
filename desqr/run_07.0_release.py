#!/usr/bin/env python
"""
Create release validation checks.
"""
import os
import yaml
import subprocess
import time
import glob

import numpy as np
import healpy

SECTIONS = ['footprint','astrometry','depth']

if __name__ == "__main__":
    from parser import Parser
    parser = Parser(description=__doc__)
    parser.add_argument('-b','--band',default=None,action='append')
    parser.add_argument('--nside',default=None,type=int)
    parser.add_argument('-r','--run',action='append',choices=SECTIONS)
    args = parser.parse_args()

    force = '-f' if args.force else ''

    config = yaml.safe_load(open(args.config))
    sections = args.run if args.run else SECTIONS
    path = os.path.dirname(os.path.abspath(__file__))

    survey = config.get('survey','des')

    if 'footprint' in sections:
        exe = os.path.join(path,'footprint.py')
        nside = 2048 if args.nside is None else args.nside
        cmd = 'python %s %s -p -n %i --survey %s'%(exe,args.config,nside,survey)
        if args.verbose: cmd += ' -v'
        if args.band: bands = args.band
        else:         bands = config['bands'] + [''.join(config['bands'])]
        cmd += ' '+' '.join(['-b '+b for b in bands])

        submit = 'csub -q %s %s'%(args.queue, cmd)
        subprocess.call(submit,shell=True)
        
    if 'depth' in sections:
        exe = os.path.join(path,'depth.py')
        nside = 1024 if args.nside is None else args.nside
        snr = 10
        mag = 'psf'
        cmd = 'python %s %s -n %i --snr %s --survey %s --mag %s'%(exe,args.config,nside,snr,survey,mag)
        if args.verbose: cmd += ' -v'

        submit = 'csub -q %s %s'%(args.queue,cmd)
        subprocess.call(submit,shell=True)

    if 'astrometry' in sections:
        """ Calculate internal/external astrometry """
        exe = os.path.join(path,'astrometry.py')
        # r-band only
        if args.band: bands = args.band
        else:         bands = ['r']

        cmd = 'python %s %s'%(exe,args.config)
        if args.verbose: cmd += ' -v'
        # Internal or external astrometry?
        cmd += ' --type internal'
        #cmd += ' --type external'

        cmd += ' '+' '.join(['-b '+b for b in bands])

        submit = 'csub -q %s %s'%(args.queue,cmd)
        subprocess.call(submit,shell=True)

    if 'photometry' in sections:
        """ Compare photometry """
        exe = os.path.join(path,'photometry.py')

        # r-band only
        bands = ['r']

        for b in bands:
            cmd = 'python %s %s -b %s'%(exe,args.config,b)
            #cmd += ' --type gaia'
            cmd += ' --type rms'
            if args.verbose: cmd += ' -v'

            submit = 'csub -q %s %s'%(args.queue,cmd)
            subprocess.call(submit,shell=True)

    
