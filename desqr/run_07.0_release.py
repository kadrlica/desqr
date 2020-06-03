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

    if 'footprint' in sections:
        exe = os.path.join(path,'footprint.py')
        survey = config.get('survey','des')
        nside = 2048 if args.nside is None else args.nside
        cmd = 'python %s %s -p -n %i --survey %s'%(exe,args.config,nside,survey)
        if args.verbose: cmd += ' -v'
        if args.band: bands = args.band
        else:         bands = config['bands'] + [''.join(config['bands'])]
        cmd += ' '+' '.join(['-b '+b for b in bands])
        if args.queue == 'local':
            print cmd
            submit = cmd
        else:
            submit = 'csub %s'%(cmd)

        subprocess.call(submit,shell=True)
        
    if 'depth' in sections:
        nside = 1024 if args.nside is None else args.nside
        exe = os.path.join(path,'depth.py')
        cmd = 'python %s %s -n %i'%(exe,args.config,nside)
        if args.verbose: cmd += ' -v'

        if args.queue == 'local':
            print cmd
            submit = cmd
        else:
            submit = 'csub %s'%(cmd)

        subprocess.call(submit,shell=True)

    if 'astrometry' in sections:
        exe = os.path.join(path,'astrometry.py')
        # r-band only
        if args.band: bands = args.band
        else:         bands = ['r']
        
        cmd = 'python %s %s'%(exe,args.config)
        cmd += ' '+' '.join(['-b '+b for b in bands])

        if args.verbose: cmd += ' -v'

        if args.queue == 'local':
            print cmd
            submit = cmd
        else:
            submit = 'csub %s'%(cmd)

        subprocess.call(submit,shell=True)


    #if 'astrometry' in sections:
    #    exe = os.path.join(path,'astrometry')
    #    cmd = 'python %s'%exe

    
