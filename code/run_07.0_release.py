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

SECTIONS = ['footprint','astrometry','depth']

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-s','--sleep',default=2,type=float)
    parser.add_argument('-q','--queue',default='local')
    parser.add_argument('--section',action='append',choices=SECTIONS)
    parser.add_argument('-v','--verbose',action='store_true')
    args = parser.parse_args()
    
    force = '-f' if args.force else ''

    config = yaml.load(open(args.config))
    sections = args.section if args.section else SECTIONS
    path = os.path.dirname(os.path.abspath(__file__))

    if 'footprint' in sections:
        exe = os.path.join(path,'footprint.py')
        cmd = 'python %s --survey %s'%(exe,config.get('survey','des'))
        if args.verbose: cmd += ' -v'
            
        if args.queue == 'local':
            print cmd
            submit = cmd
        else:
            submit = 'csub %s'%(cmd)

        subprocess.call(submit,shell=True)
        
    if 'depth' in sections:
        exe = os.path.join(path,'depth.py')
        cmd = 'python %s'%exe
        if args.verbose: cmd += ' -v'

        if args.queue == 'local':
            print cmd
            submit = cmd
        else:
            submit = 'csub %s'%(cmd)

        subprocess.call(submit,shell=True)

    if 'astrometry' in sections:
        exe = os.path.join(path,'astrometry.py')
        cmd = 'python %s %s'%(exe,args.config)
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

    
