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
    parser.add_argument('-s','--sleep',default=1,type=float)
    parser.add_argument('-q','--queue',default='condor')
    opts = parser.parse_args()
    
    force = '-f' if opts.force else ''

    config = yaml.load(open(opts.config))
    catdir = config['catdir']

    pixels = [1376,1382,1511,1960,2585]
    for pix in pixels:
        basename = 'cat_hpx_%05d.fits'%pix
        catfile = os.path.join(catdir,basename)
        logfile = 'plots/validate_%05d.log'%pix
        cmd = 'validate.py --pix %s'%pix

        if opts.queue == 'local':
            submit = cmd
        else:
            submit = 'csub -o %s %s'%(logfile,cmd)
        subprocess.call(submit,shell=True)
        if opts.queue != 'local': time.sleep(opts.sleep)
