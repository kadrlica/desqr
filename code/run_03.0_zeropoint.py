#!/usr/bin/env python
import os
import yaml
import subprocess
import time
import glob

from ugali.utils.shell import mkdir

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-s','--sleep',default=5,type=float)
    parser.add_argument('-q','--queue',default='condor')
    opts = parser.parse_args()

    config = yaml.load(open(opts.config))
    hpxdir = config['hpxdir']
    zpfile = config['zpfile']
    blfile = config['blfile']
    blacklist = '-b %s'%blfile if blfile is not None else ''
    force = '-f' if opts.force else ''

    for band in config['bands']:
        indir = os.path.join(hpxdir,band)
        infiles = sorted(glob.glob(indir+'/*.fits'))
        logdir = mkdir(os.path.join(indir,'log'))
        for infile in infiles:
            logbase = ('zp_'+os.path.basename(infile)).replace('.fits','.log')     
            logfile = os.path.join(logdir,logbase)

            cmd = 'zeropoint.py %s %s %s %s'%(infile,zpfile,blacklist,force)

            if opts.queue == 'local':
                submit = cmd
            else:
                submit = 'csub -o %s %s'%(logfile,cmd)
            subprocess.call(submit,shell=True)
            time.sleep(opts.sleep)
        if opts.queue != 'local': time.sleep(10)
            
