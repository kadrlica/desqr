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
from utils import mkdir

SECTIONS = ['footprint','stargal','color','astrometry','depth','photometry']

if __name__ == "__main__":
    from parser import Parser
    parser = Parser(description=__doc__)
    parser.set_defaults(njobs=8) 
    parser.add_argument('-b','--band',default=None,action='append')
    parser.add_argument('--nside',default=None,type=int)
    parser.add_argument('-r','--run',action='append',choices=SECTIONS)
    args = parser.parse_args()

    force = '-f' if args.force else ''

    config = yaml.safe_load(open(args.config))
    all_bands = ''.join(config['bands'])
    sections = args.run if args.run else SECTIONS
    path = os.path.dirname(os.path.abspath(__file__))
    logdir = mkdir('release/log/')

    if 'footprint' in sections:
        print("Calculating survey footprint...")
        exe = os.path.join(path,'footprint.py')

        nside = 2048 if args.nside is None else args.nside
        if args.band: bands = args.band
        else:         bands = config['bands'] + [all_bands]

        cmd = 'python %s %s -p -n %i'%(exe,args.config,nside)
        #cmd += ' --nproc 4' # number of parallel processes
        if args.verbose: cmd += ' -v'
        cmd += ' '+' '.join(['-b '+b for b in bands])

        logfile = os.path.join(logdir,'footprint.log')
        submit = 'csub -q %s -n %s -o %s %s'%(args.queue, args.njobs, logfile, cmd)
        if args.verbose: print(submit)
        subprocess.call(submit,shell=True)

    if 'stargal' in sections:
        print("Calculating star/galaxy separation...")
        exe = os.path.join(path,'footprint.py')

        nside = 1024 if args.nside is None else args.nside
        if args.band: bands = args.band
        else:         bands = config['bands'] + [all_bands]

        for cls in ['star','gal']:
            cmd = 'python %s %s -p -n %i'%(exe,args.config,nside)
            cmd += ' -c %s'%cls
            #cmd += ' --nproc 4' # number of parallel processes
            if args.verbose: cmd += ' -v'
            cmd += ' '+' '.join(['-b '+b for b in bands])

            logfile = os.path.join(logdir,'footprint_%s.log'%cls)
            submit = 'csub -q %s -n %s -o %s %s'%(args.queue, args.njobs, logfile, cmd)
            if args.verbose: print(submit)
            subprocess.call(submit,shell=True)

    if 'color' in sections:
        print("Calculating color uniformity...")
        exe = os.path.join(path,'footprint.py')

        nside = 256 if args.nside is None else args.nside
        bands = ['gr','gi']
        #bands = ['gi']

        for cls in ['star']:
            cmd = 'python %s %s -p -n %i'%(exe,args.config,nside)
            cmd += ' -c %s'%cls
            #cmd += ' --nproc 4' # number of parallel processes
            if args.verbose: cmd += ' -v'
            cmd += ' '+' '.join(['-b '+b for b in bands])

            logfile = os.path.join(logdir,'color_%s.log'%cls)
            submit = 'csub -q %s -n %s -o %s %s'%(args.queue, args.njobs, logfile, cmd)
            if args.verbose: print(submit)
            subprocess.call(submit,shell=True)
        
    if 'depth' in sections:
        print("Calculating photometric depth...")
        exe = os.path.join(path,'depth.py')
        nside = 1024 if args.nside is None else args.nside
        types = [
            ['psf',  [5,10]],
            ['auto', [5,10]],
        ]
        for t,snrs in types:
            for snr in snrs:
                cmd  = 'python %s %s -n %i'%(exe,args.config,nside)
                cmd += ' --type %s --snr %s'%(t,snr)
                #cmd += ' --nproc 4' # number of parallel processes
                if args.verbose: cmd += ' -v'

                logfile = os.path.join(logdir,'depth.log')
                submit = 'csub -q %s -n %s -o %s %s'%(args.queue, args.njobs, logfile, cmd)
                if args.verbose: print(submit)
                subprocess.call(submit,shell=True)

    if 'astrometry' in sections:
        print("Calculating astrometric scatter...")
        exe = os.path.join(path,'astrometry.py')
        
        types = [
            #['hpx_rms'  , ['g','r','i','z']],
            ['gaia_edr3', ['r']],
            #['gaia_dr2' , ['r']],
        ]

        for t,bands in types:
            if args.band: bands = args.band
            for b in bands:
                cmd = 'python %s %s -b %s --type %s'%(exe,args.config,b,t)
                cmd += ' --nproc 8' # number of parallel processes
                if args.verbose: cmd += ' -v'
                 
                logfile = os.path.join(logdir,'astrometry_%s_%s.log'%(t,b))
                submit = 'csub -q %s -n %s -o %s %s'%(args.queue, args.njobs, logfile, cmd)
                if args.verbose: print(submit)
                subprocess.call(submit,shell=True)

                if t not in ['hpx_rms']: 
                    print("External astrometry uses only one band")
                    break

    if 'photometry' in sections:
        print("Calculting photometric scatter...")
        exe = os.path.join(path,'photometry.py')

        types = [
            #['wavg_rms' , ['g','r','i','z']],
            #['hpx_rms'  , ['g','r','i','z']],
            #['des_dr2'  , ['g','r','i','z']],
            #['gaia_dr2' , ['griz']],
            ['gaia_edr3', ['griz']],
        ]

        for t,bands in types:
            if args.band: bands = args.band
            if t.startswith('gaia'): bands = ['griz']
            for b in bands:
                cmd = 'python %s %s -b %s --type %s'%(exe,args.config,b,t)
                cmd += ' --nproc 4' # number of parallel processes
                if args.verbose: cmd += ' -v'
            
                logfile = os.path.join(logdir,'photometry_%s_%s.log'%(t,b))
                submit = 'csub -q %s -n %s -o %s %s'%(args.queue,args.njobs,logfile,cmd)
                if args.verbose: print(submit)
                subprocess.call(submit,shell=True)
                    
                # Gaia uses all bands at once
                if t.startswith('gaia'): 
                    print("Gaia photometry uses all bands griz")
                    break
