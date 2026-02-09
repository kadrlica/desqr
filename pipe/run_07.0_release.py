#!/usr/bin/env python
""" Run release validation checks. """
import os
import subprocess

from utils import mkdir

SECTIONS = ['footprint','stargal','color','astrometry','depth','photometry']

if __name__ == "__main__":
    from parser import Parser
    parser = Parser(description=__doc__)
    parser.set_defaults(njobs=8) 
    parser.add_argument('--nproc', default=1, type=int,
                        help='number processes per job')
    parser.add_argument('--nside', default=None, type=int)
    #parser.add_argument('--mag', action='store_true')
    parser.add_argument('-r', '--run', action='append', choices=SECTIONS)
    args = parser.parse_args()

    config = args.config
    all_bands = ''.join(config['bands'])
    sections = args.run if args.run else SECTIONS
    path = os.path.dirname(os.path.abspath(__file__))
    logdir = mkdir('release/log')

    if 'footprint' in sections:
        print("Calculating survey footprint...")
        exe = 'footprint.py'
        #exe = 'python ' + os.path.join(path,exe)
        
        nside = 2048 if args.nside is None else args.nside
        if args.bands: bands = args.bands
        else:          bands = config['bands'] + [all_bands]

        cmd = '%s %s -p -n %i'%(exe,args.configfile,nside)
        cmd += ' --mag' # always apply magnitude selection
        cmd += ' --nproc %s'%args.nproc # number of parallel processes
        cmd += ' -v' if args.verbose else ''
        cmd += ' '+' '.join(['-b '+b for b in bands])

        logfile = os.path.join(logdir,'footprint.log')
        submit = 'csub -q %s -n %s -o %s %s'%(args.queue, args.njobs, logfile, cmd)
        if args.verbose: print(submit)
        subprocess.call(submit,shell=True)

    if 'stargal' in sections:
        print("Calculating star/galaxy separation...")
        exe = 'footprint.py'
        #exe = 'python ' + os.path.join(path,exe)
        
        nside = 1024 if args.nside is None else args.nside
        if args.bands: bands = args.bands
        else:          bands = config['bands'] + [all_bands]

        for cls in ['star','gal']:
            cmd = '%s %s -p -n %i'%(exe,args.configfile,nside)
            cmd += ' --mag' # always apply magnitude selection
            cmd += ' -e %s'%cls 
            cmd += ' --nproc %s'%args.nproc # number of parallel processes
            cmd += ' -v' if args.verbose else ''
            cmd += ' '+' '.join(['-b '+b for b in bands])

            logfile = os.path.join(logdir,'footprint_%s.log'%cls)
            submit = 'csub -q %s -n %s -o %s %s'%(args.queue, args.njobs, logfile, cmd)
            if args.verbose: print(submit)
            subprocess.call(submit,shell=True)

    if 'color' in sections:
        print("Calculating color uniformity...")
        exe = 'footprint.py'
        #exe = 'python ' + os.path.join(path,exe)
        
        nside = 256 if args.nside is None else args.nside
        bands = ['gr','gi']
        #bands = ['gi']

        for cls in ['star']:
            cmd = '%s %s -p -n %i'%(exe,args.configfile,nside)
            cmd += ' --color'
            #cmd += ' --mag' if args.mag else ''  # magnitude selection applied by color
            cmd += ' -e %s'%cls
            cmd += ' --nproc %s'%args.nproc # number of parallel processes
            cmd += ' -v' if args.verbose else ''
            cmd += ' '+' '.join(['-b '+b for b in bands])

            logfile = os.path.join(logdir,'color_%s.log'%cls)
            submit = 'csub -q %s -n %s -o %s %s'%(args.queue, args.njobs, logfile, cmd)
            if args.verbose: print(submit)
            subprocess.call(submit,shell=True)
        
    if 'depth' in sections:
        print("Calculating photometric depth...")
        exe = 'depth.py'
        #exe = 'python ' + os.path.join(path,exe)
        
        nside = 1024 if args.nside is None else args.nside
        types = [
            ['psf',     [5,10]], # sextractor psf
            ['auto',    [5,10]], # sextractor auto
            #['psf_mag', [5,10]], # fitvd psf
            #['bdf_mag', [5,10]], # fitvd bdf
        ]
        for t,snrs in types:
            for snr in snrs:
                cmd  = '%s %s -n %i'%(exe,args.configfile,nside)
                cmd += ' --type %s --snr %s'%(t,snr)
                cmd += ' --nproc %s'%args.nproc # number of parallel processes
                cmd += ' -v' if args.verbose else ''

                logfile = os.path.join(logdir,'depth.log')
                submit = 'csub -q %s -n %s -o %s %s'%(args.queue, args.njobs, logfile, cmd)
                if args.verbose: print(submit)
                subprocess.call(submit,shell=True)

    if 'astrometry' in sections:
        print("Calculating astrometric scatter...")
        exe = 'astrometry.py'
        #exe = 'python ' + os.path.join(path,exe)
        
        types = [
            #['hpx_rms'  , ['g','r','i','z']],
            #['gaia_dr2' , ['r']],
            #['gaia_edr3', ['r']],
            ['gaia_dr3', ['r']],
        ]

        for t,bands in types:
            if args.bands: bands = args.bands
            for b in bands:
                cmd = '%s %s -b %s --type %s'%(exe,args.configfile,b,t)
                cmd += ' --nproc %s'%args.nproc # number of parallel processes
                cmd += ' -v' if args.verbose else ''
                 
                logfile = os.path.join(logdir,'astrometry_%s_%s.log'%(t,b))
                submit = 'csub -q %s -n %s -o %s %s'%(args.queue, args.njobs, logfile, cmd)
                if args.verbose: print(submit)
                subprocess.call(submit,shell=True)

                if t not in ['hpx_rms']: 
                    print("External astrometry uses only one band.")
                    break

    if 'photometry' in sections:
        print("Calculting photometric scatter...")
        exe = 'photometry.py'
        #exe = 'python ' + os.path.join(path,exe)
        
        types = [
            #['wavg_rms' , ['g','r','i','z']],
            #['hpx_rms'  , ['g','r','i','z']],
            #['des_dr2'  , ['g','r','i','z']],
            #['gaia_dr2' , ['griz']],
            #['gaia_edr3', ['griz']],
            ['gaia_dr3', ['griz']],
        ]

        for t,bands in types:
            if args.bands: bands = args.bands
            if t.startswith('gaia'):
                bands = ['griz']
                print("Gaia photometry uses all bands: %s"%bands)
            for b in bands:
                cmd = '%s %s -b %s --type %s'%(exe,args.configfile,b,t)
                cmd += ' --nproc %s'%args.nproc # number of parallel processes
                cmd += ' -v' if args.verbose else ''
            
                logfile = os.path.join(logdir,'photometry_%s_%s.log'%(t,b))
                submit = 'csub -q %s -n %s -o %s %s'%(args.queue,args.njobs,logfile,cmd)
                if args.verbose: print(submit)
                subprocess.call(submit,shell=True)
