#!/usr/bin/env python
"""
Sort input catalogs into healpix pixels.
"""
import os
import time
import glob
import subprocess
from collections import OrderedDict as odict
import yaml
from functools import partial

import numpy as np

from ugali.utils.shell import mkdir

def make_paths(config, bands=None):
    """ Make lists of raw file paths """
    explist = config['explist']
    rawdir = config['rawdir']
    basename = config.get('rawbase')
    exposures = np.recfromcsv(explist)

    if not bands: bands = config['bands']
    bands = list(bands)

    paths = odict([(b,[]) for b in bands])

    for band in bands:
        dirname = os.path.join(rawdir,band)
        files = glob.glob(os.path.join(dirname,'*.fits'))

        for exp in exposures[exposures['band'] == band]:
            edict = dict(zip(exp.dtype.names,exp))
            path = os.path.join(dirname,basename.format(**edict))
            paths[band].append(path)

            
        exist = np.in1d(paths[band],files)
        if not exist.all():
            msg  = "Missing files...\n"
            msg += str(np.array(paths[band])[~exist])
            raise Exception(msg)

    return paths

if __name__ == "__main__":
    from parser import Parser
    parser = Parser(description=__doc__)
    parser.set_defaults(njobs=None)
    args = parser.parse_args()

    config = yaml.load(open(args.config))
    rawdir = config['rawdir']
    hpxdir = config['hpxdir']

    for band in config['bands']:
        outdir = mkdir(os.path.join(hpxdir,band))
        logdir = mkdir(os.path.join(outdir,'log'))
        logfile = os.path.join(logdir,'pixelize_%s.log'%band)
        start,end = config['hpxbase'].rsplit('_',1)
        outbase = '_'.join([start%dict(band=band),end])

        # Create a txt file of file paths
        paths = make_paths(config,bands=band)[band]
        pathfile = os.path.join(hpxdir,'filepaths_%s.txt'%band)
        np.savetxt(pathfile,paths,fmt='%s')

        cmd = 'pixelize.py %s %s -o %s -n %i'%(pathfile,outdir,outbase,config['nside'])

        submit = 'csub -q %s -o %s '%(args.queue, logfile)
        if args.njobs: submit += '-n %s '%args.njobs
        submit += cmd
        if args.verbose: print(submit)
        if not args.dryrun: subprocess.call(submit,shell=True)
        time.sleep(args.sleep)
