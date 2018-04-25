#!/usr/bin/env python
"""
Skim data from a file.
"""
import logging
import fitsio
import utils

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infile',
                        help="input FITS file")
    parser.add_argument('outfile',
                        help="output FITS file")
    parser.add_argument('-c','--columns',action='append')
    parser.add_argument('-s','--select',default=None)
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-v','--verbose',action='store_true')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.getLogger().setLevel(level)
    
    utils.skim_data(args.infile,args.outfile,args.columns,args.select,args.force)
