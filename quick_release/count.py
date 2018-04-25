#!/usr/bin/env python
import fitsio
import numpy as np

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infiles',nargs='+')
    parser.add_argument('--filename',action='store_true')
    parser.add_argument('--columns',action='store_true')
    parser.add_argument('--colname')
    parser.add_argument('--nan',action='store_true')
    parser.add_argument('-v','--verbose',action='store_true')
    opts = parser.parse_args()

    TOTAL = 0
    NCOLS = 0
    nfiles = len(opts.infiles)
    for i,infile in enumerate(opts.infiles):
        if i % 100 == 0: print "(%i/%i)"%(i,nfiles)
        try:
            fits = fitsio.FITS(infile)
            TOTAL += fits[1].read_header()['NAXIS2']
            if opts.filename:
                filename = fits[1].read(columns='FILENAME')
                if np.any(~np.char.startswith(filename,'D00')):
                    print infile
            if opts.columns:
                ncols = fits[1].get_ncols()
                if NCOLS == -1:
                    NCOLS = ncols
                if ncols != NCOLS: print '(%s/%s):'%(ncols,NCOLS),infile
            if opts.colname:
                if opts.colname not in fits[1].get_colnames():
                    print infile
                col = fits[1].read(columns=opts.colname)
                if  np.any(np.isnan(col)):
                    print infile
        except Exception,msg:
            print infile
            print msg

    print TOTAL
