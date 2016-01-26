#!/usr/bin/env python
import fitsio
import utils
if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infile')
    parser.add_argument('outfile')
    parser.add_argument('-c','--columns',action='append')
    parser.add_argument('-f','--force',action='store_true')
    opts = parser.parse_args()

    print "Reading %s..."%opts.infile
    data,header=fitsio.read(opts.infile,header=True,columns=opts.columns)
    print "Writing %s..."%opts.outfile
    utils.write(opts.outfile,data,header=header,force=opts.force)
