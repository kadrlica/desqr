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
    parser.add_argument('-s','--select',default=None)
    parser.add_argument('-f','--force',action='store_true')
    args = parser.parse_args()

    print("Reading %s..."%args.infile)
    data,header=fitsio.read(args.infile,header=True,columns=args.columns)

    if args.select:
        print "Applying selection: %s"%args.select
        if 'data' not in args.select:
            msg = "Invalid select statement: %s"%args.select
            raise Exception(msg)
        data = data[eval(args.select)]

    if len(data) == 0:
        print("No objects pass selection.")
    else:
        print("Writing %s..."%args.outfile)
        utils.write(args.outfile,data,header=header,force=args.force)
