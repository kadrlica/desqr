#!/usr/bin/env python
import fitsio
import numpy as np

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infiles',nargs='+')
    parser.add_argument('--input',default='CATALOG_ID')
    parser.add_argument('--output',default='QUICK_OBJECT_ID')
    parser.add_argument('-v','--verbose',action='store_true')
    opts = parser.parse_args()

    for i,infile in enumerate(opts.infiles):
        if opts.verbose:
            print infile
        f = fitsio.FITS(infile,'rw')
        h = f[1].read_header()

        keys = np.array(h.keys())
        keys = keys[np.char.startswith(keys,'TTYPE')]
        names = np.char.strip([h[k] for k in keys])

        found = False
        for n,k in zip(names,keys):
            if n == opts.input:
                f[1].write_key(k, opts.output)
                found = True
                break
        
        if not found:
            print "%s not found in %s"%(opts.input,infile)
