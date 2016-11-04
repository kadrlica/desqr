#!/usr/bin/env python
import os
import subprocess
import tempfile
import time
import sys

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', stream=sys.stdout)

import easyaccess as ea

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infiles',nargs='+')
    parser.add_argument('-t','--table',default='MyTable')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-s','--section',default='desoper')
    args = parser.parse_args()
    
    ea.set_color_term(False)
    con = ea.connect(section=args.section,quiet=True)
    cur = con.cursor()

    msg = "easyaccess v%s"%ea.__version__
    logging.info(msg)

    exists = con.check_table_exists(args.table)
    if exists and not args.force:
        msg = "Found table %s; skipping..."%args.table
        raise Exception(msg)
    elif exists and args.force:
        con.drop_table(args.table,purge=True)
     
    for i,infile in enumerate(args.infiles):
        msg = '\n' + 30*"-" + '\n'
        msg += "(%s/%s) Uploading %s..."%(i+1,len(args.infiles),infile)
        logging.info(msg)
        if i == 0:
            con.load_table(infile,args.table)
        else:
            con.append_table(infile,args.table)
        sys.stdout.flush()
        sys.stderr.flush()
        time.sleep(1)

    grant = 'grant select on %s to DES_READER'%args.table
    logging.info(msg)
    cur.execute(grant)
    
