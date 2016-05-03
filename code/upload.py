#!/usr/bin/env python
import os
import subprocess
import tempfile
import time
import sys

import easyaccess as ea
from ugali.utils.logger import logger

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infiles',nargs='+')
    parser.add_argument('-t','--table',default='MyTable')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-s','--section',default='desoper')
    opts = parser.parse_args()

    con = ea.connect(section=opts.section,quiet=True)
    cur = con.cursor()

    exists = con.check_table_exists(opts.table)
    if exists and not opts.force:
        msg = "Found table %s; skipping..."%opts.table
        raise Exception(msg)
    elif exists and opts.force:
        con.drop_table(opts.table)
     
    for i,infile in enumerate(opts.infiles):
        msg = "(%s/%s) Uploading %s..."%(i,len(opts.infiles),infile)
        logger.info(msg)
        if i == 0:
            con.load_table(infile,opts.table)
            grant = 'grant select on %s to DES_READER'%opts.table
            cur.execute(grant)
        else:
            con.append_table(infile,opts.table)
        sys.stdout.flush()
        sys.stderr.flush()
        time.sleep(1)

