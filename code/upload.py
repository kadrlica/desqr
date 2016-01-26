#!/usr/bin/env python
import os
import subprocess
import tempfile
import easyaccess as ea
import time
import sys

from ugali.utils.logger import logger
from utils import mkscratch


def upload(infile, force=False):
    basename = os.path.basename(infile)
    table = os.path.splitext(basename)[0].upper()
    
    cmd = 'easyaccess -s desoper -c '

    sql = '"'
    if opts.force: 
        sql += "DROP TABLE %s; "%table
    sql += ' load_table %s;'%infile
    sql += '"'
    cmd += sql
    subprocess.call(cmd,shell=True)

def sqlldr(infile):
    basename = os.path.basename(infile)
    ctrl = os.path.join(ldrdir,basename.replace('.fits','.ldr'))
    csv = ctrl + '.csv'
    create = ctrl + '.create.sql'
    params = (infile,table,ctrl,'CATALOG_ID')
    cmd = 'des-fits2table %s %s %s --primary-key %s --ext 1 '%params
    cmd += '--create' if i==0 else ''
        
    print cmd
    subprocess.call(cmd,shell=True)

    sqlldr = "time sqlldr <username>/<password>@leovip148.ncsa.uiuc.edu/DESOPER rows=10000 direct=true control=%s"%(ctrl)
    print sqlldr
    
    for f in [ctrl,csv,create]:
        if os.path.exists(f): os.remove(f)
    
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
        logger.info("Uploading %s..."%infile)
        if i == 0:
            con.load_table(infile,opts.table)
            grant = 'grant select on %s to DES_READER'%opts.table
            cur.execute(grant)
        else:
            con.append_table(infile,opts.table)
        sys.stdout.flush()
        sys.stderr.flush()
        time.sleep(1)

