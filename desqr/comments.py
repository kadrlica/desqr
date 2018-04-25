#!/usr/bin/env python
import os
import subprocess
import time
import sys
import logging

logging.basicConfig(stream=sys.stdout,format='%(message)s')

import easyaccess as ea
import yaml

from ugali.utils.logger import logger
from utils import mkscratch, bfields
from const import BANDS, OBJECT_ID, UNIQUE_ID
from download import create_bitmap_index

if __name__ == "__main__":
    import argparse
    description = "Add comments to table"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('-d','--dryrun',action='store_true')
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    config = yaml.load(open(args.config))
    section = config['db']
    data = yaml.load(open(config['dbtables']))

    con = ea.connect(section=section,quiet=True)
    cur = con.cursor()

    tables = data.keys()
    for table in tables:
        logging.debug(table)
        comment = data[table]['comment']
        cmd = "table %s '%s'"%(table,comment)
        logging.debug(cmd)
        con.do_add_comment(cmd)

        for key,val in data[table]['columns'].items():
            break # Don't do column comments for now
            for b in BANDS:
                column = key.format(b=b.upper())
                comment = val['comment']
                cmd = "column %s.%s '%s'"%(table,column,comment)
                logging.debug(cmd)
                if not args.dryrun:
                    con.do_add_comment(cmd)

                if column == key: break
            
        #break
