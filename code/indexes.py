#!/usr/bin/env python
import os
import subprocess
import time
import sys

import easyaccess as ea

from ugali.utils.logger import logger
from utils import mkscratch, bfields
from const import BANDS, OBJECT_ID, UNIQUE_ID
from download import create_bitmap_index


OBJECT_COLUMNS  = [OBJECT_ID]
OBJECT_COLUMNS += ['RA','DEC','HPX2048']
OBJECT_COLUMNS += bfields(["MAG_PSF","MAG_AUTO"],BANDS)
OBJECT_COLUMNS += bfields(["WAVG_MAG_PSF","WAVG_MAG_AUTO"],BANDS)
OBJECT_COLUMNS += bfields(["SPREAD_MODEL","WAVG_SPREAD_MODEL"],BANDS)

#INDEX_COLUMN = ['CATALOG_ID','FILENAME','OBJECT_NUMBER','UNIQUE_ID','TAG']
INDEX_COLUMNS = [OBJECT_ID,UNIQUE_ID,'FILENAME','OBJECT_NUMBER','TAG']

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-t','--table',default='Y2Q2_OBJECTS_V0')
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-s','--section',default='dessci')
    args = parser.parse_args()

    con = ea.connect(section=args.section,quiet=True)
    cur = con.cursor()
    
    for column in OBJECT_COLUMNS:
        version = args.table.split('_')[0]
        index = '%s_%s_BMX'%(version,column)
        query=create_bitmap_index(column,index,args.table)
        print query
        con.onecmd(query)
        #cmd = 'easyaccess -s %s -c "%s"'%(args.section,query)
        #subprocess.call(cmd,shell=True)

        sys.stdout.flush()
        sys.stderr.flush()
        time.sleep(1)


    #for column in INDEX_COLUMNS:
    #    index = 'Y2Q1IDX_%s_BMX'%column
    #    query=create_bitmap_index(column,index,'Y2Q1_OBJECTS_INDEX_V1')
    #    cmd = 'easyaccess -s desoper -c "%s"'%query
    #    subprocess.call(cmd,shell=True)
