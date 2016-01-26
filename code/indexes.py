#!/usr/bin/env python
import os
import subprocess
import time
import sys

from ugali.utils.logger import logger
from utils import mkscratch, bfields
from const import BANDS
from download import create_bitmap_index

OBJECT_COLUMNS = ['RA','DEC','HPX2048']
OBJECT_COLUMNS += bfields("WAVG_MAG_PSF",BANDS)
OBJECT_COLUMNS += bfields("WAVG_SPREAD_MODEL",BANDS)

#INDEX_COLUMN = ['CATALOG_ID','FILENAME','OBJECT_NUMBER','UNIQUE_ID','TAG']
INDEX_COLUMN = ['UNIQUE_ID','TAG']

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    opts = parser.parse_args()

    for column in INDEX_COLUMN:
        index = 'Y2Q1IDX_%s_BMX'%column
        query=create_bitmap_index(column,index,'Y2Q1_OBJECTS_INDEX_V1')
        cmd = 'easyaccess -s desoper -c "%s"'%query
        subprocess.call(cmd,shell=True)
