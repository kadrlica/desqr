#!/usr/bin/env python
"""Constants used elsewhere in the code..."""

from collections import OrderedDict as odict

ZEROSTR = '1%05d%08d'
BANDS = ['g','r','i','z','Y']
#BANDS = ['g','r']
#MINBANDS = 1
MINBANDS = 2
COLORS = odict([('g','green'),('r','red'),('i','gold'),('z','purple'),('Y','gray')])
BADMAG = 99.
BADVAL = -1
# Matches (ill-advised?) choice by DESDM 
# BADMAGERR = 1. 
BADZP = 0b10
OBJECT_ID = 'QUICK_OBJECT_ID'
#OBJECT_ID = 'COADD_OBJECTS_ID'
UNIQUE_ID = 'UNIQUE_ID'
TAGS = ['Y1A1_FINALCUT','Y2N_FIRSTCUT','Y3N_FIRSTCUT','DESGW',
        'Y2A1_FINALCUT','Y3A1_FINALCUT',
        'MAGLITES_FIRSTCUT','MAGLITES_R4','Y3A1_MAGLITES',
        'BLISS','BLISS_Y17T2','BLISS_Y18T1',
        'DELVE']
NSIDES = [2048]
EDGE = 15 # Not used yet

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    opts = parser.parse_args()
