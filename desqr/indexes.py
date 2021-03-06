#!/usr/bin/env python
import os
import subprocess
import time
import sys
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p', stream=sys.stdout)

import easyaccess as ea
import yaml

from ugali.utils.logger import logger
from utils import mkscratch, bfields
from const import BANDS, OBJECT_ID, UNIQUE_ID
from download import create_bitmap_index

QUERY = dict(
    BTX='create index {index} on {table}({column})',
    BMX='create bitmap index {index} on {table}({column})',
    PK='alter table {table} add constraint {index} primary key ({column})',
)

INDEX = dict(
    BTX='{idxname}_{column}_BTX',
    BMX='{idxname}_{column}_BMX',
    PK ='{idxname}_PK',
)

import logging


def drop_index(cursor, index):
    query = "DROP INDEX %s"%(index)
    try:
        return cursor.execute(query)
    except:
        return
    
if __name__ == "__main__":
    import argparse
    description = "Add comments to table"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config')
    parser.add_argument('-t','--table',action='append')
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
    con.drop_index = drop_index
    cursor = con.cursor()

    tables = data.keys()
    for table in tables:
        if args.table and (table not in args.table): continue

        logging.debug(table)
        idxname = data[table].get('idxname',table.split('_')[0])

        params = dict(idxname=idxname,table=table)

        for key,val in data[table]['columns'].items():
            if not val.get('index',None): continue
            for b in BANDS:
                params.update(b=b.upper())
                column = key.format(**params)
                params.update(column=column)

                idx = val['index']
                index = INDEX[idx].format(**params) if idx in INDEX else None
                params.update(index=index)

                query = QUERY[idx] if idx in QUERY else idx
                query = query.format(**params)

                logging.debug(query)


                if not args.dryrun:
                    if args.force and index is not None:
                        drop_index(cursor,index)
                    try:
                        cursor.execute(query)
                    except Exception as e:
                        print(e)

                if column == key: break

        for key,val in data[table]['indexes'].items():
            index = key.format(**params)
            params.update(index=index)
            query = val['query'].format(**params)

            logging.debug(query)
            if not args.dryrun:
                if args.force:
                    drop_index(con.cursor(),index)
                try: 
                    cursor.execute(query)
                except Exception as e:
                    print(e)
        break
