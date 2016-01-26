#!/usr/bin/env python
import os
import subprocess
import tempfile

from utils import ccdnum
from const import BANDS, TAGS

def quota():
    query="""-- Query space quota for current user
SELECT USERNAME,BYTES/1024/1024/1024 as "USED(GB)", MAX_BYTES/1024/1024/1024 as "QUOTA(GB)" from DBA_TS_QUOTAS where USERNAME = (SELECT USER FROM DUAL);
"""
    return query

def size():
    query="""-- List table sizes owned by current user
SELECT segment_name, segment_type, bytes/1024/1024/1024 GB 
FROM user_segments WHERE segment_type='TABLE' 
AND segment_name in (SELECT table_name FROM user_tables);
"""

def find_indices(table='Y1A1_OBJECTS'):
    kwargs=dict(table=table)
    query = """-- Select the indices from a give table
SELECT dba_tables.table_name, dba_indexes.index_name,
dba_ind_columns.column_name
FROM dba_tables
JOIN dba_indexes ON dba_indexes.table_name = dba_tables.table_name
JOIN dba_ind_columns ON dba_indexes.index_name =
dba_ind_columns.index_name
WHERE dba_tables.table_name='%(table)s'
ORDER BY dba_tables.table_name,dba_indexes.index_name;
"""%kwargs
    return query

def find_primary_key(table='Y1A1_OBJECTS'):
    kwargs=dict(table=table)
    query = """-- Find primary key of table
SELECT cols.table_name, cols.column_name, cols.position, cons.status, cons.owner
FROM all_constraints cons, all_cons_columns cols
WHERE cols.table_name = '%(table)s'
AND cons.constraint_type = 'P'
AND cons.constraint_name = cols.constraint_name
AND cons.owner = cols.owner
ORDER BY cols.table_name, cols.position;
"""%kwargs
    return query

def synonym(table):
    query = "create synonym Y2Q1_ZEROPOINTS_V1 for erykoff.y2n_y1a1_qslr_v6@DESSCI;"
    return query

def read_only(table):
    query = "ALTER TABLE %s READ ONLY;"%table
    return query

def primary_key():
    query = """-- Add a primary key to Y2Q1_OBJECTS
ALTER TABLE Y2Q1_OBJECTS_V1 
ADD CONSTRAINT Y2Q1_OBJECTS_PK PRIMARY KEY (CATALOG_ID);
"""
    return query

def create_bitmap_index(column,index=None,table='Y2Q1_OBJECTS_V1'):
    if index is None: index = 'Y2Q1_%s_BMX'%column
    query = "create bitmap index %(index)s on %(table)s(%(column)s);"
    return query%dict(column=column,index=index,table=table)

def y1a1_object_query(expnum=None):
    kwargs = dict(expnum=expnum)
    query="""-- Y1A1_FINALCUT single-epoch catalog download
-- fluxerr = flux * magerr / (2.5/ln(10))
SELECT CAST('D00'||ev.expnum||'_'||CAST(i.BAND AS VARCHAR(1))||'_c'||i.ccd||'_fake.fits' as VARCHAR(48)) as FILENAME,
CAST('D00'||ev.expnum as VARCHAR(9)) as UNITNAME, 
CAST(-1 AS INT) AS REQNUM, CAST(-1 AS INT) AS ATTNUM, 
CAST('Y1A1_FINALCUT' as VARCHAR(13)) as TAG,
i.EXPNUM, i.ccd as CCDNUM, CAST(i.BAND AS VARCHAR(1)) AS BAND,
ev.T_EFF, o.FWHM_WORLD, o.FLAGS,
o.OBJECT_ID as OBJECT_NUMBER, o.RA, o.DEC,
POWER(10, (o.MAG_PSF - 25.0)/-2.5) as FLUX_PSF,
POWER(10, (o.MAG_PSF - 25.0)/-2.5)*o.MAGERR_PSF/1.0857362 as FLUXERR_PSF,
o.FLUX_AUTO, o.FLUXERR_AUTO,
o.CLASS_STAR, o.SPREAD_MODEL, o.SPREADERR_MODEL
FROM y1a1_firstcut_eval ev, y1a1_objects o, y1a1_image i
WHERE ev.expnum = %(expnum)i and ev.exposureid = i.exposureid and o.imageid = i.id 
and o.mag_psf < 90 and o.mag_auto < 90 and o.flags < 4 and o.magerr_psf < 0.5; 
"""%kwargs
    return query

def y2n_object_query(expnum=None,reqnum=None,attnum=None):
    kwargs=dict(expnum=expnum,reqnum=reqnum,attnum=attnum)
    kwargs['unitname'] = 'D%(expnum)08d'%kwargs
    kwargs['filename'] = "%(unitname)s_%%_r%(reqnum)dp%(attnum)02d_%%"%kwargs
    # Robert suggests the CCDNUM can come from:
    #select o.filename, c.ccdnum from prod.se_object o, prod.catalog c where c.filename=o.filename and rownum < 10;
    # Also, it would be good to cut objects near the CCD edges.
    #... o.XWIN_IMAGE between 15 and 2048-30,YWIN_IMAGE between 15 and 4096-30
    query = """-- Y2N_FIRSTCUT single-epoch catalog download
-- magerr = 2.5/ln(10) * fluxerr/flux
SELECT CAST(o.FILENAME as VARCHAR(48)) as FILENAME, 
CAST(ev.UNITNAME AS VARCHAR(9)) as UNITNAME, 
ev.REQNUM, ev.ATTNUM, 
CAST('Y2N_FIRSTCUT' AS VARCHAR(13)) as TAG,
ev.EXPNUM, CAST(SUBSTR(o.FILENAME,14,2) AS INT) as CCDNUM,
CAST(o.BAND AS VARCHAR(1)) AS BAND, ev.T_EFF, o.FWHM_WORLD, o.FLAGS, 
o.OBJECT_NUMBER, o.RA, o.DEC,
o.FLUX_PSF, o.FLUXERR_PSF,
o.FLUX_AUTO, o.FLUXERR_AUTO,
o.CLASS_STAR, o.SPREAD_MODEL, o.SPREADERR_MODEL 
FROM prod.se_object o, prod.firstcut_eval ev
WHERE o.filename like '%(filename)s'
and ev.unitname='%(unitname)s' and ev.reqnum=%(reqnum)i and ev.attnum=%(attnum)i
and o.flux_psf > 0 and o.flux_auto > 0 and o.flags < 4
and 1.0857362*(o.fluxerr_psf/o.flux_psf) < 0.5;"""%kwargs
    return query

def y1a1_coadd_objects_query():
    query = """-- Y1A1 object download
select coadd_objects_id, ra, dec,
wavgcalib_mag_psf_g - xcorr_sfd98_g as wavg_mag_psf_g,
wavgcalib_mag_psf_r - xcorr_sfd98_r as wavg_mag_psf_r,
wavgcalib_mag_psf_i - xcorr_sfd98_i as wavg_mag_psf_i, 
wavgcalib_mag_psf_z - xcorr_sfd98_z as wavg_mag_psf_z, 
wavgcalib_mag_psf_y - xcorr_sfd98_y as wavg_mag_psf_y,
wavg_magerr_psf_g, 
wavg_magerr_psf_r, 
wavg_magerr_psf_i, 
wavg_magerr_psf_z, 
wavg_magerr_psf_y
from y1a1_coadd_objects
where abs(wavg_spread_model_r) < 0.002 
and wavgcalib_mag_psf_r - xcorr_sfd98_r < 23
and wavg_magerr_psf_r < 0.1;
"""
    return query

def y1a1_finalcut_query(band):
    kwargs = dict(band=band)
    query = """-- Download Y1A1_FINALCUT objects
SELECT OBJECT_ID, COADD_OBJECTS_ID,
RA,DEC,FLAGS,
CAST(BAND AS VARCHAR(1)) AS BAND,
MAG_PSF, MAGERR_PSF, MAG_AUTO, MAGERR_AUTO,
SPREAD_MODEL,SPREADERR_MODEL,CLASS_STAR
FROM Y1A1_OBJECTS
WHERE MAG_PSF < 90 AND MAG_AUTO < 90 
AND MAGERR_PSF < 1 AND MAGERR_AUTO < 1
AND FLAGS < 4 AND BAND like '%(band)s%%';""" % kwargs
    return query

def exposure_query(tag=None,program='survey'):
    y1a1 = y1a1_exposure_query(program)
    y2n = y2n_exposure_query(program)
    if tag is None:
        query = y1a1.split('ORDER BY')[0]
        #query += "and e.object like 'DES supernova hex SN-C1%'\n"
        query += "UNION ALL\n" 
        query += y2n.split('ORDER BY')[0]
        #query += "and e.object like 'DES supernova hex SN-C1%'\n"
        query += "ORDER BY expnum;"
        return query
    elif tag.upper() == 'Y1A1':
        return y1a1
    elif tag.upper() == 'Y2N':
        return y2n
    else:
        msg = 'Unrecognized tag: %s'%tag
        raise Exception(msg)

def y1a1_exposure_query(program='survey'):
    query = """-- Y1A1 exposures that pass firstcut eval
SELECT CAST('D00'||ev.expnum as VARCHAR(9)) as UNITNAME, 
ev.EXPNUM, TO_CHAR(-1) AS REQNUM, TO_CHAR(-1) AS ATTNUM,
e.TELRA, e.TELDEC, e.NITE, 
CAST(e.BAND AS VARCHAR(1)) AS BAND, ev.T_EFF, 'Y1A1_FINALCUT' AS TAG
from y1a1_firstcut_eval ev, y1a1_exposure e 
where e.id = ev.exposureid 
and ev.accepted = 'True' and ev.program = '%s'
ORDER BY e.expnum;"""%program
    return query    

def y2n_exposure_query(program='survey'):
    query = """-- Y2N exposures that pass firstcut eval
SELECT ev.UNITNAME, ev.EXPNUM, 
TO_CHAR(ev.REQNUM) AS REQNUM, TO_CHAR(ev.ATTNUM) AS ATTNUM, 
e.TRADEG as TELRA, e.TDECDEG as TELDEC, e.NITE, 
CAST(e.BAND AS VARCHAR(1)) AS BAND, ev.T_EFF, 'Y2N_FIRSTCUT' as TAG
from prod.firstcut_eval ev, prod.exposure e, prod.ops_proctag t
WHERE e.expnum = ev.expnum 
and ev.unitname = t.unitname and ev.reqnum = t.reqnum and ev.attnum = t.attnum
and ev.accepted = 'True' and ev.program = '%s'
ORDER BY ev.expnum;"""%program
    return query

def desgw_exposure_query():
    query = """-- DESGW exposures that pass firstcut eval
SELECT ev.UNITNAME, ev.EXPNUM, 
TO_CHAR(ev.REQNUM) AS REQNUM, TO_CHAR(ev.ATTNUM) AS ATTNUM, 
e.TRADEG as TELRA, e.TDECDEG as TELDEC, e.NITE, 
CAST(e.BAND AS VARCHAR(1)) AS BAND, ev.T_EFF, 'DESGW' as TAG
from prod.firstcut_eval ev, prod.exposure e
WHERE e.expnum = ev.expnum 
and (
(ev.expnum >= 475905 and ev.expnum <= 475986) or
(ev.expnum >= 476335 and ev.expnum <= 476353) or
(ev.expnum >= 476942 and ev.expnum <= 477025) or
(ev.expnum >= 482835 and ev.expnum <= 482849) or
(ev.expnum >= 482856 and ev.expnum <= 482924) 
) ORDER BY ev.expnum;"""
    return query

def zeropoint_query():
    query = """-- qSLR zeropoints
SELECT * from erykoff.y2n_qslr_v2 zp
ORDER BY expnum;"""

    query = """-- qSLR zeropoints
SELECT * from erykoff.y2n_y1a1_qslr_v3
ORDER BY expnum, ccdnum;"""

    query = """-- qSLR zeropoints
SELECT * from erykoff.y2n_y1a1_qslr_v4
ORDER BY expnum, ccdnum;"""

    return query

qslr_query = zeropoint_query

def gcm_query():
    query = """-- Select Y1A1 GCM zeropoints
select EXPNUM, CCD AS CCDNUM, 
CAST(BAND AS VARCHAR(1)) AS BAND, RA, DEC, 
ZEROPOINT, SIGMA_ZEROPOINT 
from Y1A1_IMAGE where ZEROPOINT is not NULL;
"""
    return query


def download(outfile,query,sqlfile=None,section='desoper',force=False):
    if os.path.exists(outfile) and not force:
        msg = "Found %s; skipping..."%outfile
        raise Exception(msg)

    if sqlfile:
        sql = open(sqlfile,'w')
        remove = False
    else:
        sql = tempfile.NamedTemporaryFile(delete=False)
        sqlfile = sql.name
        remove = True

    sql.write(query)
    sql.write('\n'+'> %s \n'%outfile)
    sql.close()
    
    cmd = 'easyaccess -s %s -l %s'%(section,sqlfile)
    print cmd
    subprocess.call(cmd,shell=True)
    if remove: os.remove(sqlfile)

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('outfile',help="Output file")
    parser.add_argument('-f','--force',action='store_true')
    parser.add_argument('-l','--loadsql',default=None)
    parser.add_argument('-e','--expnum',default=None,type=int)
    parser.add_argument('-r','--reqnum',default=None,type=int)
    parser.add_argument('-a','--attnum',default=None,type=int)
    parser.add_argument('-t','--tag',default=None,choices=TAGS)
    opts = parser.parse_args()

    section = 'desoper'
    if opts.tag == 'Y2N_FIRSTCUT':
        query = y2n_object_query(opts.expnum,opts.reqnum,opts.attnum)
    elif opts.tag == 'Y1A1_FINALCUT':
        query = y1a1_object_query(expnum=opts.expnum)
    elif opts.tag == 'DESGW':
        query = y2n_object_query(opts.expnum,opts.reqnum,opts.attnum)
    else:
        raise Exception('...')

    download(opts.outfile,query,opts.loadsql,section,opts.force)

    # Should be downloaded by default
    #ccdnum(opts.outfile,force=True)
