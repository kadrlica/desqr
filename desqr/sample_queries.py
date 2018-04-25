#!/usr/bin/env python
import subprocess
import pylab as plt
import numpy as np
import fitsio
import healpy
from os.path import join

from ugali.utils.healpix import ang2pix
from ugali.utils.shell import mkdir

def fornax_query(outfile='fornax.fits'):
    query = """--- This is a query to grab a sample of Fornax stars
select CATALOG_ID,RA,DEC,MAG_PSF_G,MAG_PSF_R,SPREAD_MODEL_R
from Y2Q1_OBJECTS_V1 sample block (1.0,1)
where RA between (40.0 - 1) and (40.0 + 1) 
and DEC between (-34.4 - 1) and (-34.4 + 1)
and MAG_PSF_G between 15 and 30
and MAG_PSF_R between 15 and 30
and abs(SPREAD_MODEL_R) < 0.01;
> %s
"""%outfile
    return query

def query(outfile,ra,dec,coord='equ',nside=128):
    query = """--- Select stars near a given location
select CATALOG_ID, RA, DEC, HPX2048,
WAVG_MAG_PSF_G, WAVG_MAG_PSF_R, WAVG_SPREAD_MODEL_R
from Y2Q1_OBJECTS_V1 
where %(position)s
and WAVG_MAG_PSF_G between 10 and 30
and WAVG_MAG_PSF_R between 10 and 30
and abs(WAVG_SPREAD_MODEL_R) < 0.002;
> %(outfile)s
"""

    if coord.lower() == 'equ':
        delta = 0.2
        position = """RA between %(ra).2f - %(delta).2f and %(ra).2f + %(delta).2f
and DEC between %(dec).2f - %(delta).2f and %(dec).2f + %(delta).2f"""%dict(ra=ra,dec=dec,delta=delta)
    elif coord.lower() == 'hpx':
        pixel = ang2pix(nside,ra,dec,nest=True)
        factor = healpy.nside2npix(2048)/healpy.nside2npix(nside)
        position = """HPX2048 > (%(pixel)i)*%(factor)i and HPX2048 < (%(pixel)i + 1)*%(factor)i"""%dict(pixel=pixel,factor=factor)
        
    return query%dict(position=position,outfile=outfile)

def messier2(outfile='messier2.fits',coord='equ'):
    ra,dec = 323.36,-0.82
    return query(outfile,ra,dec,coord,nside=128)
        
def fornax(outfile='fornax.fits',coord='equ'):
    ra,dec = 40.0,-34.4
    return query(outfile,ra,dec,coord,nside=128)


if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    opts = parser.parse_args()

    name = 'messier2'
    #outfile = 'fornax.fits'
    #sqlfile = 'fornax.sql'
    outdir = mkdir("release/examples")

    for coord in ['equ','hpx']:
        sqlfile = join(outdir,'%s_%s.sql'%(name,coord))
        outfile = join(outdir,'%s_%s.fits'%(name,coord))
        out = open(sqlfile,'w')
        out.write(messier2(outfile,coord=coord))
        out.close()

        cmd = 'easyaccess -s desoper -l %s'%sqlfile
        subprocess.call(cmd,shell=True)

        d = fitsio.read(outfile)
        
        plt.figure()
        plt.scatter(d['RA'],d['DEC'],s=2)
        plt.xlabel('RA (deg)')
        plt.ylabel('DEC (deg)')
        plt.title("%s (%s)"%(name.capitalize(),coord.upper()))
        outfile = join(outdir,'%s_%s_spatial.png'%(name,coord))
        plt.savefig(outfile,bbox_inches='tight')
         
        plt.figure()
        plt.scatter(d['WAVG_MAG_PSF_G']-d['WAVG_MAG_PSF_R'],d['WAVG_MAG_PSF_G'],s=2)
        plt.xlabel(r'$g-r$')
        plt.ylabel(r'$g$')
        plt.xlim(-1,1)
        plt.gca().invert_yaxis()
        plt.title("%s (%s)"%(name.capitalize(),coord.upper()))
        outfile = join(outdir,'%s_%s_cmd.png'%(name,coord))
        plt.savefig(outfile,bbox_inches='tight')

