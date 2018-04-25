#!/usr/bin/env python

y2n_query = "select ra,dec,flux_psf,flux_auto,x_image,y_image from prod.se_object where filename like 'D00381976_g_c%_r1203p01_red-fullcat.fits'; > tmp.fits"

y1a1_query = "select o.ra,o.dec,o.mag_psf,o.mag_auto,o.x_image,o.y_image from y1a1_objects o, y1a1_firstcut_eval ev, y1a1_image i where ev.expnum = 229272 and ev.exposureid = i.exposureid and o.imageid = i.id;"

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    opts = parser.parse_args()
