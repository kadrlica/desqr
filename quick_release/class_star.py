#!/usr/bin/env python
import numpy as np
import pylab as plt
import fitsio

from matplotlib.colors import LogNorm

from catalog import WAVGCLASS

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-p','--pix',default=2585,type=int)
    opts = parser.parse_args()

    filename = 'cat/cat_hpx_%05d.fits'%opts.pix
    data = fitsio.read(filename)
        
CLASS = ['CLASS_STAR','WAVG_CLASS_STAR','WAVG_CLASSMIN_STAR','WAVG_CLASSMAX_STAR']
CLASS = [c+'_I' for c in CLASS]

spread = 'WAVG_SPREAD_MODEL_I'
mag = 'WAVG_MAG_PSF_I'
class_star = 'CLASS_STAR_I'

mmins = [15,19,21,22]
mmaxs = [19,21,22,23]

gbins = np.linspace(0.0,0.2,50)
sbins = np.linspace(0.8,1.0,50)

good = (data[mag] < max(mmaxs))&(data[mag] > min(mmins))
d = data[good]

plt.figure()
plt.hist2d(d[mag],d[class_star],bins=50,norm=LogNorm())
plt.colorbar(label='Counts')
plt.xlabel('g (mag)')
plt.ylabel('CLASS_STAR_I')
plt.title('HPX %05d (%g < i < %g)'%(opts.pix,min(mmins),max(mmaxs)))
outfile='plots/class_star_2d_%05d.png'%opts.pix
plt.savefig(outfile,bbox_inches='tight')
"""
for mmin,mmax in zip(mmins,mmaxs):
    print '%g < i < %g'%(mmin,mmax)
    good = (data[mag] < mmax)&(data[mag] > mmin)
    d = data[good]
 
    fig,ax = plt.subplots(1,2,figsize=(12,5))
 
    for var in CLASS:
        ax[0].hist(d[var],bins=gbins,lw=2,histtype='step',label=var)
        ax[1].hist(d[var],bins=sbins,lw=2,histtype='step',label=var)
    ax[1].legend(loc='upper left',fontsize=10)
 
    fig.suptitle(r'HPX %05d (%i < i < %i)'%(opts.pix,mmin,mmax))
    outfile = 'plots/class_star_%g_%g_%05d.png'%(mmin,mmax,opts.pix)
    plt.savefig(outfile,bbox_inches='tight')
"""

mmin,mmax = 15,23
good = (data[mag] > mmin)&(data[mag] < mmax)
d = data[good]

fig,ax = plt.subplots(1,2,figsize=(12,5))
bins = [np.linspace(mmin,mmax,100),np.linspace(-0.01,0.04,100)]
ax[0].hist2d(d[mag],d[spread],bins=bins,norm=LogNorm(),alpha=0.5,cmap='gray_r')
ax[1].hist2d(d[mag],d[spread],bins=bins,norm=LogNorm(),alpha=0.5,cmap='gray_r')

for var,c in zip(CLASS,['k','r','b','g']):
    print "Contour:",var
    #star_thresh,gal_thresh = 0.9,0.1
    star_thresh,gal_thresh = np.percentile(d[var],q=[80,20])
    print "Star Threshold:",star_thresh
    print "Galaxy Threshold:",gal_thresh
    scut = d[var] > star_thresh
    gcut = d[var] < gal_thresh
    cts,yb,xb,im=ax[0].hist2d(d[mag][scut],d[spread][scut],bins=bins,visible=0)
    extent=[yb.min(),yb.max(),xb.min(),xb.max()]
    ax[0].contour(cts.T,2,extent=extent,linewidths=2,colors=c)
    cts,yb,xb,im=ax[1].hist2d(d[mag][gcut],d[spread][gcut],bins=bins,visible=0)
    extent=[yb.min(),yb.max(),xb.min(),xb.max()]
    ax[1].contour(cts.T,2,extent=extent,linewidths=2,colors=c)
    ax[1].plot(np.nan,np.nan,'-',lw=2,label=var)

ax[0].set_ylim(-0.005,0.005)
ax[1].legend(loc='upper left',fontsize=10)
ax[0].set_ylabel(spread)
ax[0].set_xlabel(mag)
ax[1].set_xlabel(mag)
fig.suptitle("SPREAD_MODEL Distribution")

outfile = 'plots/class_star_spread_%05d.png'%opts.pix
plt.savefig(outfile,bbox_inches='tight')

plt.ion()

"""
for var,c in zip(['WAVG_SPREAD_MODEL_I'],['magenta']):
    scut = np.abs(d[var]) < 0.002 
    gcut = np.abs(d[var]) > 0.002 
    cts,yb,xb,im = ax[0].hist2d(d[mag][scut],d[spread][scut],bins=50,visible=0)
    extent=[yb.min(),yb.max(),xb.min(),xb.max()]
    ax[0].contour(cts.T,2,extent=extent,lw=3,colors=c)
    ax[0].set_ylim(-0.005,0.005)
    cts,yb,xb,im = ax[1].hist2d(d[mag][gcut],d[spread][gcut],bins=50,visible=0)
    extent=[yb.min(),yb.max(),xb.min(),xb.max()]
    ax[1].contour(cts.T,2,extent=extent,lw=3,colors=c)
    ax[1].plot(np.nan,np.nan,'-',lw=2,label=var)
"""
