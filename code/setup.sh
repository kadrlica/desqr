#!/usr/bin/env bash

astrosetup

cvmfssetup
setup pyyaml
setup basemap
setup pyfits 3.3+1
setup fitsio 0.9.6+1
setup astropy 0.4.2+1
setup easyaccess 1.1.0+1
setup healpy 1.8.1+3
setup python 2.7.6+3
setup matplotlib 1.3.1+3

# localsetup
# setup fitsio 0.9.5.1rc1.cfitsio3370+5 --nolocks
# setup esutil 0.5.1+2 --nolocks
# setup astropy 0.4.2+1 --nolocks
# setup easyaccess 1.1.0+1 --nolocks
#  
# setup python 2.7.6+3 --nolocks
# setup matplotlib 1.3.1+3 --nolocks
# export TK_LIBRARY=$TK_DIR/lib/tk8.5

SOFTWARE=/home/s1/kadrlica/software
SCIRELEASE=$SOFTWARE/des-sci-release
export PATH=$PATH:$SCIRELEASE/users/kadrlica/catalog_coadd/code
export PYTHONPATH=$PYTHONPATH:$SCIRELEASE/users/kadrlica/catalog_coadd/code
export PYTHONPATH=$PYTHONPATH:$SOFTWARE/ugali/master

export PATH=$SOFTWARE/easyaccess/bin:$PATH
export PYTHONPATH=$SOFTWARE/easyaccess:$PYTHONPATH

export MAGLITES=$SOFTWARE/maglites/master
export PATH=$MAGLITES/bin:$PATH
export PYTHONPATH=$MAGLITES:$PYTHONPATH