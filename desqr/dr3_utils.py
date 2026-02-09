import scipy 
import numpy as np
from collections import OrderedDict as odict
import ugali.utils.healpix as healpix
from ugali.utils.projector import cel2gal
import healpy as hp

Y3A1 = odict([
        ('g',3.185),
        ('r',2.140),
        ('i',1.571),
        ('z',1.198),
        ('Y',1.052),
        ('y',1.052),
        ])
COEFF = Y3A1

def bdf_extended_class(data):
    """                                                                                                                                                           
    https://cdcvs.fnal.gov/redmine/projects/des-y6/wiki/Y6A2_Extended_Classifier                                                                                  
                                                                                                                                                                  
    https://github.com/des-science/science_release/blob/master/object_classification/y6/y6a2_object_classification_fitvd.py                                       
    """
    x = [-1.        ,  0.79891862,  0.90845217,  0.98558583,  1.05791208,
         1.13603715,  1.22479487,  1.33572223,  1.48983602,  1.74124395,
         2.43187589,  5.        ]
    y_1 = [0.028, 0.028, 0.008, 0.   , 0.004, 0.012, 0.012, 0.004, 0.012,
           0.024, 0.04 , 0.04 ]
    y_2 = [-0.028, -0.028, -0.04 , -0.032, -0.036, -0.032, -0.028, -0.016,
           -0.012,  0.008,  0.016,  0.016]
    y_3 = [-0.1  , -0.1  , -0.1  , -0.1  , -0.1  , -0.1  , -0.1  , -0.1  ,
           -0.012,  0.008,  0.016,  0.016]

    f_array = [scipy.interpolate.interp1d(x, y_1,bounds_error=False, fill_value=np.nan),
               scipy.interpolate.interp1d(x, y_2,bounds_error=False, fill_value=np.nan),
               scipy.interpolate.interp1d(x, y_3,bounds_error=False, fill_value=np.nan)]

    try:
        x_data = np.log10(data['BDF_S2N'])
        y_data = data['BDF_T']
    except ValueError:
        x_data = np.log10(data['bdf_s2n'])
        y_data = data['bdf_T']

    ext = np.tile(0, len(x_data))
    for f in f_array: 
        selection = (y_data > f(x_data))    
        ext += selection.astype(int)
    ext[(x_data < np.min(x)) | (x_data > np.max(x))]=-9
    ext=np.where(np.isfinite(ext), ext, -9)
    return ext

def ebv(ra,dec,ebvmap=None):
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    ebvmap = hp.read_map(ebvmap, verbose=False)
    
    glon,glat = cel2gal(ra,dec)
    ebv = healpix.get_interp_val(ebvmap,glon,glat)
    return ebv

def extinction(ebv,band):
    """
    Calculate the extinction from the E(B-V) value and the band.
    
    ebv  : The dust value
    band : The DES band (string or array)
    """
   
    band = np.repeat(band,len(ebv))
        
    bands,inv = np.unique(band,return_inverse=True)
    values = np.array([COEFF[b] for b in bands])
    return values[inv] * ebv