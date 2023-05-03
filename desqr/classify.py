#!/usr/bin/env python
"""
Star-galaxy classification module.
"""
__author__ = "Alex Drlica-Wagner"

import numpy as np
import pylab as plt
import scipy.interpolate
import pandas as pd

def hsc_concentration_class(mag_psf, mag_cmodel):
    """ This is based on the HSC PDR1 star/galaxy classification 

    Roughly inspired by Coupon et al. 2017: https://arxiv.org/abs/1705.00622

    Returns
    -------
    sel : boolean HSC classification; 0 = stars; 1 = galaxies
    """
    # Calculate the concentration parameter
    hsc_concentration_i = (mag_psf - mag_cmodel)
    sel_stars = ( hsc_concentration_i < 0.03) | \
                ((hsc_concentration_i < 0.1) & (mag_psf < 22.))

    return ~sel_stars

def hsc_color_class(mag_psf_g, mag_psf_r, mag_psf_i):
    """ Stellar locus selection based classification.

    Returns
    -------
    sel : boolean HSC classification; 0 = stars; 1 = galaxies
    """
    pass

def hsc_extended_class(extendedness):
    """ This is based on the HSC PDR2 star/galaxy classification.

    Parameters
    ----------
    data : the hsc data
    var  : the name of the extendedness variable

    Returns
    -------
    sel : boolean HSC classification; 0 = stars; 1 = galaxies
    """
    sel_stars = (extendedness == 0)
    return ~sel_stars

def wavg_interp_array_y6a2():
    x   = [-1.0        ,         -0.01,           0.01,          1.0]
    y_1 = [ 3.0 + 0.005,  0.03 + 0.005,  -0.03 + 0.005, -3.0 + 0.005]
    y_2 = [ 1.0 + 0.003,  0.01 + 0.003,  -0.01 + 0.003, -1.0 + 0.003]
    y_3 = [-0.5 + 0.001,-0.005 + 0.001,  0.005 + 0.001,  0.5 + 0.001]

    f_array = [scipy.interpolate.interp1d(x, y_1),
               scipy.interpolate.interp1d(x, y_2),
               scipy.interpolate.interp1d(x, y_3)]
    return x, f_array

def wavg_extended_class(spread_model,spread_model_error):
    """ DES DR2 extended classifier for stars and galaxies using WAVG quantities.

        0 : High-confidence stars
        1 : Likely stars
        2 : Likely galaxies
        3 : High-confidence galaxies
       -9 : Failures

    Parameters
    ----------
    spread_model      : SExtractor spread model value in a band
    spread_model_error: SExtractor spread model error in a band
    
    Returns
    -------
    extended_class    : Extended classifier output
    """
    extend_val  = ((spread_model + 3.0*spread_model_error) > 0.005)*1
    extend_val += ((spread_model + 1.0*spread_model_error) > 0.003)*1
    extend_val += ((spread_model - 0.5*spread_model_error) > 0.001)*1

    extend_val[spread_model==-1]=-9
    extend_val[(spread_model==0) & (spread_model_error==0)]=-9
    extend_val[(spread_model==1) & (spread_model_error==1)]=-9
    return(extend_val)

def coadd_interp_array_y6a2():
    x   = [-1.        ,        -0.01,         0,          0.01,         1.0 ]
    y_1 = [3.0 + 0.005, 0.03 + 0.005, 0 + 0.005, -0.03 + 0.005, -3.0 + 0.005]
    y_2 = [1.0 + 0.003, 0.01 + 0.003, 0 + 0.003, -0.01 + 0.003, -1.0 + 0.003]
    y_3 = [-1. + 0.002,-0.01 + 0.002, 0 + 0.002,  0.01 + 0.002,  1.0 + 0.002]

    f_array = [scipy.interpolate.interp1d(x, y_1),
               scipy.interpolate.interp1d(x, y_2),
               scipy.interpolate.interp1d(x, y_3)]
    return x, f_array

def coadd_extended_class(spread_model,spread_model_error):
    """ DES DR2 extended classifier for stars and galaxies using coadd quantities.

        0 : High-confidence stars
        1 : Likely stars
        2 : Likely galaxies
        3 : High-confidence galaxies
       -9 : Failures

    Parameters
    ----------
    spread_model      : SExtractor spread model value in a band
    spread_model_error: SExtractor spread model error in a band
    
    Returns
    -------
    extended_class    : Extended classifier output
    """
    extend_val  = ((spread_model + 3.0*spread_model_error) > 0.005)*1
    extend_val += ((spread_model + 1.0*spread_model_error) > 0.003)*1
    extend_val += ((spread_model - 1.0*spread_model_error) > 0.002)*1

    extend_val[spread_model==-1]=-9
    extend_val[(spread_model==0) & (spread_model_error==0)]=-9
    extend_val[(spread_model==1) & (spread_model_error==1)]=-9
    return(extend_val)


def delve_extended_class(spread_model,spread_model_error):
    """Extended classifier for stars and galaxies.
        0 : High-confidence stars
        1 : Likely stars
        2 : Likely galaxies
        3 : High-confidence galaxies
       -9 : Failures

    Parameters
    ----------
    spread_model      : SExtractor spread model value in a band
    spread_model_error: SExtractor spread model value in a band
    
    Returns
    -------
    extended_class    : Extended classifier output
    """
    extend_val=((spread_model + 3*spread_model_error ) > 0.005)*1
    extend_val+=((spread_model + spread_model_error ) > 0.003)*1
    extend_val+=((spread_model - spread_model_error ) > 0.003)*1
    extend_val[spread_model==-1]=-9
    extend_val[(spread_model==0) & (spread_model_error==0)]=-9
    extend_val[(spread_model==1) & (spread_model_error==1)]=-9
    return(extend_val)

extended_class = delve_extended_class

def bdf_interp_array_y6a1():
    x = [-1.        ,  0.79891862,  0.90845217,  0.98558583,  1.05791208,
         1.13603715,  1.22479487,  1.33572223,  1.48983602,  1.74124395,
         2.43187589,  5.        ] 
    y_1 = [0.028, 0.028, 0.008, 0.   , 0.004, 0.012, 0.012, 0.004, 0.012,
           0.024, 0.04 , 0.04 ]
    y_2 = [-0.028, -0.028, -0.04 , -0.032, -0.036, -0.032, -0.028, -0.016,
           -0.012,  0.008,  0.016,  0.016]
    y_3 = [-0.1  , -0.1  , -0.1  , -0.1  , -0.1  , -0.1  , -0.1  , -0.1  ,
           -0.012,  0.008,  0.016,  0.016]

    f_array = [scipy.interpolate.interp1d(x, y_1),
               scipy.interpolate.interp1d(x, y_2),
               scipy.interpolate.interp1d(x, y_3)]
    return x, f_array

def bdf_extended_class_y6a1(data): 
    """
    Older version of the bdf extended classifier:
    https://cdcvs.fnal.gov/redmine/projects/des-y6/wiki/Y6A2_Extended_Classifier
    """
    x, f_array = bdf_interp_array_y6a1()

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
    return np.where(np.isfinite(ext), ext, -9)

def bdf_interp_array_y6a2():
    x = [-3.       ,  0.79891862,  0.90845217,  0.98558583,  1.05791208,
         1.13603715,  1.22479487,  1.33572223,  1.48983602,  1.74124395,
         2.43187589,  6.        ] 
    y_1 = [0.028, 0.028, 0.008, 0.   , 0.004, 0.012, 0.012, 0.004, 0.012,
           0.024, 0.04 , 0.04 ]
    y_2 = [-0.028, -0.028, -0.04 , -0.032, -0.036, -0.032, -0.028, -0.016,
           -0.012,  0.008,  0.016,  0.016]
    y_3 = [-0.1  , -0.1  , -0.1  , -0.1  , -0.1  , -0.1  , -0.1  , -0.1  ,
           -0.012,  0.008,  0.016,  0.016]
    y_4 = [0.252, 0.252, 0.188, 0.14 , 0.096, 0.104, 0.052, 0.048, 0.04 ,
           0.052, 0.088, 0.088]

    f_array = [scipy.interpolate.interp1d(x, y_1),
               scipy.interpolate.interp1d(x, y_2),
               scipy.interpolate.interp1d(x, y_3),
               scipy.interpolate.interp1d(x, y_4)]
    return x, f_array


def bdf_interp_array(data):
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

def bdf_extended_class(data): 
    x, f_array = bdf_interp_array_y6a2()

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

def bdf_extended_class_y6a2(data): 
    """
    Newer version of the bdf classifier (with high-purity galaxy sample).

    https://github.com/des-science/science_release/blob/master/object_classification/y6/y6a2_object_classification_fitvd.py
    """
    x, f_array = bdf_interp_array_y6a2()

    s2n_var = 'bdf_s2n' if 'bdf_s2n' in data.dtype.names else 'BDF_S2N'
    t_var = 'bdf_T' if 'bdf_T' in data.dtype.names else 'BDF_T'

    x_data = np.log10(data[s2n_var])
    x_data = np.where(np.isfinite(x_data), x_data, x[0])
    y_data = data[t_var]

    ext = np.tile(0, len(x_data))
    for f in f_array:
        selection = (y_data > f(x_data))
        ext += selection.astype(int)
    
    # Sentinel values
    selection  = np.isclose(data[t_var], -9.999e+09) 
    selection |= np.isclose(data[s2n_var], -9.999e+09) 
    selection |= (data[s2n_var] <= 0.)
    ext[selection] = -9

    return np.where(np.isfinite(ext), ext, -9)


def bdf_interp_array_y6gold():
    f_array = []
    x = [-3.       ,  0.79891862,  0.90845217,  0.98558583,  1.05791208,
         1.13603715,  1.22479487,  1.33572223,  1.48983602,  1.74124395,
         2.43187589,  6.        ] 
    y_1 = [0.028, 0.028, 0.008, 0.   , 0.004, 0.012, 0.012, 0.004, 0.012,
           0.024, 0.04 , 0.04 ]
    y_2 = [-0.028, -0.028, -0.04 , -0.032, -0.036, -0.032, -0.028, -0.016,
           -0.012,  0.008,  0.016,  0.016]
    f_array += [scipy.interpolate.interp1d(x, y_1)]
    f_array += [scipy.interpolate.interp1d(x, y_2)]

    #x   = [-3.0,      1.3,    1.6,   1.75, 1.917, 2.083, 2.583,  2.75, 6.000]
    #y_1 = [0.0280,  0.012, 0.0170, 0.0240, 0.028, 0.032, 0.040, 0.040, 0.040]
    #y_1 = [-0.010, -0.010, 0.0150, 0.0083, 0.010, 0.012, 0.016, 0.016, 0.016]
    #y_2 = [-0.028, -0.011, -0.003, 0.0083, 0.010, 0.012, 0.016, 0.016, 0.016]
    #y_3 = [-0.010, -0.010, 0.0083, 0.0083, 0.025, 0.025, 0.058, 0.058, 0.058]

    x   = [-3.0,  0.5,  1.0,   1.5 ,  2.0,  2.5,   3.0, 3.5,   6.0 ]
    y_1 = [-0.1, -0.1, -0.1,   0.03, 0.06, 0.06,  0.06, 0.06, 0.06 ]
    y_2 = [-0.1, -0.1, -0.1,  0.007, 0.03, 0.06,  0.28, 0.28, 0.28 ]
    #y_2 = [ 0.0,  0.0,  0.0,   0.01, 0.04, 0.07, 0.08, 0.08, 0.08 ]
    y_3 = [ 0.0,  0.0,  0.0,   0.01, 0.06, 0.15, 0.15, 0.15, 0.15 ]

    f_array += [scipy.interpolate.interp1d(x, y_1)]
        
    return x, f_array

def bdf_interp_array_y6gold():
    x = [-3.       ,  0.79891862,  0.90845217,  0.98558583,  1.05791208,
         1.13603715,  1.22479487,  1.33572223,  1.48983602,  1.74124395,
         2.43187589,  6.        ] 
    y_1 = [0.028, 0.028, 0.008, 0.   , 0.004, 0.012, 0.012, 0.004, 0.012,
           0.024, 0.04 , 0.04 ]
    #y_2 = [-0.028, -0.028, -0.04 , -0.032, -0.036, -0.032, -0.028, -0.016,
    #       -0.012,  0.008,  0.016,  0.016]
    y_2 = [-0.028, -0.028, -0.028, -0.028, -0.028, -0.028, -0.028, -0.012,
           0.005 ,  0.022,  0.04 ,  0.04]
    y_3 = [-0.1  , -0.1  , -0.1  , -0.1  , -0.1  , -0.1   , -0.1  , -0.1  ,
           -0.012,  0.008,  0.016,  0.016]
    y_4 = [0.252, 0.252, 0.188, 0.14 , 0.096, 0.104, 0.052, 0.048, 0.04 ,
           0.052, 0.088, 0.088]

    f_array = [scipy.interpolate.interp1d(x, y_1),
               scipy.interpolate.interp1d(x, y_2),
               scipy.interpolate.interp1d(x, y_3),
               scipy.interpolate.interp1d(x, y_4)]
    return x, f_array


def bdf_extended_class_y6gold(data): 
    """
    Newer version of the bdf classifier (with high-purity galaxy sample).

    https://github.com/des-science/science_release/blob/master/object_classification/y6/y6a2_object_classification_fitvd.py
    """
    x, f_array = bdf_interp_array_y6gold()

    s2n_var = 'bdf_s2n' if 'bdf_s2n' in data.dtype.names else 'BDF_S2N'
    t_var = 'bdf_T' if 'bdf_T' in data.dtype.names else 'BDF_T'

    x_data = np.log10(data[s2n_var])
    x_data = np.where(np.isfinite(x_data), x_data, x[0])
    y_data = data[t_var]

    ext = np.tile(0, len(x_data))
    for f in f_array:
        selection = (y_data > f(x_data))
        ext += selection.astype(int)
    
    # Sentinel values
    selection  = np.isclose(data[t_var], -9.999e+09) 
    selection |= np.isclose(data[s2n_var], -9.999e+09) 
    selection |= (data[s2n_var] <= 0.)
    ext[selection] = -9

    return np.where(np.isfinite(ext), ext, -9)


bdf_extended_class = bdf_extended_class_y6a2


def bdf_concentration_class(data):
    """ This is based on the HSC PDR1 star/galaxy classification 

    Roughly inspired by Coupon et al. 2017: https://arxiv.org/abs/1705.00622

    Returns
    -------
    sel : boolean HSC classification; 0 = stars; 1 = galaxies
    """
    # Calculate the concentration parameter
    concentration_i = data['PSF_MAG_APER_8_I'] - data['BDF_MAG_I']
    sel_stars = ( concentration_i < 0.04) | \
                ((concentration_i < 0.1) & (data['PSF_MAG_APER_8_I'] < 21.))

    return ~sel_stars


def wavg_extended_class(spread_model,spread_model_error):
    """ DES DR2 extended classifier for stars and galaxies using WAVG quantities.

        0 : High-confidence stars
        1 : Likely stars
        2 : Likely galaxies
        3 : High-confidence galaxies
       -9 : Failures

    Parameters
    ----------
    spread_model      : SExtractor spread model value in a band
    spread_model_error: SExtractor spread model error in a band
    
    Returns
    -------
    extended_class    : Extended classifier output
    """
    extend_val  = ((spread_model + 3.0*spread_model_error) > 0.005)*1
    extend_val += ((spread_model + 1.0*spread_model_error) > 0.003)*1
    extend_val += ((spread_model - 0.5*spread_model_error) > 0.001)*1

    extend_val[spread_model==-1]=-9
    extend_val[(spread_model==0) & (spread_model_error==0)]=-9
    extend_val[(spread_model==1) & (spread_model_error==1)]=-9
    return(extend_val)

def create_feature_matrix(datasets, names=None):
    PARAMS = ['SPREAD_MODEL_I','SPREADERR_MODEL_I', 
              'WAVG_SPREAD_MODEL_I', 'WAVG_SPREADERR_MODEL_I', 
              'BDF_T', 'BDF_S2N']
    if names is None: names = PARAMS

    X = np.zeros((0,len(names)),dtype=float)

    datasets = np.atleast_1d(datasets)
    for data in datasets:
        #data   = dataset['data']
        #data   = recfn.rec_append_fields(data,['LOG_BDF_S2N'],[np.log10(data['BDF_S2N'])])
        #data   = recfn.rec_append_fields(data,['CONCENTRATION_I'],
        #                                 [data['PSF_MAG_APER_8_I'] - data['BDF_MAG_I']])
     
        x = np.zeros((len(data),len(names)),dtype=float)
        for i,n in enumerate(names):
            x[:,i] = data[n]
        X = np.concatenate([X, x])

    return X

def create_feature_matrix(data, names=None):
    PARAMS = ['SPREAD_MODEL_I','SPREADERR_MODEL_I', 
              'WAVG_SPREAD_MODEL_I', 'WAVG_SPREADERR_MODEL_I', 
              'BDF_T', 'BDF_S2N']
    if names is None: names = PARAMS

    if isinstance(data,pd.DataFrame):
        X = data[names]
    else:
        X = pd.DataFrame(data[names])

    return X

def create_truth_matrix(datasets):
    y = np.concatenate([data for data in datasets])
    return y

def combo_random_forest_y6gold(filename,data):
    import joblib
    pipe = joblib.load(filename)
    X = create_feature_matrix(data,names=pipe[-1].feature_names_)
    return pipe.predict_proba(X)[:,1]

def y6gold_xgboost(filename,data):
    import joblib
    pipe = joblib.load(filename)
    X = create_feature_matrix(data,names=pipe[-1].feature_names_)
    return pipe.predict_proba(X)[:,1]



