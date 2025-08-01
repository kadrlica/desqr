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

def wavg_extended_class_y6a2(spread_model,spread_model_error):
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


def wavg_extended_class_delve(spread_model,spread_model_error):
    """ DELVE DR3 extended classifier for stars and galaxies using WAVG quantities.

        0 : High-confidence stars
        1 : Likely stars
        2 : Likely galaxies
        3 : High-confidence galaxies
       -9 : Failures

    Parameters
    ----------
    spread_model      : SExtractor wavg spread model value in a band
    spread_model_error: SExtractor wavg spread model error in a band
    
    Returns
    -------
    extended_class    : Extended classifier output
    """
    extend_val  = ((spread_model + 3.0*spread_model_error) > 0.005)*1
    extend_val += ((spread_model + 1.0*spread_model_error) > 0.003)*1
    extend_val += ((spread_model - 0.5*spread_model_error) > 0.001)*1

    extend_val[spread_model==-99]=-9
    extend_val[spreaderr_model==-99]=-9
    return(extend_val)

wavg_extended_class_y6gold = wavg_extended_class_y6a2
wavg_extended_class = wavg_extended_class_y6a2

def coadd_interp_array_y6a2():
    x   = [-1.        ,        -0.01,         0,          0.01,         1.0 ]
    y_1 = [3.0 + 0.005, 0.03 + 0.005, 0 + 0.005, -0.03 + 0.005, -3.0 + 0.005]
    y_2 = [1.0 + 0.003, 0.01 + 0.003, 0 + 0.003, -0.01 + 0.003, -1.0 + 0.003]
    y_3 = [-1. + 0.002,-0.01 + 0.002, 0 + 0.002,  0.01 + 0.002,  1.0 + 0.002]

    f_array = [scipy.interpolate.interp1d(x, y_1),
               scipy.interpolate.interp1d(x, y_2),
               scipy.interpolate.interp1d(x, y_3)]
    return x, f_array

def coadd_extended_class_y6a2(spread_model,spread_model_error):
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

coadd_extended_class_y6gold = coadd_extended_class_y6a2
coadd_extended_class = coadd_extended_class_y6a2

def delve_extended_class(spread_model,spread_model_error):
    """Extended classifier for stars and galaxies. Used for DELVE DR1 & DR2

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
    """ BDF interpolation array. Includes new high completeness stellar sample.
    """

    x = [-3.     ,  0.79891862,  0.90845217,  0.98558583,  1.05791208,
         1.13603715,  1.22479487,  1.33572223,  1.48983602,  1.74124395,
         2.43187589,  6.     ] 
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

        0 : High-confidence stars
        1 : Likely stars
        2 : Likely galaxies
        3 : High-confidence galaxies
        4 : High-purity galaxies
       -9 : Failures

    """
    x, f_array = bdf_interp_array_y6gold()

    data = pd.DataFrame(data)

    s2n_var = 'bdf_s2n' if 'bdf_s2n' in data.columns else 'BDF_S2N'
    t_var = 'bdf_T' if 'bdf_T' in data.columns else 'BDF_T'

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
    selection |= (data[s2n_var] <= 0.0)
    selection |= (x_data > x[-1])
    ext[selection] = -9

    return np.where(np.isfinite(ext), ext, -9)


def bdf_extended_class_dr3gold(data): 
    """
    Star/galaxy classification criteria for DELVE DR3 based on DES Y6 Gold. More robust to spurious measurements.

        0 : High-confidence stars
        1 : Likely stars
        2 : Likely galaxies
        3 : High-confidence galaxies
        4 : High-purity galaxies
       -9 : Failures

    """
    # Use the gold interpolations
    x, f_array = bdf_interp_array_y6gold()

    data = pd.DataFrame(data)
    s2n_var = 'bdf_s2n' if 'bdf_s2n' in data.columns else 'BDF_S2N'
    t_var = 'bdf_T' if 'bdf_T' in data.columns else 'BDF_T'

    x_data = np.log10(data[s2n_var])
    y_data = data[t_var]

    # Sentinel values for spurious measurements
    sentinel  = np.isclose(data[t_var], -9.999e+09) 
    sentinel |= np.isclose(data[s2n_var], -9.999e+09) 
    sentinel |= (data[s2n_var] <= 0.0)
    sentinel |= (x_data < x[0]) | (x_data > x[-1])
    sentinel |= ~np.isfinite(x_data)
    
    # Restrict values to interpolation range (output sentinel)
    x_data = np.clip(x_data, x[0], x[-1])
    x_data[np.isnan(x_data)] = x[0]

    ext = np.tile(0, len(x_data))
    for f in f_array:
        selection = (y_data > f(x_data))
        ext += selection.astype(int)
    
    ext[sentinel] = -9

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



def y6gold_ext_mash(ext_fitvd, ext_wavg, ext_coadd):
    """ Create MASH of fitvd, wavg, and coadd classifiers """
    ext_mash = -9*np.ones(len(ext_fitvd),dtype=int)
    sel_coadd = ext_coadd > -9
    ext_mash[sel_coadd] = ext_coadd[sel_coadd]
    sel_wavg = ext_wavg > -9
    ext_mash[sel_wavg] = ext_wavg[sel_wavg]
    sel_fitvd = ext_fitvd > -9
    ext_mash[sel_fitvd] = ext_fitvd[sel_fitvd]
    return ext_mash


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

def y6gold_sklearn(filename,data):
    import joblib
    pipe = joblib.load(filename)
    X = create_feature_matrix(data,names=pipe[-1].feature_names_)
    return pipe.predict_proba(X)[:,1]

def xgb_extended_class_y6gold(xgb_pred): 
    """ XGBoost extended classification.

        0 : High-confidence stars
        1 : Likely stars
        2 : Likely galaxies
        3 : High-confidence galaxies
        4 : Ultra-pure galaxies
       -9 : Failures
    
    Parameters
    ----------
    xgb_pred : XGBoost classifier prediction pseudo-probability

    Returns
    -------
    eclass : Extended classification
    """
    #eclass_cuts = [0.8, 0.5, 0.1, 0.05]
    eclass_cuts = [0.865, 0.11, 0.045, 0.015]
    eclass = np.zeros(len(xgb_pred), dtype=int)

    for cut in eclass_cuts:
        eclass += xgb_pred < cut
    
    eclass[xgb_pred < 0] = -9
    return eclass
