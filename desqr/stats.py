#!/usr/bin/env python
"""
Statistical functions
"""
__author__ = "Alex Drlica-Wagner"

import numpy as np

def confusion_matrix(class_pos,class_neg,true_pos,true_neg,mag=None,bins=None):
    """
    Define a confusion matrix for statistical classification. For more
    information on this object, see the wikipedia page here:
    https://en.wikipedia.org/wiki/Confusion_matrix

    Parameters:
    -----------
    class_pos : Boolean array of objects classified as positive
    class_neg : Boolean array of objects classified as negative
    true_pos  : Boolean array of objects that are true positive
    true_neg  : Boolean array of objects that are true negative
    mag       : A parameter over which the calculate statistics
    bins      : Bins of mag to use to calculate statistics

    Returns:
    --------
    matrix    : A confusion matrix dict (see reference for keys)
    """
    # See https://en.wikipedia.org/wiki/Sensitivity_and_specificity
    # and https://en.wikipedia.org/wiki/Confusion_matrix

    scalar=False
    if mag is None or np.isscalar(class_pos):
        scalar = True
    if mag is None:
        mag = np.ones_like(class_pos,dtype=float)
    if bins is None:
        bins = [np.min(mag),np.max(mag)]

    # True Positive (TP)
    TP,junk = np.histogram(mag,bins,weights=(class_pos&true_pos).astype(float))
    # True Negative (TN)
    TN,junk = np.histogram(mag,bins,weights=(class_neg&true_neg).astype(float))
    # False Positive (FP)
    FP,junk = np.histogram(mag,bins,weights=(class_pos&true_neg).astype(float))
    # False Negative (FN)
    FN,junk = np.histogram(mag,bins,weights=(class_neg&true_pos).astype(float))

    # These are generally called:
    # P = true_npos, N = true_nneg
    class_npos,j = np.histogram(mag,bins,weights=class_pos.astype(float))
    class_nneg,j = np.histogram(mag,bins,weights=class_neg.astype(float))
    true_npos,j  = np.histogram(mag,bins,weights=true_pos.astype(float) )
    true_nneg,j  = np.histogram(mag,bins,weights=true_neg.astype(float) )

    # True Positive Rate (TPR) aka Sensitivity aka Efficiency aka Completeness
    TPR = TP/true_npos
    # False negative rate (FNR) aka Miss rate
    FNR = FN/true_npos
    # True negative rate (TNR) aka Specificity (SPC)
    TNR = TN/true_nneg
    # False positive rate (FPR) aka Fall-out
    FPR = FP/true_nneg

    # False Discovery Rate (FDR) aka contamination
    FDR = FP/class_npos
    # Positive predictive value (PPV) aka precision aka purity
    PPV = TP/class_npos
    # False omission rate (FOR)
    FOR = FN/class_nneg
    # Negative predictive value (NPV)
    NPV = TN/class_nneg

    #Accuracy (ACC)
    ACC = (TP + TN)/(true_npos + true_nneg)
    # Prevalence (PRV)
    PRV = class_npos/(true_npos + true_nneg)

    # Positive likelihood ratio (PLR)
    PLR = TPR/FPR
    # Negative likelihood ration (NLR)
    NLR = FNR/TNR

    # Return a dictionary of all values
    ret = dict(TP=TP,TN=TN,FP=FP,FN=FN)
    ret.update(TPR=TPR, FNR=FNR, TNR=TNR, FPR=FPR)
    ret.update(PPV=PPV, FDR=FDR, FOR=FOR, NPV=NPV)
    ret.update(ACC=ACC,PRV=PRV,PLR=PLR,NLR=NLR)
    ret.update(P=true_npos,N=true_nneg)
    # Also define lowercase aliases
    for key,val in list(ret.items()):
        ret[key.lower()] = val

    ret.update(true_npos=true_npos,true_nneg=true_nneg)
    ret.update(class_npos=class_npos,class_nneg=class_nneg)
    ret.update(mag=mag, bins=bins)

    if scalar:
        for key,val in ret.items():
            ret[key] = np.asscalar(val)

    return ret
