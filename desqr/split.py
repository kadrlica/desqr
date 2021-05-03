"""
Split incorrectly matched catalog objects.

Adapted from Eric Neilsen's split_qcat:
/data/des60.b/data/neilsen/split_qcat/split_qcat.py
/data/des60.b/data/neilsen/split_qcat/split_qcat_original.py
"""
__author__ = "Eric Neilsen"

from collections import OrderedDict as odict
import pandas as pd
import numpy as np

__author__ = "Eric Neilsen"

SPLIT_FLAG = odict([
        ('single'         , 0   ), 
        ('doubled'        , 2**0),
        ('isolated double', 2**1), 
        ('multi'          , 2**2),
        ('split'          , 2**3),
        ('failed split'   , 2**4),
        ('merged'         , 2**5),
])

def split_doubles(epochs, kmeans_iterations=3, objid='QUICK_OBJECT_ID'):
    """Split a DataFrame of doubled objects into sub-objects

    The DataFrame is assumed to be prepared to have exactly two
    detections per object per exposure.

    The DataFrame requires the following columns:
      RA, DEC
    and a MultiIndex row index with the following levels:
      QUICK_OBJECT_ID, EXPNUM

    If SUB_ID or SPLIT_FLAG are present, they will be overwritten in the
    return.

    Parameters
    ----------
    epochs: the pandas DataFrame with the objects to split
    kmeans_iterations: iterations to refine sub-object locations

    Returns
    -------
    df: a copy of the pandas DataFrame supplied with two additonal columns:
        SUB_ID -- identifying which subobject each row is identified with
        SPLIT_FLAG -- a string containing "split" if the split was successful, 
          or "failed split" otherwise
    """
    # epochs assumed to be cleaned to objects and epochs with two
    # instances of each object in each epoch

    initial_epochs_cols = epochs.columns.values

    epochs.loc[:,'RA_RAD'] = np.radians(epochs.RA)
    epochs.loc[:,'DEC_RAD'] = np.radians(epochs.DEC)
    
    # Get the reference EXPNUM for each object
    obj_expnums = pd.Series(epochs.index.get_level_values('EXPNUM'),
                            index=epochs.index.get_level_values(objid))
    ref_expnums = obj_expnums.groupby(level=0).min()

    # Get coordinates for objects in reference image
    ref_coords = epochs.loc[zip(ref_expnums.index, ref_expnums),['RA_RAD','DEC_RAD']]
    ref_coords.columns = ['REF_RA_RAD', 'REF_DEC_RAD']
    
    # Split each reference pair of objects
    ref_coords.loc[:,'SUB_ID'] = ref_coords.groupby(level=(0,1)).cumcount()
    ref_coords.reset_index(inplace=True)
    ref_coords.drop(['EXPNUM'], axis=1, inplace=True)
    ref_coords.set_index([objid, 'SUB_ID'], inplace=True)

    for kmeans_iteration in range(kmeans_iterations):
        if kmeans_iteration > 0:
            # if this is not the first iteration,
            # revise reference coordinates based on the previous iteration
            ref_coord_input = obj_offsets.query("SPLIT_FLAG=='split'")
            ref_coord_input.reset_index(inplace=True)
            ref_coord_input.set_index([objid, 'SUB_ID'], inplace=True)

            # because of the spherical coordinate wrap-around, we cannot just
            # average RA and DEC to refine our reference center.
            # So, take the mean in cartesian coordinates.
            #
            # This will not work if the points are far away, but in this
            # application they should all be very close, so it should be fine.
            #
            # Note that below, Xprime = X/cos(DEC_RAD), Yprime=Y/cos(DEC_RAD)
            #
            # This will fail to calculate dec optimally exactly at the poles
            # 
            ref_coord_input.eval("Xprime=sin(RA_RAD)", inplace=True)
            ref_coord_input.eval("Yprime=cos(RA_RAD)", inplace=True)
            ref_coord_input.eval("Z=sin(DEC_RAD)", inplace=True)
            
            ref_coords = ref_coord_input.groupby(level=(0,1)).agg({'Xprime': 'mean',
                                                                   'Yprime': 'mean',
                                                                   'Z': 'mean'})
            ref_coords.eval("REF_DEC_RAD=arcsin(Z)", inplace=True)
            ref_coords.eval("REF_RA_RAD=arctan2(Xprime,Yprime)", inplace=True)
            ref_coords = ref_coords.loc[:, ['REF_RA_RAD', 'REF_DEC_RAD']]
            
        # Get distances of all detections to each set of reference coordinates
        epochs = epochs.reset_index()
        epochs.drop(['SUB_ID'], axis=1, inplace=True)
        ref_coords = ref_coords.reset_index()
        obj_offsets = epochs.merge(ref_coords, on=objid)
        # haversine formula
        obj_offsets.eval(
            'HAVANG = sin((DEC_RAD-REF_DEC_RAD)/2)**2 + cos(DEC_RAD)*cos(REF_DEC_RAD)*(sin((RA_RAD-REF_RA_RAD)/2)**2)',
            inplace=True)
            
        # The object is closest where cos(angle) is largest
        # so get the largest cos(angle) for each detection,
        # and select only object/subobject pairs with that angle
        obj_offsets.set_index([objid, 'EXPNUM', 'RA', 'DEC'], inplace=True)
        min_obj_offsets = obj_offsets.groupby(level=(0,1,2,3)).agg({'HAVANG': 'min'})
        min_obj_offsets.columns = ['MIN_HAVANG']
        obj_offsets = obj_offsets.merge(min_obj_offsets, left_index=True, right_index=True)
        obj_offsets = obj_offsets.query('MIN_HAVANG == HAVANG')

        # Flag detections of objects where both detections in the
        # same exposure have been idenitified with the same subobject
        obj_subobject_counts = obj_offsets.groupby(level=(0,1,2,3)).agg({'HAVANG': 'count'})
        obj_subobject_counts.columns = ['COUNTS']
        failed_split = obj_subobject_counts.query('COUNTS!=1').index

        obj_offsets.loc[:,'SPLIT_FLAG'] = 'split'
        if len(failed_split)>0:
            obj_offsets.loc[failed_split, ['SPLIT_FLAG']] == 'failed split'

        epochs = obj_offsets.reset_index().set_index([objid, 'EXPNUM'])
        epochs.drop(['REF_RA_RAD', 'REF_DEC_RAD', 'MIN_HAVANG'], axis=1, inplace=True)

    epochs = epochs.loc[:, initial_epochs_cols]

    return epochs
    
def split_qcat(qcat, objid='QUICK_OBJECT_ID'):
    """Split a DataFrame of detected objects into subobjects when necessary

    The DataFrame requires the following columns:
      QUICK_OBJECT_ID, EXPNUM, RA, DEC

    It creates two additional columns:
      - SUB_ID: a subobject-id identifying which sub-object any given row 
                corresponds to
      - SPLIT_FLAG: a string representing the state of the split:
        'single': a single, unsplit detection
        'isolated double': an object doubled in this exposure, but none of the others
            No attempt is made to split these; other instances are classified as 
            'single'
        'multi': an object that has more than two detections in an exposure
            No attempt is made to split these.
        'split': a successfully split pair
        'failed split': a pair that was not split successfully
        'merged': detected once in this exposure, but pairs in others

    Parameters
    ----------
      qcat: the pandas DataFrame with the objects to split

    Returns
    -------
    df: a copy of the pandas DataFrame supplied with two additonal columns:
        SUB_ID -- identifying which subobject each row is identified with
        SPLIT_FLAG -- a string containing "split" if the split was successful, 
          or "failed split" otherwise
    """

    epochs = qcat.copy(deep=True)
    epochs.loc[:,'orig_row_num'] = np.arange(len(epochs))
    epochs.reset_index(inplace=True)
    epochs.set_index([objid, 'EXPNUM'], inplace=True)
    epochs.sort_index(inplace=True)

    # Objects are considered single until set otherwise
    epochs.loc[:,'SPLIT_FLAG'] = 'single'
    epochs.loc[:,'SUB_ID'] = 0

    # Count the detections of each object in each epoch
    detections_by_epoch = epochs[['RA']].groupby(level=(0,1)).count()
    detections_by_epoch.columns = ['COUNTS']

    # Initialize SPLIT_FLAG for doubled objects
    # This will by overwritten later for successfully split objects
    doubled_in_epoch = detections_by_epoch.query('COUNTS==2').index
    if len(doubled_in_epoch) > 0:
        epochs.loc[doubled_in_epoch, 'SPLIT_FLAG'] = 'doubled'

    # Set SPLIT_FLAG for multiply detected objects (never to be updated)
    many_in_epoch = detections_by_epoch.query('COUNTS>2').index
    if len(many_in_epoch) > 0:
        epochs.loc[many_in_epoch, 'SPLIT_FLAG'] = 'multi'

    #
    # Mark object/epochs split in only one epoch
    #

    # count numbers of epochs in which exposures are doubled
    times_doubled = detections_by_epoch.loc[doubled_in_epoch, :].groupby(level=0).count()
    times_doubled.columns = ['TIMES_DOUBLED']

    # get index for object/epochs split only once, and mark SPLIT_FLAG
    objs_doubled_once = times_doubled.query('TIMES_DOUBLED<2').index
    single_doubles = detections_by_epoch.loc[pd.IndexSlice[objs_doubled_once, :], :].query('COUNTS==2').index
    if len(single_doubles)>0:
        epochs.loc[single_doubles, 'SPLIT_FLAG'] = 'isolated double'
    
    # Mark detections with one detection in an exposure, but doubled more than one other
    objs_doubled = times_doubled.query('TIMES_DOUBLED>1').index
    merged = epochs.loc[pd.IndexSlice[objs_doubled, :], :].query("SPLIT_FLAG=='single'").index
    if len(merged)>0:
        epochs.loc[merged, 'SPLIT_FLAG'] = 'merged'

    to_split = epochs.query("SPLIT_FLAG=='doubled'").copy()
    split = split_doubles(to_split,objid=objid)

    unsplit = epochs.query("SPLIT_FLAG!='doubled'")
    epochs = pd.concat([split, unsplit])
    epochs.sort_values(by='orig_row_num', inplace=True)
    epochs.drop(['orig_row_num'], axis = 1, inplace = True, errors = 'ignore')
    epochs.reset_index(inplace=True)
    if None not in qcat.index.names:
        epochs.set_index(qcat.index.names, inplace=True)
    
    return epochs
