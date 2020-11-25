"""
This file is part of the openPMD-viewer.

It defines a set of helper data and functions which
are used by the other files.

Copyright 2015-2016, openPMD-viewer contributors
Authors: Remi Lehe, Axel Huebl
License: 3-Clause-BSD-LBNL
"""
import numpy as np


def get_data(series, record_component, i_slice=None, pos_slice=None,
             output_type=np.float64):
    """
    Extract the data from a (possibly constant) dataset
    Slice the data according to the parameters i_slice and pos_slice

    Parameters:
    -----------
    series: openpmd_api.Series
        An open, readable openPMD-api series object

    record_component: an openPMD.Record_Component

    pos_slice: int or list of int, optional
        Slice direction(s).
        When None, no slicing is performed

    i_slice: int or list of int, optional
       Indices of slices to be taken.

    output_type: a numpy type
       The type to which the returned array should be converted

    Returns:
    --------
    An np.ndarray (non-constant dataset) or a single double (constant dataset)
    """
    # For back-compatibility: Convert pos_slice and i_slice to
    # single-element lists if they are not lists (e.g. float
    # and int respectively).
    if pos_slice is not None and not isinstance(pos_slice, list):
        pos_slice = [pos_slice]
    if i_slice is not None and not isinstance(i_slice, list):
        i_slice = [i_slice]

    if pos_slice is None:
        data = record_component[()]
    else:
        # Get largest element of pos_slice
        max_pos = max(pos_slice)
        # Create list of indices list_index of type
        # [:, :, :, ...] where Ellipsis starts at max_pos + 1
        list_index = [np.s_[:]] * (max_pos + 2)
        list_index[max_pos + 1] = np.s_[...]
        # Fill list_index with elements of i_slice
        for count, dir_index in enumerate(pos_slice):
            list_index[dir_index] = i_slice[count]
        # Convert list_index into a tuple
        tuple_index = tuple(list_index)
        # Slice dset according to tuple_index
        data = record_component[tuple_index]

    series.flush()

    # Convert to the right type
    if data.dtype != output_type:
        data = data.astype( output_type )
    # Scale by the conversion factor
    if output_type in [ np.float64, np.float32, np.float16 ]:
        if record_component.unit_SI != 1.0:
            data *= record_component.unit_SI

    return data


def join_infile_path(*paths):
    """
    Join path components using '/' as separator.
    This method is defined as an alternative to os.path.join, which uses '\\'
    as separator in Windows environments and is therefore not valid to navigate
    within data files.

    Parameters:
    -----------
    *paths: all strings with path components to join

    Returns:
    --------
    A string with the complete path using '/' as separator.
    """
    # Join path components
    path = '/'.join(paths)
    # Correct double slashes, if any is present
    path = path.replace('//', '/')

    return path
