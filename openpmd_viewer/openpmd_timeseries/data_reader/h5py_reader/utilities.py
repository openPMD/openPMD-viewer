"""
This file is part of the openPMD-viewer.

It defines a set of helper data and functions which
are used by the other files.

Copyright 2015-2016, openPMD-viewer contributors
Authors: Remi Lehe, Axel Huebl
License: 3-Clause-BSD-LBNL
"""
import os
import h5py
import numpy as np


def list_files(path_to_dir):
    """
    Return a list of the hdf5 files in this directory,
    and a list of the corresponding iterations

    Parameter
    ---------
    path_to_dir : string
        The path to the directory where the hdf5 files are.

    Returns
    -------
    A tuple with:
    - an array of integers which correspond to the iteration of each file
    - a dictionary that matches iterations to the corresponding filename
    """
    # group based encoding?
    is_single_file = os.path.isfile(path_to_dir)

    if is_single_file:
        all_files = [path_to_dir]
    else:
        # Find all the files in the provided directory
        all_files = os.listdir(path_to_dir)

    # Select the hdf5 files, and fill dictionary of correspondence
    # between iterations and files
    iteration_to_file = {}
    for filename in all_files:
        # Use only the name that end with .h5 or .hdf5
        if filename.endswith('.h5') or filename.endswith('.hdf5'):
            if is_single_file:
                full_name = filename
            else:
                full_name = os.path.join(
                    os.path.abspath(path_to_dir), filename)
            # extract all iterations from hdf5 file
            f = h5py.File(full_name, 'r')
            iterations = list(f['/data'].keys())
            f.close()
            # Add iterations to dictionary
            for key_iteration in iterations:
                iteration_to_file[ int(key_iteration) ] = full_name

    # Extract iterations and sort them
    iterations = np.array( sorted( list( iteration_to_file.keys() ) ) )

    return iterations, iteration_to_file


def is_scalar_record(record):
    """
    Determine whether a record is a scalar record or a vector record

    Parameter
    ---------
    record: an h5py Dataset or an h5py Group

    Return
    ------
    A boolean indicating whether the record is scalar
    """
    scalar = False
    if 'value' in record.attrs:
        scalar = True
    elif isinstance(record, h5py.Dataset):
        scalar = True

    return(scalar)


def get_data(dset, i_slice=None, pos_slice=None, output_type=None):
    """
    Extract the data from a (possibly constant) dataset
    Slice the data according to the parameters i_slice and pos_slice

    Parameters:
    -----------
    dset: an h5py.Dataset or h5py.Group (when constant)
        The object from which the data is extracted

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
    # Case of a constant dataset
    if isinstance(dset, h5py.Group):
        shape = dset.attrs['shape']
        # Restrict the shape if slicing is enabled
        if pos_slice is not None:
            shape = [ x for index, x in enumerate(shape) if
                      index not in pos_slice ]
        # Create the corresponding dataset
        data = dset.attrs['value'] * np.ones(shape)

    # Case of a non-constant dataset
    elif isinstance(dset, h5py.Dataset):
        if pos_slice is None:
            data = dset[...]
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
            data = dset[tuple_index]

    # Convert to the right type
    if (output_type is not None) and (data.dtype != output_type):
        data = data.astype( output_type )
    # Scale by the conversion factor
    if np.issubdtype(data.dtype, np.floating) or \
        np.issubdtype(data.dtype, np.complexfloating):
        if dset.attrs['unitSI'] != 1.0:
            data *= dset.attrs['unitSI']

    return(data)


def get_shape(dset):
    """
    Extract the shape of a (possibly constant) dataset

    Parameters:
    -----------
    dset: an h5py.Dataset or h5py.Group (when constant)
        The object whose shape is extracted

    Returns:
    --------
    A tuple corresponding to the shape
    """
    # Case of a constant dataset
    if isinstance(dset, h5py.Group):
        shape = dset.attrs['shape']
    # Case of a non-constant dataset
    elif isinstance(dset, h5py.Dataset):
        shape = dset.shape

    return(shape)


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
