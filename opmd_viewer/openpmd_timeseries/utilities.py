"""
This file is part of the openPMD-viewer.

It defines a number of helper functions that are used in main.py

Copyright 2015-2016, openPMD-viewer contributors
Authors: Remi Lehe
License: 3-Clause-BSD-LBNL
"""

import os
import re
import numpy as np
from .data_reader.particle_reader import read_species_data


def list_h5_files(path_to_dir):
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
    - a list of strings which correspond to the absolute path of each file
    - a list of integers which correspond to the iteration of each file
    """
    # Find all the files in the provided directory
    all_files = os.listdir(path_to_dir)

    # Select the hdf5 files
    iters_and_names = []
    for filename in all_files:
        # Use only the name that end with .h5 or .hdf5
        if filename[-3:] == '.h5' or filename[-5:] == '.hdf5':
            # Extract the iteration, using regular expressions (regex)
            regex_match = re.search('(\d+).h[df]*5', filename)
            if regex_match is None:
                print('Ill-formated HDF5 file: %s\n File names should end with'
                      ' the iteration number, followed by ".h5"' % filename)
            else:
                iteration = int(regex_match.groups()[-1])
                full_name = os.path.join(
                    os.path.abspath(path_to_dir), filename)
                # Create list of tuples (which can be sorted together)
                iters_and_names.append((iteration, full_name))

    # Sort the list of tuples according to the iteration
    iters_and_names.sort()
    # Extract the list of filenames and iterations
    filenames = [name for (it, name) in iters_and_names]
    iterations = [it for (it, name) in iters_and_names]

    return(filenames, iterations)


def apply_selection(file_handle, data_list, select, species, extensions):
    """
    Select the elements of each particle quantities in data_list,
    based on the selection rules in `select`

    Parameters
    ----------
    file_handle: h5py.File object
        The HDF5 file from which to extract data

    data_list: list of 1darrays
        A list of arrays with one element per macroparticle, that represent
        different particle quantities

    select: dict
        A dictionary of rules to select the particles
        'x' : [-4., 10.]   (Particles having x between -4 and 10 microns)
        'ux' : [-0.1, 0.1] (Particles having ux between -0.1 and 0.1 mc)
        'uz' : [5., None]  (Particles with uz above 5 mc)

    species: string
       Name of the species being requested

    extensions: list of strings
        The extensions that the current OpenPMDTimeSeries complies with

    Returns
    -------
    A list of 1darrays that correspond to data_list, but were only the
    macroparticles that meet the selection rules are kept
    """
    # Create the array that determines whether the particle
    # should be selected or not.
    Ntot = len(data_list[0])
    select_array = np.ones(Ntot, dtype='bool')

    # Loop through the selection rules, and aggregate results in select_array
    for quantity in select.keys():
        q = read_species_data(file_handle, species, quantity, extensions)
        # Check lower bound
        if select[quantity][0] is not None:
            select_array = np.logical_and(
                select_array,
                q > select[quantity][0])
        # Check upper bound
        if select[quantity][1] is not None:
            select_array = np.logical_and(
                select_array,
                q < select[quantity][1])

    # Use select_array to reduce each quantity
    for i in range(len(data_list)):
        if len(data_list[i]) > 1:  # Do not apply selection on scalar records
            data_list[i] = data_list[i][select_array]

    return(data_list)
