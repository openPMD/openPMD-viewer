"""
This file is part of the openPMD-viewer.

It defines a number of helper functions that are used in main.py

Copyright 2015-2017, openPMD-viewer contributors
Authors: Remi Lehe, Richard Pausch
License: 3-Clause-BSD-LBNL
"""

import os
import numpy as np
import h5py
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
    - an array of integers which correspond to the iteration of each file
    """
    # Find all the files in the provided directory
    all_files = os.listdir(path_to_dir)

    # Select the hdf5 files
    iters_and_names = []
    for filename in all_files:
        # Use only the name that end with .h5 or .hdf5
        if filename[-3:] == '.h5' or filename[-5:] == '.hdf5':
            full_name = os.path.join(
                os.path.abspath(path_to_dir), filename)
            # extract all iterations from hdf5 file
            f = h5py.File(full_name, 'r')
            iterations = list(f['/data'].keys())
            f.close()
            # for each found iteration create list of tuples
            # (which can be sorted together)
            for key_iteration in iterations:
                iters_and_names.append((int(key_iteration), full_name))

    # Sort the list of tuples according to the iteration
    iters_and_names.sort()
    # Extract the list of filenames and iterations
    filenames = [name for (it, name) in iters_and_names]
    iterations = np.array([it for (it, name) in iters_and_names])

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


def fit_bins_to_grid( hist_size, grid_size, grid_range ):
    """
    Given a tentative number of bins `hist_size` for a histogram over
    the range `grid_range`, return a modified number of bins `hist_size`
    and a modified range `hist_range` so that the spacing of the histogram
    bins is an integer multiple (or integer divisor) of the grid spacing.

    Parameters:
    ----------
    hist_size: integer
        The number of bins in the histogram along the considered direction

    grid_size: integer
        The number of cells in the grid

    grid_range: list of floats (in meters)
        The extent of the grid

    Returns:
    --------
    hist_size: integer
        The new number of bins

    hist_range: list of floats (in microns)
        The new range of the histogram
    """
    # The new histogram range is the same as the grid range
    hist_range = grid_range

    # Calculate histogram tentative spacing, and grid spacing
    hist_spacing = ( hist_range[1] - hist_range[0] ) * 1. / hist_size
    grid_spacing = ( grid_range[1] - grid_range[0] ) * 1. / grid_size

    # Modify the histogram spacing, so that either:
    if hist_spacing >= grid_spacing:
        # - The histogram spacing is an integer multiple of the grid spacing
        hist_spacing = int( hist_spacing / grid_spacing ) * grid_spacing
    else:
        # - The histogram spacing is an integer divisor of the grid spacing
        hist_spacing = grid_spacing / int( grid_spacing / hist_spacing )

    # Get the corresponding new number of bins, and the new range
    hist_size = int( ( hist_range[1] - hist_range[0] ) / hist_spacing )
    hist_range[1] = hist_range[0] + hist_size * hist_spacing

    # Convert the range to microns (since this is how particle positions
    # are returned in the openPMD-viewer)
    hist_range = [ 1.e6 * hist_range[0], 1.e6 * hist_range[1] ]

    return( hist_size, hist_range )


def combine_cylindrical_components( Fr, Ft, theta, coord, info ):
    """
    Calculate the catesian field Fx or Fy,
    from the cylindrical components Fr and Ft.

    Parameters:
    -----------
    Fr, Ft: 3darrays or 2Darrays (depending on whether `theta` is None)
        Contains the value of the fields
    theta: float or None
        Indicates the angle of the plane in which Fr and Ft where taken
    coord: string
        Either 'x' or 'y' ; indicates which component to calculate
    info: FieldMetaInformation object
        Contains info on the coordinate system
    """
    if theta is not None:
        # Fr and Fr are 2Darrays
        assert (Fr.ndim == 2) and (Ft.ndim == 2)

        if coord == 'x':
            F = np.cos(theta) * Fr - np.sin(theta) * Ft
        elif coord == 'y':
            F = np.sin(theta) * Fr + np.cos(theta) * Ft
        # Revert the sign below the axis
        F[: int(F.shape[0] / 2)] *= -1

    else:
        # Fr, Ft are 3Darrays, info corresponds to Cartesian data
        assert (Fr.ndim == 3) and (Ft.ndim == 3)

        # Calculate cos(theta) and sin(theta) in the transverse Cartesian plane
        # while avoiding divisions by 0
        r = np.sqrt( info.x[:,np.newaxis]**2 + info.y[np.newaxis,:]**2 )
        inv_r = 1./np.where( r!=0, r, 1. )
        # The value `1.`` is a placeholder in the above (to avoid division by 0)
        # The lines below replace this placeholder value.
        cos = np.where( r!=0, info.x[:,np.newaxis]*inv_r, 1. )
        sin = np.where( r!=0, info.y[np.newaxis,:]*inv_r, 0. )
        if coord == 'x':
            F = cos[:,:,np.newaxis] * Fr - sin[:,:,np.newaxis] * Ft
        elif coord == 'y':
            F = sin[:,:,np.newaxis] * Fr + cos[:,:,np.newaxis] * Ft

    return F


def construct_3d_from_circ( F3d, Fcirc, x_array, y_array, modes,
    nx, ny, nz, nr, nmodes, inv_dr, rmax ):
    """
    Reconstruct the field from a quasi-cylindrical simulation (`Fcirc`), as
    a 3D cartesian array (`F3d`).
    """
    for ix in range(nx):
        x = x_array[ix]
        for iy in range(ny):
            y = y_array[iy]
            r = np.sqrt( x**2 + y**2 )
            ir = nr - 1 - int( (rmax - r) * inv_dr + 0.5 )
            # Handle out-of-bounds
            if ir < 0:
                ir = 0
            if ir >= nr:
                ir = nr-1
            # Loop over all modes and recontruct data
            if r == 0:
                expItheta = 1. + 0.j
            else:
                expItheta = (x+1.j*y)/r
            for im in range(nmodes):
                mode = modes[im]
                if mode==0:
                    F3d[ix, iy, :] += Fcirc[0, ir, :]
                else:
                    cos = (expItheta**mode).real
                    sin = (expItheta**mode).imag
                    F3d[ix, iy, :] += Fcirc[2*mode-1, ir, :]*cos \
                                    + Fcirc[2*mode, ir, :]*sin
