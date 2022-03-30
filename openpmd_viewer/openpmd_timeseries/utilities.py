"""
This file is part of the openPMD-viewer.

It defines a number of helper functions that are used in main.py

Copyright 2015-2017, openPMD-viewer contributors
Authors: Remi Lehe, Richard Pausch
License: 3-Clause-BSD-LBNL
"""

import copy
import math
import numpy as np
from .numba_wrapper import jit

def sanitize_slicing(slice_across, slice_relative_position):
    """
    Return standardized format for `slice_across` and `slice_relative_position`:
    - either `slice_across` and `slice_relative_position` are both `None` (no slicing)
    - or `slice_across` and `slice_relative_position` are both lists,
    with the same number of elements

    Parameters
    ----------
    slice_relative_position : float, or list of float, or None

    slice_across : str, or list of str, or None
       Direction(s) across which the data should be sliced
    """
    # Skip None and empty lists
    if slice_across is None or slice_across == []:
        return None, None

    # Convert to lists
    if not isinstance(slice_across, list):
        slice_across = [slice_across]
    if slice_relative_position is None:
        slice_relative_position = [0]*len(slice_across)
    if not isinstance(slice_relative_position, list):
        slice_relative_position = [slice_relative_position]
    # Check that the length are matching
    if len(slice_across) != len(slice_relative_position):
        raise ValueError(
            'The argument `slice_relative_position` is erroneous: \nIt should have'
            'the same number of elements as `slice_across`.')

    # Return a copy. This is because the rest of the `openPMD-viewer` code
    # sometimes modifies the objects returned by `sanitize_slicing`.
    # Using a copy avoids directly modifying objects that the user may pass
    # to this function (and live outside of openPMD-viewer, e.g. directly in
    # a user's notebook)
    return copy.copy(slice_across), copy.copy(slice_relative_position)

def apply_selection(iteration, data_reader, data_list,
                    select, species, extensions):
    """
    Select the elements of each particle quantities in data_list,
    based on the selection rules in `select`

    Parameters
    ----------
    iteration: int
        The iteration at which to apply the selection

    data_reader: a DataReader object
        Contains the method that read particle data

    data_list: list of 1darrays
        A list of arrays with one element per macroparticle, that represent
        different particle quantities

    select: dict
        A dictionary of rules to select the particles
        'x' : [-4., 10.]   (Particles having x between -4 and 10)
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
        q = data_reader.read_species_data(
            iteration, species, quantity, extensions)
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


def try_array( L ):
    """
    Attempt to convert L to a single array.
    """
    try:
        # Stack the arrays
        return np.stack( L, axis=0 )
    except ValueError:
        # Do not stack
        return L


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

    grid_range: list of floats (in)
        The extent of the grid

    Returns:
    --------
    hist_size: integer
        The new number of bins

    hist_range: list of floats
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
        if coord == 'x':
            F = np.cos(theta) * Fr - np.sin(theta) * Ft
        elif coord == 'y':
            F = np.sin(theta) * Fr + np.cos(theta) * Ft
        # Revert the sign below the axis
        if info.axes[0] == 'r':
            F[ : int(F.shape[0]/2) ] *= -1
        elif (F.ndim == 2) and (info.axes[1] == 'r'):
            F[ : , : int(F.shape[1]/2) ] *= -1
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

@jit
def histogram_cic_1d( q1, w, nbins, bins_start, bins_end ):
    """
    Return an 1D histogram of the values in `q1` weighted by `w`,
    consisting of `nbins` evenly-spaced bins between `bins_start`
    and `bins_end`. Contribution to each bins is determined by the
    CIC weighting scheme (i.e. linear weights).
    """
    # Define various scalars
    bin_spacing = (bins_end-bins_start)/nbins
    inv_spacing = 1./bin_spacing
    n_ptcl = len(w)

    # Allocate array for histogrammed data
    hist_data = np.zeros( nbins, dtype=np.float64 )

    # Go through particle array and bin the data
    for i in range(n_ptcl):
        # Calculate the index of lower bin to which this particle contributes
        q1_cell = (q1[i] - bins_start) * inv_spacing
        i_low_bin = int( math.floor( q1_cell ) )
        # Calculate corresponding CIC shape and deposit the weight
        S_low = 1. - (q1_cell - i_low_bin)
        if (i_low_bin >= 0) and (i_low_bin < nbins):
            hist_data[ i_low_bin ] += w[i] * S_low
        if (i_low_bin + 1 >= 0) and (i_low_bin + 1 < nbins):
            hist_data[ i_low_bin + 1 ] += w[i] * (1. - S_low)

    return( hist_data )


@jit
def histogram_cic_2d( q1, q2, w,
    nbins_1, bins_start_1, bins_end_1,
    nbins_2, bins_start_2, bins_end_2 ):
    """
    Return an 2D histogram of the values in `q1` and `q2` weighted by `w`,
    consisting of `nbins_1` bins in the first dimension and `nbins_2` bins
    in the second dimension.
    Contribution to each bins is determined by the
    CIC weighting scheme (i.e. linear weights).
    """
    # Define various scalars
    bin_spacing_1 = (bins_end_1-bins_start_1)/nbins_1
    inv_spacing_1 = 1./bin_spacing_1
    bin_spacing_2 = (bins_end_2-bins_start_2)/nbins_2
    inv_spacing_2 = 1./bin_spacing_2
    n_ptcl = len(w)

    # Allocate array for histogrammed data
    hist_data = np.zeros( (nbins_1, nbins_2), dtype=np.float64 )

    # Go through particle array and bin the data
    for i in range(n_ptcl):

        # Calculate the index of lower bin to which this particle contributes
        q1_cell = (q1[i] - bins_start_1) * inv_spacing_1
        q2_cell = (q2[i] - bins_start_2) * inv_spacing_2
        i1_low_bin = int( math.floor( q1_cell ) )
        i2_low_bin = int( math.floor( q2_cell ) )

        # Calculate corresponding CIC shape and deposit the weight
        S1_low = 1. - (q1_cell - i1_low_bin)
        S2_low = 1. - (q2_cell - i2_low_bin)
        if (i1_low_bin >= 0) and (i1_low_bin < nbins_1):
            if (i2_low_bin >= 0) and (i2_low_bin < nbins_2):
                hist_data[ i1_low_bin, i2_low_bin ] += w[i]*S1_low*S2_low
            if (i2_low_bin+1 >= 0) and (i2_low_bin+1 < nbins_2):
                hist_data[ i1_low_bin, i2_low_bin+1 ] += w[i]*S1_low*(1.-S2_low)
        if (i1_low_bin+1 >= 0) and (i1_low_bin+1 < nbins_1):
            if (i2_low_bin >= 0) and (i2_low_bin < nbins_2):
                hist_data[ i1_low_bin+1, i2_low_bin ] += w[i]*(1.-S1_low)*S2_low
            if (i2_low_bin+1 >= 0) and (i2_low_bin+1 < nbins_2):
                hist_data[ i1_low_bin+1, i2_low_bin+1 ] += w[i]*(1.-S1_low)*(1.-S2_low)

    return( hist_data )


@jit
def construct_3d_from_circ( F3d, Fcirc, x_array, y_array, modes,
    nx, ny, nz, nr, nmodes, inv_dr, rmax, rz_switch = False ):
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

            # Calculate linear projection from ir and ir-1
            if ir>0:
                s0 = ir + 0.5 - r* inv_dr
                s1 = 1. - s0
                if not rz_switch:
                    Fcirc_proj = s1*Fcirc[:, ir, :] + s0*Fcirc[:, ir-1, :]
                else:
                    Fcirc_proj = s1*Fcirc[:, :, ir] + s0*Fcirc[:, :, ir-1]
            else:
                if not rz_switch:
                    Fcirc_proj = Fcirc[:, ir, :]
                else:
                    Fcirc_proj = Fcirc[:, :, ir]

            # Loop over all modes and recontruct data
            if r == 0:
                expItheta = 1. + 0.j
            else:
                expItheta = (x+1.j*y)/r

            for im in range(nmodes):
                mode = modes[im]
                if mode==0:
                    F3d[ix, iy, :] += Fcirc_proj[0, :]
                else:
                    cos = (expItheta**mode).real
                    sin = (expItheta**mode).imag
                    F3d[ix, iy, :] += Fcirc_proj[2*mode-1,:]*cos + \
                        Fcirc_proj[2*mode,:]*sin
