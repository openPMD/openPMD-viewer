"""
This file is part of the openPMD-viewer.

It defines functions that can read the fields from an HDF5 file.

Copyright 2015-2016, openPMD-viewer contributors
Author: Remi Lehe
License: 3-Clause-BSD-LBNL
"""

import h5py
import numpy as np
from .utilities import get_shape, get_data, join_infile_path
from openpmd_viewer.openpmd_timeseries.field_metainfo import FieldMetaInformation
from openpmd_viewer.openpmd_timeseries.utilities import construct_3d_from_circ


def read_field_cartesian( filename, iteration, field, coord, axis_labels,
                          slice_relative_position, slice_across ):
    """
    Extract a given field from an HDF5 file in the openPMD format,
    when the geometry is cartesian (1d, 2d or 3d).

    Parameters
    ----------
    filename : string
       The absolute path to the HDF5 file

    iteration : int
        The iteration at which to obtain the data

    field : string, optional
       Which field to extract

    coord : string, optional
       Which component of the field to extract

    axis_labels: list of strings
       The name of the dimensions of the array (e.g. ['x', 'y', 'z'])

    slice_across : list of str or None
       Direction(s) across which the data should be sliced
       Elements can be:
         - 1d: 'z'
         - 2d: 'x' and/or 'z'
         - 3d: 'x' and/or 'y' and/or 'z'
       Returned array is reduced by 1 dimension per slicing.

    slice_relative_position : list of float or None
       Number(s) between -1 and 1 that indicate where to slice the data,
       along the directions in `slice_across`
       -1 : lower edge of the simulation box
       0 : middle of the simulation box
       1 : upper edge of the simulation box

    Returns
    -------
    A tuple with
       F : a ndarray containing the required field
       info : a FieldMetaInformation object
       (contains information about the grid; see the corresponding docstring)
    """
    # Open the HDF5 file
    dfile = h5py.File( filename, 'r' )
    # Extract the dataset and and corresponding group
    if coord is None:
        field_path = field
    else:
        field_path = join_infile_path( field, coord )
    group, dset = find_dataset( dfile, iteration, field_path )

    # Dimensions of the grid
    shape = list( get_shape( dset ) )
    grid_spacing = list( group.attrs['gridSpacing'] )
    global_offset = list( group.attrs['gridGlobalOffset'] )

    # Slice selection
    if slice_across is not None:
        # Get the integer that correspond to the slicing direction
        list_slicing_index = []
        list_i_cell = []
        for count, slice_across_item in enumerate(slice_across):
            slicing_index = axis_labels.index(slice_across_item)
            list_slicing_index.append(slicing_index)
            # Number of cells along the slicing direction
            n_cells = shape[ slicing_index ]
            # Index of the slice (prevent stepping out of the array)
            i_cell = int( 0.5 * (slice_relative_position[count] + 1.) * n_cells )
            i_cell = max( i_cell, 0 )
            i_cell = min( i_cell, n_cells - 1)
            list_i_cell.append(i_cell)

        # Remove metainformation relative to the slicing index
        # Successive pops starting from last coordinate to slice
        shape = [ x for index, x in enumerate(shape)
                  if index not in list_slicing_index ]
        grid_spacing = [ x for index, x in enumerate(grid_spacing)
                         if index not in list_slicing_index ]
        global_offset = [ x for index, x in enumerate(global_offset)
                          if index not in list_slicing_index ]
        axis_labels = [ x for index, x in enumerate(axis_labels)
                         if index not in list_slicing_index ]

        axes = { i: axis_labels[i] for i in range(len(axis_labels)) }
        # Extract data
        F = get_data( dset, list_i_cell, list_slicing_index )
        info = FieldMetaInformation( axes, shape, grid_spacing, global_offset,
                group.attrs['gridUnitSI'], dset.attrs['position'] )
    else:
        F = get_data( dset )
        axes = { i: axis_labels[i] for i in range(len(axis_labels)) }
        info = FieldMetaInformation( axes, F.shape,
            group.attrs['gridSpacing'], group.attrs['gridGlobalOffset'],
            group.attrs['gridUnitSI'], dset.attrs['position'] )

    # Close the file
    dfile.close()
    return( F, info )


def read_field_circ( filename, iteration, field, coord,
                     slice_relative_position, slice_across, m=0, theta=0.,
                     max_resolution_3d=None ):
    """
    Extract a given field from an HDF5 file in the openPMD format,
    when the geometry is thetaMode

    Parameters
    ----------
    filename : string
       The absolute path to the HDF5 file

    iteration : int
        The iteration at which to obtain the data

    field : string, optional
       Which field to extract

    coord : string, optional
       Which component of the field to extract

    m : int or string, optional
       The azimuthal mode to be extracted

    theta : float or None
       Angle of the plane of observation with respect to the x axis
       If `theta` is not None, then this function returns a 2D array
       corresponding to the plane of observation given by `theta` ;
       otherwise it returns a full 3D Cartesian array

    slice_across : list of str or None
       Direction(s) across which the data should be sliced
       Elements can be 'r' and/or 'z'
       Returned array is reduced by 1 dimension per slicing.

    slice_relative_position : list of float or None
       Number(s) between -1 and 1 that indicate where to slice the data,
       along the directions in `slice_across`
       -1 : lower edge of the simulation box
       0 : middle of the simulation box
       1 : upper edge of the simulation box

    max_resolution_3d : list of int or None
        Maximum resolution that the 3D reconstruction of the field (when
        `theta` is None) can have. The list should contain two values,
        e.g. `[200, 100]`, indicating the maximum longitudinal and transverse
        resolution, respectively. This is useful for performance reasons,
        particularly for 3D visualization.

    Returns
    -------
    A tuple with
       F : a 3darray or 2darray containing the required field,
           depending on whether `theta` is None or not
       info : a FieldMetaInformation object
       (contains information about the grid; see the corresponding docstring)
    """
    # Open the HDF5 file
    dfile = h5py.File( filename, 'r' )
    # Extract the dataset and and corresponding group
    if coord is None:
        field_path = field
    else:
        field_path = join_infile_path( field, coord )
    group, dset = find_dataset( dfile, iteration, field_path )

    # Extract the metainformation
    coord_labels = {ii: coord.decode() for (ii,coord) in
                                enumerate(group.attrs['axisLabels'])}
    if coord_labels[0] == 'r':
        rz_switch = False  # fastest varying index is z
        Nm, Nr, Nz = get_shape( dset )
        N_pair = (Nr, Nz)
    else:
        rz_switch = True  # fastest varying index is r
        Nm, Nz, Nr = get_shape( dset )
        N_pair = (Nz, Nr)
    info = FieldMetaInformation( coord_labels, N_pair,
        group.attrs['gridSpacing'], group.attrs['gridGlobalOffset'],
        group.attrs['gridUnitSI'], dset.attrs['position'], thetaMode=True )

    # Convert to a 3D Cartesian array if theta is None
    if theta is None:

        # Get cylindrical info
        rmax = info.rmax
        inv_dr = 1./info.dr
        Fcirc = get_data( dset )  # (Extracts all modes)
        if m == 'all':
            modes = [ mode for mode in range(0, int(Nm / 2) + 1) ]
        else:
            modes = [ m ]
        modes = np.array( modes, dtype='int' )
        nmodes = len(modes)

        # If necessary, reduce resolution of 3D reconstruction
        if max_resolution_3d is not None:
            max_res_lon, max_res_transv = max_resolution_3d
            if Nz > max_res_lon:
                # Calculate excess of elements along z
                excess_z = int(np.round(Nz/max_res_lon))
                # Preserve only one every excess_z elements
                if not rz_switch:
                    Fcirc = Fcirc[:, :, ::excess_z]
                else:
                    Fcirc = Fcirc[:, ::excess_z, :]
                # Update info accordingly
                info.z = info.z[::excess_z]
                info.dz = info.z[1] - info.z[0]
            if Nr > max_res_transv/2:
                # Calculate excess of elements along r
                excess_r = int(np.round(Nr/(max_res_transv/2)))
                # Preserve only one every excess_r elements
                if not rz_switch:
                    Fcirc = Fcirc[:, ::excess_r, :]
                else:
                    Fcirc = Fcirc[:, :, ::excess_r]
                # Update info and necessary parameters accordingly
                info.r = info.r[::excess_r]
                info.dr = info.r[1] - info.r[0]
                inv_dr = 1./info.dr
                # Update Nr after reducing radial resolution.
                if not rz_switch:
                    Nr = Fcirc.shape[1]
                else:
                    Nr = Fcirc.shape[2]

        # Convert cylindrical data to Cartesian data
        info._convert_cylindrical_to_3Dcartesian()
        nx, ny, nz = len(info.x), len(info.y), len(info.z)
        F_total = np.zeros( (nx, ny, nz) )
        construct_3d_from_circ( F_total, Fcirc, info.x, info.y, modes,
            nx, ny, nz, Nr, nmodes, inv_dr, rmax, rz_switch=rz_switch )

    else:

        # Extract the modes and recombine them properly
        if not rz_switch:
            F_total = np.zeros( (2 * Nr, Nz ) )
        else:
            F_total = np.zeros( (Nz, 2 * Nr ) )
        if m == 'all':
            # Sum of all the modes
            # - Prepare the multiplier arrays
            mult_above_axis = [1]
            mult_below_axis = [1]
            for mode in range(1, int(Nm / 2) + 1):
                cos = np.cos( mode * theta )
                sin = np.sin( mode * theta )
                mult_above_axis += [cos, sin]
                mult_below_axis += [ (-1) ** mode * cos, (-1) ** mode * sin ]
            mult_above_axis = np.array( mult_above_axis )
            mult_below_axis = np.array( mult_below_axis )
            # - Sum the modes
            F = get_data( dset )  # (Extracts all modes)
            if not rz_switch:
                F_total[Nr:, :] = np.tensordot( mult_above_axis,
                                                F, axes=(0, 0) )[:, :]
                F_total[:Nr, :] = np.tensordot( mult_below_axis,
                                                F, axes=(0, 0) )[::-1, :]
            else:
                F_total[:, Nr:] = np.tensordot( mult_above_axis,
                                                F, axes=(0, 0) )[:, :]
                F_total[:, :Nr] = np.tensordot( mult_below_axis,
                                                F, axes=(0, 0) )[:, ::-1]
        elif m == 0:
            # Extract mode 0
            F = get_data( dset, 0, 0 )
            if not rz_switch:
                F_total[Nr:, :] = F[:, :]
                F_total[:Nr, :] = F[::-1, :]
            else:
                F_total[:, Nr:] = F[:, :]
                F_total[:, :Nr] = F[:, ::-1]

        else:
            # Extract higher mode
            cos = np.cos( m * theta )
            sin = np.sin( m * theta )
            F_cos = get_data( dset, 2 * m - 1, 0 )
            F_sin = get_data( dset, 2 * m, 0 )
            F = cos * F_cos + sin * F_sin
            if not rz_switch:
                F_total[Nr:, :] = F[:, :]
                F_total[:Nr, :] = (-1) ** m * F[::-1, :]
            else:
                F_total[:, Nr:] = F[:, :]
                F_total[:, :Nr] = (-1) ** m * F[:, ::-1]

    # Perform slicing if needed
    if slice_across is not None:
        # Slice field and clear metadata
        inverted_axes_dict = {info.axes[key]: key for key in info.axes.keys()}
        for count, slice_across_item in enumerate(slice_across):
            slicing_index = inverted_axes_dict[slice_across_item]
            coord_array = getattr( info, slice_across_item )
            # Number of cells along the slicing direction
            n_cells = len(coord_array)
            # Index of the slice (prevent stepping out of the array)
            i_cell = int( 0.5 * (slice_relative_position[count] + 1.) * n_cells )
            i_cell = max( i_cell, 0 )
            i_cell = min( i_cell, n_cells - 1)
            F_total = np.take( F_total, [i_cell], axis=slicing_index )
        F_total = np.squeeze(F_total)
        # Remove the sliced labels from the FieldMetaInformation
        for slice_across_item in slice_across:
            info._remove_axis(slice_across_item)

    # Close the file
    dfile.close()

    return( F_total, info )


def find_dataset( dfile, iteration, field_path ):
    """
    Extract the dataset that corresponds to field_path,
    and the corresponding group

    (In the case of scalar records, the group and the dataset are identical.
    In the case of vector records, the group contains all the components
    and the dataset corresponds to one given component.)

    Parameters
    ----------
    dfile: an h5Py.File object
       The file from which to extract the dataset

    field_path : string
       The relative path to the requested field, from the openPMD meshes path
       (e.g. 'rho', 'E/r', 'B/x')

    Returns
    -------
    A tuple with:
    - an h5py.Group object
    - an h5py.Dataset object
    """
    # Find the meshes path
    base_path = '/data/{0}'.format( iteration )
    relative_meshes_path = dfile.attrs["meshesPath"].decode()

    # Get the proper dataset
    full_field_path = join_infile_path(
        base_path, relative_meshes_path, field_path )
    dset = dfile[ full_field_path ]
    # Get the proper group
    group_path = field_path.split('/')[0]
    full_group_path = join_infile_path(
        base_path, relative_meshes_path, group_path )
    group = dfile[ full_group_path ]

    return( group, dset )


def get_grid_parameters( filename, iteration, avail_fields, metadata ):
    """
    Return the parameters of the spatial grid (grid size and grid range)
    in two dictionaries

    Parameters:
    -----------
    filename : string
       The absolute path to the HDF5 file

    iteration : int
        The iteration at which to obtain the data

    avail_fields: list
       A list of the available fields
       e.g. ['B', 'E', 'rho']

    metadata: dictionary
      A dictionary whose keys are the fields of `avail_fields` and
      whose values are dictionaries that contain metadata (e.g. geometry)

    Returns:
    --------
    A tuple with `grid_size_dict` and `grid_range_dict`
    Both objects are dictionaries, with their keys being the labels of the axis
    of the grid (e.g. 'x', 'y', 'z')
    The values of `grid_size_dict` are the number of gridpoints along each axis
    The values of `grid_range_dict` are lists of two floats, which correspond
    to the min and max of the grid, along each axis.
    """
    # Open the HDF5 file
    dfile = h5py.File( filename, 'r' )
    # Pick field with the highest dimensionality ('3d'>'thetaMode'>'2d')
    # (This function is for the purpose of histogramming the particles;
    # in this case, the highest dimensionality ensures that more particle
    # quantities can be properly histogrammed.)
    geometry_ranking = {'1dcartesian': 0, '2dcartesian': 1,
                        'thetaMode': 2, '3dcartesian': 3}
    fields_ranking = [ geometry_ranking[ metadata[field]['geometry'] ]
                        for field in avail_fields ]
    index_best_field = fields_ranking.index( max(fields_ranking) )
    field_name = avail_fields[ index_best_field ]

    # Get the corresponding field data
    group, dset = find_dataset( dfile, iteration, field_name )
    if metadata[field_name]['type'] == 'vector':
        # For field vector, extract the first coordinate, to get the dataset
        first_coord = next(iter(group.keys()))
        dset = group[first_coord]

    # Extract relevant quantities
    labels = group.attrs['axisLabels']
    grid_spacing = group.attrs['gridSpacing'] * group.attrs['gridUnitSI']
    grid_offset = group.attrs['gridGlobalOffset'] * group.attrs['gridUnitSI']
    grid_size = dset.shape
    if metadata[field_name]['geometry'] == 'thetaMode':
        # In thetaMode: skip the first number of dset.shape, as this
        # corresponds to the number of modes
        grid_size = dset.shape[1:]

    # Build the dictionaries grid_size_dict and grid_range_dict
    grid_size_dict = {}
    grid_range_dict = {}
    for i in range(len(labels)):
        coord = labels[i].decode()
        grid_size_dict[coord] = grid_size[i]
        grid_range_dict[coord] = \
            [ grid_offset[i], grid_offset[i] + grid_size[i] * grid_spacing[i] ]
    # Close the file
    dfile.close()
    return( grid_size_dict, grid_range_dict )
