"""
This file is part of the openPMD-viewer.

It defines functions that can read the fields from an HDF5 file.

Copyright 2020, openPMD-viewer contributors
Author: Axel Huebl
License: 3-Clause-BSD-LBNL
"""

import numpy as np
from .utilities import get_data
from openpmd_viewer.openpmd_timeseries.field_metainfo import FieldMetaInformation
from openpmd_viewer.openpmd_timeseries.utilities import construct_3d_from_circ


def read_field_cartesian( series, iteration, field_name, component_name,
                          axis_labels, slice_relative_position, slice_across ):
    """
    Extract a given field from a file in the openPMD format,
    when the geometry is cartesian (1d, 2d or 3d).

    Parameters
    ----------
    series: openpmd_api.Series
        An open, readable openPMD-api series object

    iteration: integer
        Iteration from which parameters should be extracted

    field_name : string, optional
       Which field to extract

    component_name : string, optional
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
    it = series.iterations[iteration]

    # Extract the dataset and and corresponding group
    field = it.meshes[field_name]
    if field.scalar:
        component = next(field.items())[1]
    else:
        component = field[component_name]

    # Dimensions of the grid
    shape = component.shape
    # FIXME here and in h5py reader, we need to invert the order on 'F'
    grid_spacing = field.grid_spacing
    global_offset = field.grid_global_offset
    grid_unit_SI = field.grid_unit_SI
    grid_position = component.position

    # Slice selection
    #   TODO put in general utilities
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
        F = get_data( series, component, list_i_cell, list_slicing_index )
        info = FieldMetaInformation( axes, shape, grid_spacing, global_offset,
                grid_unit_SI, grid_position )
    else:
        F = get_data( series, component )
        axes = { i: axis_labels[i] for i in range(len(axis_labels)) }
        info = FieldMetaInformation( axes, F.shape,
            grid_spacing, global_offset,
            grid_unit_SI, grid_position )

    return F, info


def read_field_circ( series, iteration, field_name, component_name,
                     slice_relative_position, slice_across, m=0, theta=0.,
                     max_resolution_3d=None ):
    """
    Extract a given field from a file in the openPMD format,
    when the geometry is thetaMode

    Parameters
    ----------
    series: openpmd_api.Series
        An open, readable openPMD-api series object

    iteration: integer
        Iteration from which parameters should be extracted

    field_name : string, optional
       Which field to extract

    component_name : string, optional
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
    it = series.iterations[iteration]

    # Extract the dataset and and corresponding group
    field = it.meshes[field_name]
    if field.scalar:
        component = next(field.items())[1]
    else:
        component = field[component_name]

    # Extract the metainformation
    #   FIXME here and in h5py reader, we need to invert the order on 'F' for
    #         grid spacing/offset/position

    coord_labels = {ii: coord for (ii, coord) in enumerate(field.axis_labels)}
    if coord_labels[0] == 'r':
        rz_switch = False
        Nm, Nr, Nz = component.shape
        N_pair = (Nr, Nz)
    else:
        rz_switch = True
        Nm, Nz, Nr = component.shape
        N_pair = (Nz, Nr)

    # Nm, Nr, Nz = component.shape
    info = FieldMetaInformation( coord_labels, N_pair,
        field.grid_spacing, field.grid_global_offset,
        field.grid_unit_SI, component.position, thetaMode=True )

    # Convert to a 3D Cartesian array if theta is None
    if theta is None:

        # Get cylindrical info
        rmax = info.rmax
        inv_dr = 1./info.dr
        Fcirc = get_data( series, component )  # (Extracts all modes)
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
            nx, ny, nz, Nr, nmodes, inv_dr, rmax, rz_switch=rz_switch)

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
            F = get_data( series, component )  # (Extracts all modes)
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
            F = get_data( series, component, 0, 0 )
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
            F_cos = get_data( series, component, 2 * m - 1, 0 )
            F_sin = get_data( series, component, 2 * m, 0 )
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

    return F_total, info


# FIXME this looks like it can be generalized from already read meta-data
def get_grid_parameters( series, iteration, avail_fields, metadata ):
    """
    Return the parameters of the spatial grid (grid size and grid range)
    in two dictionaries

    Parameters:
    -----------
    series: openpmd_api.Series
        An open, readable openPMD-api series object

    iteration: integer
        Iteration from which parameters should be extracted

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
    it = series.iterations[iteration]

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

    # Extract the dataset and and corresponding group
    # For field vector, extract the first component, to get the dataset
    field = it.meshes[field_name]
    component = next(field.items())[1]

    # Extract relevant quantities
    #   FIXME here and in h5py reader, we need to invert the order on 'F' for
    #         grid spacing/offset/position
    labels = field.axis_labels
    grid_spacing = np.array(field.grid_spacing) * field.grid_unit_SI
    grid_offset = np.array(field.grid_global_offset) * field.grid_unit_SI
    grid_size = component.shape
    if metadata[field_name]['geometry'] == 'thetaMode':
        # In thetaMode: skip the first number of dset.shape, as this
        # corresponds to the number of modes
        grid_size = component.shape[1:]

    # Build the dictionaries grid_size_dict and grid_range_dict
    grid_size_dict = {}
    grid_range_dict = {}
    for i in range(len(labels)):
        coord = labels[i]
        grid_size_dict[coord] = grid_size[i]
        grid_range_dict[coord] = \
            [ grid_offset[i], grid_offset[i] + grid_size[i] * grid_spacing[i] ]

    return grid_size_dict, grid_range_dict
