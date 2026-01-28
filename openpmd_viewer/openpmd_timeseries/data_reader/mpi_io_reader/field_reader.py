"""
MPI-enabled parallel versions of io_reader field_reader functions.

Copyright 2020, openPMD-viewer contributors
License: 3-Clause-BSD-LBNL
"""
import numpy as np
from .utilities import _read_field_portion, _read_field_circ_portion
from ..io_reader.utilities import get_data, chunk_to_slice
from ...data_order import RZorder, order_error_msg
from openpmd_viewer.openpmd_timeseries.field_metainfo import FieldMetaInformation
from openpmd_viewer.openpmd_timeseries.utilities import construct_3d_from_circ


def read_field_cartesian(series, iteration, field_name, component_name,
                        axis_labels, slice_relative_position, slice_across, comm):
    """
    Extract a given field from a file in the openPMD format,
    when the geometry is cartesian (1d, 2d or 3d).
    
    Workload is divided among ranks, each rank collects its share,
    and results are gathered to rank 0.

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
    
    comm : MPI communicator
        MPI communicator for parallel operations

    Returns
    -------
    A tuple with
       F : a ndarray containing the required field
       info : a FieldMetaInformation object
       (contains information about the grid; see the corresponding docstring)
    """
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    it = series.iterations[iteration]
    
    # Extract the dataset and corresponding group
    field = it.meshes[field_name]
    if field.scalar:
        component = next(field.items())[1]
    else:
        component = field[component_name]
    
    # Dimensions of the grid
    full_shape = component.shape
    
    # Handle dataOrder 'F' (Fortran order)
    if field.get_attribute('dataOrder') == 'F':
        grid_spacing = field.grid_spacing[::-1]
        global_offset = field.grid_global_offset[::-1]
        grid_position = component.position[::-1]
    else:
        grid_spacing = field.grid_spacing
        global_offset = field.grid_global_offset
        grid_position = component.position
    
    grid_unit_SI = field.grid_unit_SI
    time = (it.time + field.time_offset) * it.time_unit_SI
    
    field_attrs = {a: field.get_attribute(a) for a in field.attributes}
    component_attrs = {a: component.get_attribute(a) for a in component.attributes}
    
    # Determine final shape after slicing
    if slice_across is not None:
        # Get the integer that correspond to the slicing direction
        list_slicing_index = []
        list_i_cell = []
        for count, slice_across_item in enumerate(slice_across):
            slicing_index = axis_labels.index(slice_across_item)
            list_slicing_index.append(slicing_index)
            # Number of cells along the slicing direction
            n_cells = full_shape[slicing_index]
            # Index of the slice (prevent stepping out of the array)
            i_cell = int(0.5 * (slice_relative_position[count] + 1.) * n_cells)
            i_cell = max(i_cell, 0)
            i_cell = min(i_cell, n_cells - 1)
            list_i_cell.append(i_cell)
        
        # Calculate final shape after slicing
        final_shape = [x for index, x in enumerate(full_shape)
                      if index not in list_slicing_index]
        final_grid_spacing = [x for index, x in enumerate(grid_spacing)
                             if index not in list_slicing_index]
        final_global_offset = [x for index, x in enumerate(global_offset)
                              if index not in list_slicing_index]
        final_axis_labels = [x for index, x in enumerate(axis_labels)
                            if index not in list_slicing_index]
    else:
        final_shape = list(full_shape)
        final_grid_spacing = list(grid_spacing)
        final_global_offset = list(global_offset)
        final_axis_labels = list(axis_labels)
        list_slicing_index = None
        list_i_cell = None
    
    # Divide workload along the first dimension
    if len(final_shape) == 0:
        # Scalar data - no parallelization needed
        if rank == 0:
            if slice_across is not None:
                F = get_data(series, component, list_i_cell, list_slicing_index)
            else:
                F = get_data(series, component)
        else:
            F = None
        F = comm.bcast(F, root=0)
    else:
        # Divide work along the first dimension
        dim0_size = final_shape[0]
        chunk_size = dim0_size // size
        remainder = dim0_size % size
        
        # Calculate start and end indices for this rank
        start_idx = rank * chunk_size + min(rank, remainder)
        if rank < remainder:
            end_idx = start_idx + chunk_size + 1
        else:
            end_idx = start_idx + chunk_size
        
        # Each rank reads only its portion using a helper function
        F = _read_field_portion(series, component, slice_across, list_slicing_index,
                               list_i_cell, final_shape, start_idx, end_idx)
        
        # Gather all chunks to rank 0
        gathered_data = comm.gather(F, root=0)
        
        if rank == 0:
            # Concatenate all chunks along the first dimension
            F = np.concatenate(gathered_data, axis=0)
        else:
            F = None
            # Rank 0 will broadcast the final result
        F = comm.bcast(F, root=0)
    
    # Create FieldMetaInformation (same for all ranks)
    axes = {i: final_axis_labels[i] for i in range(len(final_axis_labels))}
    info = FieldMetaInformation(axes, F.shape, final_grid_spacing, final_global_offset,
            grid_unit_SI, grid_position,
            time, iteration, field_attrs=field_attrs,
            component_attrs=component_attrs)
    
    return F, info


def read_field_circ(series, iteration, field_name, component_name,
                   slice_relative_position, slice_across, comm,
                   m=0, theta=0., max_resolution_3d=None):
    """
    Extract a given field from a file in the openPMD format,
    when the geometry is thetaMode.
    
    Workload is divided among ranks, each rank collects its share,
    and results are gathered to rank 0.

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
        e.g. `[200, 100]`, indicating the maximum longitudinal and
        transverse resolution, respectively. This is useful for
        performance reasons, particularly for 3D visualization.
    
    comm : MPI communicator
        MPI communicator for parallel operations

    Returns
    -------
    A tuple with
       F : a 3darray or 2darray containing the required field,
           depending on whether `theta` is None or not
       info : a FieldMetaInformation object
       (contains information about the grid; see the corresponding docstring)
    """
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    it = series.iterations[iteration]
    
    # Extract the dataset and corresponding group
    field = it.meshes[field_name]
    if field.scalar:
        component = next(field.items())[1]
    else:
        component = field[component_name]
    
    field_attrs = {a: field.get_attribute(a) for a in field.attributes}
    component_attrs = {a: component.get_attribute(a) for a in component.attributes}
    
    # Extract the metainformation
    coord_labels = {ii: coord for (ii, coord) in enumerate(field.axis_labels)}
    
    if coord_labels[0] == 'r':
        coord_order = RZorder.mrz
        Nm, Nr, Nz = component.shape
        N_pair = (Nr, Nz)
    elif coord_labels[1] == 'r':
        Nm, Nz, Nr = component.shape
        N_pair = (Nz, Nr)
        coord_order = RZorder.mzr
    else:
        raise Exception(order_error_msg)
    
    time = (it.time + field.time_offset) * it.time_unit_SI
    
    info = FieldMetaInformation(coord_labels, N_pair,
        field.grid_spacing, field.grid_global_offset,
        field.grid_unit_SI, component.position, time, iteration,
        thetaMode=True, field_attrs=field_attrs,
        component_attrs=component_attrs)
    
    # Convert to a 3D Cartesian array if theta is None
    if theta is None:
        # Get cylindrical info (will be updated if max_resolution_3d is used)
        rmax = info.rmax
        inv_dr = 1./info.dr
        
        # Divide workload: each rank reads portion of modes
        # For simplicity, divide along the first spatial dimension (r or z)
        # Determine which dimension to divide along
        if coord_order is RZorder.mrz:
            divide_dim = 1  # Divide along r dimension
            dim_size = Nr
        else:  # RZorder.mzr
            divide_dim = 1  # Divide along z dimension  
            dim_size = Nz
        
        chunk_size = dim_size // size
        remainder = dim_size % size
        start_idx = rank * chunk_size + min(rank, remainder)
        if rank < remainder:
            end_idx = start_idx + chunk_size + 1
        else:
            end_idx = start_idx + chunk_size
        
        # Each rank reads its portion of the data using helper function
        if coord_order is RZorder.mrz:
            Fcirc = _read_field_circ_portion(series, component, start_idx, end_idx, 
                                             coord_order, dim=1)
        else:  # RZorder.mzr
            Fcirc = _read_field_circ_portion(series, component, start_idx, end_idx,
                                             coord_order, dim=2)
        
        # Handle max_resolution_3d if needed (on full data before division)
        if max_resolution_3d is not None:
            max_res_lon, max_res_transv = max_resolution_3d
            if Nz > max_res_lon:
                excess_z = int(np.round(Nz/max_res_lon))
                if coord_order is RZorder.mrz:
                    Fcirc = Fcirc[:, :, ::excess_z]
                else:  # RZorder.mzr
                    Fcirc = Fcirc[:, ::excess_z, :]
                # Update info (only rank 0 needs to do this, then broadcast)
                if rank == 0:
                    info.z = info.z[::excess_z]
                    info.dz = info.z[1] - info.z[0]
            if Nr > max_res_transv/2:
                excess_r = int(np.round(Nr/(max_res_transv/2)))
                if coord_order is RZorder.mrz:
                    Fcirc = Fcirc[:, ::excess_r, :]
                    Nr_local = Fcirc.shape[1]
                else:  # RZorder.mzr
                    Fcirc = Fcirc[:, :, ::excess_r]
                    Nr_local = Fcirc.shape[2]
                # Update info (only rank 0 needs to do this, then broadcast)
                if rank == 0:
                    info.r = info.r[::excess_r]
                    info.dr = info.r[1] - info.r[0]
                    inv_dr = 1./info.dr
                    rmax = info.rmax  # Update rmax after downsampling r
                    Nr = Fcirc.shape[1] if coord_order is RZorder.mrz else Fcirc.shape[2]
        
        # Broadcast updated info from rank 0
        if max_resolution_3d is not None:
            if rank == 0:
                info_to_bcast = {
                    'z': getattr(info, 'z', None),
                    'dz': getattr(info, 'dz', None),
                    'r': getattr(info, 'r', None),
                    'dr': getattr(info, 'dr', None),
                    'inv_dr': inv_dr,
                    'rmax': rmax,
                    'Nr': Nr
                }
            else:
                info_to_bcast = None
            info_to_bcast = comm.bcast(info_to_bcast, root=0)
            if rank != 0:
                if info_to_bcast['z'] is not None:
                    info.z = info_to_bcast['z']
                    info.dz = info_to_bcast['dz']
                if info_to_bcast['r'] is not None:
                    info.r = info_to_bcast['r']
                    info.dr = info_to_bcast['dr']
                    inv_dr = info_to_bcast['inv_dr']
                    rmax = info_to_bcast['rmax']
                    Nr = info_to_bcast['Nr']
        
        # Determine modes to extract
        if m == 'all':
            modes = [mode for mode in range(0, int(Nm / 2) + 1)]
        else:
            modes = [m]
        modes = np.array(modes, dtype='int')
        nmodes = len(modes)
        
        # Convert cylindrical data to Cartesian data
        # First, gather all Fcirc portions to rank 0 to reconstruct full Fcirc
        gathered_Fcirc = comm.gather(Fcirc, root=0)
        
        if rank == 0:
            # Concatenate along the divided dimension to reconstruct full Fcirc
            if coord_order is RZorder.mrz:
                Fcirc_full = np.concatenate(gathered_Fcirc, axis=1)
            else:  # RZorder.mzr
                Fcirc_full = np.concatenate(gathered_Fcirc, axis=2)
        else:
            Fcirc_full = None
        
        # Convert cylindrical data to Cartesian data
        info._convert_cylindrical_to_3Dcartesian()
        nx, ny, nz = len(info.x), len(info.y), len(info.z)
        
        # Rank 0 performs the full 3D conversion
        if rank == 0:
            F_total = np.zeros((nx, ny, nz), dtype=component.dtype)
            construct_3d_from_circ(F_total, Fcirc_full, info.x, info.y, modes,
                                   nx, ny, nz, Nr, nmodes, inv_dr, rmax, coord_order)
        else:
            F_total = None
        
        # Broadcast final result to all ranks
        F_total = comm.bcast(F_total, root=0)
        
    else:
        # theta is not None - 2D projection
        # Divide workload along the first spatial dimension (r or z)
        if coord_order is RZorder.mrz:
            divide_dim = 1  # Divide along r dimension
            dim_size = Nr
        else:  # RZorder.mzr
            divide_dim = 1  # Divide along z dimension
            dim_size = Nz
        
        chunk_size = dim_size // size
        remainder = dim_size % size
        start_idx = rank * chunk_size + min(rank, remainder)
        if rank < remainder:
            end_idx = start_idx + chunk_size + 1
        else:
            end_idx = start_idx + chunk_size
        
        # Extract the modes and recombine them properly
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
            # - Sum the modes (read full 2D data, then slice locally)
            F = get_data( series, component )  # (Extracts all modes)
            # Slice locally
            if coord_order is RZorder.mrz:
                F_local = F[:, start_idx:end_idx, :]
                F_total_local = np.zeros( (2 * Nr, Nz ), dtype=F.dtype )
                F_total_local[Nr:, :] = np.tensordot( mult_above_axis,
                                                    F_local, axes=(0, 0) )[:, :]
                F_total_local[:Nr, :] = np.tensordot( mult_below_axis,
                                                    F_local, axes=(0, 0) )[::-1, :]
                # Extract only our portion
                F_total_local = F_total_local[start_idx:end_idx, :]
            elif coord_order is RZorder.mzr:
                F_local = F[:, :, start_idx:end_idx]
                F_total_local = np.zeros( (Nz, 2 * Nr ), dtype=F.dtype )
                F_total_local[:, Nr:] = np.tensordot( mult_above_axis,
                                                    F_local, axes=(0, 0) )[:, :]
                F_total_local[:, :Nr] = np.tensordot( mult_below_axis,
                                                    F_local, axes=(0, 0) )[:, ::-1]
                # Extract only our portion
                F_total_local = F_total_local[:, start_idx:end_idx]
        elif m == 0:
            # Extract mode 0
            F = get_data( series, component, 0, 0 )
            if coord_order is RZorder.mrz:
                F_local = F[start_idx:end_idx, :]
                F_total_local = np.zeros( (2 * Nr, Nz ), dtype=F.dtype )
                F_total_local[Nr:, :] = F_local[:, :]
                F_total_local[:Nr, :] = F_local[::-1, :]
                F_total_local = F_total_local[start_idx:end_idx, :]
            elif coord_order is RZorder.mzr:
                F_local = F[:, start_idx:end_idx]
                F_total_local = np.zeros( (Nz, 2 * Nr ), dtype=F.dtype )
                F_total_local[:, Nr:] = F_local[:, :]
                F_total_local[:, :Nr] = F_local[:, ::-1]
                F_total_local = F_total_local[:, start_idx:end_idx]
        else:
            # Extract higher mode
            cos = np.cos( m * theta )
            sin = np.sin( m * theta )
            F_cos = get_data( series, component, 2 * m - 1, 0 )
            F_sin = get_data( series, component, 2 * m, 0 )
            if coord_order is RZorder.mrz:
                F_cos_local = F_cos[start_idx:end_idx, :]
                F_sin_local = F_sin[start_idx:end_idx, :]
                F_local = cos * F_cos_local + sin * F_sin_local
                F_total_local = np.zeros( (2 * Nr, Nz ), dtype=F_local.dtype )
                F_total_local[Nr:, :] = F_local[:, :]
                F_total_local[:Nr, :] = (-1) ** m * F_local[::-1, :]
                F_total_local = F_total_local[start_idx:end_idx, :]
            elif coord_order is RZorder.mzr:
                F_cos_local = F_cos[:, start_idx:end_idx]
                F_sin_local = F_sin[:, start_idx:end_idx]
                F_local = cos * F_cos_local + sin * F_sin_local
                F_total_local = np.zeros( (Nz, 2 * Nr ), dtype=F_local.dtype )
                F_total_local[:, Nr:] = F_local[:, :]
                F_total_local[:, :Nr] = (-1) ** m * F_local[:, ::-1]
                F_total_local = F_total_local[:, start_idx:end_idx]
        
        # Gather all portions to rank 0
        gathered_data = comm.gather(F_total_local, root=0)
        
        if rank == 0:
            # Concatenate along the divided dimension
            if coord_order is RZorder.mrz:
                F_total = np.concatenate(gathered_data, axis=0)
            else:  # RZorder.mzr
                F_total = np.concatenate(gathered_data, axis=1)
        else:
            F_total = None
        
        # Broadcast final result
        F_total = comm.bcast(F_total, root=0)
    
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


def get_grid_parameters(series, iteration, avail_fields, metadata, comm):
    """
    Return the parameters of the spatial grid (grid size and grid range)
    in two dictionaries.
    
    Only rank 0 performs the grid parameter extraction, and broadcasts
    the results to all other ranks.

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
    
    comm : MPI communicator
        MPI communicator for parallel operations

    Returns:
    --------
    A tuple with `grid_size_dict` and `grid_range_dict`
    Both objects are dictionaries, with their keys being the labels of the axis
    of the grid (e.g. 'x', 'y', 'z')
    The values of `grid_size_dict` are the number of gridpoints along each axis
    The values of `grid_range_dict` are lists of two floats, which correspond
    to the min and max of the grid, along each axis.
    """
    from ..io_reader.field_reader import get_grid_parameters as _get_grid_parameters
    
    rank = comm.Get_rank()
    
    if rank == 0:
        grid_size_dict, grid_range_dict = _get_grid_parameters(
            series, iteration, avail_fields, metadata)
    else:
        grid_size_dict = None
        grid_range_dict = None
    
    # Broadcast results from rank 0 to all ranks
    grid_size_dict = comm.bcast(grid_size_dict, root=0)
    grid_range_dict = comm.bcast(grid_range_dict, root=0)
    
    return grid_size_dict, grid_range_dict
