"""
MPI-enabled parallel versions of io_reader utilities.

Copyright 2020, openPMD-viewer contributors
License: 3-Clause-BSD-LBNL
"""
import numpy as np
from ..io_reader.utilities import get_data, chunk_to_slice
from ...data_order import RZorder


def _read_field_portion(series, component, slice_across, list_slicing_index,
                       list_i_cell, final_shape, start_idx, end_idx):
    """
    Helper function to read only a portion of the field data.
    Reads only chunks that overlap with the specified range along the first dimension.
    
    Parameters
    ----------
    series : openpmd_api.Series
        The openPMD series object
    component : openPMD Record_Component
        The component to read from
    slice_across : list or None
        Slicing directions
    list_slicing_index : list or None
        Indices of slicing directions
    list_i_cell : list or None
        Cell indices for slicing
    final_shape : list
        Final shape after slicing
    start_idx : int
        Start index along first dimension for this rank
    end_idx : int
        End index along first dimension for this rank
    
    Returns
    -------
    numpy array with the portion of data for this rank
    """
    # Determine the shape we need to allocate
    portion_shape = [end_idx - start_idx] + final_shape[1:]
    
    # Get NaN value for masking
    NaN_value = np.nan if (np.issubdtype(component.dtype, np.floating) or 
                           np.issubdtype(component.dtype, np.complexfloating)) else 0
    
    # Allocate array for this rank's portion
    data = np.full(portion_shape, NaN_value, dtype=component.dtype)
    
    # Get available chunks
    chunks = component.available_chunks()
    
    # Build the slice for reading from the full dataset
    if slice_across is not None:
        # For sliced data, first read the full slice, then extract portion
        F_full = get_data(series, component, list_i_cell, list_slicing_index)
        return F_full[start_idx:end_idx]
    else:
        # For non-sliced data, read only chunks that overlap with our portion
        for chunk in chunks:
            chunk_slice = chunk_to_slice(chunk)
            
            # Check if chunk overlaps with our portion along first dimension
            chunk_start = chunk_slice[0].start
            chunk_end = chunk_slice[0].stop
            
            # Calculate overlap
            overlap_start = max(chunk_start, start_idx)
            overlap_end = min(chunk_end, end_idx)
            
            if overlap_start < overlap_end:
                # This chunk overlaps with our portion
                # Read the overlapping region
                read_slice_list = []
                for i, cs in enumerate(chunk_slice):
                    if i == 0:
                        # First dimension: only the overlap region
                        read_slice_list.append(slice(overlap_start, overlap_end))
                    else:
                        # Other dimensions: full chunk extent
                        read_slice_list.append(cs)
                read_slice = tuple(read_slice_list)
                
                # Read data
                x = component[read_slice]
                series.flush()
                
                # Place in our portion array
                # Calculate target slice: offset by start_idx for first dimension
                target_slice_list = [slice(overlap_start - start_idx, overlap_end - start_idx)]
                # For remaining dimensions, use absolute coordinates
                # (our portion array has full extent in those dimensions)
                for i in range(1, len(chunk_slice)):
                    chunk_dim_start = chunk_slice[i].start
                    chunk_dim_end = chunk_slice[i].stop
                    target_slice_list.append(slice(chunk_dim_start, chunk_dim_end))
                target_slice = tuple(target_slice_list)
                
                data[target_slice] = x
    
    # Scale by unit_SI if needed
    if component.unit_SI != 1.0:
        if np.issubdtype(data.dtype, np.floating) or \
           np.issubdtype(data.dtype, np.complexfloating):
            data *= component.unit_SI
        else:
            data = data * component.unit_SI
    
    return data


def _read_field_circ_portion(series, component, start_idx, end_idx, coord_order, dim):
    """
    Helper function specifically for 3D cylindrical field data (theta is None case).
    
    Parameters
    ----------
    series : openpmd_api.Series
        The openPMD series object
    component : openPMD Record_Component
        The component to read from
    start_idx : int
        Start index along the divided dimension for this rank
    end_idx : int
        End index along the divided dimension for this rank
    coord_order : RZorder
        Coordinate order (mrz or mzr)
    dim : int
        Dimension to divide along (1 for r in mrz, 2 for z in mzr)
    
    Returns
    -------
    numpy array with the portion of cylindrical data for this rank
    """
    # Get full shape
    full_shape = component.shape
    Nm = full_shape[0]
    
    # Determine local shape
    if coord_order is RZorder.mrz:
        # Shape: (Nm, local_size, Nz)
        local_size = end_idx - start_idx
        local_shape = (Nm, local_size, full_shape[2])
    else:  # RZorder.mzr
        # Shape: (Nm, Nz, local_size)
        local_size = end_idx - start_idx
        local_shape = (Nm, full_shape[1], local_size)
    
    # Get NaN value for masking
    NaN_value = np.nan if (np.issubdtype(component.dtype, np.floating) or 
                           np.issubdtype(component.dtype, np.complexfloating)) else 0
    
    # Allocate array for this rank's portion
    data = np.full(local_shape, NaN_value, dtype=component.dtype)
    
    # Get available chunks
    chunks = component.available_chunks()
    
    # Read only chunks that overlap with our portion along the divided dimension
    for chunk in chunks:
        chunk_slice = chunk_to_slice(chunk)
        
        # Check overlap along the divided dimension
        chunk_start = chunk_slice[dim].start
        chunk_end = chunk_slice[dim].stop
        
        overlap_start = max(chunk_start, start_idx)
        overlap_end = min(chunk_end, end_idx)
        
        if overlap_start < overlap_end:
            # This chunk overlaps with our portion
            # Build read slice
            read_slice_list = []
            for i, cs in enumerate(chunk_slice):
                if i == dim:
                    # Divided dimension: only the overlap region
                    read_slice_list.append(slice(overlap_start, overlap_end))
                else:
                    # Other dimensions: full chunk extent
                    read_slice_list.append(cs)
            read_slice = tuple(read_slice_list)
            
            # Read data
            x = component[read_slice]
            series.flush()
            
            # Place in our portion array
            target_slice_list = []
            for i in range(len(chunk_slice)):
                if i == dim:
                    # Offset by start_idx for divided dimension
                    target_slice_list.append(slice(overlap_start - start_idx, 
                                                  overlap_end - start_idx))
                else:
                    # Use absolute coordinates for other dimensions
                    chunk_dim_start = chunk_slice[i].start
                    chunk_dim_end = chunk_slice[i].stop
                    target_slice_list.append(slice(chunk_dim_start, chunk_dim_end))
            target_slice = tuple(target_slice_list)
            
            data[target_slice] = x
    
    # Scale by unit_SI if needed
    if component.unit_SI != 1.0:
        if np.issubdtype(data.dtype, np.floating) or \
           np.issubdtype(data.dtype, np.complexfloating):
            data *= component.unit_SI
        else:
            data = data * component.unit_SI
    
    return data


def _read_species_portion(series, component, start_idx, end_idx, output_type=np.float64):
    """
    Helper function to read a portion of particle/species data.
    
    Parameters
    ----------
    series : openpmd_api.Series
        The openPMD series object
    component : openPMD Record_Component
        The component to read from
    start_idx : int
        Start particle index for this rank
    end_idx : int
        End particle index for this rank
    output_type : numpy dtype, optional
        Output data type (default: np.float64, use np.uint64 for 'id' field)
    
    Returns
    -------
    numpy array with the portion of particle data for this rank
    """
    # Get NaN value for masking
    NaN_value = np.nan if (np.issubdtype(output_type, np.floating) or 
                           np.issubdtype(output_type, np.complexfloating)) else 0
    
    # Allocate array for this rank's portion
    portion_size = end_idx - start_idx
    data = np.full(portion_size, NaN_value, dtype=output_type)
    
    # Get available chunks
    chunks = component.available_chunks()
    
    # Read only chunks that overlap with our portion
    for chunk in chunks:
        chunk_slice = chunk_to_slice(chunk)
        
        # For 1D particle data, check overlap along first dimension
        chunk_start = chunk_slice[0].start
        chunk_end = chunk_slice[0].stop
        
        overlap_start = max(chunk_start, start_idx)
        overlap_end = min(chunk_end, end_idx)
        
        if overlap_start < overlap_end:
            # This chunk overlaps with our portion
            read_slice = (slice(overlap_start, overlap_end),)
            
            # Read data
            x = component[read_slice]
            series.flush()
            
            # Place in our portion array
            target_start = overlap_start - start_idx
            target_end = overlap_end - start_idx
            data[target_start:target_end] = x
    
    # Scale by unit_SI if needed
    if component.unit_SI != 1.0:
        if np.issubdtype(data.dtype, np.floating) or \
           np.issubdtype(data.dtype, np.complexfloating):
            data *= component.unit_SI
        else:
            data = data * component.unit_SI
    
    return data
