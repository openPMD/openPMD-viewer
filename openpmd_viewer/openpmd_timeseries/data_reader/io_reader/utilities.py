"""
This file is part of the openPMD-viewer.

It defines a set of helper data and functions which
are used by the other files.

Copyright 2015-2016, openPMD-viewer contributors
Authors: Remi Lehe, Axel Huebl
License: 3-Clause-BSD-LBNL
"""
import numpy as np


def chunk_to_slice(chunk):
    """
    Convert an openPMD_api.ChunkInfo to slice
    """
    stops = [a + b for a, b in zip(chunk.offset, chunk.extent)]
    indices_per_dim = zip(chunk.offset, stops)
    index_tuple = map(lambda s: slice(s[0], s[1], None), indices_per_dim)
    return tuple(index_tuple)


def get_data(series, record_component, i_slice=None, pos_slice=None,
             output_type=None):
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

    # ADIOS2: Actual chunks, all other: one chunk
    chunks = record_component.available_chunks()

    # read whole data set
    if pos_slice is None:
        # mask invalid regions with NaN
        #   note: full_like triggers a full read, thus we avoid it #340
        data = np.full(record_component.shape, np.nan, record_component.dtype)
        for chunk in chunks:
            chunk_slice = chunk_to_slice(chunk)

            # skip empty slices
            # https://github.com/ornladios/ADIOS2
            volume = 1
            for csl in chunk_slice:
                volume *= csl.stop - csl.start
            if volume == 0:
                continue

            # read only valid region
            x = record_component[chunk_slice]
            series.flush()
            data[chunk_slice] = x
    # slice: read only part of the data set
    else:
        full_shape = record_component.shape

        slice_shape = list(full_shape)      # copy
        pos_slice_sorted = pos_slice.copy() # copy for in-place sort
        pos_slice_sorted.sort(reverse=True)
        for dir_index in pos_slice_sorted:  # remove indices in list
            del slice_shape[dir_index]

        # mask invalid regions with NaN
        data = np.full(slice_shape, np.nan, dtype=record_component.dtype)

        # build requested ND slice with respect to full data
        s = []
        for d in range(len(full_shape)):
            if d in pos_slice:
                s.append(i_slice[pos_slice.index(d)]) # one index in such directions
            else: # all indices in other direction
                s.append(slice(None, None, None))
        s = tuple(s)

        # now we check which chunks contribute to the slice
        for chunk in chunks:
            skip_this_chunk = False
            s_valid = list(s)  # same as s but reduced to valid regions in chunk
            s_target = []  # starts and stops in sliced array
            chunk_slice = chunk_to_slice(chunk)

            # skip empty slices
            # https://github.com/ornladios/ADIOS2
            volume = 1
            for csl in chunk_slice:
                volume *= csl.stop - csl.start
            if volume == 0:
                continue

            # read only valid region
            for d, slice_d in enumerate(s):
                start = chunk_slice[d].start
                stop = chunk_slice[d].stop
                if isinstance(slice_d, int):
                    # Nothing to do for s_target (dimension sliced out)
                    # Nothing to do for s_valid (dimension index is set)
                    if slice_d < start or slice_d >= stop:
                        # chunk not in slice line/plane
                        skip_this_chunk = True
                else:
                    if slice_d.start is None or slice_d.start < start:
                        s_valid[d] = slice(start, s_valid[d].stop)
                    if slice_d.stop is None or slice_d.stop > stop:
                        s_valid[d] = slice(s_valid[d].start, stop)
                    s_target.append(slice(start, stop))

            s_valid = tuple(s_valid)
            s_target = tuple(s_target)

            # read
            if not skip_this_chunk:
                x = record_component[s_valid]
                series.flush()
                data[s_target] = x

    # Convert to the right type
    if (output_type is not None) and (data.dtype != output_type):
        data = data.astype( output_type )
    # Scale by the conversion factor
    if np.issubdtype(data.dtype, np.floating) or \
        np.issubdtype(data.dtype, np.complexfloating):
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
