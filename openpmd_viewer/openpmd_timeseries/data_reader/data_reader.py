"""
This file is part of the openPMD-viewer

It routes the calls to the data reader to either
the h5py data reader, or to openpmd-api.

Copyright 2020, openPMD-viewer contributors
Authors: Remi Lehe
License: 3-Clause-BSD-LBNL
"""
import numpy as np
import os
import re

available_backends = []

try:
    import openpmd_api as io
    from . import io_reader
    available_backends.append('openpmd-api')
except ImportError:
    pass

try:
    from . import h5py_reader
    available_backends.append('h5py')
except ImportError:
    pass

if len(available_backends) == 0:
    raise ImportError('No openPMD backend found.\n'
        'Please install either `h5py` or `openpmd-api`:\n'
        'e.g. with `pip install h5py` or `pip install openpmd-api`')

class DataReader( object ):
    """
    Class that performs various type of access the openPMD file.

    The methods of this class are agnostic of the actual backend package
    used in order to access the openPMD file (e.g. h5py or openpmd-api).
    The backend that is used in practice depends on which package is
    available on the current environment.
    """

    def __init__(self, backend):
        """
        Initialize the DataReader class.
        """
        self.backend = backend

        # Point to the correct reader module
        if self.backend == 'h5py':
            self.iteration_to_file = {}
        elif self.backend == 'openpmd-api':
            pass
        else:
            raise RuntimeError('Unknown backend: %s' % self.backend)

    def list_iterations(self, path_to_dir):
        """
        Return a list of the iterations that correspond to the files
        in this directory. (The correspondance between iterations and
        files is stored internally.)

        Parameter
        ---------
        path_to_dir : string
            The path to the directory where the hdf5 files are.

        Returns
        -------
        an array of integers which correspond to the iteration of each file
        (in sorted order)
        """
        if self.backend == 'h5py':
            iterations, iteration_to_file = \
                h5py_reader.list_files( path_to_dir )
            # Store dictionary of correspondence between iteration and file
            self.iteration_to_file = iteration_to_file
            if len(iterations) == 0:
                raise RuntimeError(
                    "Found no valid files in directory {0}.\n"
                    "Please check that this is the path to the openPMD files."
                    "Valid files must have the extension '.h5' if you "
                    "use the `h5py` backend. For ADIOS '.bp' and other files, "
                    "please install the `openpmd-api` package."
                    .format(path_to_dir))
        elif self.backend == 'openpmd-api':
            # guess file ending from first file in directory
            first_file_name = None

            is_single_file = os.path.isfile(path_to_dir)
            if is_single_file:
                first_file_name = path_to_dir
            else:
                for file_name in os.listdir( path_to_dir ):
                    if file_name.split(os.extsep)[-1] in io.file_extensions:
                        first_file_name = file_name
            if first_file_name is None:
                raise RuntimeError(
                    "Found no valid files in directory {0}.\n"
                    "Please check that this is the path to the openPMD files."
                    "(valid files must have one of the following extensions: {1})"
                    .format(path_to_dir, io.file_extensions))

            if is_single_file:
                file_path = path_to_dir
                series_name = file_path
            else:
                # match last occurance of integers and replace with %T wildcards
                # examples: data00000100.h5 diag4_00000500.h5 io12.0.bp
                #           te42st.1234.yolo.json scan7_run14_data123.h5
                file_path = re.sub(r'(\d+)(\.(?!\d).+$)', r'%T\2', first_file_name)
                series_name = os.path.join( path_to_dir, file_path)

            self.series = io.Series(
                series_name,
                io.Access.read_only )
            iterations = np.array( self.series.iterations )

        return iterations

    def read_openPMD_params(self, iteration, extract_parameters=True):
        """
        Extract the time and some openPMD parameters from a file

        Parameter
        ---------
        iteration: int
            The iteration at which the parameters should be extracted

        extract_parameters: bool, optional
            Whether to extract all parameters or only the time
            (Function execution is faster when extract_parameters is False)

        Returns
        -------
        A tuple with:
        - A float corresponding to the time of this iteration in SI units
        - A dictionary containing several parameters, such as the geometry, etc
         When extract_parameters is False, the second argument returned is None
        """
        if self.backend == 'h5py':
            filename = self.iteration_to_file[iteration]
            return h5py_reader.read_openPMD_params(
                    filename, iteration, extract_parameters)

        elif self.backend == 'openpmd-api':
            return io_reader.read_openPMD_params(
                    self.series, iteration, extract_parameters)

    def read_field_cartesian( self, iteration, field, coord, axis_labels,
                          slice_relative_position, slice_across ):
        """
        Extract a given field from an openPMD file in the openPMD format,
        when the geometry is cartesian (1d, 2d or 3d).

        Parameters
        ----------
        iteration : int
           The iteration at which to extract the fields

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
        if self.backend == 'h5py':
            filename = self.iteration_to_file[iteration]
            return h5py_reader.read_field_cartesian(
                filename, iteration, field, coord, axis_labels,
                slice_relative_position, slice_across )
        elif self.backend == 'openpmd-api':
            return io_reader.read_field_cartesian(
                self.series, iteration, field, coord, axis_labels,
                slice_relative_position, slice_across )

    def read_field_circ( self, iteration, field, coord, slice_relative_position,
                        slice_across, m=0, theta=0., max_resolution_3d=None ):
        """
        Extract a given field from an openPMD file in the openPMD format,
        when the geometry is thetaMode

        Parameters
        ----------
        iteration : int
           The iteration at which to extract the fields

        field : string, optional
           Which field to extract
           Either 'rho', 'E', 'B' or 'J'

        coord : string, optional
           Which component of the field to extract
           Either 'r', 't' or 'z'

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

        Returns
        -------
        A tuple with
           F : a 3darray or 2darray containing the required field,
               depending on whether `theta` is None or not
           info : a FieldMetaInformation object
           (contains information about the grid; see the corresponding docstring)
        """
        if self.backend == 'h5py':
            filename = self.iteration_to_file[iteration]
            return h5py_reader.read_field_circ(
                filename, iteration, field, coord, slice_relative_position,
                slice_across, m, theta, max_resolution_3d )
        elif self.backend == 'openpmd-api':
            return io_reader.read_field_circ(
                self.series, iteration, field, coord, slice_relative_position,
                slice_across, m, theta, max_resolution_3d )

    def read_species_data( self, iteration, species, record_comp, extensions):
        """
        Extract a given species' record_comp

        Parameters
        ----------
        iteration: int
            The iteration at which to extract the species data

        species: string
            The name of the species to extract (in the openPMD file)

        record_comp: string
            The record component to extract
            Either 'x', 'y', 'z', 'ux', 'uy', 'uz', or 'w'

        extensions: list of strings
            The extensions that the current OpenPMDTimeSeries complies with
        """
        if self.backend == 'h5py':
            filename = self.iteration_to_file[iteration]
            return h5py_reader.read_species_data(
                    filename, iteration, species, record_comp, extensions )
        elif self.backend == 'openpmd-api':
            return io_reader.read_species_data(
                    self.series, iteration, species, record_comp, extensions )

    def get_grid_parameters(self, iteration, avail_fields, metadata ):
        """
        Return the parameters of the spatial grid (grid size and grid range)
        in two dictionaries

        Parameters:
        -----------
        iteration: int
            The iteration at which to extract the parameters

        avail_fields: list
           A list of the available fields
           e.g. ['B', 'E', 'rho']

        metadata: dictionary
          A dictionary whose keys are the fields of `avail_fields` and
          whose values are dictionaries that contain metadata (e.g. geometry)

        Returns:
        --------
        A tuple with `grid_size_dict` and `grid_range_dict`
        Both objects are dictionaries, with their keys being the labels of
        the axis of the grid (e.g. 'x', 'y', 'z')
        The values of `grid_size_dict` are the number of gridpoints along
        each axis.
        The values of `grid_range_dict` are lists of two floats, which
        correspond to the min and max of the grid, along each axis.
        """
        if self.backend == 'h5py':
            filename = self.iteration_to_file[iteration]
            return h5py_reader.get_grid_parameters(
                filename, iteration, avail_fields, metadata )
        elif self.backend == 'openpmd-api':
            return io_reader.get_grid_parameters(
                self.series, iteration, avail_fields, metadata )
