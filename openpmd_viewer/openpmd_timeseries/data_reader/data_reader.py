"""
This file is part of the openPMD-viewer

It routes the calls to the data reader to either
the h5py data reader, or to openpmd-api.

Copyright 2020, openPMD-viewer contributors
Authors: Remi Lehe
License: 3-Clause-BSD-LBNL
"""
from . import h5py_reader


class DataReader( object ):
    """
    TODO
    """

    def __init__(self, backend='h5py'):
        """
        TODO
        """
        self.backend = backend

        # Point to the correct reader module
        if self.backend == 'h5py':
            self.iteration_to_file = {}

        else:
            raise RuntimeError('Unknown backend: %s' % backend)

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
                    filename, extract_parameters)

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
                filename, field, coord, axis_labels,
                slice_relative_position, slice_across )

    def read_field_circ( self, iteration, field, coord, slice_relative_position,
                        slice_across, m=0, theta=0. ):
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
                filename, field, coord, slice_relative_position,
                slice_across, m, theta )

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
                    filename, species, record_comp, extensions )

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
                filename, avail_fields, metadata )
