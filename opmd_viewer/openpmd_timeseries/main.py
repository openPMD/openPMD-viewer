"""
This file is part of the openPMD-viewer.

It defines the main OpenPMDTimeSeries class.

Copyright 2015-2016, openPMD-viewer contributors
Authors: Remi Lehe, Axel Huebl
License: 3-Clause-BSD-LBNL
"""

import os
import re
import numpy as np
import h5py as h5
from .plotter import Plotter
from .data_reader.params_reader import read_openPMD_params
from .data_reader.particle_reader import read_species_data
from .data_reader.field_reader import read_field_2d, \
    read_field_circ, read_field_3d


# Check wether the interactive interface can be loaded
try:
    # If yes, use the InteractiveViewer as a parent class
    from .interactive import InteractiveViewer
    parent_class = InteractiveViewer
except ImportError:
    # Otherwise, use the default parent class
    print('[opmd_viewer] Failed to import the interactive interface.\n'
          '(Make sure that ipywidgets is installed.)\n'
          'The opmd_viewer API is nonetheless working.')
    parent_class = object


# Define a custom Exception
class OpenPMDException(Exception):
    "Exception raised for invalid use of the openPMD-viewer API"
    pass


# Define the OpenPMDTimeSeries class and have it inherit
# from the parent class defined above
class OpenPMDTimeSeries(parent_class):
    """
    Main class for the exploration of an openPMD timeseries

    For more details, see the docstring of the following methods:
    - get_field
    - get_particle
    - slider
    """

    def __init__(self, path_to_dir, check_all_files=True):
        """
        Initialize an openPMD time series

        More precisely, scan the directory and extract the openPMD files,
        as well as some useful openPMD parameters

        Parameters
        ----------
        path_to_dir: string
            The path to the directory where the openPMD files are.
            For the moment, only HDF5 files are supported. There should be
            one file per iteration, and the name of the files should end
            with the iteration number, followed by '.h5' (e.g. data0005000.h5)

        check_all_files: bool, optional
            Check that all the files in the timeseries are consistent
            (i.e. that they contain the same fields and particles,
            with the same metadata)
            For fast access to the files, this can be changed to False.
        """
        # Extract the files and the iterations
        self.h5_files, self.iterations = list_h5_files(path_to_dir)

        # Check that there are HDF5 files in this directory
        if len(self.h5_files) == 0:
            print("Error: Found no HDF5 files in the specified directory.\n"
                  "Please check that this is the path to the HDF5 files.")
            return(None)

        # Go through the files of the series, extract the time
        # and a few parameters.
        N_files = len(self.h5_files)
        self.t = np.zeros(N_files)

        # - Extract parameters from the first file
        t, params0 = read_openPMD_params(self.h5_files[0])
        self.t[0] = t
        self.avail_fields = params0['avail_fields']
        self.extensions = params0['extensions']
        if self.avail_fields is not None:
            self.geometry = params0['geometry']
            self.axis_labels = params0['axis_labels']
            self.avail_circ_modes = params0['avail_circ_modes']
        self.avail_species = params0['avail_species']
        self.avail_record_components = \
            params0['avail_record_components']

        # - Extract the time for each file and, if requested, check
        #   that the other files have the same parameters
        for k in range(1, N_files):
            t, params = read_openPMD_params(self.h5_files[k], check_all_files)
            self.t[k] = t
            if check_all_files:
                for key in params0.keys():
                    if params != params0:
                        print("Warning: File %s has different openPMD "
                              "parameters than the rest of the time series."
                              % self.h5_files[k])
                        break

        # - Set the current iteration and time
        self.current_i = 0
        self.current_t = self.t[0]
        # - Find the min and the max of the time
        self.tmin = self.t.min()
        self.tmax = self.t.max()

        # - Initialize a plotter object, which holds information about the time
        self.plotter = Plotter(self.t, self.iterations)

    def get_particle(self, var_list=None, species=None, t=None,
                     iteration=None, select=None, output=True,
                     plot=False, nbins=150, **kw):
        """
        Extract a list of particle variables
        from an HDF5 file in the openPMD format.

        In the case of positions, the result is returned in microns

        Plot the histogram of the returned quantity.
        If two quantities are requested by the user, this plots
        a 2d histogram of these quantities.

        Parameters
        ----------
        var_list : list of string, optional
            A list of the particle variables to extract. If var_list is not
            provided, the available particle quantities are printed

        species: string
            A string indicating the name of the species

        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        output : bool, optional
           Whether to return the requested quantity

        select: dict, optional
            Either None or a dictionary of rules
            to select the particles, of the form
            'x' : [-4., 10.]   (Particles having x between -4 and 10 microns)
            'ux' : [-0.1, 0.1] (Particles having ux between -0.1 and 0.1 mc)
            'uz' : [5., None]  (Particles with uz above 5 mc)

        plot : bool, optional
           Whether to plot the requested quantity
           Plotting support is only available when requesting one or two
           quantities (i.e. when var_list is of length 1 or 2)

        nbins : int, optional
           Number of bins for the histograms

        **kw : dict, otional
           Additional options to be passed to matplotlib's
           hist or hist2d.

        Returns
        -------
        A list of 1darray corresponding to the data requested in `var_list`
        (one 1darray per element of 'var_list', returned in the same order)
        """
        # Check that the species and quantity required are present
        if self.avail_species is None:
            raise OpenPMDException('No particle data in this time series')
        if species not in self.avail_species:
            species_list = '\n - '.join(self.avail_species)
            raise OpenPMDException(
                "The argument `species` is missing or erroneous.\n"
                "The available species are: \n - %s\nPlease set the "
                "argument `species` accordingly." % species_list)

        # Check the list of variables
        valid_var_list = True
        if not isinstance(var_list, list):
            valid_var_list = False
        else:
            for quantity in var_list:
                if quantity not in self.avail_record_components[species]:
                    valid_var_list = False
        if not valid_var_list:
            quantity_list = '\n - '.join(
                self.avail_record_components[species])
            raise OpenPMDException(
                "The argument `var_list` is missing or erroneous.\n"
                "It should be a list of strings representing species record "
                "components.\n The available quantities for species '%s' are:"
                "\n - %s\nPlease set the argument `var_list` "
                "accordingly." % (species, quantity_list) )

        # Check the selection quantities
        if select is not None:
            valid_select_list = True
            if not isinstance(select, dict):
                valid_select_list = False
            else:
                for quantity in select.keys():
                    if (quantity in self.avail_record_components[species]) \
                            is False:
                        valid_select_list = False
            if not valid_select_list:
                quantity_list = '\n - '.join(
                    self.avail_record_components[species])
                raise OpenPMDException(
                    "The argument `select` is erroneous.\n"
                    "It should be a dictionary whose keys represent particle "
                    "quantities.\n The available quantities are: "
                    "\n - %s\nPlease set the argument `select` "
                    "accordingly." % quantity_list)

        # Find the output that corresponds to the requested time/iteration
        # (Modifies self.current_i and self.current_t)
        self._find_output(t, iteration)
        # Get the corresponding file name & open the file
        file_name = self.h5_files[self.current_i]
        file_handle = h5.File(file_name, 'r')

        # Extract the list of particle quantities
        data_list = []
        for quantity in var_list:
            data_list.append(read_species_data(
                file_handle, species, quantity, self.extensions))
        # Apply selection if needed
        if select is not None:
            data_list = apply_selection(
                file_handle, data_list, select, species, self.extensions)

        # Plotting
        if plot:

            # Extract the weights, if they are available
            if 'w' in self.avail_record_components[species]:
                w = read_species_data(
                    file_handle, species, 'w', self.extensions)
                if select is not None:
                    w, = apply_selection(
                        file_handle, [w], select, species, self.extensions)
            # Otherwise consider that all particles have a weight of 1
            else:
                w = np.ones_like(data_list[0])

            # - In the case of only one quantity
            if len(data_list) == 1:
                # Do the plotting
                self.plotter.hist1d(data_list[0], w, var_list[0], species,
                                    self.current_i, nbins, **kw)
            # - In the case of two quantities
            elif len(data_list) == 2:
                # Do the plotting
                self.plotter.hist2d(data_list[0], data_list[1], w,
                                    var_list[0], var_list[1], species,
                                    self.current_i, nbins, **kw)
        # Close the file
        file_handle.close()

        # Output
        if output:
            return(data_list)

    def get_field(self, field=None, coord=None, t=None, iteration=None,
                  m='all', theta=0., slicing=0., slicing_dir='y',
                  output=True, plot=False, **kw):
        """
        Extract a given field from an HDF5 file in the openPMD format.

        Parameters
        ----------

        field : string, optional
           Which field to extract
           Either 'rho', 'E', 'B' or 'J'

        coord : string, optional
           Which component of the field to extract
           Either 'r', 't' or 'z'

        m : int or str, optional
           Only used for thetaMode geometry
           Either 'all' (for the sum of all the modes)
           or an integer (for the selection of a particular mode)

        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        theta : float, optional
           Only used for thetaMode geometry
           The angle of the plane of observation, with respect to the x axis

        slicing : float, optional
           Only used for 3dcartesian geometry
           A number between -1 and 1 that indicates where to slice the data,
           along the direction `slicing_dir`
           -1 : lower edge of the simulation box
           0 : middle of the simulation box
           1 : upper edge of the simulation box
           If slicing is None, the full 3D grid is returned.

        slicing_dir : str, optional
           Only used for 3dcartesian geometry
           The direction along which to slice the data
           Either 'x', 'y' or 'z'

        output : bool, optional
           Whether to return the requested quantity

        plot : bool, optional
           Whether to plot the requested quantity

        **kw : dict, otional
           Additional options to be passed to matplotlib's imshow.

        Returns
        -------
        A tuple with
           F : a 2darray containing the required field
           info : a FieldMetaInformation object
           (see the corresponding docstring)
        """
        # Check that the field required is present
        if self.avail_fields is None:
            raise OpenPMDException('No field data in this time series')
        # Check the field type
        if field not in self.avail_fields:
            field_list = '\n - '.join(self.avail_fields)
            raise OpenPMDException(
                "The `field` argument is missing or erroneous.\n"
                "The available fields are: \n - %s\nPlease set the `field` "
                "argument accordingly." % field_list)
        # Check the coordinate (for vector fields)
        if self.avail_fields[field] == 'vector':
            available_coord = ['x', 'y', 'z']
            if self.geometry == 'thetaMode':
                available_coord += ['r', 't']
            if coord not in available_coord:
                coord_list = '\n - '.join(available_coord)
                raise OpenPMDException(
                    "The field %s is a vector field, \nbut the `coord` "
                    "argument is missing or erroneous.\nThe available "
                    "coordinates are: \n - %s\nPlease set the `coord` "
                    "argument accordingly." % (field, coord_list))
        # Check the mode (for thetaMode)
        if self.geometry == "thetaMode":
            if str(m) not in self.avail_circ_modes:
                mode_list = '\n - '.join(self.avail_circ_modes)
                raise OpenPMDException(
                    "The requested mode '%s' is not available.\n"
                    "The available modes are: \n - %s" % (m, mode_list))

        # Find the output that corresponds to the requested time/iteration
        # (Modifies self.current_i and self.current_t)
        self._find_output(t, iteration)
        # Get the corresponding filename
        filename = self.h5_files[self.current_i]

        # Find the proper path for vector or scalar fields
        if self.avail_fields[field] == 'scalar':
            field_path = field
            field_label = field
        elif self.avail_fields[field] == 'vector':
            field_path = os.path.join(field, coord)
            field_label = field + coord

        # Get the field data
        # - For 2D
        if self.geometry == "2dcartesian":
            F, info = read_field_2d(filename, field_path, self.axis_labels)
        # - For 3D
        elif self.geometry == "3dcartesian":
            F, info = read_field_3d(
                filename, field_path, self.axis_labels, slicing, slicing_dir)
        # - For thetaMode
        elif self.geometry == "thetaMode":
            if (coord in ['x', 'y']) and \
                    (self.avail_fields[field] == 'vector'):
                # For Cartesian components, combine r and t components
                Fr, info = read_field_circ(filename, field + '/r', m, theta)
                Ft, info = read_field_circ(filename, field + '/t', m, theta)
                if coord == 'x':
                    F = np.cos(theta) * Fr - np.sin(theta) * Ft
                elif coord == 'y':
                    F = np.sin(theta) * Fr + np.cos(theta) * Ft
                # Revert the sign below the axis
                F[: int(F.shape[0] / 2)] *= -1
            else:
                # For cylindrical or scalar components, no special treatment
                F, info = read_field_circ(filename, field_path, m, theta)

        # Plot the resulting field
        # Deactivate plotting when there is no slice selection
        if (self.geometry == "3dcartesian") and (slicing is None):
            plot = False
        if plot:
            self.plotter.show_field(F, info, slicing_dir, m,
                                    field_label, self.geometry,
                                    self.current_i, **kw)

        # Return the result
        return(F, info)

    def _find_output(self, t, iteration):
        """
        Find the output that correspond to the requested `t` or `iteration`
        Modify self.current_i accordingly.

        Parameter
        ---------
        t : float (in seconds)
            Time requested

        iteration : int
            Iteration requested
        """
        # Check the arguments
        if (t is not None) and (iteration is not None):
            raise OpenPMDException(
                "Please pass either a time (`t`) \nor an "
                "iteration (`iteration`), but not both.")
        # If a time is requested
        elif (t is not None):
            # Make sur the time requested does not exceed the allowed bounds
            if t < self.tmin:
                self.current_i = 0
            elif t > self.tmax:
                self.current_i = len(self.t) - 1
            # Find the last existing output
            else:
                self.current_i = self.t[self.t <= t].argmax()
        # If an iteration is requested
        elif (iteration is not None):
            if (iteration in self.iterations):
                self.current_i = self.iterations.index(iteration)
            else:
                iter_list = '\n - '.join([str(it) for it in self.iterations])
                print("The requested iteration '%s' is not available.\nThe "
                      "available iterations are: \n - %s\nThe first iteration "
                      "is used instead." % (iteration, iter_list))
                self.current_i = 0
        else:
            pass  # self.current_i retains its previous value

        # Register the value in the object
        self.current_t = self.t[self.current_i]


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
