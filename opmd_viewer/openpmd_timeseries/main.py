"""
This file is part of the openPMD-viewer.

It defines the main OpenPMDTimeSeries class.

Copyright 2015-2016, openPMD-viewer contributors
Authors: Remi Lehe, Axel Huebl
License: 3-Clause-BSD-LBNL
"""

import os
import numpy as np
import h5py as h5
from .utilities import list_h5_files, apply_selection, fit_bins_to_grid
from .plotter import Plotter
from .particle_tracker import ParticleTracker
from .data_reader.params_reader import read_openPMD_params
from .data_reader.particle_reader import read_species_data
from .data_reader.field_reader import read_field_1d, read_field_2d, \
    read_field_circ, read_field_3d, get_grid_parameters
from .interactive import InteractiveViewer


# Define a custom Exception
class OpenPMDException(Exception):
    "Exception raised for invalid use of the openPMD-viewer API"
    pass


class OpenPMDTimeSeries(InteractiveViewer):
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
        self.extensions = params0['extensions']
        self.avail_fields = params0['avail_fields']
        if self.avail_fields is not None:
            self.fields_metadata = params0['fields_metadata']
            self.avail_geom = set( self.fields_metadata[field]['geometry']
                                for field in self.avail_fields )
        # Extract information of the particles
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
        self._current_i = 0
        self.current_iteration = self.iterations[0]
        self.current_t = self.t[0]
        # - Find the min and the max of the time
        self.tmin = self.t.min()
        self.tmax = self.t.max()

        # - Initialize a plotter object, which holds information about the time
        self.plotter = Plotter(self.t, self.iterations)

    def get_particle(self, var_list=None, species=None, t=None, iteration=None,
            select=None, output=True, plot=False, nbins=150,
            plot_range=[[None, None], [None, None]],
            use_field_mesh=True, histogram_deposition='cic', **kw):
        """
        Extract a list of particle variables
        from an HDF5 file in the openPMD format.

        Plot the histogram of the returned quantity.
        If two quantities are requested by the user, this plots
        a 2d histogram of these quantities.

        In the case of positions, the result is returned in microns
        In the case of momenta, the result is returned as:
        - unitless momentum (i.e. gamma*beta) for particles with non-zero mass
        - in kg.m.s^-1 for particles with zero mass

        Parameters
        ----------
        var_list : list of string, optional
            A list of the particle variables to extract. If var_list is not
            provided, the available particle quantities are printed

        species: string
            A string indicating the name of the species
            This is optional if there is only one species

        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        output : bool, optional
           Whether to return the requested quantity

        select: dict or ParticleTracker object, optional
            - If `select` is a dictionary:
            then it lists a set of rules to select the particles, of the form
            'x' : [-4., 10.]   (Particles having x between -4 and 10 microns)
            'ux' : [-0.1, 0.1] (Particles having ux between -0.1 and 0.1 mc)
            'uz' : [5., None]  (Particles with uz above 5 mc)
            - If `select` is a ParticleTracker object:
            then it returns particles that have been selected at another
            iteration ; see the docstring of `ParticleTracker` for more info.

        plot : bool, optional
           Whether to plot the requested quantity
           Plotting support is only available when requesting one or two
           quantities (i.e. when var_list is of length 1 or 2)

        nbins : int, optional
           (Only used when `plot` is True)
           Number of bins for the histograms

        plot_range : list of lists
           A list containing 2 lists of 2 elements each
           Indicates the values between which to perform the histogram,
           along the 1st axis (first list) and 2nd axis (second list)
           Default: the range is automatically determined

        use_field_mesh: bool, optional
           (Only used when `plot` is True)
           Whether to use the information of the spatial mesh (whenever
           possible) in order to choose the parameters of the histograms.
           More precisely, when this is True:
           - The extent of the histogram (along any spatial dimension) is
             automatically chosen to be roughly the extent of the spatial mesh.
           - The number of bins (along any spatial dimension) is slightly
             modified (from the value `nbins` provided by the user) so that
             the spacing of the histogram is an integer multiple of the grid
             spacing. This avoids artifacts in the plot, whenever particles
             are regularly spaced in each cell of the spatial mesh.

        histogram_deposition : string
            Either `ngp` (Nearest Grid Point) or `cic` (Cloud-In-Cell)
            When plotting the particle histogram, this determines how
            particles affects neighboring bins.
            `cic` (which is the default) leads to smoother results than `ngp`.

        **kw : dict, otional
           Additional options to be passed to matplotlib's
           hist or hist2d.

        Returns
        -------
        A list of 1darray corresponding to the data requested in `var_list`
        (one 1darray per element of 'var_list', returned in the same order)
        """
        # Check that the species required are present
        if self.avail_species is None:
            raise OpenPMDException('No particle data in this time series')
        # If there is only one species, infer that the user asks for that one
        if species is None and len(self.avail_species) == 1:
            species = self.avail_species[0]
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

        # Check the format of the particle selection
        if select is None or isinstance(select, ParticleTracker):
            pass
        elif isinstance(select, dict):
            # Dictionary: Check that all selection quantities are available
            valid_select_list = True
            for quantity in select.keys():
                if not (quantity in self.avail_record_components[species]):
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
        else:
            raise OpenPMDException("The argument `select` is erroneous.\n"
            "It should be either a dictionary or a ParticleTracker object.")

        # Find the output that corresponds to the requested time/iteration
        # (Modifies self._current_i, self.current_iteration and self.current_t)
        self._find_output(t, iteration)
        # Get the corresponding file name & open the file
        file_name = self.h5_files[self._current_i]
        file_handle = h5.File(file_name, 'r')

        # Extract the list of particle quantities
        data_list = []
        for quantity in var_list:
            data_list.append(read_species_data(
                file_handle, species, quantity, self.extensions))
        # Apply selection if needed
        if isinstance( select, dict ):
            data_list = apply_selection(
                file_handle, data_list, select, species, self.extensions)
        elif isinstance( select, ParticleTracker ):
            data_list = select.extract_tracked_particles(
                file_handle, data_list, species, self.extensions )

        # Plotting
        if plot and len(var_list) in [1, 2]:

            # Extract the weights, if they are available
            if 'w' in self.avail_record_components[species]:
                w = read_species_data(
                    file_handle, species, 'w', self.extensions)
                if isinstance( select, dict ):
                    w, = apply_selection(
                        file_handle, [w], select, species, self.extensions)
                elif isinstance( select, ParticleTracker ):
                    w, = select.extract_tracked_particles(
                        file_handle, [w], species, self.extensions )
            # Otherwise consider that all particles have a weight of 1
            else:
                w = np.ones_like(data_list[0])

            # Determine the size of the histogram bins
            # - First pick default values
            hist_range = []
            for i_data in range(len(data_list)):
                data = data_list[i_data]
                # Check if the user specified a value
                if (plot_range[i_data][0] is not None) and \
                        (plot_range[i_data][1] is not None):
                    hist_range.append( plot_range[i_data] )
                # Else use min and max of data
                elif len(data) != 0:
                    hist_range.append( [ data.min(), data.max() ] )
                else:
                    hist_range.append( [ -1., 1. ] )
            hist_bins = [ nbins for i_data in range(len(data_list)) ]
            # - Then, if required by the user, modify this values by
            #   fitting them to the spatial grid
            if use_field_mesh and self.avail_fields is not None:
                # Extract the grid resolution
                grid_size_dict, grid_range_dict = get_grid_parameters(
                    file_handle, self.avail_fields, self.fields_metadata )
                # For each direction, modify the number of bins, so that
                # the resolution is a multiple of the grid resolution
                for i_var in range(len(var_list)):
                    var = var_list[i_var]
                    if var in grid_size_dict.keys():
                        # Check that the user indeed allowed this dimension
                        # to be determined automatically
                        if (plot_range[i_var][0] is None) or \
                                (plot_range[i_var][1] is None):
                            hist_bins[i_var], hist_range[i_var] = \
                                fit_bins_to_grid(hist_bins[i_var],
                                grid_size_dict[var], grid_range_dict[var] )

            # - In the case of only one quantity
            if len(data_list) == 1:
                # Do the plotting
                self.plotter.hist1d(data_list[0], w, var_list[0], species,
                        self._current_i, hist_bins[0], hist_range[0],
                        deposition=histogram_deposition, **kw)
            # - In the case of two quantities
            elif len(data_list) == 2:
                # Do the plotting
                self.plotter.hist2d(data_list[0], data_list[1], w,
                    var_list[0], var_list[1], species,
                    self._current_i, hist_bins, hist_range,
                    deposition=histogram_deposition, **kw)
        # Close the file
        file_handle.close()

        # Output
        if output:
            return(data_list)

    def get_field(self, field=None, coord=None, t=None, iteration=None,
                  m='all', theta=0., slicing=0., slicing_dir='y',
                  output=True, plot=False,
                  plot_range=[[None, None], [None, None]], **kw):
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

        plot_range : list of lists
           A list containing 2 lists of 2 elements each
           Indicates the values between which to clip the plot,
           along the 1st axis (first list) and 2nd axis (second list)
           Default: plots the full extent of the simulation box

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
        if self.fields_metadata[field]['type'] == 'vector':
            available_coord = ['x', 'y', 'z']
            if self.fields_metadata[field]['geometry'] == 'thetaMode':
                available_coord += ['r', 't']
            if coord not in available_coord:
                coord_list = '\n - '.join(available_coord)
                raise OpenPMDException(
                    "The field %s is a vector field, \nbut the `coord` "
                    "argument is missing or erroneous.\nThe available "
                    "coordinates are: \n - %s\nPlease set the `coord` "
                    "argument accordingly." % (field, coord_list))
        # Check the mode (for thetaMode)
        if self.fields_metadata[field]['geometry'] == "thetaMode":
            avail_circ_modes = self.fields_metadata[field]['avail_circ_modes']
            if str(m) not in avail_circ_modes:
                mode_list = '\n - '.join(avail_circ_modes)
                raise OpenPMDException(
                    "The requested mode '%s' is not available.\n"
                    "The available modes are: \n - %s" % (m, mode_list))

        # Find the output that corresponds to the requested time/iteration
        # (Modifies self._current_i, self.current_iteration and self.current_t)
        self._find_output(t, iteration)
        # Get the corresponding filename
        filename = self.h5_files[self._current_i]

        # Find the proper path for vector or scalar fields
        if self.fields_metadata[field]['type'] == 'scalar':
            field_path = field
            field_label = field
        elif self.fields_metadata[field]['type'] == 'vector':
            field_path = os.path.join(field, coord)
            field_label = field + coord

        # Get the field data
        geometry = self.fields_metadata[field]['geometry']
        axis_labels = self.fields_metadata[field]['axis_labels']
        # - For 1D
        if geometry == "1dcartesian":
            F, info = read_field_1d(filename, field_path, axis_labels)
        # - For 2D
        elif geometry == "2dcartesian":
            F, info = read_field_2d(filename, field_path, axis_labels)
        # - For 3D
        elif geometry == "3dcartesian":
            F, info = read_field_3d(
                filename, field_path, axis_labels, slicing, slicing_dir)
        # - For thetaMode
        elif geometry == "thetaMode":
            if (coord in ['x', 'y']) and \
                    (self.fields_metadata[field]['type'] == 'vector'):
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
        if (geometry == "3dcartesian") and (slicing is None):
            plot = False
        if plot:
            if geometry == "1dcartesian":
                self.plotter.show_field_1d(F, info, field_label,
                self._current_i, plot_range=plot_range, **kw)
            else:
                self.plotter.show_field_2d(F, info, slicing_dir, m,
                        field_label, geometry, self._current_i,
                        plot_range=plot_range, **kw)

        # Return the result
        return(F, info)

    def _find_output(self, t, iteration):
        """
        Find the output that correspond to the requested `t` or `iteration`
        Modify self._current_i accordingly.

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
                self._current_i = 0
            elif t > self.tmax:
                self._current_i = len(self.t) - 1
            # Find the last existing output
            else:
                self._current_i = self.t[self.t <= t].argmax()
        # If an iteration is requested
        elif (iteration is not None):
            if (iteration in self.iterations):
                # Get the index that corresponds to this iteration
                self._current_i = abs(iteration - self.iterations).argmin()
            else:
                iter_list = '\n - '.join([str(it) for it in self.iterations])
                print("The requested iteration '%s' is not available.\nThe "
                      "available iterations are: \n - %s\nThe first iteration "
                      "is used instead." % (iteration, iter_list))
                self._current_i = 0
        else:
            pass  # self._current_i retains its previous value

        # Register the value in the object
        self.current_t = self.t[self._current_i]
        self.current_iteration = self.iterations[self._current_i]
