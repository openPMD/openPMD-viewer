"""
This file is part of the openPMD-viewer.

This file contains diagnostics relevant for laser-plasma acceleration

Copyright 2015-2016, openPMD-viewer contributors
Authors: Soeren Jalas, Remi Lehe
License: 3-Clause-BSD-LBNL
"""

# Class that inherits from OpenPMDTimeSeries, and implements
# some standard diagnostics (emittance, etc.)
from openpmd_viewer import OpenPMDTimeSeries, FieldMetaInformation
import numpy as np
import scipy.constants as const
from scipy.optimize import curve_fit
from openpmd_viewer.openpmd_timeseries.utilities import sanitize_slicing
from openpmd_viewer.openpmd_timeseries.plotter import check_matplotlib
from scipy.signal import hilbert
try:
    import matplotlib.pyplot as plt
except ImportError:
    # Any error will be caught later by `check_matplotlib`
    pass


class LpaDiagnostics( OpenPMDTimeSeries ):

    def __init__( self, path_to_dir, check_all_files=True, backend=None ):
        """
        Initialize an OpenPMD time series with various methods to diagnose the
        data

        Parameter
        ---------
        path_to_dir : string
            The path to the directory where the openPMD files are.
            For the moment, only HDF5 files are supported. There should be
            one file per iteration, and the name of the files should end
            with the iteration number, followed by '.h5' (e.g. data0005000.h5)

        check_all_files: bool, optional
            Check that all the files in the timeseries are consistent
            (i.e. that they contain the same fields and particles,
            with the same metadata)
            For fast access to the files, this can be changed to False.

        backend: string
            Backend to be used for data reading. Can be `openpmd-api`
            or `h5py`. If not provided will use `openpmd-api` if available
            and `h5py` otherwise.
        """
        OpenPMDTimeSeries.__init__( self, path_to_dir,
                                    check_all_files=check_all_files, backend=backend )

    def get_energy_spread( self, t=None, iteration=None, species=None,
                        select=None, center='mean', width='std', property='energy' ):
        """
        Calculate the central energy and energy spread according to the
        particle weights

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        species : str
            Particle species to use for calculations

        select : dict, optional
            Either None or a dictionary of rules
            to select the particles, of the form
            'x' : [-4., 10.]   (Particles having x between -4 and 10 meters)
            'z' : [0, 100] (Particles having z between 0 and 100 meters)

        center : str
            Method to find the central energy of the particle distribution
            Should be one of

            - 'mean'
            - 'median'

        width : str
            Method to find the energy spread of the particle distribution
            Should be one of

            - 'std'
            - 'mad'

        property : str
            Unit of energy. Should be one of

            - 'energy' returns energy in MeV
            - 'gamma' returns Lorentz factor

        Returns
        -------
        A tuple of floats with:
        - central energy
        - energy spread
        Returns NaN if particle selection is empty
        """
        # Get particle data
        ux, uy, uz, w, m = self.get_particle(
            var_list=['ux', 'uy', 'uz', 'w', 'mass'], select=select,
            species=species, t=t, iteration=iteration)
        if len(w) == 0:
            # Return NaN if no particles are found
            return np.nan, np.nan
        # Calculate Lorentz factor and energy for all particles
        gamma = np.sqrt(1 + ux ** 2 + uy ** 2 + uz ** 2)
        if property == 'energy':
            prop = (gamma - 1) * m * const.c ** 2 / const.e * 1e-6

        elif property == 'gamma':
            prop = gamma
        else:
            raise ValueError('Invalid output property: %s'%property)

        # Calculate weighted center
        if center == 'mean':
            p_center = w_ave(prop, w)
        elif center == 'median':
            p_center = w_median(prop, w)
        else:
            raise ValueError('Invalid center property: %s' % center)

        # Calculate weighted width
        if width == 'std':
            p_width = w_std(prop, w)
        elif width == 'mad':
            p_width = w_mad(prop, w)
        else:
            raise ValueError('Invalid width property: %s' % width)

        return p_center, p_width

    def get_mean_gamma( self, t=None, iteration=None, species=None,
                        select=None ):
        """
        Calculate the mean gamma and standard deviation according to the
        particle weights

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        species : string
            Particle species to use for calculations

        select : dict, optional
            Either None or a dictionary of rules
            to select the particles, of the form
            'x' : [-4., 10.]   (Particles having x between -4 and 10 meters)
            'z' : [0, 100] (Particles having z between 0 and 100 meters)

        Returns
        -------
        A tuple of floats with:
        - mean weighted gamma
        - weighted standard deviation of gamma
        Returns NaN if particle selection is empty
        """
        # Get particle data
        ux, uy, uz, w = self.get_particle(
            var_list=['ux', 'uy', 'uz', 'w'], select=select,
            species=species, t=t, iteration=iteration )
        if len(w) == 0:
            # Return NaN if no particles are found
            return np.nan, np.nan
        # Calculate Lorentz factor for all particles
        gamma = np.sqrt(1 + ux ** 2 + uy ** 2 + uz ** 2)
        # Calculate weighted mean and average
        try:
            # Calculate mean gamma for selected particles
            mean_gamma = np.average(gamma, weights=w)
        except ZeroDivisionError:
            # If selection is empty or all particles have weight zero,
            # return NaN
            mean_gamma = np.nan
        std_gamma = w_std(gamma, w)
        # Return the result
        return mean_gamma, std_gamma

    def get_sigma_gamma_slice(self, dz, t=None, iteration=None, species=None,
                              select=None, plot=False, **kw):
        """
        Calculate the standard deviation of gamma for particles in z-slices of
        width dz

        Parameters
        ----------
        dz : float (in micrometers)
            Width of slices in which to calculate sigma gamma

        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        species : string
            Particle species to use for calculations

        select : dict, optional
            Either None or a dictionary of rules
            to select the particles, of the form
            'x' : [-4., 10.]   (Particles having x between -4 and 10 meters)
            'z' : [0, 100] (Particles having z between 0 and 100 meters)

        plot : bool, optional
           Whether to plot the requested quantity

        **kw : dict, otional
           Additional options to be passed to matplotlib's `plot` method

        Returns
        -------
        A tuple of arrays:
        - Sigma gamma in each slice
        - Central z position of each slice

        """
        z, uz, ux, uy, w = self.get_particle(t=t, species=species,
            select=select, var_list=['z', 'uz', 'ux', 'uy', 'w'],
            iteration=iteration)
        # Calculate gamma of each particle
        gamma = np.sqrt(1 + (uz ** 2 + ux ** 2 + uy ** 2))
        z0 = min(z)
        zend = max(z)
        N = int((zend - z0) / dz)
        spreads = np.zeros(N + 1)
        z_pos = np.linspace(z0, zend, N + 1)
        zi = z0 + dz / 2.
        i = 0
        # Iterate over slices and calculate sigma gamma
        while zi < zend:
            z_filter = (z > zi - dz / 2.) & (z < zi + dz / 2.)
            spreads[i] = w_std(gamma[z_filter], w[z_filter])
            zi += dz
            i += 1
        # Plot the result if needed
        if plot:
            check_matplotlib()
            iteration = self.iterations[ self._current_i ]
            time_s = self.t[ self._current_i ]
            plt.plot(z_pos, spreads, **kw)
            plt.title("Slice energy spread at %.2e s   (iteration %d)"
                % (time_s, iteration), fontsize=self.plotter.fontsize)
            plt.xlabel('$z \;(m)$', fontsize=self.plotter.fontsize)
            plt.ylabel('$\sigma_\gamma (\Delta_z=%s m)$' % dz,
                       fontsize=self.plotter.fontsize)
        return(spreads, z_pos)

    def get_charge( self, t=None, iteration=None, species=None, select=None ):
        """
        Calculate the charge of the selcted particles.

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        species : string
            Particle species to use for calculations

        select : dict, optional
            Either None or a dictionary of rules
            to select the particles, of the form
            'x' : [-4., 10.]   (Particles having x between -4 and 10 meters)
            'z' : [0, 100] (Particles having z between 0 and 100 meters)

        Returns
        -------
        A float with the electric charge of the selected particles in Coulomb
        """
        # Get particle data
        w, q = self.get_particle( var_list=['w', 'charge'], species=species,
            select=select, t=t, iteration=iteration )
        # Calculate charge
        charge = np.sum(w * q)
        # Return the result
        return( charge )

    def get_divergence( self, t=None, iteration=None, species=None,
                        select=None ):
        """
        Calculate the divergence of the selected particles.

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        species : string
            Particle species to use for calculations

        select : dict, optional
            Either None or a dictionary of rules
            to select the particles, of the form
            'x' : [-4., 10.]   (Particles having x between -4 and 10 meters)
            'z' : [0, 100] (Particles having z between 0 and 100 meters)

        Returns
        -------
        A tuple with:
        - divergence in x plane in rad
        - divergence in y plane in rad
        Returns NaN if particle selection is empty
        """
        # Get particle data
        ux, uy, uz, w = self.get_particle( var_list=['ux', 'uy', 'uz', 'w'],
                                           t=t, iteration=iteration,
                                           species=species, select=select )
        if len(w) == 0:
            # Return NaN if no particles are found
            return np.nan, np.nan
        # Calculate divergence
        div_x = w_std( np.arctan2(ux, uz), w )
        div_y = w_std( np.arctan2(uy, uz), w )
        # Return the result
        return div_x, div_y

    def get_emittance(self, t=None, iteration=None, species=None,
                      select=None, kind='normalized', description='projected',
                      nslices=0, beam_length=None):
        """
        Calculate the RMS emittance.
        (See K Floetmann: Some basic features of beam emittance. PRSTAB 2003)

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        species : string
            Particle species to use for calculations

        select : dict, optional
            Either None or a dictionary of rules to select the particles,
            of the form:
            'x' : [-4., 10.]   (Particles having x between -4 and 10 meters)
            'z' : [0, 100] (Particles having z between 0 and 100 meters).

        kind : string, optional
            Kind of emittance to be computed. Can be 'normalized' or 'trace'.

        description : string, optional
            Type of emittance to be computed. Available options:
               - 'projected' : projected emittance
               - 'all-slices' : emittance within slices taken along the z
                              direction
               - 'slice-averaged' : slice emittance averaged over all slices.

        nslices : integer, optional
            Number of slices to compute slice emittance. Required if
            description='slice-average' or 'all-slices'.

        beam_length : float (in meters), optional
            Beam length, used to calculate slice positions when nslices>1.
            By default, it is 4 times the standard deviation in z.

        Returns
        -------
        If description='projected' or 'slice-averaged':
            - beam emittance in the x plane (pi m rad)
            - beam emittance in the y plane (pi m rad)
        If description='all-slices':
            - A 1d array with beam emittance in the x plane
              (pi m rad) for each slice
            - A 1d array with beam emittance in the y plane
              (pi m rad) for each slice
            - A 1d array with number of electrons in each slice
            - A 1d array with slice centers
        """
        if kind not in ['normalized', 'trace']:
            raise ValueError('Argument `kind` not recognized.')
        if description not in ['projected', 'slice-averaged', 'all-slices']:
            raise ValueError('Argument `description` not recognized.')
        # Wheter to compute slice emittance
        do_slice_emittance = ( description in ['slice-averaged',
                                               'all-slices'] )
        if do_slice_emittance and not nslices > 0:
            raise ValueError(
                'nslices must be given if `description`=' + description + '.')
        # Get particle data
        x, y, z, ux, uy, uz, w = self.get_particle(
            var_list=['x', 'y', 'z', 'ux', 'uy', 'uz', 'w'], t=t,
            iteration=iteration, species=species, select=select )
        # Normalized or trace-space emittance
        if kind == 'normalized':
            ux = ux
            uy = uy
        if kind == 'trace':
            ux = ux / uz
            uy = uy / uz

        if do_slice_emittance:
            # Get slice locations
            zavg = w_ave(z, w)
            z = z - zavg
            if beam_length is None:
                std = w_std(z, w)
                beam_length = 4 * std
            bins = np.linspace( -beam_length / 2, beam_length / 2, nslices )
            binwidth = .5 * bins[1] - .5 * bins[0]
            slice_centers = bins + .5 * binwidth
            # Initialize slice emittance arrays
            emit_slice_x = np.zeros(len(bins[:-1]))
            emit_slice_y = np.zeros(len(bins[:-1]))
            slice_weights = np.zeros(len(bins[:-1]))
            # Loop over slices
            for count, leftedge in enumerate(bins[:-1]):
                # Get emittance in this slice
                current_slice = (np.abs(z - slice_centers[count]) <= binwidth)
                slice_weights[count] = np.sum(w[current_slice])
                if slice_weights[count] > 0:
                    emit_x, emit_y = emittance_from_coord(
                        x[current_slice], y[current_slice],
                        ux[current_slice], uy[current_slice],
                        w[current_slice])
                    emit_slice_x[count] = emit_x
                    emit_slice_y[count] = emit_y
            if description == 'all-slices':
                return (emit_slice_x, emit_slice_y,
                    slice_weights, slice_centers)
            else:
                emit_x = w_ave(emit_slice_x, slice_weights)
                emit_y = w_ave(emit_slice_y, slice_weights)
                return emit_x, emit_y
        else:
            return emittance_from_coord(x, y, ux, uy, w)

    def get_current( self, t=None, iteration=None, species=None, select=None,
                     bins=100, plot=False, **kw ):
        """
        Calculate the electric current along the z-axis for selected particles.

        Parameters
        ----------
         t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user

        species : string
            Particle species to use for calculations

        select : dict, optional
            Either None or a dictionary of rules
            to select the particles, of the form
            'x' : [-4., 10.]   (Particles having x between -4 and 10 meters)
            'z' : [0, 100] (Particles having z between 0 and 100 meters)

        bins : int, optional
            Number of bins along the z-axis in which to calculate the current

        plot : bool, optional
           Whether to plot the requested quantity

        **kw : dict, otional
           Additional options to be passed to matplotlib's `plot` method
        Returns
        -------
        A tuple of arrays containig
        - The current in each bin in Ampere
        - A FieldMetaInformation object
          (see object's docstring for more details)

        """
        # Get particle data
        z, uz, uy, ux, w, q = self.get_particle(
            var_list=['z', 'uz', 'uy', 'ux', 'w', 'charge'],
            t=t, iteration=iteration,
            species=species, select=select )
        # Length to be seperated in bins

        if w.size > 0:
            min_z = np.min(z)
            len_z = np.max(z) - min_z
            # Calculate Lorentz factor for all particles
            gamma = np.sqrt(1 + ux ** 2 + uy ** 2 + uz ** 2)
            # Calculate particle velocities
            vz = uz / gamma * const.c
            # Length to be seperated in bins
            len_z = np.max(z) - np.min(z)
            vzq_sum, _ = np.histogram(z, bins=bins, weights=(vz * w * q))
            # Calculate the current in each bin
            current = np.abs(vzq_sum * bins / (len_z))
        else:
            current = np.zeros(bins)
            len_z = 0
            min_z = 0
        # Info object with central position of the bins
        info = FieldMetaInformation( {0: 'z'}, current.shape,
            grid_spacing=(len_z / bins, ), grid_unitSI=1,
            global_offset=(min_z + len_z / bins / 2,), position=(0,))
        # Plot the result if needed
        if plot:
            check_matplotlib()
            iteration = self.iterations[ self._current_i ]
            time_s = self.t[ self._current_i ]
            plt.plot( info.z, current, **kw)
            plt.title("Current at %.2e s   (iteration %d)"
                % (time_s, iteration ), fontsize=self.plotter.fontsize)
            plt.xlabel('$z \;(m)$', fontsize=self.plotter.fontsize)
            plt.ylabel('$I \;(A)$', fontsize=self.plotter.fontsize)
        # Return the current and bin centers
        return(current, info)

    def get_laser_envelope( self, t=None, iteration=None, pol=None,
                            laser_propagation='z',
                            m='all', theta=0, slice_across=None,
                            slice_relative_position=None, plot=False,
                            plot_range=[[None, None], [None, None]], **kw ):
        """
        Calculate a laser field by filtering out high frequencies. Can either
        return the envelope slice-wise or a full 2D envelope.

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        pol : string
            Polarization of the field. Options are 'x', 'y', 'z'

        laser_propagation : string, optional
            Coordinate along which laser field propagates.
            Default is 'z'.

        m : int or str, optional
           Only used for thetaMode geometry
           Either 'all' (for the sum of all the modes)
           or an integer (for the selection of a particular mode)

        theta : float, optional
           Only used for thetaMode geometry
           The angle of the plane of observation, with respect to the x axis

        slice_across : str or list of str, optional
           Direction(s) across which the data should be sliced
           + In cartesian geometry, elements can be:
               - 1d: 'z'
               - 2d: 'x' and/or 'z'
               - 3d: 'x' and/or 'y' and/or 'z'
           + In cylindrical geometry, elements can be 'r' and/or 'z'
           Returned array is reduced by 1 dimension per slicing.
           If slice_across is None, the full grid is returned.
           Default is None.

        slice_relative_position : float or list of float, optional
           Number(s) between -1 and 1 that indicate where to slice the data,
           along the directions in `slice_across`
           -1 : lower edge of the simulation box
           0 : middle of the simulation box
           1 : upper edge of the simulation box
           Default is 0.

        plot : bool, optional
           Whether to plot the requested quantity

        plot_range : list of lists
           A list containing 2 lists of 2 elements each
           Indicates the values between which to clip the plot,
           along the 1st axis (first list) and 2nd axis (second list)
           Default: plots the full extent of the simulation box

        **kw : dict, otional
           Additional options to be passed to matplotlib's `plot`(1D) or
           `imshow` (2D) method

        Returns
        -------
        A tuple with:
        - Envelope data (1D or 2D array)
        - A FieldMetaInformation object
        """
        # Check if polarization has been entered
        if pol not in ['x', 'y', 'z']:
            raise ValueError('The `pol` argument is missing or erroneous.')

        # Prevent slicing across `laser_propagation`, when extracting the raw electric field
        # (`laser_propagation` axis is needed for calculation of envelope)
        # but record whether the user asked for slicing across `laser_propagation`,
        # and whether a corresponding `slice_relative_position` coordinate
        # along `laser_propagation` was given, so as to perform this slicing later in this function.
        slicing_coord_laser = None
        if slice_across is not None:
            slice_across, slice_relative_position = \
                sanitize_slicing(slice_across, slice_relative_position)
            if laser_propagation in slice_across:
                index_slicing_coord_laser = slice_across.index(laser_propagation)
                slice_across.pop(index_slicing_coord_laser)
                slicing_coord_laser = slice_relative_position.pop(index_slicing_coord_laser)
        # Get field data, and perform Hilbert transform
        field, info = self.get_field( t=t, iteration=iteration, field='E',
                              coord=pol, theta=theta, m=m,
                              slice_across=slice_across,
                              slice_relative_position=slice_relative_position )
        inverted_axes_dict = {info.axes[key]: key for key in info.axes.keys()}
        e_complx = hilbert(field, axis=inverted_axes_dict[laser_propagation])
        envelope = np.abs(e_complx)
        # If the user asked for slicing along `laser_propagation`, do it now
        if slicing_coord_laser is not None:
            slicing_index = inverted_axes_dict[laser_propagation]
            coord_array = getattr( info, laser_propagation )
            # Number of cells along the slicing direction
            n_cells = len(coord_array)
            # Index of the slice (prevent stepping out of the array)
            i_cell = int( 0.5 * (slicing_coord_laser + 1.) * n_cells )
            i_cell = max( i_cell, 0 )
            i_cell = min( i_cell, n_cells - 1)
            envelope = np.take( envelope, [i_cell], axis=slicing_index )
            envelope = np.squeeze(envelope)
            # Remove the sliced labels from the FieldMetaInformation
            info._remove_axis(laser_propagation)

        # Plot the result if needed
        if plot:
            geometry = self.fields_metadata['E']['geometry']
            field_label = 'E%s (envelope)' %pol
            if envelope.ndim == 1:
                self.plotter.show_field_1d(envelope, info, field_label,
                self._current_i, plot_range=plot_range, **kw)
            elif envelope.ndim == 2:
                self.plotter.show_field_2d(envelope, info, slice_across, m,
                    field_label, geometry, self._current_i,
                    plot_range=plot_range, **kw)
        # Return the result
        return( envelope, info )

    def get_main_frequency( self, t=None, iteration=None, pol=None, m='all',
                            method='max'):
        """
        Calculate the angular frequency of a laser pulse.

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        pol : string
            Polarization of the field. Options are 'x', 'y'

        m : int or str, optional
           Only used for thetaMode geometry
           Either 'all' (for the sum of all the modes)
           or an integer (for the selection of a particular mode)

        method : string, optional
            Method which is used to calculate the frequency of the pulse
            'fit' : Fit a Gaussian curve to find central frequency
            'max' : Take frequency with highest intensity in the spectrum

        Returns
        -------
        A float with mean angular frequency
        """
        # Extract the spectrum
        spectrum, info = self.get_spectrum( t, iteration, pol, m )

        # Calculate the main frequency
        # Use the maximum
        i_max = np.argmax( spectrum )
        omega0 = info.omega[i_max]
        if method == 'max':
            return( omega0 )
        # Gaussian fit
        elif method == 'fit':
            # Guess start values for fit
            f0 = omega0
            fmax = np.amax( spectrum )
            fsigma = info.omegamax - info.omegamin
            params, _ = curve_fit( gaussian_profile, info.omega,
                                   spectrum, p0=[f0, fmax, fsigma])
            return( params[0] )
        else:
            raise ValueError('Unknown method: {:s}'.format(method))

    def get_spectrum( self, t=None, iteration=None, pol=None,
                      m='all', plot=False, **kw ):
        """
        Return the spectrum of the laser
        (Absolute value of the Fourier transform of the fields.)

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        pol : string
            Polarization of the field. Options are 'x', 'y'

        m : int or str, optional
           Only used for thetaMode geometry
           Either 'all' (for the sum of all the modes)
           or an integer (for the selection of a particular mode)

        plot: bool, optional
           Whether to plot the data

        **kw : dict, otional
           Additional options to be passed to matplotlib's `plot` method

        Returns
        -------
        A tuple with:
            - The 1D spectrum on axis
            - A FieldMetaInformation object
        """
        # Check if polarization has been entered
        if pol not in ['x', 'y']:
            raise ValueError('The `pol` argument is missing or erroneous.')
        # Get a lineout along the 'z' axis,
        slice_across = self._get_slicing_for_longitudinal_lineout()

        # Get field data
        field1d, info = self.get_field( t=t, iteration=iteration, field='E',
                                coord=pol, m=m, slice_across=slice_across )
        # FFT of 1d data
        dt = (info.z[1] - info.z[0]) / const.c  # Integration step for the FFT
        fft_field = np.fft.fft(field1d) * dt
        # Take half of the data (positive frequencies only)
        spectrum = abs( fft_field[ : int( len(fft_field) / 2 ) ] )
        # Create a FieldMetaInformation object
        T = (info.zmax - info.zmin) / const.c
        spect_info = FieldMetaInformation( {0: 'omega'}, spectrum.shape,
            grid_spacing=( 2 * np.pi / T, ), grid_unitSI=1,
            global_offset=(0,), position=(0,))

        # Plot the field if required
        if plot:
            check_matplotlib()
            iteration = self.iterations[ self._current_i ]
            time_s = self.t[ self._current_i ]
            plt.plot( spect_info.omega, spectrum, **kw )
            plt.xlabel('$\omega \; (rad.s^{-1})$',
                       fontsize=self.plotter.fontsize )
            plt.ylabel('Spectrum', fontsize=self.plotter.fontsize )
            plt.title("Spectrum at %.2e s   (iteration %d)"
                % (time_s, iteration ), fontsize=self.plotter.fontsize)
        return( spectrum, spect_info )

    def get_a0( self, t=None, iteration=None, pol=None ):
        """
        Gives the laser strength a0 given by a0 = Emax * e / (me * c * omega)

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        pol : string
            Polarization of the field. Options are 'x', 'y'

        Returns
        -------
        Float with normalized vector potential a0
        """
        # Get a lineout along the 'z' axis,
        slice_across = self._get_slicing_for_longitudinal_lineout()

        # Get the peak field from field envelope
        Emax = np.amax(self.get_laser_envelope(t=t, iteration=iteration,
                                       pol=pol, slice_across=slice_across)[0])
        # Get mean frequency
        omega = self.get_main_frequency(t=t, iteration=iteration, pol=pol)
        # Calculate a0
        a0 = Emax * const.e / (const.m_e * const.c * omega)
        return( a0 )

    def get_ctau( self, t=None, iteration=None, pol=None, method='fit' ):
        """
        Calculate the length of a (gaussian) laser pulse. Here 'length' means
        the 'longitudinal waist' (i.e sqrt(2) * sigma_z).

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        pol : string
            Polarization of the field. Options are 'x', 'y'

        method : str, optional
           The method which is used to compute ctau
           'fit': Gaussian fit of the longitudinal profile
           'rms': RMS radius, weighted by the longitudinal profile
           ('rms' tends to give more weight to the "wings" of the pulse)

        Returns
        -------
        Float with ctau in meters
        """
        # Get a lineout along the 'z' axis,
        slice_across = self._get_slicing_for_longitudinal_lineout()

        # Get the field envelope
        E, info = self.get_laser_envelope(t=t, iteration=iteration,
                                            pol=pol, slice_across=slice_across)
        # Calculate ctau with RMS value
        ctau = np.sqrt(2) * w_std(info.z, E)
        if method == 'rms':
            return( ctau )

        # Calculate ctau from Gaussian fit
        elif method == 'fit':
            # Start value for E0 and x0 in fit
            E0 = np.amax( E )
            z0 = info.z[np.argmax( E )]
            # Perform the fit
            params, _ = curve_fit( gaussian_profile, info.z,
                                   E, p0=[ z0, E0, ctau ])
            return( params[2] )

        else:
            raise ValueError('Unknown method: {:s}'.format(method))

    def get_laser_waist( self, t=None, iteration=None, pol=None, theta=0,
                         method='fit' ):
        """
        Calculate the waist of a (gaussian) laser pulse. ( sqrt(2) * sigma_r)

        In 3D, this function takes a slice across `y`, and thus computes the
        waist in the `x-z` plane.

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        pol : string
            Polarization of the field. Options are 'x', 'y'

        theta : float, optional
           Only used for thetaMode geometry
           The angle of the plane of observation, with respect to the x axis

        method : str, optional
           The method which is used to compute the waist
           'fit': Gaussian fit of the transverse profile
           'rms': RMS radius, weighted by the transverse profile
           ('rms' tends to give more weight to the "wings" of the pulse)

        Returns
        -------
        Float with laser waist in meters
        """
        # In 3D, slice across 'y' by default
        geometry = self.fields_metadata['E']['geometry']
        if geometry == '3dcartesian':
            slice_across = 'y'
        else:
            slice_across = None

        # Get the field envelope (as 2D array)
        field, info = self.get_laser_envelope(t=t, iteration=iteration,
                         pol=pol, slice_across=slice_across, theta=theta)
        assert field.ndim == 2
        # Find the indices of the maximum field, and
        # pick the corresponding transverse slice
        itrans_max, iz_max = np.unravel_index(
            np.argmax( field ), field.shape )
        trans_slice = field[ :, iz_max ]
        # Get transverse positons
        trans_pos = getattr(info, info.axes[0])

        # Compute waist with RMS value
        # (serves as initial guess when method=='fit')
        w0 = np.sqrt(2) * w_std(trans_pos, trans_slice)
        if method == 'rms':
            return( w0 )

        # Compute waist with Gaussian fit
        elif method == 'fit':
            # Get initial guess for the amplitude
            E0 = field[ itrans_max, iz_max ]
            # Assume that the pulse is centered
            x0 = 0
            # Perform the fit
            params, _ = curve_fit( gaussian_profile, trans_pos,
                                   trans_slice, p0=[x0, E0, w0 ])
            return( params[2] )

        else:
            raise ValueError('Unknown method: {:s}'.format(method))

    def get_spectrogram( self, t=None, iteration=None, pol=None,
                          plot=False, **kw ):
        """
        Calculates the spectrogram of a laserpulse, by the FROG method.

        Mathematically:
        $$ s(\omega, \tau) = | \int_{-\infty}^{\infty} E(t) |E(t-\tau)|^2
            \exp( -i\omega t) dt |^2 $$
        See Trebino, R: Frequency Resolved Optical Gating: The measurements of
        Ultrashort Laser Pulses: year 2000: formula 5.2

        The time is centered around the laser pulse.

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        pol : string
            Polarization of the laser field. Options are 'x', 'y'

        plot: bool, optional
            Whether to plot the spectrogram

        **kw : dict, otional
           Additional options to be passed to matplotlib's `imshow` method

        Returns
        -------
        - A 2d array with spectrogram
        - info : a FieldMetaInformation object
           (see the corresponding docstring)
        """
        # Get a lineout along the 'z' axis
        slice_across = self._get_slicing_for_longitudinal_lineout()

        # Get the field envelope
        env, _ = self.get_laser_envelope(t=t, iteration=iteration,
                                    pol=pol, slice_across=slice_across)
        # Get the field
        E, info = self.get_field( t=t, iteration=iteration, field='E',
                                    coord=pol, slice_across=slice_across)
        Nz = len(E)
        # Get time domain of the data
        tmin = info.zmin / const.c
        tmax = info.zmax / const.c
        T = tmax - tmin
        dt = T / Nz
        # Normalize the Envelope
        env /= np.sqrt(np.trapz(env ** 2, dx=dt))
        # Allocate array for the gating function and the spectrogran
        E_shift = np.zeros_like(E)
        spectrogram = np.zeros((2 * Nz, Nz))
        # Loop over the time variable of the spectrogram
        for i in range( Nz * 2):
            itau = i % Nz
            # Shift the E field and fill the rest with zeros
            if i < Nz:
                E_shift[:itau] = env[ Nz - itau: Nz]
                E_shift[itau:] = 0
            else:
                E_shift[itau:] = env[: Nz - itau]
                E_shift[:itau] = 0
            EE = E * E_shift ** 2
            fft_EE = np.fft.fft(EE)
            spectrogram[i, :] = np.abs(fft_EE) ** 2
        # Rotate and flip array to have input form of imshow
        spectrogram = np.flipud(np.rot90(spectrogram[:, int(Nz / 2):]))
        # Find the time at which the wigner transform is the highest
        maxi, maxj = np.unravel_index(spectrogram.argmax(), spectrogram.shape)
        tmin = -(T - T / spectrogram.shape[1] * maxj)
        info = FieldMetaInformation( {0: 'omega', 1: 't'}, spectrogram.shape,
            grid_spacing=( 2 * np.pi / T, dt / 2. ), grid_unitSI=1,
            global_offset=(0, tmin), position=(0, 0))

        # Plot the result if needed
        if plot:
            check_matplotlib()
            iteration = self.iterations[ self._current_i ]
            time_s = self.t[ self._current_i ]
            plt.imshow( spectrogram, extent=info.imshow_extent, aspect='auto',
                        **kw)
            plt.title("Spectrogram at %.2e s   (iteration %d)"
                % (time_s, iteration ), fontsize=self.plotter.fontsize)
            plt.xlabel('$t \;(s)$', fontsize=self.plotter.fontsize )
            plt.ylabel('$\omega \;(rad.s^{-1})$',
                       fontsize=self.plotter.fontsize )
        return( spectrogram, info )


    def _get_slicing_for_longitudinal_lineout(self):
        """
        Return the `slice_across` argument which results in a 1D slice
        along `z`, for the current geometry.
        """
        geometry = self.fields_metadata['E']['geometry']
        if geometry == "2dcartesian":
            return 'x'
        elif geometry == "3dcartesian":
            return ['x', 'y']
        elif geometry == "thetaMode":
            return 'r'
        else:
            raise ValueError('Unknown geometry: %s' %geometry)


def w_ave( a, weights ):
    """
    Calculate the weighted average of array `a`

    Parameters
    ----------
    a : 1d array
        Calculate the weighted average for these a.

    weights : 1d array
        An array of weights for the values in a.

    Returns
    -------
    Float with the weighted average
    Returns nan if input array is empty
    """
    # Check if input contains data
    if not np.any(weights) and not np.any(a):
        # If input is empty return NaN
        return np.nan
    else:
        # Calculate the weighted average
        average = np.average(a, weights=weights)
        return( average )


def w_std( a, weights ):
    """
    Calculate the weighted standard deviation.

    Parameters
    ----------
    a : array_like
        Calculate the weighted standard deviation for these a.

    weights : array_like
        An array of weights for the values in a.

    Returns
    -------
    Float with the weighted standard deviation.
    Returns nan if input array is empty
    """
    # Check if input contains data
    if not np.any(weights) and not np.any(a):
        # If input is empty return NaN
        return np.nan
    else:
        # Calculate the weighted standard deviation
        average = np.average(a, weights=weights)
        variance = np.average((a - average) ** 2, weights=weights)
        return( np.sqrt(variance) )

def w_median(a, weights):
    """
    Compute the weighted median of a 1D numpy array.
    Parameters
    ----------
    a : ndarray
        Input array (one dimension).
    weights : ndarray
        Array with the weights of the same size of `data`.
    Returns
    -------
    median : float
        The output value.
    """
    quantile = .5
    if not isinstance(a, np.matrix):
        a = np.asarray(a)
    if not isinstance(weights, np.matrix):
        weights = np.asarray(weights)
    if a.shape != weights.shape:
        raise TypeError("the length of data and weights must be the same")
    ind_sorted = np.argsort(a)
    sorted_data = a[ind_sorted]
    sorted_weights = weights[ind_sorted]

    Sn = np.cumsum(sorted_weights)
    # Center and normalize the cumsum (i.e. divide by the total sum)
    Pn = (Sn - 0.5 * sorted_weights) / Sn[-1]
    # Get the value of the weighted median
    return np.interp(quantile, Pn, sorted_data)

def w_mad(a, w):
    """
    Compute the weighted median absolute deviation of a 1D numpy array.
    Parameters
    ----------
    a : ndarray
        Input array (one dimension).
    weights : ndarray
        Array with the weights of the same size of `data`.
    Returns
    -------
    mad : float
        The output value.
    """
    med = w_median(a, w)
    mad = w_median(np.abs(a - med), w)
    return mad

def gaussian_profile( x, x0, E0, w0 ):
    """
    Returns a Gaussian profile with amplitude E0 and waist w0.
    (Used in order to fit the transverse laser profile and find the waist.)

    Parameters
    ----------
    x: 1darray of floats
        An array of transverse positions (in meters)

    x0: float
        Position of the peak of the profile

    E0: float
        The amplitude at the peak of the profile

    w0: float
        The waist of the profile

    Returns
    -------
    A 1darray of floats, of the same length as x
    """
    return( E0 * np.exp( -(x - x0) ** 2 / w0 ** 2 ) )


def emittance_from_coord(x, y, ux, uy, w):
    """
    Calculate emittance from arrays of particle coordinates.

    Parameters
    ----------
    x : arrays of floats
        x position of particles
    y : arrays of floats
        y position of particles
    ux : arrays of floats
        ux normalized momentum of particles
    uy : arrays of floats
        uy normalized momentum of particles
    w : arrays of floats
        Particle weights

    Returns
    -------
    emit_x : float
        emittance in the x direction (m*rad)
    emit_y : float
        emittance in the y direction (m*rad)
    """
    xm = w_ave( x, w )
    xsq = w_ave( (x - xm) ** 2, w )
    ym = w_ave( y, w )
    ysq = w_ave( (y - ym) ** 2, w )
    uxm = w_ave( ux, w )
    uxsq = w_ave( (ux - uxm) ** 2, w )
    uym = w_ave( uy, w )
    uysq = w_ave( (uy - uym) ** 2, w )
    xux = w_ave( (x - xm) * (ux - uxm), w )
    yuy = w_ave( (y - ym) * (uy - uym), w )
    emit_x = ( abs(xsq * uxsq - xux ** 2) )**.5
    emit_y = ( abs(ysq * uysq - yuy ** 2) )**.5
    return emit_x, emit_y
