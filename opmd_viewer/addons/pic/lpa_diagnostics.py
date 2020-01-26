"""
This file is part of the openPMD-viewer.

This file contains diagnostics relevant for laser-plasma acceleration

Copyright 2015-2016, openPMD-viewer contributors
Authors: Soeren Jalas, Remi Lehe
License: 3-Clause-BSD-LBNL
"""

# Class that inherits from OpenPMDTimeSeries, and implements
# some standard diagnostics (emittance, etc.)
from opmd_viewer import OpenPMDTimeSeries, FieldMetaInformation
import numpy as np
import scipy.constants as const
from scipy.optimize import curve_fit
from opmd_viewer.openpmd_timeseries.plotter import check_matplotlib
from scipy.signal import hilbert
try:
    import matplotlib.pyplot as plt
except ImportError:
    # Any error will be caught later by `check_matplotlib`
    pass


class LpaDiagnostics( OpenPMDTimeSeries ):

    def __init__( self, path_to_dir, check_all_files=True ):
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
        """
        OpenPMDTimeSeries.__init__( self, path_to_dir,
                                    check_all_files=check_all_files )

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
            'x' : [-4., 10.]   (Particles having x between -4 and 10 microns)
            'z' : [0, 100] (Particles having x between 0 and 100 microns)

        Returns
        -------
        A tuple of floats with:
        - mean weighted gamma
        - weighted standard deviation of gamma
        """
        # Get particle data
        ux, uy, uz, w = self.get_particle(
            var_list=['ux', 'uy', 'uz', 'w'], select=select,
            species=species, t=t, iteration=iteration )
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
        return( mean_gamma, std_gamma )

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
            'x' : [-4., 10.]   (Particles having x between -4 and 10 microns)
            'z' : [0, 100] (Particles having x between 0 and 100 microns)

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
            time_fs = 1.e15 * self.t[ self._current_i ]
            plt.plot(z_pos, spreads, **kw)
            plt.title("Slice energy spread at %.1f fs   (iteration %d)"
                % (time_fs, iteration), fontsize=self.plotter.fontsize)
            plt.xlabel('$z \;(\mu m)$', fontsize=self.plotter.fontsize)
            plt.ylabel('$\sigma_\gamma (\Delta_z=%s\mu m)$' % dz,
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
            'x' : [-4., 10.]   (Particles having x between -4 and 10 microns)
            'z' : [0, 100] (Particles having x between 0 and 100 microns)

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
            'x' : [-4., 10.]   (Particles having x between -4 and 10 microns)
            'z' : [0, 100] (Particles having x between 0 and 100 microns)

        Returns
        -------
        A tuple with:
        - divergence in x plane in rad
        - divergence in y plane in rad
        """
        # Get particle data
        ux, uy, uz, w = self.get_particle( var_list=['ux', 'uy', 'uz', 'w'],
                                           t=t, iteration=iteration,
                                           species=species, select=select )
        # Calculate divergence
        div_x = w_std( np.arctan2(ux, uz), w )
        div_y = w_std( np.arctan2(uy, uz), w )
        # Return the result
        return( div_x, div_y )

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
            'x' : [-4., 10.]   (Particles having x between -4 and 10 microns);
            'z' : [0, 100] (Particles having x between 0 and 100 microns).

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
        x *= 1.e-6
        y *= 1.e-6
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
            'x' : [-4., 10.]   (Particles having x between -4 and 10 microns)
            'z' : [0, 100] (Particles having x between 0 and 100 microns)

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
            current = np.abs(vzq_sum * bins / (len_z * 1.e-6))
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
            time_fs = 1.e15 * self.t[ self._current_i ]
            plt.plot( info.z, current, **kw)
            plt.title("Current at %.1f fs   (iteration %d)"
                % (time_fs, iteration ), fontsize=self.plotter.fontsize)
            plt.xlabel('$z \;(\mu m)$', fontsize=self.plotter.fontsize)
            plt.ylabel('$I \;(A)$', fontsize=self.plotter.fontsize)
        # Return the current and bin centers
        return(current, info)

    def get_laser_envelope( self, t=None, iteration=None, pol=None, m='all',
                            index='center', theta=0, slicing_dir='y', plot=False, **kw ):
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
            Polarization of the field. Options are 'x', 'y'

        m : int or str, optional
           Only used for thetaMode geometry
           Either 'all' (for the sum of all the modes)
           or an integer (for the selection of a particular mode)

        index : int or str, optional
            Transversal index of the slice from which to calculate the envelope
            Default is 'center', using the center slice.
            Use 'all' to calculate a full 2D envelope

        theta : float, optional
           Only used for thetaMode geometry
           The angle of the plane of observation, with respect to the x axis

        slicing_dir : str, optional
           Only used for 3dcartesian geometry
           The direction along which to slice the data
           Either 'x', 'y'

        plot : bool, optional
           Whether to plot the requested quantity

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
        if pol is None:
            raise ValueError('The `pol` argument is missing or erroneous.')
        # Get field data
        field = self.get_field( t=t, iteration=iteration, field='E',
                                coord=pol, theta=theta, m=m,
                                slicing_dir=slicing_dir )
        info = field[1]
        if index == 'all':
            # Filter the full 2D array
            e_complx = hilbert(field[0], axis=1)
        elif index == 'center':
            # Filter the central slice (1D array)
            field_slice = field[0][int( field[0].shape[0] / 2), :]
            e_complx = hilbert(field_slice)
        else:
            # Filter the requested slice (2D array)
            field_slice = field[0][index, :]
            e_complx = hilbert(field_slice)
        envelope = np.abs(e_complx)
        # Restrict the metainformation to 1d if needed
        if index != 'all':
            info.restrict_to_1Daxis( info.axes[1] )

        # Plot the result if needed
        if plot:
            check_matplotlib()
            iteration = self.iterations[ self._current_i ]
            time_fs = 1.e15 * self.t[ self._current_i ]
            if index != 'all':
                plt.plot( 1.e6 * info.z, envelope, **kw)
                plt.ylabel('$E_%s \;(V/m)$' % pol,
                           fontsize=self.plotter.fontsize)
            else:
                plt.imshow( envelope, extent=1.e6 * info.imshow_extent,
                            aspect='auto', **kw)
                plt.colorbar()
                plt.ylabel('$%s \;(\mu m)$' % pol,
                            fontsize=self.plotter.fontsize)
            plt.title("Laser envelope at %.1f fs   (iteration %d)"
                % (time_fs, iteration ), fontsize=self.plotter.fontsize)
            plt.xlabel('$z \;(\mu m)$', fontsize=self.plotter.fontsize)
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
        if pol == 'x':
            slicing_dir = 'y'
            theta = 0
        else:
            slicing_dir = 'x'
            theta = np.pi / 2.

        # Get field data
        field, info = self.get_field( t=t, iteration=iteration, field='E',
                                coord=pol, theta=theta, m=m,
                                slicing_dir=slicing_dir )
        # Get central field lineout
        field1d = field[ int( field.shape[0] / 2 ), :]
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
            time_fs = 1.e15 * self.t[ self._current_i ]
            plt.plot( spect_info.omega, spectrum, **kw )
            plt.xlabel('$\omega \; (rad.s^{-1})$',
                       fontsize=self.plotter.fontsize )
            plt.ylabel('Spectrum', fontsize=self.plotter.fontsize )
            plt.title("Spectrum at %.1f fs   (iteration %d)"
                % (time_fs, iteration ), fontsize=self.plotter.fontsize)
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
        if pol not in ['x', 'y']:
            raise ValueError('The `pol` argument is missing or erroneous.')

        if pol == 'x':
            slicing_dir = 'y'
            theta = 0
        else:
            slicing_dir = 'x'
            theta = np.pi / 2.
        # Get the peak field from field envelope
        Emax = np.amax(self.get_laser_envelope(t=t, iteration=iteration,
                                               pol=pol, theta=theta,
                                               slicing_dir=slicing_dir)[0])
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
        if pol not in ['x', 'y']:
            raise ValueError('The `pol` argument is missing or erroneous.')
        if pol == 'x':
            slicing_dir = 'y'
            theta = 0
        else:
            slicing_dir = 'x'
            theta = np.pi / 2.
        # Get the field envelope
        E, info = self.get_laser_envelope(t=t, iteration=iteration,
                                            pol=pol, theta=theta,
                                            slicing_dir=slicing_dir)
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
                         slicing_dir='y', method='fit' ):
        """
        Calculate the waist of a (gaussian) laser pulse. ( sqrt(2) * sigma_r)

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

        slicing_dir : str, optional
           Only used for 3dcartesian geometry
           The direction along which to slice the data
           Either 'x', 'y'

        method : str, optional
           The method which is used to compute the waist
           'fit': Gaussian fit of the transverse profile
           'rms': RMS radius, weighted by the transverse profile
           ('rms' tends to give more weight to the "wings" of the pulse)

        Returns
        -------
        Float with laser waist in meters
        """
        # Get the field envelope
        field, info = self.get_laser_envelope(t=t, iteration=iteration,
                                                pol=pol, index='all',
                                                slicing_dir=slicing_dir,
                                                theta=theta)
        # Find the indices of the maximum field, and
        # pick the corresponding transverse slice
        itrans_max, iz_max = np.unravel_index(
            np.argmax( field ), dims=field.shape )
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

    def get_spectrogram( self, t=None, iteration=None, pol=None, theta=0,
                          slicing_dir='y', plot=False, **kw ):
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
        # Get the field envelope
        env, _ = self.get_laser_envelope(t=t, iteration=iteration, pol=pol)
        # Get the field
        E, info = self.get_field( t=t, iteration=iteration, field='E',
                                    coord=pol, theta=theta,
                                    slicing_dir=slicing_dir )
        # Get central slice
        E = E[ int(E.shape[0] / 2), :]
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
            time_fs = 1.e15 * self.t[ self._current_i ]
            plt.imshow( spectrogram, extent=info.imshow_extent, aspect='auto',
                        **kw)
            plt.title("Spectrogram at %.1f fs   (iteration %d)"
                % (time_fs, iteration ), fontsize=self.plotter.fontsize)
            plt.xlabel('$t \;(s)$', fontsize=self.plotter.fontsize )
            plt.ylabel('$\omega \;(rad.s^{-1})$',
                       fontsize=self.plotter.fontsize )
        return( spectrogram, info )

    def get_z_phase( self, t=None, iteration=None,
                      m='all', plot2D=False, plot1D=False, start_ratio=1.e-2, wavelength_guess = 2.e-5, **kw ):
        """
        Return the z position of phase points of E_z.

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        m : int or str, optional
           Only used for thetaMode geometry
           Either 'all' (for the sum of all the modes)
           or an integer (for the selection of a particular mode)

        plot2D: bool, optional
           Whether to plot the 2D field in pseudocolor map

        plot1D: bool, optional
           Whether to plot the axial field lineout

        start_ratio: float, optional
           The start Ez ratio respect to the maximum value

        wavelength_guess: float, optional
           The guess of plasma wavelength, in unit of metre

        **kw : dict, otional
           Additional options to be passed to matplotlib's `plot` method

        Returns
        -------
        A array with z position of 5 points:
            - The point with the first E_z exceeding a threshold
            - The first maximum of E_z
            - The first zero-crossing of E_z
            - The first minimum of E_z
            - The second zero-crossing of E_z
            - The second maximum of E_z
        """
        # Get field data
        field, info = self.get_field( t=t, iteration=iteration, field='E',
                                coord='z', theta=0, m=m,
                                slicing_dir='y', plot=plot2D )
        # Get central field lineout
        field1d = field[ int( field.shape[0] / 2 ), :]
        xaxis = getattr( info, 'z' )

        # Initialization
        z_array = np.full(6, np.nan)
        # Looking for start point
        try: z_array[0], ind_start = start_point(field1d, xaxis, np.max(np.abs(field1d))*start_ratio)
        except RuntimeError: return z_array

        wavelength_range = int(wavelength_guess/(xaxis[1]-xaxis[0]))
        # Looking for 1st maximum ind_max
        try: ind_max = next_local_max(field1d, start_index=ind_start, local_range = wavelength_range//2)
        except RuntimeError: return z_array
        z_array[1] = xaxis[ind_max]

        # Look for ind_start again
        z_array[0], ind_start = start_point(field1d, xaxis, field1d[ind_max]*start_ratio)

        # Looking for 1st minimum ind_min
        try: ind_min = next_local_min(field1d, start_index=ind_max, local_range = wavelength_range)
        except RuntimeError: return z_array
        z_array[3] = xaxis[ind_min]

        # Looking for 1st zero point
        z_array[2] = zero_point(field1d, ind_min, ind_max, xaxis)

        # Looking for 2nd maximum ind_max2
        try: ind_max2 = next_local_max(field1d, start_index=ind_min, local_range = int(wavelength_range*0.6))
        except RuntimeError: return z_array
        z_array[5] = xaxis[ind_max2]

        # Looking for 2nd zero point
        z_array[4] = zero_point(field1d, ind_max2, ind_min, xaxis)
        # Plot the field if required
        if plot1D:
            check_matplotlib()
            iteration = self.iterations[ self._current_i ]
            time_fs = 1.e15 * self.t[ self._current_i ]
            plt.figure()
            plt.plot( 1.e6*xaxis, field1d, **kw )
            plt.xlabel('$z \; [\mu m]$',
                       fontsize=self.plotter.fontsize )
            plt.ylabel('$E_z$', fontsize=self.plotter.fontsize )
            plt.title("$E_z$ lineout at %.1f fs   (iteration %d)"
                % (time_fs, iteration ), fontsize=self.plotter.fontsize)
            plt.plot(1.e6*z_array, np.array([0., field1d[ind_max], 0., field1d[ind_min], 0., field1d[ind_max2]]),'ro',fillstyle='none')
        plt.show()
        return z_array


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

def start_point(F, x, threshold = 1.e-2):
    '''Find the start point based on the index when F abslute value > threshold.

    Parameters
    ----------
    F : 1D array of floats
        Field array. If F is not 1D, raise an exception.
    x : 1D array of float
        The axis of points
    threshold : float

    Return
    ----------
    pos_start : float
                The position of start point.
    ind_start: int
               The right most index when F abslute value > threshold.'''

    if F.ndim!=1: raise RuntimeError('F shold be 1 dimensional!')
    for i in range(len(F)-1, 0, -1):
        if np.absolute(F[i])>threshold:
            #return  x[i]-F[i]*(x[i]-x[i-1])/(F[i]-F[i-1]), i
            return  x[i], i
    raise RuntimeError('Cannot find start point! You can try to reduce threshold.')

def next_local_max(F, start_index=None, local_range = None):
    '''
    Find the next local maximum point for a 1D data F from the start_index to start_index-local_range.

    Parameters
    ----------
    F : 1D array of floats
        Field array. If F is not 1D, raise an exception.
    start_index : int
                  The upper index in F to search for. If start_index is None, start_index = len(F).
    local_range : int
                  The length to search for. This function will search the range from start_index-local_range to start_index. In general, local_range should be approximately plasma wavelength/dz. If local_range is None, local_range = len(F)//2

    Return
    ----------
    A integer, the index of local maximum
    '''
    if F.ndim!=1: raise RuntimeError('F shold be 1 dimensional!')
    if start_index is None: start_index = len(F)
    if local_range is None: local_range = len(F)//2
    till_index = start_index - local_range
    if till_index<0: till_index = 0
    return_ind=np.argmax(F[till_index:start_index])+till_index
    if return_ind<1: raise RuntimeError('Local maximum not found!')
    return return_ind

def next_local_min(F, start_index=None, local_range = None):
    '''
    Find the next local minimum point for a 1D data F from the start_index to start_index-local_range.

    Parameters
    ----------
    F : 1D array of floats
        Field array. If F is not 1D, raise an exception.
    start_index : int
                  The upper index in F to search for. If start_index is None, start_index = len(F).
    local_range : int
                  The length to search for. This function will search the range from start_index-local_range to start_index. In general, local_range should be approximately plasma wavelength/dz. If local_range is None, local_range = len(F)//2

    Return
    ----------
    A integer, the index of local minimum
    '''
    return next_local_max(-F, start_index=start_index, local_range = local_range)

def zero_point(F, start, stop, xaxis):
    '''
    Find the zero point for a 1D data F from index of start to stop.

    Parameters
    ----------
    F : 1D array of floats
        Field array. If F is not 1D, raise an exception.
    start : int
    stop : int
    xaxis : 1D array of float
            The axis of points

    Return
    ----------
    A float, the location of zero point
    '''
    if F.ndim!=1: raise RuntimeError('F shold be 1 dimensional!')
    if np.sign(F[start])*np.sign(F[stop])>0:
        #raise RuntimeError('F[start] and F[stop] have the same sign! There may be no zero point between them.')
        return np.nan
    ind1 = start
    ind2 = stop
    while ind2-ind1 > 1:
        ind_mid = (ind1 + ind2) // 2
        if np.sign(F[ind1])*np.sign(F[ind_mid]) > 0: ind1 = ind_mid
        else: ind2 = ind_mid
    # Do linear interpolation between ind1 and ind2 for the zero point location
    x1 = xaxis[ind1]
    x2 = xaxis[ind2]
    y1 = F[ind1]
    y2 = F[ind2]
    return (x1*y2-x2*y1)/(y2-y1)
