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
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
from scipy.optimize import curve_fit


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
        std_gamma = wstd(gamma, w)
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
            spreads[i] = wstd(gamma[z_filter], w[z_filter])
            zi += dz
            i += 1
        # Plot the result if needed
        if plot:
            iteration = self.iterations[ self.current_i ]
            time_fs = 1.e15 * self.t[ self.current_i ]
            plt.plot( z_pos, spreads, **kw)
            plt.title("Slice energy spread at %.1f fs   (iteration %d)"
                % (time_fs, iteration ), fontsize=self.plotter.fontsize)
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
        # Calculate diveregence
        div_x = wstd( np.arctan2(ux, uz), w )
        div_y = wstd( np.arctan2(uy, uz), w )
        # Return the result
        return( div_x, div_y )

    def get_emittance( self, t=None, iteration=None, species=None,
                       select=None ):
        """
        Calculate the normalized RMS emittance.
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
            Either None or a dictionary of rules
            to select the particles, of the form
            'x' : [-4., 10.]   (Particles having x between -4 and 10 microns)
            'z' : [0, 100] (Particles having x between 0 and 100 microns)

        Returns
        -------
        A tuple with :
        - normalized beam emittance in the x plane (pi m rad)
        - normalized beam emittance in the y plane (pi m rad)
        """
        # Get particle data
        x, y, ux, uy, w = self.get_particle(
            var_list=['x', 'y', 'ux', 'uy', 'w'],
            t=t, iteration=iteration,
            species=species, select=select )
        # Calculate the necessary RMS values
        x *= 1.e-6
        y *= 1.e-6
        xsq = np.average( x ** 2, weights=w )
        ysq = np.average( y ** 2, weights=w )
        uxsq = np.average( ux ** 2, weights=w )
        uysq = np.average( uy ** 2, weights=w )
        xux = np.average( x * ux, weights=w )
        yuy = np.average( y * uy, weights=w )
        # Calculate the beam emittances
        emit_x = np.sqrt( xsq * uxsq - xux ** 2 )
        emit_y = np.sqrt( ysq * uysq - yuy ** 2 )
        # Return the results
        return( emit_x, emit_y )

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
        # Calculate Lorentz factor for all particles
        gamma = np.sqrt(1 + ux ** 2 + uy ** 2 + uz ** 2)
        # Calculate particle velocities
        vz = uz / gamma * const.c
        # Length to be seperated in bins
        len_z = np.max(z) - np.min(z)
        vzq_sum, _ = np.histogram(z, bins=bins, weights=(vz * w * q))
        # Calculete the current in each bin
        current = np.abs(vzq_sum * bins / (len_z * 1.e-6))
        # Info object with central position of the bins
        info = FieldMetaInformation( {0: 'z'}, current.shape,
            grid_spacing=(len_z / bins, ), grid_unitSI=1,
            global_offset=(np.min(z) + len_z / bins / 2,), position=(0,))
        # Plot the result if needed
        if plot:
            iteration = self.iterations[ self.current_i ]
            time_fs = 1.e15 * self.t[ self.current_i ]
            plt.plot( info.z, current, **kw)
            plt.title("Current at %.1f fs   (iteration %d)"
                % (time_fs, iteration ), fontsize=self.plotter.fontsize)
            plt.xlabel('$z \;(\mu m)$', fontsize=self.plotter.fontsize)
            plt.ylabel('$I \;(A)$', fontsize=self.plotter.fontsize)
        # Return the current and bin centers
        return(current, info)

    def get_laser_envelope( self, t=None, iteration=None, pol=None, m='all',
                            freq_filter=40, index='center', theta=0,
                            slicing_dir='y', plot=False, **kw ):
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

        freq_filter : float, optional
            Range of frequencies in percent which to filter: Frequencies higher
            than freq_filter/100 times the dominant frequencies will be
            filtered out

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
        if index == 'center':
            # Get central slice
            field_slice = field[0][int( field[0].shape[0] / 2), :]
            # Calculate inverse FFT of filtered FFT array
            envelope = self._fft_filter(field_slice, freq_filter)
        elif index == 'all':
            envelope = np.array([ self._fft_filter(field[0][i, :], freq_filter)
                                  for i in range(field[0].shape[0]) ])
        else:
            field_slice = field[0][index, :]
            # Calculate inverse FFT of filtered FFT array
            envelope = self._fft_filter(field_slice, freq_filter)

        # Restrict the metainformation to 1d if needed
        if index != 'all':
            info.restrict_to_1Daxis( info.axes[1] )

        # Plot the result if needed
        if plot:
            iteration = self.iterations[ self.current_i ]
            time_fs = 1.e15 * self.t[ self.current_i ]
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

    def _fft_filter(self, field, freq_filter):
        """
        Filters out high frequencies in input data. Frequencies higher than
        freq_filter / 100 times the dominant frequency will be filtered.

        Parameters
        ----------
        field : 1D array
            Array with input data in time/space domain

        freq_filter : float
            Frequency range in percent around the dominant frequency which will
            not be filtered out

        Returns
        -------
        A 1D array with filtered input data in time/space domain
        """
        # Number of sample points
        N = field.size
        # Fourier transform of the field slice
        fft_field_slice = np.fft.fft(field)
        fft_freqs = np.fft.fftfreq(N)
        # Find central frequency
        central_freq_i = np.argmax(np.abs(fft_field_slice[:int(N / 2)]))
        central_freq = fft_freqs[central_freq_i]
        # Filter frequencies higher than central_freq * freq_filter/100
        filter_bound = central_freq * freq_filter / 100.
        # Find index from where to filter
        filter_i = np.argmin(np.abs(filter_bound - fft_freqs))
        filter_freq_range_i = central_freq_i - filter_i
        # Write filtered FFT array
        filtered_fft = np.zeros_like( field, dtype=np.complex )
        filtered_fft[int(N / 2) - filter_freq_range_i:
                     int(N / 2) + filter_freq_range_i] \
        = fft_field_slice[central_freq_i - filter_freq_range_i:
                          central_freq_i + filter_freq_range_i]
        # Calculate inverse FFT of filtered FFT array
        envelope = np.abs(np.fft.ifft(np.fft.fftshift(2 * filtered_fft)))

        # Return the result
        return( envelope )

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
            iteration = self.iterations[ self.current_i ]
            time_fs = 1.e15 * self.t[ self.current_i ]
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
        ctau = np.sqrt(2) * wstd(info.z, E)
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
        w0 = np.sqrt(2) * wstd(trans_pos, trans_slice)
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
            iteration = self.iterations[ self.current_i ]
            time_fs = 1.e15 * self.t[ self.current_i ]
            plt.imshow( spectrogram, extent=info.imshow_extent, aspect='auto',
                        **kw)
            plt.title("Spectrogram at %.1f fs   (iteration %d)"
                % (time_fs, iteration ), fontsize=self.plotter.fontsize)
            plt.xlabel('$t \;(s)$', fontsize=self.plotter.fontsize )
            plt.ylabel('$\omega \;(rad.s^{-1})$',
                       fontsize=self.plotter.fontsize )

        return( spectrogram, info )


def wstd( a, weights ):
    """
    Calcualte the weighted standard deviation.

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
