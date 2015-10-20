# Class that inherits from OpenPMDTimeSeries, and implements
# some standard diagnostics (emittance, etc.)
from opmd_viewer import OpenPMDTimeSeries
import numpy as np
import scipy.constants as const


class LpaDiagnostics( OpenPMDTimeSeries ):

    def __init__( self, path_to_dir ):
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
        """
        super(LpaDiagnostics, self).__init__( path_to_dir )

    def get_rms_gamma( self, t=None, iteration=None, species=None,
                       select=None ):
        """
        Calculate the rms gamma and standard deviation according to the
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
        # Find the output that corresponds to the requested time/iteration
        # (Modifies self.current_i and self.current_t)
        self._find_output( t, iteration )
        # Get particle data
        ux, uy, uz, w = self.get_particle(
                         var_list=['ux', 'uy', 'uz', 'w'],
                         species=species, t=t, iteration=iteration )
        # Calculate Lorentz factor for all particles
        gamma = np.sqrt(1 + ux ** 2 + uy ** 2 + uz ** 2)
        # Calculate weighted mean and average
        try:
            # Calculate mean_gamma for selected particles
            mean_gamma = np.sqrt(np.average(gamma ** 2, weights=w))
        except ZeroDivisionError:
            # If selection is empty or all particles have weight zero,
            # return NaN
            mean_gamma = np.nan
        std_gamma = wstd(gamma, w)
        # Return the result
        return( mean_gamma, std_gamma )

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
        # Find the output that corresponds to the requested time/iteration
        # (Modifies self.current_i and self.current_t)
        self._find_output( t, iteration )
        # Get particle data
        w, q = self.get_particle( var_list=['w', 'charge'], species=species,
                                 t=t, iteration=iteration )
        # Calculate charge
        charge = np.sum(w) * q
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
        # Find the output that corresponds to the requested time/iteration
        # (Modifies self.current_i and self.current_t)
        self._find_output( t, iteration )
        # Get particle data
        ux, uy, uz, w = self.get_particle( var_list=['ux', 'uy', 'uz', 'w'],
                                           t=t, iteration=iteration,
                                           species=species )
        # Calculate diveregence
        div_x = wstd( ux / uz, w )
        div_y = wstd( uy / uz, w )
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
        - normalized beam get_emittance in the y plane (pi m rad)
        """
        # Find the output that corresponds to the requested time/iteration
        # (Modifies self.current_i and self.current_t)
        self._find_output( t, iteration )
        # Get particle data
        x, y, ux, uy, w = self.get_particle(
                                    var_list=['x', 'y', 'ux', 'uy', 'w'],
                                    t=t, iteration=iteration,
                                    species=species )
        # Calculate the necessary RMS values
        x *= 1.e-6
        y *= 1.e-6
        xsq = np.average( x ** 2, weights=w )
        ysq = np.average( y ** 2, weights=w )
        uxsq = np.average( ux ** 2, weights=w )
        uysq = np.average( uy ** 2, weights=w )
        xpx = np.average( x * ux, weights=w )
        ypy = np.average( y * uy, weights=w )
        # Calculate the beam emittances
        emit_x = np.sqrt( xsq * uxsq - xpx ** 2 )
        emit_y = np.sqrt( ysq * uysq - ypy ** 2 )
        # Return the results
        return( emit_x, emit_y )

    def get_current( self, t=None, iteration=None, species=None, select=None,
                     bins=100 ):
        """
        Calculate the electric current along the z-axis for selected particles.

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

        bins : float
            Number of bins along the z-axis in which to calculate the current

        Returns
        -------
        A tuple of arrays containig
        - The current in each bin in Ampere
        - The z positions of the bin edges
        """
        # Find the output that corresponds to the requested time/iteration
        # (Modifies self.current_i and self.current_t)
        self._find_output( t, iteration )
        # Get particle data
        z, uz, uy, ux, w, q = self.get_particle(
                               var_list=['z', 'uz', 'uy', 'ux', 'w', 'charge'],
                               t=t, iteration=iteration,
                               species=species )
        # Calculate Lorentz factor for all particles
        gamma = np.sqrt(1 + ux ** 2 + uy ** 2 + uz ** 2)
        # Length to be seperated in bins
        len_z = np.max(z) - np.min(z)
        # Define bins
        bins_list = np.linspace(np.min(z), np.max(z), bins)
        # Get bin index of each particle
        particle_bin = np.digitize(z, bins_list)
        # Calculate sum of vz in each bin
        vz_sum = np.zeros_like(bins_list)
        for i in range(bins):
            vz_sum[i] = np.sum(uz[particle_bin == i] / gamma[particle_bin == i]
                               * const.c * w[particle_bin == i] )
        # Calculete the current in each bin
        current = vz_sum * np.abs(q) * bins / (len_z * 1.e-6)

        # Return the current and bin edges
        return(current, bins_list)

    def get_laser_envelope( self, t=None, iteration=None, pol=None, m='all',
                            freq_filter=40, index='center' ):
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

        pol : float
            Polarization angel of the field relative to the x plane. Can be
            freely chosen between 0 and Pi/2 in thetaMode. For carthesian
            coordinates it as to be either 0 or Pi/2 (x or y plane) 

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

        Returns
        -------
        A tuple with:
        - Envelope data (1D or 2D array)
        - A FieldMetaInformation object (note that this contains obsolete
            information about the transversal coordinate in 1D case)
        """
        # Check if polarization has been entered
        if pol is None:
            raise ValueError('The `pol` argument is missing or erroneous.')
        # Find the output that corresponds to the requested time/iteration
        # (Modifies self.current_i and self.current_t)
        self._find_output( t, iteration )
        # Test that polarization is in x or y plane for 3D cart. coordinates
        if self.geometry == "3dcartesian" or self.geometry == "2dcartesian":
            if pol == 0:
                slicing_dir = 'y'
                coord = 'x'
                pol = None
            elif pol == np.pi/2:
                slicing_dir = 'x'
                coord = 'y'
                pol = None
            else:
                raise ValueError('Only polarization in the x or y plane is '
                                 'supported for carthesian coordinates')
        else:
            # Extract radial if in thetaMode
            coord = 'r'
            slicing_dir = None
        # Get field data
        field = self.get_field( t=t, iteration=iteration, field='E',
                                coord=coord, theta=pol, m=m,
                                slicing_dir=slicing_dir )
        extent = field[1]
        if index == 'center':
            # Get central slice
            field_slice = field[0][int( field[0].shape[0] / 2), :]
            # Calculate inverse FFT of filtered FFT array
            envelope = self._fft_filter(field_slice, freq_filter)
        elif index == 'all':
            envelope = [ self._fft_filter(field[0][i, :], freq_filter)
                         for i in range(field[0].shape[0]) ]
        else:
            field_slice = field[0][index, :]
            # Calculate inverse FFT of filtered FFT array
            envelope = self._fft_filter(field_slice, freq_filter)
        # Return the result
        return( envelope, extent )

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
        central_freq_i = np.argmax(fft_field_slice[:N / 2])
        central_freq = fft_freqs[central_freq_i]
        # Filter frequencies higher than central_freq * freq_filter/100
        filter_bound = central_freq * freq_filter / 100.
        # Find index from where to filter
        filter_i = np.argmin(np.abs(filter_bound - fft_freqs))
        filter_freq_range_i = central_freq_i - filter_i
        # Write filtered FFT array
        filtered_fft = np.zeros_like( field, dtype=np.complex )
        filtered_fft[N / 2 - filter_freq_range_i:
                     N / 2 + filter_freq_range_i] \
        = fft_field_slice[central_freq_i - filter_freq_range_i:
                          central_freq_i + filter_freq_range_i]
        # Calculate inverse FFT of filtered FFT array
        envelope = np.abs(np.fft.ifft(np.fft.fftshift(2 * filtered_fft)))

        # Return the result
        return( envelope )

    def get_mean_frequency( self, t=None, iteration=None, pol=None, m='all'):
        """
        Calculate the rms angular frequency of a laser pulse.

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        pol : float
            Polarization angel of the field relative to the x plane. Can be
            freely chosen between 0 and Pi/2 in thetaMode. For carthesian
            coordinates it as to be either 0 or Pi/2 (x or y plane) 

        m : int or str, optional
           Only used for thetaMode geometry
           Either 'all' (for the sum of all the modes)
           or an integer (for the selection of a particular mode)

        Returns
        -------
        A float with rms angular frequency
        """
        # Check if polarization has been entered
        if pol is None:
            raise ValueError('The `pol` argument is missing or erroneous.')
        # Find the output that corresponds to the requested time/iteration
        # (Modifies self.current_i and self.current_t)
        self._find_output( t, iteration )
        # Test that polarization is in x or y plane for 3D cart. coordinates
        if self.geometry == "3dcartesian" or self.geometry == "2dcartesian":
            if pol == 0:
                slicing_dir = 'y'
                coord = 'x'
                pol = None
            elif pol == np.pi/2:
                slicing_dir = 'x'
                coord = 'y'
                pol = None
            else:
                raise ValueError('Only polarization in the x or y plane is '
                                 'supported for carthesian coordinates')
        else:
            # Extract radial if in thetaMode
            coord = 'r'
            slicing_dir = None
        # Get field data
        field = self.get_field( t=t, iteration=iteration, field='E',
                                coord=coord, theta=pol, m=m,
                                slicing_dir=slicing_dir )
        z_length = field[1].zmax - field[1].zmin
        # Get central field lineout
        field1d = field[0][field[0].shape[0]/2, :]
        # FFT of 1d data
        fft_field = np.fft.fft(field1d)
        # Corresponding angular frequency
        frq = np.fft.fftfreq(field1d.size, z_length /
                             field1d.size * 1 / const.c) * 2 * np.pi
        # Calculate the RMS of the frequencies
        rms = np.sqrt(np.average(frq[:frq.size/2]**2, weights=np.abs(
                      fft_field[:frq.size/2])))
        return( rms )

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

        pol : float
            Polarization angel of the field relative to the x plane. Can be
            freely chosen between 0 and Pi/2 in thetaMode. For carthesian
            coordinates it as to be either 0 or Pi/2 (x or y plane)

        Returns
        -------
        Float with normalized vector potential a0
        """
        # Get the peak field from field envelope
        Emax = np.amax(self.get_laser_envelope(t=t, iteration=iteration,
                                                 pol=pol)[0])
        # Get mean frequency
        omega = self.get_mean_frequency(t=t, iteration=iteration, pol=pol)
        # Calculate a0
        a0 = Emax * const.e / (const.m_e * const.c * omega)
        return( a0 )

    def get_ctau( self, t=None, iteration=None, pol=None ):
        """
        Calculate the length of a (gaussian) laser pulse. ( sqrt(2) * sigma_z)

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        pol : float
            Polarization angel of the field relative to the x plane. Can be
            freely chosen between 0 and Pi/2 in thetaMode. For carthesian
            coordinates it as to be either 0 or Pi/2 (x or y plane)

        Returns
        -------
        Float with ctau in meters
        """
        # Get the field envelope
        E, extent = self.get_laser_envelope(t=t, iteration=iteration,
                                            pol=pol)
        # Calculate standard deviation
        sigma = wstd(extent.z, E)
        # Return ctau = sqrt(2) * sigma
        return( np.sqrt(2) * sigma )

    def get_laser_waist( self, t=None, iteration=None, pol=None ):
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

        pol : float
            Polarization angel of the field relative to the x plane. Can be
            freely chosen between 0 and Pi/2 in thetaMode. For carthesian
            coordinates it as to be either 0 or Pi/2 (x or y plane)

        Returns
        -------
        Float with laser waist in meters
        """
        # Get the field envelope
        field, extent = self.get_laser_envelope(t=t, iteration=iteration,
                                                pol=pol, index='all')

        # Find the maximum of the envelope along the transverse axis
        trans_max = np.amax(field, axis=1)
        # Get transverse positons
        trans_pos = getattr(extent, extent.axes[0])
        # Calculate standard deviation
        sigma_r = wstd(trans_pos, trans_max)
        # Return the laser waist = sqrt(2) * sigma_r
        return(np.sqrt(2) * sigma_r)

    def wigner_transform( self, t=None, iteration=None, pol=None ):
        """
        Calculates the wigner transformation of a laserpulse.
        See Trebino, R: Frequency Resolved Optical Gating: The measurements of
        Ultrashort Laser Pulses: year 2000: formula 5.2
        Resulting wigner transform is not taken to power of two.

        Parameters
        ----------
        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        pol : float
            Polarization angel of the field relative to the x plane. Can be
            freely chosen between 0 and Pi/2 in thetaMode. For carthesian
            coordinates it as to be either 0 or Pi/2 (x or y plane)

        Returns
        -------
        - A 2d array with wigner transform
        - extents (time centered around the laser pulse in s and the
          angular frequency)
        """
        # Test that polarization is in x or y plane for 3D cart. coordinates
        if self.geometry == "3dcartesian" or self.geometry == "2dcartesian":
            if pol == 0:
                slicing_dir = 'y'
                coord = 'x'
            elif pol == np.pi/2:
                slicing_dir = 'x'
                coord = 'y'
            else:
                raise ValueError('Only polarization in the x or y plane is '
                                 'supported for carthesian coordinates')
        else:
            # Extract radial if in thetaMode
            coord = 'r'
            slicing_dir = None
        # Get the field envelope
        env, _ = self.get_laser_envelope(t=t, iteration=iteration,
                                         pol=pol)
        # Get the field
        E, extent = self.get_field( t=t, iteration=iteration, field='E',
                                    coord=coord, theta=pol,
                                    slicing_dir=slicing_dir )
        # Get central slice
        E = E[E.shape[0] / 2, :]
        # Get time domain of the data
        tmin = extent.zmin / const.c
        tmax = extent.zmax / const.c
        T = tmax - tmin
        dt = T / E.size
        # Normalize the Envelope
        env /= np.sqrt(np.trapz(env ** 2, dx=dt))
        # Allocate array for shifted E field and wigner transform
        E_shift = np.zeros_like(E)
        wigner = np.zeros((2 * E.size, E.size))
        # Construct shifted E field
        for i in range(E.size * 2):
            itau = i % E.size
            if i < E.size:
                E_shift[:itau] = env[E.size - itau: E.size]
                E_shift[itau:] = 0
            else:
                E_shift[itau:] = env[: E.size - itau]
                E_shift[:itau] = 0
            EE = E * E_shift ** 2
            fftwigner = np.fft.fft(EE)
            wigner[i, :] = np.abs(fftwigner) ** 2
        # Rotate and flip array to have input form of imshow
        wigner = np.flipud(np.rot90(wigner[:, E.size / 2:]))
        # Calculate the axis range
        maxi, maxj = np.unravel_index(wigner.argmax(), wigner.shape)
        xmin = -(T - T / wigner.shape[1] * maxj)
        xmax = -(0 - T / wigner.shape[1] * maxj)
        # Return the wigner transform an the extents
        return( wigner, [xmin, xmax, 0, np.pi / dt] )

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
        variance = np.average((a-average)**2, weights=weights)
        return( np.sqrt(variance) )
