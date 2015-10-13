# Class that inherits from OpenPMDTimeSeries, and implements
# some standard diagnostics (emittance, etc.)
from opmd_viewer import OpenPMDTimeSeries
import numpy as np
import scipy.constants as const


class LpaDiagnostics( OpenPMDTimeSeries ):

    def __init__( self, path_to_dir ):
        """
        Initialize an OpenPMD time series

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

    def get_charge( self, t=None, iteration=None, species=None, select=None,
                    q=const.e ):
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

        q : float, optional
            Charge of `species` (default e)

        Returns
        -------
        A float with the electric charge of the selected particles
        """
        # Find the output that corresponds to the requested time/iteration
        # (Modifies self.current_i and self.current_t)
        self._find_output( t, iteration )
        # Get particle data
        w = self.get_particle( var_list=['w'], species=species, t=t,
                               iteration=iteration )
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
        - divergence in x plane
        - divergence in y plane
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
        - normalized beam emittance in the x plane
        - normalized beam get_emittance in the y plane
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
                     bins=100, q=const.e ):
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
        - The current in each bin
        - The z positions of the bin edges
        """
        # Find the output that corresponds to the requested time/iteration
        # (Modifies self.current_i and self.current_t)
        self._find_output( t, iteration )
        # Get particle data
        z, uz, uy, ux, w = self.get_particle(
                                     var_list=['z', 'uz', 'uy', 'ux', 'w'],
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
        current = vz_sum * q * bins / (len_z * 1.e-6)

        # Return the current and bin edges
        return(current, bins_list)




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
