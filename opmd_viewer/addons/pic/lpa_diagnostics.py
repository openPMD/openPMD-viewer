# Class that inherits from OpenPMDTimeSeries, and implements
# some standard diagnostics (emittance, etc.)
# For the moment the class is empty
from opmd_viewer import OpenPMDTimeSeries
import numpy as np
import scipy.constants as const

class LpaDiagnostics( OpenPMDTimeSeries ):

    def __init__(self, path_to_dir):
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

    def get_mean_gamma(self, t=None, iteration=None, species=None,
                       select=None):
        """
        Calculate the mean energy and standard deviation according to their
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
            'gamma' : [5., None]  (Particles with gamma above 5)

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
        x, y, z, ux, uy, uz , w= self.get_particle(
                                    var_list=['x','y','z','ux','uy','uz','w'], 
                                    species=species, t=t, iteration=iteration )
        # Calculate Lorentz factor for all particles
        gamma = np.sqrt(1 + ux ** 2 + uy ** 2 + uz ** 2)
        # Quantity dictionary for particle selection
        quantities = {'x': x, 'y': y, 'z': z, 'ux': ux, 'uy': uy,
                       'uz': uz, 'w': w, 'gamma': gamma} 
        # Create selection array
        select_array = apply_selection(select, quantities, np.size(w))
        # Trim gamma and w arrays to selction
        gamma = gamma[select_array]
        w = w[select_array]

        if not np.any(w):
            # If selection results in deselection of all particles, return nan
            mean_gamma = np.nan
            std_gamma = np.nan
        else:
            # Calculate weighted mean and average
            mean_gamma = np.ma.average(gamma, weights=w)
            std_gamma = wstd(gamma, w)

        # Return the result
        return [mean_gamma, std_gamma]

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
            'gamma' : [5., None]  (Particles with gamma above 5)

        q : float, optional
            Charge of `species` (default e)
        """
        # Find the output that corresponds to the requested time/iteration
        # (Modifies self.current_i and self.current_t)
        self._find_output( t, iteration )
        # Get particle data
        x, y, z, ux, uy, uz, w = self.get_particle(
                                    var_list=['x','y','z','ux','uy','uz','w'], 
                                    species=species, t=t, iteration=iteration )
        # Calculate Lorentz factor for all particles
        gamma = np.sqrt(1 + ux ** 2 + uy ** 2 + uz ** 2)
        # Quantity dictionary for particle selection
        quantities = {'x': x, 'y': y, 'z': z, 'ux': ux, 'uy': uy,
                       'uz': uz, 'w': w, 'gamma': gamma} 
        # Create selection array
        select_array = apply_selection(select, quantities, np.size(w))
        # Trim gamma and w arrays to selction
        w = w[select_array]
        # Calculate charge
        charge = np.sum(w) * q
        # Return the result
        return charge

def apply_selection(select, quantities, N):
    """
        Apply the rules of self.select to determine which
        particles should be written

        Parameters
        ----------
        select : dict
            Dictionary of rules
            to select the particles, of the form
            'x' : [-4., 10.]   (Particles having x between -4 and 10 microns)
            'z' : [0, 100] (Particles having x between 0 and 100 microns)
            'gamma' : [5., None]  (Particles with gamma above 5)

        quantities : dict
            Dictionary of particle quantities with corresponding data

        N : float
            Length of array to apply rules on

        Returns
        -------
        A 1d array of the same shape as that particle array
        containing True for the particles that satify all
        the rules of self.select
        """
    # Initialize an array filled with True
    select_array = np.ones( N, dtype='bool' )

    # Apply the rules successively
    if select is not None :
        # Go through the quantities on which a rule applies
        for quantity in select.keys() :
            quantity_array = quantities[quantity]
            # Lower bound
            if select[quantity][0] is not None :
                select_array = np.logical_and(
                    quantity_array > select[quantity][0],
                    select_array )
            # Upper bound
            if select[quantity][1] is not None :
                select_array = np.logical_and(
                    quantity_array < select[quantity][1],
                    select_array )
    return select_array

def wstd(a, weights):
    """
    Calcualted the weighted standard deviation.

    Parameters
    ----------
    a : array_like
        Calculate the weighted standard deviation for these a.

    weights : array_like
        An array of weights for the values in a.

    Returns
    -------
    Float with the weighted standard deviation
    """
    average = np.average(a, weights=weights)
    variance = np.average((a-average)**2, weights=weights)
    return np.sqrt(variance)
