"""
This file is part of the OpenPMD viewer.

It defines the main OpenPMDTimeSeries class.
"""
import os
import re
import h5py
import numpy as np
from .plotter import Plotter
from .data_reader.params_reader import read_openPMD_params
from .data_reader.particle_reader import read_particle
from .data_reader.field_reader import read_field

# Check wether the interactive interface can be loaded
try:
    # If yes, use the InteractiveViewer as a parent class
    from .interactive import InteractiveViewer
    parent_class = InteractiveViewer
except ImportError:
    # Otherwise, use the default parent class
    print('[opmd_viewer] Failed to import the interactive interface.\n'
        '(Make sure that ipywidgets and IPython.display are installed.)\n'
        'The opmd_viewer API is nonetheless working.')
    parent_class = object

# Define the OpenPMDTimeSeries class and have it inherit
# from the parent class defined above

class OpenPMDTimeSeries(parent_class) :
    """
    Main class for the exploration of an openPMD timeseries

    For more details, see the docstring of the following methods:
    - get_field
    - get_particle
    - slider
    """

    def __init__( self, path_to_dir ) :
        """
        Initialize an openPMD time series
        
        More precisely, scan the directory and extract the openPMD files,
        as well as some useful openPMD parameters

        Parameter
        ---------
        path_to_dir : string
            The path to the directory where the openPMD files are.
            For the moment, only HDF5 files are supported. There should be
            one file per iteration, and the name of the files should end
            with the iteration number, followed by '.h5' (e.g. data0005000.h5)
        """
        # Extract the files and the iterations
        self.h5_files, self.iterations = list_h5_files( path_to_dir )

        # Check that there are HDF5 files in this directory
        if len(self.h5_files) == 0:
            print("Error: Found no HDF5 files in the specified directory.\n"
                "Please check that this is the path to the HDF5 files.")
            return(None)
            
        # Go through the files of the series, extract the time
        # and a few parameters.
        N_files = len(self.h5_files) 
        self.t = np.zeros( N_files )

        # - Extract parameters from the first file
        t, params0 = read_openPMD_params( self.h5_files[0] )
        self.t[0] = t
        self.avail_fields = params0['avail_fields']
        self.extension = params0['extension']
        if self.avail_fields is not None:
            self.geometry = params0['geometry']
            self.avail_circ_modes = params0['avail_circ_modes']
        self.avail_species = params0['avail_species']
        if self.avail_species is not None:
            self.avail_ptcl_quantities = params0['avail_ptcl_quantities']

        # - Check that the other files have the same parameters
        for k in range( 1, N_files ):
            t, params = read_openPMD_params( self.h5_files[k] )
            self.t[k] = t
            for key in params0.keys():
                if params != params0:
                    print(
                        "Warning: File %s has different openPMD parameters "
                        "than the rest of the time series." %self.h5_files[k])
            
        # - Set the current iteration and time
        self.current_i = 0
        self.current_t = self.t[0]
        # - Find the min and the max of the time
        self.tmin = self.t.min()
        self.tmax = self.t.max()

        # - Initialize a plotter object, which holds information about the time
        self.plotter = Plotter( self.t, self.iterations )

    def get_particle( self, t, quantity1='z', quantity2=None,
                   species='electrons', output=True, plot=False,
                   nbins=150, **kw ) :
        """
        Extract one (or two) given particle quantity
        from an HDF5 file in the OpenPMD format.

        In the case of positions, the result is returned
        in microns

        Plot the histogram of the returned quantity.
        If two quantities are requested by the user, this plots
        a 2d histogram of these quantities.

        Parameters
        ----------
        t : float (in seconds)
            Time at which to plot the file

        quantity1 : string, optional
           Which quantity to extract
           Either 'x', 'y', 'z', 'ux', 'uy', 'uz', or 'w'
           Default : 'z'

        quantity2 : string, optional
           Which second quantity to extract
           Either 'x', 'y', 'z', 'ux', 'uy', 'uz', or 'w'
           Default : no second quantity

        output : bool, optional
           Whether to return the requested quantity
           
        plot : bool, optional
           Whether to plot the requested quantity

        nbins : int, optional
           Number of bins for the histograms

        **kw : dict, otional
           Additional options to be passed to matplotlib's
           hist or hist2d.

        Returns
        -------
        A 1darray if only one quantity is requested by the user.
        A tuple of two 1darrays if two quantities are requested.
        """
        # Check that the species and quantity required are present
        if self.avail_species is None:
            print('No particle data in this time series')
            return(None)
        if (species in self.avail_species)==False:
            species_list = '\n - '.join( self.avail_species )
            print("The requested species '%s' is not available.\nThe "
                "available species are: \n - %s" %(species, species_list))
            return(None)
        if (quantity1 in self.avail_ptcl_quantities)==False or \
            (quantity2 in self.avail_ptcl_quantities + [None, 'None'])==False:
            quantity_list = '\n - '.join( self.avail_ptcl_quantities )
            print("One of the requested quantities is not available.\nThe "
                  "available quantities are: \n - %s" %quantity_list )
            return(None)

        # The requested time may not correspond exactly to an available
        # iteration. Therefore, find the last available iteration before
        # this time (modify self.current_i accordingly)
        self._find_last_output(t)
        # Get the corresponding filename
        filename = self.h5_files[ self.current_i ]
        
        # In the case of only one quantity
        if quantity2 is None or quantity2=='None' :
            # Extract from file
            q1 = read_particle( filename, species, quantity1 )
            # Plot
            if plot :
                # Extract weights for the histogram
                w = read_particle( filename, species, 'w')
                # Do the plotting
                self.plotter.hist1d( q1, w, quantity1, self.current_i,
                                     nbins, **kw )
            # Output
            if output :
                return(q1)

        # In the case of two quantities
        else :
            # Extract from file
            q1 = read_particle( filename, species, quantity1 )
            q2 = read_particle( filename, species, quantity2 )
            # Plot
            if plot :
                # Extract weights for the histogram
                w = read_particle( filename, species, 'w')
                # Do the plotting
                self.plotter.hist2d( q1, q2, w, quantity1, quantity2,
                                     self.current_i, nbins, **kw )
            # Output
            if output :
                return(q1, q2)


    def get_field(self, t, field='E', coord='z', m='all', theta=0., slicing=0.,
                  slicing_dir='y', output=True, plot=False, **kw ) :
        """
        Extract a given field from an HDF5 file in the OpenPMD format.

        Parameters
        ----------
        t : float (in seconds)
            Time at which to plot the file

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
           extent : a 1darray with 4 elements, containing the extent
        """
        # Check that the field required is present
        if self.avail_fields is None:
            print('No field data in this time series')
            return(None)
        # Check the field type
        if (field in self.avail_fields)==False:
            field_list = '\n - '.join( self.avail_fields )
            print("The requested field '%s' is not available.\nThe "
                "available fields are: \n - %s" %(field, field_list))
            return(None)
        # Check the coordinate (for vector fields)
        if self.avail_fields[field]=='vector':
            coord_available = False
            if coord in ['x', 'y', 'z']:
                coord_available = True
            elif self.geometry=='thetaMode' and (coord in ['r', 't']):
                coord_available = True
            if coord_available==False:
                print("The requested coordinate '%s' is not available." %coord)
                return(None)
        # Check the mode (for thetaMode)
        if self.geometry=="thetaMode":
            if (str(m) in self.avail_circ_modes) == False:
                mode_list = '\n - '.join(self.avail_circ_modes)
                print("The requested mode '%s' is not available.\nThe "
                    "available modes are: \n - %s" %(m, mode_list))
                return(None)
                
        # The requested time may not correspond exactly to an available
        # iteration. Therefore, find the last available iteration before
        # this time (modify self.current_i accordingly)
        self._find_last_output(t)
        # Get the corresponding filename
        filename = self.h5_files[ self.current_i ]

        # Find the proper path for vector or scalar fields
        if self.avail_fields[field] == 'scalar':
            field_path = field
            field_label = field
        elif self.avail_fields[field] == 'vector':
            field_path = os.path.join( field, coord )
            field_label = field + coord

        # Get the field data
        if self.geometry == "thetaMode":
            F, extent = read_field_circ( filename, field_path, m, theta )
        elif self.geometry == "2dcartesian":
            F, extent = read_field_2d( filename, field_path )
        elif self.geometry == "3dcartesian":
            F, extent = read_field_3d( filename, field_path,
                                       slicing, slicing_dir)

        # Plot the resulting field
        # Deactivate plotting when there is no slice selection
        if (self.geometry=="3dcartesian") and (slicing is None):
            plot = False
        if plot==True:
            self.plotter.show_field( F, extent, slicing_dir, m,
                        field_label, self.geometry, self.current_i, **kw )

        # Return the result
        return( F, extent )


    def _find_last_output(self, t) :
        """
        Find the last file that was output before t
        and store the corresponding value in self.current_t
        and self.current_i

        Parameter
        ---------
        t : float (in seconds)
            Time requested
        """
        # Make sur the time requested does not exceed
        # the allowed bounds
        if t < self.tmin :
            i = 0
            print('Reached first iteration')
        elif t > self.tmax :
            i = len(self.t) -1
            print('Reached last iteration')
        # Find the last output
        else :
            i = self.t[ self.t <= t ].argmax()

        # Register the value in the object
        self.current_i = i
        self.current_t = self.t[i]


def list_h5_files( path_to_dir ) :
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
    all_files = os.listdir( path_to_dir )

    # Select the hdf5 files
    iters_and_names = []
    for filename in all_files :
        # Use only the name that end with .h5 or .hdf5
        if filename[-3:] == '.h5' or filename[-5:] == '.hdf5':
            # Extract the iteration, using regular expressions (regex)
            regex_match = re.search('(\d+).h[df]*5', filename)
            if regex_match is None:
                print('Ill-formated HDF5 file: %s\n File names should end '
                      'with the iteration number, followed by ".h5"' %filename)
            else:
                iteration = int( regex_match.groups()[-1] )
                full_name = os.path.join(
                    os.path.abspath( path_to_dir ), filename)
                # Create list of tuples (which can be sorted together)
                iters_and_names.append( (iteration, full_name) )
                
    # Sort the list of tuples according to the iteration
    iters_and_names.sort()
    # Extract the list of filenames and iterations
    filenames = [ name for (iteration, name) in iters_and_names ]
    iterations = [ iteration for (iteration, name) in iters_and_names ]

    return( filenames, iterations )

