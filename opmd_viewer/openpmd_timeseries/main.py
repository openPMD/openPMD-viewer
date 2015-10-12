"""
This file is part of the OpenPMD viewer.

It defines the main OpenPMDTimeSeries class.
"""
import os
import re
import numpy as np
from .plotter import Plotter
from .data_reader.params_reader import read_openPMD_params
from .data_reader.particle_reader import read_particle
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

    def get_particle( self, var_list=None, species=None, t=None,
                iteration=None, output=True, plot=False, nbins=150, **kw ) :
        """
        Extract a list of particle variables
        from an HDF5 file in the OpenPMD format.

        In the case of positions, the result is returned in microns

        Plot the histogram of the returned quantity.
        If two quantities are requested by the user, this plots
        a 2d histogram of these quantities.

        Parameters
        ----------
        var_list : list of string, optional
           A list of the particle variables to extract. If var_list is not
           provided, the available particle quantities are printed

        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.
           
        output : bool, optional
           Whether to return the requested quantity
           
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
            print('No particle data in this time series')
            return(None)
        if (species in self.avail_species)==False:
            species_list = '\n - '.join( self.avail_species )
            print("The argument `species` is missing or erroneous.\nThe "
                "available species are: \n - %s" %species_list )
            print("Please set the argument `species` accordingly.") 
            return(None)

        # Check the list of variables
        valid_var_list = True
        if type(var_list) != list:
            valid_var_list = False
        else:
            for quantity in var_list:
                if (quantity in self.avail_ptcl_quantities) == False:
                    valid_var_list = False
        if valid_var_list == False:
            quantity_list = '\n - '.join( self.avail_ptcl_quantities )
            print("The argument `var_list` is missing or erroneous.\nIt "
                  "should be a list of strings representing particle "
                  "quantities.\n The available quantities are: "
                  "\n - %s" %quantity_list )
            print("Please set the argument `var_list` accordingly.")
            return(None)

        # Find the output that corresponds to the requested time/iteration
        # (Modifies self.current_i and self.current_t)
        self._find_output( t, iteration )
        # Get the corresponding filename
        filename = self.h5_files[ self.current_i ]
        
        # Extract the list of particle quantities
        data_list = []
        for quantity in var_list:
            data_list.append( read_particle( filename, species, quantity ) )

        # Plotting
        # - In the case of only one quantity
        if len(data_list) == 1:
            if plot :
                # Extract weights for the histogram
                w = read_particle( filename, species, 'w')
                # Do the plotting
                self.plotter.hist1d( data_list[0], w, var_list[0],
                                     self.current_i, nbins, **kw )

        # - In the case of two quantities
        elif len(data_list) == 2:
            # Plot
            if plot :
                # Extract weights for the histogram
                w = read_particle( filename, species, 'w')
                # Do the plotting
                self.plotter.hist2d( data_list[0], data_list[1], w,
                                     var_list[0], var_list[1],
                                     self.current_i, nbins, **kw )
        # Output
        if output :
            return( data_list )


    def get_field(self, field=None, coord=None, t=None, iteration=None,
                  m='all', theta=0., slicing=0., slicing_dir='y', 
                  output=True, plot=False, **kw ) :
        """
        Extract a given field from an HDF5 file in the OpenPMD format.

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
           extent : a 1darray with 4 elements, containing the extent
        """
        # Check that the field required is present
        if self.avail_fields is None:
            print('No field data in this time series')
            return(None,None)
        # Check the field type
        if (field in self.avail_fields)==False:
            field_list = '\n - '.join( self.avail_fields )
            print("The `field` argument is missing or erroneous.\nThe "
                "available fields are: \n - %s" %(field_list))
            print("Please set the `field` argument accordingly.")
            return(None,None)
        # Check the coordinate (for vector fields)
        if self.avail_fields[field]=='vector':
            available_coord = ['x', 'y', 'z']
            if self.geometry=='thetaMode':
                available_coord += ['r', 't']
            if (coord in available_coord) == False:
                coord_list = '\n - '.join( available_coord )
                print("The field %s is a vector field, but the `coord` "
                      "argument is missing or erroneous.\nThe available "
                      "coordinates are: \n - %s" %(field, coord_list) ) 
                print("Please set the `coord` argument accordingly.")
                return(None,None)
        # Check the mode (for thetaMode)
        if self.geometry=="thetaMode":
            if (str(m) in self.avail_circ_modes) == False:
                mode_list = '\n - '.join(self.avail_circ_modes)
                print("The requested mode '%s' is not available.\nThe "
                    "available modes are: \n - %s" %(m, mode_list))
                return(None,None)
        
        # Find the output that corresponds to the requested time/iteration
        # (Modifies self.current_i and self.current_t)
        self._find_output( t, iteration )
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
        # - For 2D
        if self.geometry == "2dcartesian":
            F, extent = read_field_2d( filename, field_path )
        # - For 3D
        elif self.geometry == "3dcartesian":
            F, extent = read_field_3d(
                filename, field_path, slicing, slicing_dir)
        # - For thetaMode
        elif self.geometry == "thetaMode":
            if (coord in ['x', 'y']) and (self.avail_fields[field]=='vector'):
                # For Cartesian components, combine r and t components
                Fr, extent = read_field_circ( filename, field+'/r', m, theta )
                Ft, extent = read_field_circ( filename, field+'/t', m, theta )
                if coord == 'x':
                    F = np.cos(theta)*Fr - np.sin(theta)*Ft
                elif coord == 'y':
                    F = np.sin(theta)*Fr + np.cos(theta)*Ft
                # Revert the sign below the axis
                F[:len(F)/2] *= -1
            else:
                # For cylindrical or scalar components, no special treatment
                F, extent = read_field_circ( filename, field_path, m, theta )

        # Plot the resulting field
        # Deactivate plotting when there is no slice selection
        if (self.geometry=="3dcartesian") and (slicing is None):
            plot = False
        if plot==True:
            self.plotter.show_field( F, extent, slicing_dir, m,
                        field_label, self.geometry, self.current_i, **kw )

        # Return the result
        return( F, extent )


    def _find_output(self, t, iteration ) :
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
            raise ValueError("Please pass either a time (`t`) or an "
                             "iteration (`iteration`), but not both.")
        # If a time is requested
        elif (t is not None):
            # Make sur the time requested does not exceed the allowed bounds
            if t < self.tmin :
                self.current_i = 0
                print('Reached first iteration')
            elif t > self.tmax :
                self.current_i = len(self.t) -1
                print('Reached last iteration')
            # Find the last existing output
            else :
                self.current_i = self.t[ self.t <= t ].argmax()
        # If an iteration is requested 
        elif (iteration is not None):
            if (iteration in self.iterations):
                self.current_i = self.iterations.index(iteration)
            else:
                iter_list = '\n - '.join([ str(it) for it in self.iterations])
                print("The requested iteration '%s' is not available.\nThe "
                    "available modes are: \n - %s" %(iteration, iter_list))
                print("The first iteration is used instead.")
                self.current_i = 0
        else:
            pass # self.current_i retains its previous value

        # Register the value in the object
        self.current_t = self.t[ self.current_i ]

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
    filenames = [ name for (it, name) in iters_and_names ]
    iterations = [ it for (it, name) in iters_and_names ]

    return( filenames, iterations )

