"""
This file is part of the OpenPMD viewer.

It defines the main OpenPMDTimeSeries class.
"""
import os
from math import *
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

# Check wether the interactive interface can be loaded
try:
    from .interactive import InteractiveViewer, mode_dict
    parent_class = InteractiveViewer
except ImportError:
    print(
        '[opmd_viewer] Failed to import the interactive interface.\n'
        '(Make sure that ipywidgets and IPython.display are installed.)\n'
        'The opmd_viewer API is nonetheless working.')
    # Default parent class
    parent_class = object

# General dictionaries
mode_dict = { '0' : '0',
            '1 (real)' : '1', '1 (imag.)' : '2',
            '1' : '1 (real)', '2' : '1 (imag.)' }
slice_dict = { 'x':0, 'y':1, 'z':2 }
    
class OpenPMDTimeSeries(parent_class) :
    """
    Class that contains
    - the name of the files of the simulation
    - the list of iterations
    - the list of the corresponding times
    """

    def __init__( self, path_to_dir ) :
        """
        Initialize an OpenPMD time series

        Parameter
        ---------
        path_to_dir : string
            The absoluter path to the directory where the hdf5 files are.
        """
        # Extract the files and the iterations
        self.h5_files = list_h5_files( path_to_dir )
        self.iterations = np.array(
            [ int(filename[-10:-3]) for filename in self.h5_files ])

        # Check that there are HDF5 files in this directory
        if len(self.h5_files) == 0 :
            print("Error : Found no HDF5 files in the specified directory. \n"
                  "Please check that this is the full path to the HDF5 files.")
            return(None)

        # Go through the files of the series, check them and extract the time
        N_files = len(self.h5_files) 
        self.t = np.zeros( N_files )
        for k in range( N_files ) :
            # Open the file, and do a version check
            f = h5py.File( self.h5_files[k], 'r')
            version = f.attrs['openPMD'].decode()
            extension = f.attrs['openPMDextension']
            if version[:2] != '1.':
                raise RuntimeError(
                    "File {:s} is not supported: Invalid openPMD "
                    "version: {:s})".format( self.h5_files[k], version) )
            if extension != 1:
                raise RuntimeError(
                    "File {:s} is not supported: Not formatted with "
                    "the PIC extension".format(self.h5_files[k]) )
            # Find the base path, extract the time
            bpath = f[ f.attrs["basePath"] ]
            self.t[k] = bpath.attrs["time"] * bpath.attrs["timeUnitSI"]
            f.close()
        self.current_i = 0
        self.current_t = self.t[0]
        # Find the min and max of time
        self.tmin = self.t.min()
        self.tmax = self.t.max()
        
        # Find out whether fields and particles are present or not
        # If the fields are present, extract their geometry
        f = h5py.File( self.h5_files[0] )
        # - Find the base path object, and the path to the data
        bpath = f[ f.attrs["basePath"] ]   
        field_path = f.attrs['meshesPath'].decode().strip('/')
        particle_path = f.attrs['particlesPath'].decode().strip('/')
        # - Check for the fields
        if field_path in bpath.keys() :
            self.has_fields = True
            # Pick the first field and inspect its geometry
            first_field_path = next( iter(bpath[ field_path ]) )
            first_field = bpath[ os.path.join(field_path, first_field_path) ]
            self.geometry = first_field.attrs['geometry']
            if self.geometry == "cartesian":
                # Check if this a 2d or 3d Cartesian timeseries
                dim = len( first_field.attrs['axisLabels'] )
                if dim == 2:
                    self.geometry = "2dcartesian"
                elif dim==3:
                    self.geometry = "3dcartesian"
        else :
            self.has_fields = False
        # - Check for particles
        if particle_path in bpath.keys() :
            self.has_particles = True
        else :
            self.has_particles = False
        f.close()

        # Default potting parameters
        self.fontsize = 18
        
    def get_particle( self, t, quantity1='z', quantity2=None,
                   species='electrons', output=True, plot=False,
                   nbins=50, cmap='Blues', vmin=None, vmax=None, **kw ) :
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
        A tuple of 1darrays if two quantities are requested.
        """
        # Check that there is particle data
        if self.has_particles == False:
            print('No particle data in this time series')
            return()
        
        # Find the index of the time array that correponds to this time
        self.find_last_output(t)
        # Get the corresponding filename
        filename = self.h5_files[ self.current_i ]
        iteration = self.iterations[ self.current_i ]
        time_fs = 1.e15*self.current_t

        # In the case of only one quantity
        if quantity2 is None or quantity2=='None' :
            # Extract from file
            q1 = get_particle( filename, species, quantity1 )
            # Plot
            if plot :
                # Extract weights for the histogram
                w = get_particle( filename, species, 'w')
                # Do the plotting
                plt.hist(q1, bins=nbins, weights=w, **kw )
                plt.xlabel(quantity1, fontsize=self.fontsize)
                plt.title("t =  %.0f fs    (iteration %d)" \
                    %(time_fs, iteration), fontsize=self.fontsize)
            # Output
            if output :
                return(q1)

        # In the case of two quantities
        else :
            # Extract from file
            q1 = get_particle( filename, species, quantity1 )
            q2 = get_particle( filename, species, quantity2 )
            # Plot
            if plot :
                # Extract weights for the histogram
                w = get_particle( filename, species, 'w')
                # Do the plotting
                plt.hist2d(q1, q2, bins=nbins, cmap=cmap,
                    vmin=vmin, vmax=vmax, weights=w, **kw )
                plt.colorbar()

                plt.xlabel(quantity1, fontsize=self.fontsize)
                plt.ylabel(quantity2, fontsize=self.fontsize)
                plt.title("t =  %.1f fs   (iteration %d)"  \
                    %(time_fs, iteration ), fontsize=self.fontsize )
                    
            # Output
            if output :
                return(q1, q2)


    def get_field(self, t, field='E', coord='z', m=1, slicing=0.,
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

        m : int, optional
           Only used for thetaMode geometry
           0 : mode 0
           1 : real part of mode 1
           2 : imaginary part of mode 1

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
           z : a 1d array containing the values of the positions
           t : a float containing the time (in fs)
        """
        # Check that there is field data
        if self.has_fields == False:
            print('No field data in this time series')
            return()
        
        # Find the index of the time array that correponds to this time
        self.find_last_output(t)
        # Get the corresponding filename
        filename = self.h5_files[ self.current_i ]
        iteration = self.iterations[ self.current_i ]
        time_fs = 1.e15*self.current_t

        # Get the correct quantity
        if field == 'rho' :
            quantity = 'rho'
        else :
            quantity = '%s%s' %(field, coord)

        # Get the title and labels            
        if plot==True:
            # Cylindrical geometry
            if self.geometry == "thetaMode" :
                mode = mode_dict[str(m)]
                plt.title("%s in the mode %s at %.1f fs   (iteration %d)" \
                            %(quantity, mode, time_fs, iteration ),
                            fontsize=self.fontsize)
                plt.xlabel('$z \;(\mu m)$', fontsize=self.fontsize )
                plt.ylabel('$r \;(\mu m)$', fontsize=self.fontsize )
            # 2D Cartesian geometry
            elif self.geometry =="2dcartesian" :
                plt.title("%s at %.1f fs   (iteration %d)" \
                    %(quantity, time_fs, iteration ), fontsize=self.fontsize)
                plt.xlabel('$z \;(\mu m)$', fontsize=self.fontsize )
                plt.ylabel('$x \;(\mu m)$', fontsize=self.fontsize )
            # 3D Cartesian geometry
            elif self.geometry=="3dcartesian":
                plt.title("%s sliced across %s at %.1f fs  (iteration %d)" \
                    %(quantity, slicing_dir, time_fs, iteration ),
                    fontsize=self.fontsize)
                if slicing_dir=='x':
                    plt.xlabel('$z \;(\mu m)$', fontsize=self.fontsize )
                    plt.ylabel('$y \;(\mu m)$', fontsize=self.fontsize )
                elif slicing_dir=='y':
                    plt.xlabel('$z \;(\mu m)$', fontsize=self.fontsize )
                    plt.ylabel('$x \;(\mu m)$', fontsize=self.fontsize )
                elif slicing_dir=='z':
                    plt.xlabel('$y \;(\mu m)$', fontsize=self.fontsize )
                    plt.ylabel('$x \;(\mu m)$', fontsize=self.fontsize )

        # Get the fields
        return( get_field( filename, field, coord, m, slicing, slicing_dir,
                           output, plot, geometry=self.geometry, **kw ) )

    def find_last_output(self, t) :
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
    Return a list of the hdf5 files in this directory

    Parameter
    ---------
    path_to_dir : string
        The absoluter path to the directory where the hdf5 files are.

    Returns
    -------
    A list of strings, which correspond to the absolute path of each file.

    """
    # Find all the files in the provided directory
    all_files = os.listdir( path_to_dir )

    # Select the hdf5 files
    h5_files = []
    for filename in all_files :
        if filename[-3:] == '.h5' :
            h5_files.append( os.path.join( path_to_dir, filename) )

    # Sort them
    h5_files.sort()

    return( h5_files )


def get_particle( filename, species, quantity ) :
    """
    Extract a given particle quantity
    
    In the case of positions, the result
    is returned in microns
    
    Parameters
    ----------
    filename : string
        The name of the file from which to extract data
    
    species : string
        The name of the species to extract (in the OpenPMD file)

    quantity : string
        The quantity to extract
        Either 'x', 'y', 'z', 'ux', 'uy', 'uz', or 'w'

    """
    # Translate the quantity to the OpenPMD format
    dict_quantity = { 'x' : 'position/x',
                      'y' : 'position/y',
                      'z' : 'position/z',
                      'ux' : 'momentum/x',
                      'uy' : 'momentum/y',
                      'uz' : 'momentum/z',
                      'w' : 'weighting'}
    
    opmd_quantity = dict_quantity[quantity]

    # Open the HDF5 file
    dfile = h5py.File( filename, 'r' )
    base_path =  dfile.attrs['basePath'].decode()
    particles_path = dfile.attrs['particlesPath'].decode()

    # Find the right dataset
    species_grp =  dfile[ os.path.join( base_path, particles_path, species ) ]
    data = get_data( species_grp[ opmd_quantity ] )

    # - Return positions in microns, with an offset
    if quantity in ['x', 'y', 'z']:
        offset = get_data( species_grp[ 'positionOffset/%s' %quantity ] )
        return( 1.e6 * (data + offset) )
    # - Return momentum in normalized units
    elif quantity in ['ux', 'uy', 'uz' ]: 
        norm_factor = 1./( get_data( species_grp['mass'] ) * constants.c )
        return( data * norm_factor )
    # - Return the other quantities unchanged
    else :
        return( data )

def get_field( filename, field='E', coord='r', m=1, slicing=0., slicing_dir='y',
               output=True, plot=True, geometry="thetaMode", **kw ) :
    """
    Extract a given field from an HDF5 file in the OpenPMD format.
    
    Parameters
    ----------
    filename : string
       The absolute path to the HDF5 file
       
    field : string, optional
        Which field to extract
        Either 'rho', 'E', 'B' or 'J'

    coord : string, optional
       Which component of the field to extract
       Either 'x', 'y' or 'z' in '2dcartesian' or '3dcartesian' geometry
       Either 'r', 't' or 'z' in 'thetaMode' geometry

    m : int, optional
       Component to be extracted, in Circ mode
       0 : mode 0
       1 : real part of mode 1
       2 : imaginary part of mode 1

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

    geometry : string, optional
       OpenPMD designation of the geometry
       Either "2dcartesian", "3dcartesian" or "thetaMode"
       
    **kw : dict
       Additional options to be passed to matplotlib's imshow.

    Returns
    -------
    A tuple with
       F : a 2darray containing the required field
       extent : a 1darray with 4 elements, containing the extent
       z : a 1d array containing the values of the positions
       t : a float containing the time (in fs)
    """
    # Open the HDF5 file
    dfile = h5py.File( filename, 'r' )
    base_path = dfile.attrs["basePath"].decode()
    relative_fields_path = dfile.attrs["meshesPath"].decode()

    # Get the group of data and the corresponding information
    field_path = os.path.join( base_path, relative_fields_path, field )
    group = dfile[ field_path ]
    # Check the geometry
    if geometry=="thetaMode" :
        coords = ['r', 't', 'z']
        if coord==None : coord='r'
    elif geometry in ["2dcartesian", "3dcartesian"] :
        coords = ['x', 'y', 'z']
        if coord==None : coord='x'

    # Get the proper dataset
    if field == 'rho' :
        dset = dfile[ field_path ]
    elif field in ['E', 'B', 'J'] :
        if coord in coords :
            dset = dfile[ os.path.join( field_path, coord ) ]
        else :
            raise ValueError('Invalid `coord` : %s' %coord)
    else :
        raise ValueError('Invalid `field` : %s' %field)

    # Extract the data in cylindrical
    if geometry=="thetaMode" :
        F = np.array( dset[m,:,:] ) * dset.attrs["unitSI"]
        # Extract the extend
        Nr, Nz = F.shape
        dr, dz = group.attrs['gridSpacing']
        rmin, zmin = group.attrs['gridGlobalOffset']
        extent = np.array([ zmin-0.5*dz, zmin+0.5*dz+dz*Nz, 0., (Nr+1)*dr ])
        z = zmin + dz*np.arange(Nz)
    
    # Extract the data in 2D Cartesian
    elif geometry=="2dcartesian" :
        F = np.array( dset[:,:] ) * dset.attrs["unitSI"]
        # Extract the extend
        Nx, Nz = F.shape
        dx, dz = group.attrs['gridSpacing']
        xmin, zmin = group.attrs['gridGlobalOffset']
        extent = np.array([ zmin-0.5*dz, zmin+0.5*dz+dz*Nz,
                            xmin-0.5*dx, xmin+0.5*dx+dx*Nx ])
        z = zmin + dz*np.arange(Nz)
    
    # Extract the data in 3D Cartesian
    elif geometry=="3dcartesian" :
        # Dimensions of the grid
        Nx, Ny, Nz = dset.shape
        dx, dy, dz = group.attrs['gridSpacing']
        xmin, ymin, zmin = group.attrs['gridGlobalOffset']
        # Slice selection
        if slicing is not None:
            # Number of cells along the slicing direction
            n_cells = dset.shape[ slice_dict[slicing_dir] ]
            # Index of the slice (prevent stepping out of the array)
            i_cell = int( 0.5*(slicing+1.)*n_cells )
            i_cell = max( i_cell, 0 )
            i_cell = min( i_cell, n_cells-1)
            # Extraction of the data
            if slicing_dir=='x':
                F = np.array( dset[i_cell,:,:] ) * dset.attrs["unitSI"]
                extent = np.array([ zmin-0.5*dz, zmin+0.5*dz+dz*Nz,
                        xmin-0.5*dx, xmin+0.5*dx+dx*Nx ])
            elif slicing_dir=='y':
                F = np.array( dset[:,i_cell,:] ) * dset.attrs["unitSI"]
                extent = np.array([ zmin-0.5*dz, zmin+0.5*dz+dz*Nz,
                        ymin-0.5*dy, ymin+0.5*dy+dy*Ny ])
            elif slicing_dir=='z':
                F = np.array( dset[:,:,i_cell] ) * dset.attrs["unitSI"]
                extent = np.array([ ymin-0.5*dy, ymin+0.5*dy+dy*Ny,
                        xmin-0.5*dx, xmin+0.5*dx+dx*Nx ])
        else:
            # No slice selection, deactivate plotting
            plot = False
            F = np.array( dset[:,:,:] * dset.attrs["unitSI"] )
            extent = np.array([ zmin-0.5*dz, zmin+0.5*dz+dz*Nz,
                        ymin-0.5*dy, ymin+0.5*dy+dy*Ny,
                        xmin-0.5*dx, xmin+0.5*dx+dx*Nx ])

        # For the moment, do not return a particular axis (to be changed)
        z = None
        
    # Extract the time
    t = dfile[ base_path ].attrs['time']

    # Plot the result
    if plot :
        plt.imshow( F, extent=1.e6*extent, origin='lower',
            interpolation='nearest', aspect='auto', **kw )
        plt.colorbar()

    if output :
        return( F, extent, z, t )

def get_data( dset ) :
    """
    Extract the data from a (possibly constant) dataset

    Parameters:
    -----------
    dset: an h5py.Dataset or h5py.Group (when constant)
        The object from which the data is extracted

    Returns:
    --------
    An np.ndarray (non-constant dataset) or a single double (constant dataset)
    """
    # Case of a constant dataset
    if type(dset) is h5py.Group:
        data = dset.attrs['value']
    # Case of a non-constant dataset
    elif type(dset) is h5py.Dataset:
        data = dset[...]

    # Scale by the conversion factor
    data = data * dset.attrs['unitSI']

    return(data)

def get_shape( dset ) :
    """
    Extract the shape of a (possibly constant) dataset

    Parameters:
    -----------
    dset: an h5py.Dataset or h5py.Group (when constant)
        The object whose shape is extracted

    Returns:
    --------
    A tuple corresponding to the shape
    """
    # Case of a constant dataset
    if type(dset) is h5py.Group:
        shape = dset.attrs['shape']
    # Case of a non-constant dataset
    elif type(dset) is h5py.Dataset:
        shape = dset.shape

    return(shape)
    

