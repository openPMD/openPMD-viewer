"""
This file is part of the OpenPMD viewer.

It defines a function that reads particle data from an openPMD file
"""
import os
import h5py
from scipy import constants
import numpy as np
from .utilities import get_data

def read_particle( filename, species, quantity ) :
    """
    Extract a given particle quantity
    
    In the case of positions, the result is returned in microns
    
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
        data = 1.e6 * (data + offset)
    # - Return momentum in normalized units
    elif quantity in ['ux', 'uy', 'uz' ]: 
        norm_factor = 1./( get_data( species_grp['mass'] ) * constants.c )
        data = data * norm_factor

    # Close the HDF5 file and return the data
    dfile.close()
    return( data )
