"""
This file is part of the openPMD viewer.

It defines functions that can read the fields from an HDF5 file.
"""
import os
import h5py
import numpy as np
from .utilities import slice_dict, get_shape, get_data

def read_field( filename, field_path, m=0, slicing=0.,
               slicing_dir='y', geometry="thetaMode" ) :
    """
    Extract a given field from an HDF5 file in the OpenPMD format.
    
    Parameters
    ----------
    filename : string
       The absolute path to the HDF5 file
       
    field_path : string
       The relative path to the requested field, from the openPMD meshes path
       (e.g. 'rho', 'E/r', 'B/x')
    
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

    Returns
    -------
    A tuple with
       F : a 2darray containing the required field
       extent : a 1darray with 4 elements, containing the extent
    """
    # Open the HDF5 file
    dfile = h5py.File( filename, 'r' )
    base_path = dfile.attrs["basePath"].decode()
    relative_meshes_path = dfile.attrs["meshesPath"].decode()

    # Get the proper dataset
    full_field_path = os.path.join(base_path, relative_meshes_path, field_path)
    dset = dfile[ full_field_path ]
    # Get the proper group
    group_path = field_path.split('/')[0]
    full_group_path = os.path.join(base_path, relative_meshes_path, group_path)
    group = dfile[ full_group_path ]

    # Extract the data in cylindrical
    if geometry=="thetaMode" :
        F = get_data( dset, m, 0 )
        # Extract the extent
        Nr, Nz = F.shape
        dr, dz = group.attrs['gridSpacing']
        rmin, zmin = group.attrs['gridGlobalOffset']
        extent = np.array([ zmin-0.5*dz, zmin+0.5*dz+dz*Nz, 0., (Nr+1)*dr ])
        z = zmin + dz*np.arange(Nz)
    
    # Extract the data in 2D Cartesian
    elif geometry=="2dcartesian" :
        F = get_data( dset )
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
        Nx, Ny, Nz = get_shape( dset )
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
                F = get_data( dset, i_cell, 0 )
                extent = np.array([ zmin-0.5*dz, zmin+0.5*dz+dz*Nz,
                        xmin-0.5*dx, xmin+0.5*dx+dx*Nx ])
            elif slicing_dir=='y':
                F = get_data( dset, i_cell, 1 )
                extent = np.array([ zmin-0.5*dz, zmin+0.5*dz+dz*Nz,
                        ymin-0.5*dy, ymin+0.5*dy+dy*Ny ])
            elif slicing_dir=='z':
                F = get_data( dset, i_cell, 2 )
                extent = np.array([ ymin-0.5*dy, ymin+0.5*dy+dy*Ny,
                        xmin-0.5*dx, xmin+0.5*dx+dx*Nx ])
        else:
            F = get_data( dset )
            extent = np.array([ zmin-0.5*dz, zmin+0.5*dz+dz*Nz,
                        ymin-0.5*dy, ymin+0.5*dy+dy*Ny,
                        xmin-0.5*dx, xmin+0.5*dx+dx*Nx ])

    return( F, extent )
