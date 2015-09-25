"""
This file is part of the openPMD viewer.

It defines functions that can read the fields from an HDF5 file.
"""
import os
import h5py
import numpy as np
from .utilities import slice_dict, get_shape, get_data

def read_field_2d( filename, field_path ):
    """
    Extract a given field from an HDF5 file in the OpenPMD format,
    when the geometry is 2d cartesian.
    
    Parameters
    ----------
    filename : string
       The absolute path to the HDF5 file
       
    field_path : string
       The relative path to the requested field, from the openPMD meshes path
       (e.g. 'rho', 'E/r', 'B/x')
    """
    # Open the HDF5 file
    dfile = h5py.File( filename, 'r' )
    # Extract the dataset and and corresponding group
    group, dset = find_dataset( dfile, field_path )
    
    # Extract the data in 2D Cartesian
    F = get_data( dset )
    # Extract the extent
    Nx, Nz = F.shape
    dx, dz = group.attrs['gridSpacing']
    xmin, zmin = group.attrs['gridGlobalOffset']
    extent = np.array([ zmin-0.5*dz, zmin+0.5*dz+dz*Nz,
            xmin-0.5*dx, xmin+0.5*dx+dx*Nx ])
    
    # Close the file
    dfile.close()
    return( F, extent )    

def read_field_circ( filename, field_path, m=0, theta=0. ) :
    """
    Extract a given field from an HDF5 file in the OpenPMD format,
    when the geometry is 2d cartesian.
    
    Parameters
    ----------
    filename : string
       The absolute path to the HDF5 file
       
    field_path : string
       The relative path to the requested field, from the openPMD meshes path
       (e.g. 'rho', 'E/r', 'B/x')

    m : int or string, optional
       The azimuthal mode to be extracted

    theta : float, optional
       Angle of the plane of observation with respect to the x axis

    Returns
    -------
    A tuple with
       F : a 2darray containing the required field
       extent : a 1darray with 4 elements, containing the extent
    """
    # Open the HDF5 file
    dfile = h5py.File( filename, 'r' )
    # Extract the dataset and and corresponding group
    group, dset = find_dataset( dfile, field_path )

    # Extract the modes and recombine them properly
    Nm, Nr, Nz = get_shape( dset )
    F_total = np.zeros( (2*Nr, Nz ) )
    if m=='all':
        # Sum of all the modes
        # - Prepare the multiplier arrays
        mult_above_axis = [1]
        mult_below_axis = [1]
        for mode in range(1,int(Nm/2)+1):
            cos = np.cos( mode*theta )
            sin = np.sin( mode*theta )
            mult_above_axis += [cos, sin]
            mult_below_axis += [ (-1)**mode*cos, (-1)**mode*sin ]
        mult_above_axis = np.array( mult_above_axis )
        mult_below_axis = np.array( mult_below_axis )
        # - Sum the modes
        F = get_data( dset ) # (Extracts all modes)
        F_total[Nr:,:] = np.tensordot( mult_above_axis, F, axes=(0,0) )[:,:]
        F_total[:Nr,:] = np.tensordot( mult_below_axis, F, axes=(0,0) )[::-1,:]
    elif m==0:
        # Extract mode 0
        F = get_data( dset, 0, 0 )
        F_total[Nr:,:] = F[:,:]
        F_total[:Nr,:] = F[::-1,:]
    else:
        # Extract higher mode
        cos = np.cos( m*theta )
        sin = np.sin( m*theta )
        F_cos = get_data( dset, 2*m-1, 0 )
        F_sin = get_data( dset, 2*m, 0 )
        F = cos*F_cos + sin*F_sin
        F_total[Nr:,:] = F[:,:]
        F_total[:Nr,:] = (-1)**m * F[::-1,:] 
    # Extract the extent
    dr, dz = group.attrs['gridSpacing']
    rmin, zmin = group.attrs['gridGlobalOffset']
    extent = np.array([ zmin-0.5*dz, zmin+0.5*dz+dz*Nz,
                        -(Nr+1)*dr, (Nr+1)*dr ])

    # Close the file
    dfile.close()
    return( F_total, extent )


def read_field_3d( filename, field_path, slicing=0., slicing_dir='y' ) :
    """
    Extract a given field from an HDF5 file in the OpenPMD format,
    when the geometry is 3d cartesian.
    
    Parameters
    ----------
    filename : string
       The absolute path to the HDF5 file
       
    field_path : string
       The relative path to the requested field, from the openPMD meshes path
       (e.g. 'rho', 'E/r', 'B/x')

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

    Returns
    -------
    A tuple with
       F : a 2darray containing the required field
       extent : a 1darray with 4 elements, containing the extent
    """
    # Open the HDF5 file
    dfile = h5py.File( filename, 'r' )
    # Extract the dataset and and corresponding group
    group, dset = find_dataset( dfile, field_path )
    
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

    # Close the file
    dfile.close()
    return( F, extent )

        
def find_dataset( dfile, field_path ):
    """
    Extract the dataset that corresponds to field_path,
    and the corresponding group

    (In the case of scalar records, the group and the dataset are identical.
    In the case of vector records, the group contains all the components
    and the dataset corresponds to one given component.)

    Parameters
    ----------
    dfile: an h5Py.File object
       The file from which to extract the dataset
     
    field_path : string
       The relative path to the requested field, from the openPMD meshes path
       (e.g. 'rho', 'E/r', 'B/x')

    Returns
    -------
    A tuple with:
    - an h5py.Group object
    - an h5py.Dataset object
    """
    # Find the meshes path
    base_path = dfile.attrs["basePath"].decode()
    relative_meshes_path = dfile.attrs["meshesPath"].decode()

    # Get the proper dataset
    full_field_path = os.path.join(base_path, relative_meshes_path, field_path)
    dset = dfile[ full_field_path ]
    # Get the proper group
    group_path = field_path.split('/')[0]
    full_group_path = os.path.join(base_path, relative_meshes_path, group_path)
    group = dfile[ full_group_path ]

    return( group, dset )
