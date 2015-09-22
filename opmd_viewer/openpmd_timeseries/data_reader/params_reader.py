"""
This file is part of the openPMD viewer.

It defines a function that can read standard parameters from an openPMD file.
"""
import os
import h5py

def read_openPMD_params( filename ):
    """
    Extract the time and some openPMD parameters from a file

    Parameter
    ---------
    filename: string
        The path to the file from which parameters should be extracted

    Returns
    -------
    A tuple with:
    - A float corresponding to the time of this iteration in SI units
    - A dictionary containing several parameters, such as the geometry, etc.
    """
    params = {}
            
    # Open the file, and do a version check
    f = h5py.File( filename, 'r')
    version = f.attrs['openPMD'].decode()
    if version[:2] != '1.':
        raise ValueError(
            "File %s is not supported: Invalid openPMD version: "
            "%s)" %( filename, version) )
    params['extension'] = f.attrs['openPMDextension']

    # Find the base path object, and extract the time
    bpath = f[ f.attrs["basePath"] ]
    t = bpath.attrs["time"] * bpath.attrs["timeUnitSI"]

    # Find out whether fields are present and extract their geometry
    field_path = f.attrs['meshesPath'].decode().strip('/')
    if field_path in bpath.keys():
        params['has_fields'] = True
        # Pick the first field and inspect its geometry
        first_field_path = next( iter(bpath[ field_path ]) )
        first_field = bpath[ os.path.join(field_path, first_field_path) ]
        params['geometry'] = first_field.attrs['geometry']
        if params['geometry'] == "cartesian":
            # Check if this a 2d or 3d Cartesian timeseries
            dim = len( first_field.attrs['axisLabels'] )
            if dim == 2:
                params['geometry'] = "2dcartesian"
            elif dim==3:
                params['geometry'] = "3dcartesian"
    else :
        params['has_fields'] = False

    # Find out whether particles are present
    particle_path = f.attrs['particlesPath'].decode().strip('/')
    if particle_path in bpath.keys():
        params['has_particles'] = True
    else :
        params['has_particles'] = False
    
    # Close the file and return the parameters
    f.close()
    return( t, params )
