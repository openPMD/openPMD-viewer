"""
This file is part of the openPMD-viewer.

It defines a function that can read standard parameters from an openPMD file.

Copyright 2015-2016, openPMD-viewer contributors
Authors: Remi Lehe, Axel Huebl
License: 3-Clause-BSD-LBNL
"""

import h5py
import numpy as np
from .utilities import is_scalar_record, get_shape, get_bpath, join_infile_path


def read_openPMD_params(filename, extract_parameters=True):
    """
    Extract the time and some openPMD parameters from a file

    Parameter
    ---------
    filename: string
        The path to the file from which parameters should be extracted

    extract_parameters: bool, optional
        Whether to extract all parameters or only the time
        (Function execution is faster when extract_parameters is False)

    Returns
    -------
    A tuple with:
    - A float corresponding to the time of this iteration in SI units
    - A dictionary containing several parameters, such as the geometry, etc.
      When extract_parameters is False, the second argument returned is None.
    """
    # Open the file, and do a version check
    f = h5py.File(filename, 'r')
    version = f.attrs['openPMD'].decode()
    if version[:2] != '1.':
        raise ValueError(
            "File %s is not supported: Invalid openPMD version: "
            "%s)" % (filename, version))

    # Find the base path object, and extract the time
    bpath = f[get_bpath(f)]
    t = bpath.attrs["time"] * bpath.attrs["timeUnitSI"]

    # If the user did not request more parameters, close file and exit
    if not extract_parameters:
        f.close()
        return(t, None)

    # Otherwise, extract the rest of the parameters
    params = {}

    # Find out supported openPMD extensions claimed by this file
    # note: a file might implement multiple extensions
    known_extensions = {'ED-PIC': np.uint32(1)}
    bitmask_all_extensions = f.attrs['openPMDextension']
    params['extensions'] = []
    for extension, bitmask in known_extensions.items():
        # This uses a bitmask to identify activated extensions
        # efficiently in static programming languages via
        # a single attribute and a binary AND (&) operation.
        # Standard: https://git.io/vwnMw
        # Bitmasks: https://en.wikipedia.org/wiki/Mask_%28computing%29
        if bitmask_all_extensions & bitmask == bitmask:
            params['extensions'].append(extension)

    # Find out whether fields are present and extract their metadata
    fields_available = False
    if ('meshesPath' in f.attrs):        # Check for openPMD 1.1 files
        meshes_path = f.attrs['meshesPath'].decode().strip('/')
        if meshes_path in bpath.keys():  # Check for openPMD 1.0 files
            fields_available = True
    if fields_available:
        params['avail_fields'] = []
        params['fields_metadata'] = {}

        # Loop through the available fields
        for field_name in bpath[meshes_path].keys():
            field = bpath[join_infile_path(meshes_path, field_name)]
            metadata = {}
            metadata['geometry'] = field.attrs['geometry'].decode()
            metadata['axis_labels'] = [ coord.decode() for coord in
                                        field.attrs['axisLabels'] ]
            # Swap the order of the labels if the code that wrote the HDF5 file
            # was Fortran order (i.e. reverse order with respect to Python)
            if field.attrs['dataOrder'].decode() == 'F':
                metadata['axis_labels'] = metadata['axis_labels'][::-1]
            # Check whether the field is a vector or a scalar
            if is_scalar_record(field):
                metadata['type'] = 'scalar'
            else:
                metadata['type'] = 'vector'
            # Check the number of modes
            if metadata['geometry'] == "thetaMode":
                if is_scalar_record(field):
                    Nm, _, _ = get_shape(field)
                else:
                    coord = list(field.keys())[0]
                    Nm, _, _ = get_shape(field[coord])
                metadata['avail_circ_modes'] = ['all'] + \
                    [str(m) for m in range(int(Nm / 2) + 1)]
            # Check if this a 1d, 2d or 3d Cartesian
            elif metadata['geometry'] == "cartesian":
                dim = len(metadata['axis_labels'])
                if dim == 1:
                    metadata['geometry'] = "1dcartesian"
                elif dim == 2:
                    metadata['geometry'] = "2dcartesian"
                elif dim == 3:
                    metadata['geometry'] = "3dcartesian"
                metadata['avail_circ_modes'] = []

            params['avail_fields'].append( field_name )
            params['fields_metadata'][field_name] = metadata

    else:
        params['avail_fields'] = None

    # Find out whether particles are present, and if yes of which species
    particles_available = False
    if ('particlesPath' in f.attrs):        # Check for openPMD 1.1 files
        particle_path = f.attrs['particlesPath'].decode().strip('/')
        if particle_path in bpath.keys():   # Check for openPMD 1.0 files
            # Check that there is at least one species
            if len(bpath[particle_path].keys()) > 0:
                particles_available = True
    if particles_available:
        # Particles are present ; extract the species
        params['avail_species'] = []
        for species_name in bpath[particle_path].keys():
            params['avail_species'].append(species_name)
        # dictionary with list of record components for each species
        record_components = {}
        # Go through all species
        for species_name in iter(params['avail_species']):
            species = bpath[join_infile_path(particle_path, species_name)]
            record_components[species_name] = []

            # Go through all the particle records of this species
            for record_name in species.keys():
                # Skip the particlePatches, which are not used here.
                if record_name == 'particlePatches':
                    continue
                record = species[record_name]
                if is_scalar_record(record):
                    # Add the name of the scalar record
                    record_components[species_name]. \
                        append(record_name)
                else:
                    # Add each component of the vector record
                    for coord in record.keys():
                        record_components[species_name]. \
                            append(join_infile_path(record_name, coord))
            # Simplify the name of some standard openPMD records
            record_components[species_name] = \
                simplify_record(record_components[species_name])
        params['avail_record_components'] = record_components
        # deprecated
        first_species_name = next(iter(params['avail_species']))
        params['avail_ptcl_quantities'] = \
            record_components[first_species_name]
    else:
        # Particles are absent
        params['avail_species'] = None
        params['avail_record_components'] = None
        # deprecated
        params['avail_ptcl_quantities'] = None

    # Close the file and return the parameters
    f.close()
    return(t, params)


def simplify_record(record_comps):
    """
    Replace the names of some standard record by shorter names

    Parameter
    ---------
    record_comps: a list of strings
        A list of available particle record components

    Returns
    -------
    A list with shorter names, where applicable
    """
    # Replace the names of the positions
    if ('position/x' in record_comps) and ('positionOffset/x' in record_comps):
        record_comps.remove('position/x')
        record_comps.remove('positionOffset/x')
        record_comps.append('x')
    if ('position/y' in record_comps) and ('positionOffset/y' in record_comps):
        record_comps.remove('position/y')
        record_comps.remove('positionOffset/y')
        record_comps.append('y')
    if ('position/z' in record_comps) and ('positionOffset/z' in record_comps):
        record_comps.remove('position/z')
        record_comps.remove('positionOffset/z')
        record_comps.append('z')

    # Replace the names of the momenta
    if 'momentum/x' in record_comps:
        record_comps.remove('momentum/x')
        record_comps.append('ux')
    if 'momentum/y' in record_comps:
        record_comps.remove('momentum/y')
        record_comps.append('uy')
    if 'momentum/z' in record_comps:
        record_comps.remove('momentum/z')
        record_comps.append('uz')

    # Replace the name for 'weights'
    if 'weighting' in record_comps:
        record_comps.remove('weighting')
        record_comps.append('w')

    return(record_comps)
