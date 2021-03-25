"""
This file is part of the openPMD-viewer.

It defines a function that can read standard parameters from an openPMD file.

Copyright 2020, openPMD-viewer contributors
Authors: Axel Huebl
License: 3-Clause-BSD-LBNL
"""

import numpy as np
from .utilities import join_infile_path


def read_openPMD_params(series, iteration, extract_parameters=True):
    """
    Extract the time and some openPMD parameters from a file

    Parameter
    ---------
    series: openpmd_api.Series
        An open, readable openPMD-api series object

    iteration: integer
        Iteration from which parameters should be extracted

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
    it = series.iterations[iteration]

    # extract the time
    if callable(it.time):  # prior to openPMD-api 0.13.0
        t = it.time() * it.time_unit_SI()
    else:
        t = it.time * it.time_unit_SI

    # If the user did not request more parameters, close file and exit
    if not extract_parameters:
        return t, None

    # Otherwise, extract the rest of the parameters
    params = {}

    # Find out supported openPMD extensions claimed by this file
    # note: a file might implement multiple extensions
    known_extensions = {'ED-PIC': np.uint32(1)}
    bitmask_all_extensions = series.openPMD_extension
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
    fields_available = len(it.meshes) > 0
    if fields_available:
        params['avail_fields'] = []
        params['fields_metadata'] = {}

        # Loop through the available fields
        for field_name, field in it.meshes.items():
            metadata = {}
            metadata['geometry'] = field.get_attribute('geometry')
            metadata['axis_labels'] = field.axis_labels

            # Swap the order of the labels if the code that wrote the HDF5 file
            # was Fortran order (i.e. reverse order with respect to Python)
            if field.get_attribute('dataOrder') == 'F':
                metadata['axis_labels'] = metadata['axis_labels'][::-1]
            # Check whether the field is a vector or a scalar
            if field.scalar:
                metadata['type'] = 'scalar'
                components = []
            else:
                metadata['type'] = 'vector'
                components = [comp for comp, _ in field.items()]
            # Register available components
            metadata['avail_components'] = components
            # Check the number of modes
            if metadata['geometry'] == "thetaMode":
                # simply check first record component
                field_component = next(field.items())[1]
                Nm = field_component.shape[0]
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
    particles_available = len(it.particles) > 0
    if particles_available:
        # Particles are present ; extract the species
        params['avail_species'] = list(it.particles)
        # dictionary with list of record components for each species
        record_components = {}
        # Go through all species
        # TODO: I did this more elegant in ParaView... for later
        for species_name, species in it.particles.items():
            record_components[species_name] = []

            # Go through all the particle records of this species
            for record_name, record in species.items():
                # Skip the particlePatches, which are not used here.
                # API should filter this...
                # if record_name == 'particlePatches':
                #     continue
                if record.scalar:
                    # Add the name of the scalar record
                    record_components[species_name]. \
                        append(record_name)
                else:
                    # Add each component of the vector record
                    for compo_name in list(record):
                        record_components[species_name]. \
                            append(join_infile_path(record_name, compo_name))
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

    return t, params


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

    return record_comps
