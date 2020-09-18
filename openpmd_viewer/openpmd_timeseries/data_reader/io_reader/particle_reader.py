"""
This file is part of the openPMD-viewer.

It defines a function that reads a species record component (data & meta)
from an openPMD file

Copyright 2020, openPMD-viewer contributors
Authors: Axel Huebl
License: 3-Clause-BSD-LBNL
"""

import numpy as np
from scipy import constants
from .utilities import get_data


def read_species_data(series, iteration, species_name, component_name,
                      extensions):
    """
    Extract a given species' record_comp

    Parameters
    ----------
    series: openpmd_api.Series
        An open, readable openPMD-api series object

    iteration: integer
        Iteration from which parameters should be extracted

    species_name: string
        The name of the species to extract (in the openPMD file)

    component_name: string
        The record component to extract
        Either 'x', 'y', 'z', 'ux', 'uy', 'uz', or 'w'

    extensions: list of strings
        The extensions that the current OpenPMDTimeSeries complies with
    """
    it = series.iterations[iteration]

    # Translate the record component to the openPMD format
    dict_record_comp = {'x': ['position', 'x'],
                        'y': ['position', 'y'],
                        'z': ['position', 'z'],
                        'ux': ['momentum', 'x'],
                        'uy': ['momentum', 'y'],
                        'uz': ['momentum', 'z'],
                        'w': ['weighting', None]}
    
    if component_name in dict_record_comp:
        ompd_record_name, ompd_record_comp_name = \
            dict_record_comp[component_name]
    elif component_name.find('/') != -1:
        ompd_record_name, ompd_record_comp_name = \
            component_name.split('/')
    else:
        ompd_record_name = component_name
        ompd_record_comp_name = None

    # Extract the right dataset
    species = it.particles[species_name]
    record = species[ompd_record_name]
    if record.scalar:
        component = next(record.items())[1]
    else:
        component = record[ompd_record_comp_name]

    if ompd_record_name == 'id':
        output_type = np.uint64
    else:
        output_type = np.float64
    data = get_data( series, component, output_type=output_type )

    # For ED-PIC: if the data is weighted for a full macroparticle,
    # divide by the weight with the proper power
    # (Skip this if the current record component is the weight itself)
    if 'ED-PIC' in extensions and ompd_record_name != 'weighting':
        macro_weighted = record.get_attribute('macroWeighted')
        weighting_power = record.get_attribute('weightingPower')
        if (macro_weighted == 1) and (weighting_power != 0):
            w_component = next(species['weighting'].items())[1]
            w = get_data( w_component )
            data *= w ** (-weighting_power)

    # - Return positions, with an offset
    if component_name in ['x', 'y', 'z']:
        offset = get_data(series, species['positionOffset'][component_name])
        data += offset
    # - Return momentum in normalized units
    elif component_name in ['ux', 'uy', 'uz' ]:
        mass_component = next(species['mass'].items())[1]
        m = get_data(series, mass_component)
        # Normalize only if the particle mass is non-zero
        if np.all( m != 0 ):
            norm_factor = 1. / (m * constants.c)
            data *= norm_factor


    # Return the data
    return data
