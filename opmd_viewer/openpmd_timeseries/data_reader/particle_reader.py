"""
This file is part of the openPMD viewer.

It defines a function that reads a species record component (data & meta)
from an openPMD file

__authors__ = "Remi Lehe, Axel Huebl"
__copyright__ = "Copyright 2015-2016, openPMD viewer contributors"
__license__ = "3-Clause-BSD-LBNL"
"""

import os
from scipy import constants
from .utilities import get_data, get_bpath


def read_species_data(file_handle, species, record_comp, extensions):
    """
    Extract a given species' record_comp

    In the case of positions, the result is returned in microns

    Parameters
    ----------
    file_handle: h5py.File object
        The HDF5 file from which to extract data

    species: string
        The name of the species to extract (in the OpenPMD file)

    record_comp: string
        The record component to extract
        Either 'x', 'y', 'z', 'ux', 'uy', 'uz', or 'w'

    extensions: list of strings
        The extensions that the current openPMDTimeseries complies with
    """
    # Translate the record component to the OpenPMD format
    dict_record_comp = {'x': 'position/x',
                        'y': 'position/y',
                        'z': 'position/z',
                        'ux': 'momentum/x',
                        'uy': 'momentum/y',
                        'uz': 'momentum/z',
                        'w': 'weighting'}
    if record_comp in dict_record_comp:
        opmd_record_comp = dict_record_comp[record_comp]
    else:
        opmd_record_comp = record_comp

    # Open the HDF5 file
    base_path = get_bpath(file_handle)
    particles_path = file_handle.attrs['particlesPath'].decode()

    # Extract the right dataset
    species_grp = file_handle[
        os.path.join(base_path, particles_path, species) ]
    data = get_data(species_grp[ opmd_record_comp ])

    # For ED-PIC: if the data is weighted for a full macroparticle,
    # divide by the weight with the proper power
    # (Skip this if the current record component is the weight itself)
    if 'ED-PIC' in extensions and opmd_record_comp != 'weighting':
        opmd_record = opmd_record_comp.split('/')[0]
        record_dset = species_grp[ opmd_record ]
        macro_weighted = record_dset.attrs['macroWeighted']
        weighting_power = record_dset.attrs['weightingPower']
        if (macro_weighted == 1) and (weighting_power != 0):
            w = get_data( species_grp[ 'weighting' ] )
            data *= w ** (-weighting_power)

    # - Return positions in microns, with an offset
    if record_comp in ['x', 'y', 'z']:
        offset = get_data(species_grp['positionOffset/%s' % record_comp])
        data += offset
        data *= 1.e6
    # - Return momentum in normalized units
    elif record_comp in ['ux', 'uy', 'uz' ]:
        norm_factor = 1. / (get_data(species_grp['mass']) * constants.c)
        data *= norm_factor

    # Return the data
    return(data)
