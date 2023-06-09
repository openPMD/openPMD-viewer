"""
This file is part of the openPMD-viewer.

It defines a function that reads a species record component (data & meta)
from an openPMD file

Copyright 2015-2016, openPMD-viewer contributors
Authors: Remi Lehe, Axel Huebl
License: 3-Clause-BSD-LBNL
"""

import h5py
import numpy as np
from scipy import constants
from .utilities import get_data, join_infile_path


def read_species_data(filename, iteration, species, record_comp, extensions):
    """
    Extract a given species' record_comp

    Parameters
    ----------
    filename : string
       The absolute path to the HDF5 file

    iteration : int
        The iteration at which to obtain the data

    species: string
        The name of the species to extract (in the openPMD file)

    record_comp: string
        The record component to extract
        Either 'x', 'y', 'z', 'ux', 'uy', 'uz', or 'w'

    extensions: list of strings
        The extensions that the current OpenPMDTimeSeries complies with
    """
    # Open the HDF5 file
    dfile = h5py.File( filename, 'r' )
    # Translate the record component to the openPMD format
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
    base_path = '/data/{0}'.format( iteration )
    particles_path = dfile.attrs['particlesPath'].decode()

    # Extract the right dataset
    species_grp = dfile[
        join_infile_path(base_path, particles_path, species) ]
    if opmd_record_comp == 'id':
        output_type = np.uint64
    else:
        output_type = np.float64
    data = get_data( species_grp[ opmd_record_comp ], output_type=output_type )

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

    # - Return positions, with an offset
    if record_comp in ['x', 'y', 'z']:
        offset = get_data(species_grp['positionOffset/%s' % record_comp])
        data += offset
    # - Return momentum in normalized units
    elif record_comp in ['ux', 'uy', 'uz' ]:
        m = get_data(species_grp['mass'])
        # Normalize only if the particle mass is non-zero
        if np.all( m != 0 ):
            norm_factor = 1. / (m * constants.c)
            data *= norm_factor

    # Close the file
    dfile.close()
    # Return the data
    return(data)
