"""
MPI-enabled parallel versions of io_reader particle_reader functions.

Copyright 2020, openPMD-viewer contributors
License: 3-Clause-BSD-LBNL
"""
import numpy as np
from scipy import constants
from .utilities import _read_species_portion
from ..io_reader.utilities import get_data


def read_species_data(series, iteration, species_name, component_name,
                     extensions, comm):
    """
    Extract a given species' record_comp.
    
    Workload is divided among ranks, each rank collects its share,
    and results are gathered to rank 0.

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
        Either 'x', 'y', 'z', 'r', 'ux', 'uy', 'uz', 'ur', or 'w'

    extensions: list of strings
        The extensions that the current OpenPMDTimeSeries complies with
    
    comm : MPI communicator
        MPI communicator for parallel operations
    """
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    it = series.iterations[iteration]
    
    # Translate the record component to the openPMD format
    dict_record_comp = {'x': ['position', 'x'],
                        'y': ['position', 'y'],
                        'z': ['position', 'z'],
                        'r': ['position', 'r'],
                        'ux': ['momentum', 'x'],
                        'uy': ['momentum', 'y'],
                        'uz': ['momentum', 'z'],
                        'ur': ['momentum', 'r'],
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
    
    species_obj = it.particles[species_name]
    record = species_obj[ompd_record_name]
    if record.scalar:
        component = next(record.items())[1]
    else:
        component = record[ompd_record_comp_name]
    
    if rank == 0:
        shape = component.shape
        if len(shape) > 0:
            total_particles = shape[0]
        else:
            # Scalar value
            total_particles = 1
    else:
        total_particles = None
    
    total_particles = comm.bcast(total_particles, root=0)
    
    if total_particles > 1:
        chunk_size = total_particles // size
        remainder = total_particles % size
        
        start_idx = rank * chunk_size + min(rank, remainder)
        if rank < remainder:
            end_idx = start_idx + chunk_size + 1
        else:
            end_idx = start_idx + chunk_size
    else:
        # Scalar or single particle
        start_idx = 0
        end_idx = total_particles
    
    if ompd_record_name == 'id':
        output_type = np.uint64
    else:
        output_type = np.float64
    
    if total_particles > 1:
        data = _read_species_portion(series, component, start_idx, end_idx, output_type)
    else:
        if rank == 0:
            if ompd_record_name == 'id':
                output_type = np.uint64
            else:
                output_type = np.float64
            data = get_data(series, component, output_type=output_type)
        else:
            data = None
        data = comm.bcast(data, root=0)
        if component_name in ['x', 'y', 'z', 'r']:
            if rank == 0:
                offset = get_data(series, species_obj['positionOffset'][component_name])
            else:
                offset = None
            offset = comm.bcast(offset, root=0)
            data += offset
        elif component_name in ['ux', 'uy', 'uz', 'ur']:
            if rank == 0:
                mass_component = next(species_obj['mass'].items())[1]
                m = get_data(series, mass_component)
            else:
                m = None
            m = comm.bcast(m, root=0)
            if np.all(m != 0):
                norm_factor = 1. / (m * constants.c)
                data *= norm_factor
        return data
    
    if 'ED-PIC' in extensions and ompd_record_name != 'weighting':
        macro_weighted = record.get_attribute('macroWeighted')
        weighting_power = record.get_attribute('weightingPower')
        if (macro_weighted == 1) and (weighting_power != 0):
            w_component = next(species_obj['weighting'].items())[1]
            w = _read_species_portion(series, w_component, start_idx, end_idx, np.float64)
            data *= w ** (-weighting_power)
    
    if component_name in ['x', 'y', 'z', 'r']:
        if rank == 0:
            offset_full = get_data(series, species_obj['positionOffset'][component_name])
        else:
            offset_full = None
        offset_full = comm.bcast(offset_full, root=0)
        
        if np.isscalar(offset_full) or (hasattr(offset_full, 'shape') and offset_full.shape == ()):
            data += offset_full
        else:
            offset_portion = offset_full[start_idx:end_idx]
            data += offset_portion
    
    elif component_name in ['ux', 'uy', 'uz', 'ur']:
        mass_component = next(species_obj['mass'].items())[1]
        
        if rank == 0:
            mass_shape = mass_component.shape
            if len(mass_shape) > 0 and mass_shape[0] > 1:
                mass_is_per_particle = True
            else:
                mass_is_per_particle = False
        else:
            mass_is_per_particle = None
        
        mass_is_per_particle = comm.bcast(mass_is_per_particle, root=0)
        
        if mass_is_per_particle:
            m = _read_species_portion(series, mass_component, start_idx, end_idx, np.float64)
            if np.all(m != 0):
                norm_factor = 1. / (m * constants.c)
                data *= norm_factor
        else:
            if rank == 0:
                m = get_data(series, mass_component)
            else:
                m = None
            m = comm.bcast(m, root=0)
            if np.all(m != 0):
                norm_factor = 1. / (m * constants.c)
                data *= norm_factor
    
    gathered_data = comm.gather(data, root=0)
    
    if rank == 0:
        data = np.concatenate(gathered_data)
    else:
        data = None
    
    data = comm.bcast(data, root=0)
    
    return data
