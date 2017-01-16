"""
This file is part of the openPMD-viewer.

It defines the ParticleTracker class.

Copyright 2015-2016, openPMD-viewer contributors
Authors: Remi Lehe
License: 3-Clause-BSD-LBNL
"""
import numpy as np
from .data_reader.particle_reader import read_species_data
try:
    import numba
    numba_available = True
except ImportError:
    numba_available = False


class ParticleTracker( object ):
    """
    TODO: Finish docstring, explain usage
Say this requires that the id are output

    """

    def __init__(self, ts, species=None, t=None,
                iteration=None, select=None, preserve_particle_index=False):
        """
        TODO: Finish docstring + explain preserve_particle_index

        Parameters
        ----------
        ts: an OpenPMDTimeSeries object
            Contains the data on the particles

        species: string
            A string indicating the name of the species
            This is optional if there is only one species

        t : float (in seconds), optional
            Time at which to obtain the data (if this does not correspond to
            an available file, the last file before `t` will be used)
            Either `t` or `iteration` should be given by the user.

        iteration : int
            The iteration at which to obtain the data
            Either `t` or `iteration` should be given by the user.

        select: dict, optional
            Either None or a dictionary of rules
            to select the particles, of the form
            'x' : [-4., 10.]  (Particles having x between -4 and 10 microns)
            'ux' : [-0.1, 0.1] (Particles having ux between -0.1 and 0.1 mc)
            'uz' : [5., None]  (Particles with uz above 5 mc)
        """
        # Extract the particle id and sort them
        self.selected_pid, = ts.get_particle(['id'], species=species,
                                select=select, t=t, iteration=iteration)
        self.selected_pid.sort()

        # Register a few metadata
        self.N_selected = len( self.selected_pid )
        self.species = species
        self.preserve_particle_index = preserve_particle_index

    def extract_tracked_particles( self, file_handle, data_list,
                                    species, extensions ):
        """
        Select the elements of each particle quantities in data_list,
        so as to only return those that correspond to the tracked particles

        Parameters
        ----------
        file_handle: h5py.File object
            The HDF5 file from which to extract data

        data_list: list of 1darrays
            A list of arrays with one element per macroparticle, that represent
            different particle quantities

        species: string
            Name of the species being requested

        extensions: list of strings
            The extensions that the current OpenPMDTimeSeries complies with

        Returns
        -------
        A list of 1darrays that correspond to data_list, but were only the
        macroparticles that are tracked are kept
        # TODO: Say more about the NaNs
        """
        # Extract the particle id, and get the extraction indices
        pid = read_species_data(file_handle, species, 'id', extensions)
        selected_indices = self.get_extraction_indices( pid )

        # For each particle quantity, select only the tracked particles
        for i in range(len(data_list)):
            if len(data_list[i]) > 1:  # Do not apply selection on scalars
                data_list[i] = self.extract_quantity(
                    data_list[i], selected_indices )

        return( data_list )

    def extract_quantity( self, q, selected_indices ):
        """
        TODO: Do docstring
        """
        selected_q = q[ selected_indices ]
        if self.preserve_particle_index:
            if q.dtype in [ np.float64, np.float32 ]:
                # Fill the position of absent particles by NaNs
                selected_q = np.where( selected_indices == -1,
                                        np.nan, selected_q)
            else:
                # The only non-float quantity in openPMD-viewer is  particle id
                selected_q = self.selected_pid
        return( selected_q )

    def get_extraction_indices( self, pid ):
        """
        TODO: Do docstring
        """
        # TODO: Explain the algorithmic approach: how to avoid N x N_selected
        # Sort the pid, and keep track of the original index
        # at which each pid was
        original_indices = pid.argsort()
        sorted_pid = pid[ original_indices ]

        # Extract only the indices for which sorted_pid is one of pid
        # in self.sselected_pid (i.e. which correpond to one
        # of the original particles)
        selected_indices = np.empty( self.N_selected, dtype=np.int64 )
        N_extracted = extract_indices( original_indices, selected_indices,
            sorted_pid, self.selected_pid, self.preserve_particle_index)

        # If there are less particles then self.N_selected
        # (i.e. not all the pid in self.selected_pid were in sorted_pid)
        if N_extracted < self.N_selected:
            if self.preserve_particle_index:
                # Finish filling the array and indicate that absent particles
                selected_indices[N_extracted:] = -1
            else:
                # Resize the array
                selected_indices = \
                    selected_indices[:N_extracted]

        return( selected_indices )


def extract_indices( original_indices, selected_indices,
                        pid, selected_pid, preserve_particle_index ):
    # TODO: Explain algorithm: both set of pids are sorted
    i = 0
    i_select = 0
    i_fill = 0
    N = len(pid)
    N_selected = len(selected_pid)

    # Go through both sorted arrays (pid and selected_pid) and match them.
    # i.e. whenever the same number appears in both arrays,
    # record the corresponding original index in selected_indices
    while i < N and i_select < N_selected:

        if pid[i] < selected_pid[i_select]:
            i += 1
        elif pid[i] == selected_pid[i_select]:
            selected_indices[i_fill] = original_indices[i]
            i_fill += 1
            i_select += 1
        elif pid[i] > selected_pid[i_select]:
            i_select += 1
            if preserve_particle_index:
                # Fill the index, to indicate that the particle is absent
                selected_indices[i_fill] = -1
                i_fill += 1

    return( i_fill )

if numba_available:
    # Compile the Python function, to avoid bottleneck
    # (all other operations are numpy operations)
    extract_indices = numba.jit( extract_indices, nopython=True )
