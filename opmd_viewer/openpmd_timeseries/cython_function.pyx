import numpy as np
cimport numpy as np
from cpython cimport bool

def extract_indices_cython(
    np.ndarray[np.int64_t, ndim=1] original_indices,
    np.ndarray[np.int64_t, ndim=1] selected_indices,
    np.ndarray[np.uint64_t, ndim=1] pid,
    np.ndarray[np.uint64_t, ndim=1] selected_pid,
    bool preserve_particle_index ):
    """
    Go through the sorted arrays `pid` and `selected_pid`, and record
    the indices (of the array `pid`) where they match, by storing them
    in the array `selected_indices` (this array is thus modified in-place)

    Return the number of elements that were filled in `selected_indices`
    """
    cdef unsigned int i = 0
    cdef unsigned int i_select = 0
    cdef unsigned int i_fill = 0
    cdef unsigned int N = pid.shape[0]
    cdef unsigned int N_selected = selected_pid.shape[0]

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
