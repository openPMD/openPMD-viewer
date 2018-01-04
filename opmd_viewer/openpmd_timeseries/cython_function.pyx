import numpy as np
cimport numpy as np
from cpython cimport bool
cimport cython
from libc.math cimport floor

@cython.boundscheck(False)
@cython.wraparound(False)
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

@cython.boundscheck(False)
@cython.wraparound(False)
def histogram_cic_1d(
    np.ndarray[np.float64_t, ndim=1] q1,
    np.ndarray[np.float64_t, ndim=1] w,
    int nbins, double bins_start, double bins_end ):
    """
    # TODO
    """
    # Define various scalars
    cdef double bin_spacing = (bins_end-bins_start)/nbins
    cdef double inv_spacing = 1./bin_spacing
    cdef int n_ptcl = len(q1)
    cdef int i_low_bin
    cdef double q1_cell
    cdef double S_low

    # Allocate array for histogrammed data
    hist_data = np.zeros( nbins, dtype=np.float64 )

    # Go through particle array and bin the data
    for i in xrange(n_ptcl):
        # Calculate the index of lower bin to which this particle contributes
        q1_cell = (q1[i] - bins_start) * inv_spacing
        i_low_bin = <int> floor( q1_cell )
        # Calculate corresponding CIC shape and deposit the weight
        S_low = 1. - (q1_cell - i_low_bin)
        if (i_low_bin >= 0) and (i_low_bin < nbins):
            hist_data[ i_low_bin ] += w[i] * S_low
        if (i_low_bin + 1 >= 0) and (i_low_bin + 1 < nbins):
            hist_data[ i_low_bin + 1 ] += w[i] * (1. - S_low)

    return( hist_data )
