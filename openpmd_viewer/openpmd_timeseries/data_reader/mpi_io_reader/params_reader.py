"""
MPI-enabled parallel versions of io_reader params_reader functions.

Copyright 2020, openPMD-viewer contributors
License: 3-Clause-BSD-LBNL
"""
from ..io_reader.params_reader import read_openPMD_params as _read_openPMD_params


def read_openPMD_params(series, iteration, comm, extract_parameters=True):
    """
    Extract the time and some openPMD parameters from a file.
    
    Only rank 0 performs the parameter reading operation, and broadcasts
    the results to all other ranks.

    Parameter
    ---------
    series: openpmd_api.Series
        An open, readable openPMD-api series object

    iteration: integer
        Iteration from which parameters should be extracted
    
    comm : MPI communicator
        MPI communicator for parallel operations

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
    rank = comm.Get_rank()
    
    if rank == 0:
        time, params = _read_openPMD_params(series, iteration, extract_parameters)
    else:
        time = None
        params = None
    
    # Broadcast time and params from rank 0 to all ranks
    time = comm.bcast(time, root=0)
    params = comm.bcast(params, root=0)
    
    return time, params
