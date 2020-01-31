"""
This file is part of the openPMD-viewer.

It defines a wrapper around numba.

Copyright 2019, openPMD-viewer contributors
Author: Remi Lehe
License: 3-Clause-BSD-LBNL
"""
import warnings

try:
    # Import jit decorator from numba
    import numba
    numba_installed = True
    jit = numba.njit(cache=True)

except ImportError:
    numba_installed = False
    # Create dummy decorator: warns about installing numba when calling
    # the decorated function.
    def jit(f):
        def decorated_f(*args, **kwargs):
            warnings.warn(
                '\nOne of the functions called by openPMD-viewer ' +\
                '(%s)\n' %f.__name__ +\
                'could have been faster if `numba` had been installed.\n' +\
                'Please consider installing `numba` (e.g. `pip install numba`)')
            return f(*args, **kwargs)
        return decorated_f
