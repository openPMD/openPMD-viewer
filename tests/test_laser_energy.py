"""
This test file is part of the openPMD-viewer.

It checks that the function `get_electromagnetic_energy` works correctly

Usage:
This file is meant to be run from the root directory of openPMD-viewer,
by any of the following commands
$ python tests/test_laser_energy.py
$ py.test
$ python -m pytest tests

Copyright 2024, openPMD-viewer contributors
Authors: Remi Lehe
License: 3-Clause-BSD-LBNL
"""
from openpmd_viewer.addons import LpaDiagnostics
from scipy.constants import c
import numpy as np

# Download required datasets
import os
def download_if_absent( dataset_name ):
    "Function that downloads and decompress a chosen dataset"
    if os.path.exists( dataset_name ) is False:
        import wget, tarfile
        tar_name = "%s.tar.gz" %dataset_name
        url = "https://github.com/openPMD/openPMD-example-datasets/raw/draft/%s" %tar_name
        wget.download(url, tar_name)
        with tarfile.open( tar_name ) as tar_file:
            tar_file.extractall()
        os.remove( tar_name )
download_if_absent( 'example-thetaMode' )

def test_laser_energy():
    """
    Check that the function `get_electromagnetic_energy` gives the same result
    as the engineering formula for a Gaussian pulse

    E = 2.7e-5*(a0*w0/lambd)**2*(tau[fs])
    """
    ts = LpaDiagnostics('./example-thetaMode/hdf5')
    iteration = 300

    # Evaluate the laser energy using the engineering formula
    a0 = ts.get_a0(iteration=iteration, pol='y')
    w0 = ts.get_laser_waist(iteration=iteration, pol='y')
    tau_fs = ts.get_ctau(iteration=iteration, pol='y') / c * 1e15
    omega = ts.get_main_frequency(iteration=iteration, pol='y')
    lambd = 2*np.pi*c/omega
    E_eng = 2.7e-5*(a0*w0/lambd)**2*tau_fs

    # Evaluate the laser energy using the function
    E_func = ts.get_electromagnetic_energy(iteration=iteration)

    # Compare the two results
    assert np.isclose(E_eng, E_func, rtol=0.04)
