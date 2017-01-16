"""
openPMD-viewer

Usage
-----
See the class OpenPMDTimeSeries to open a set of openPMD files
"""
# Make the OpenPMDTimeSeries object accessible from outside the package
from .openpmd_timeseries import OpenPMDTimeSeries, FieldMetaInformation, \
    ParticleTracker

# Define the version number
from .__version__ import __version__
__all__ = ['OpenPMDTimeSeries', 'FieldMetaInformation',
           'ParticleTracker', '__version__']
