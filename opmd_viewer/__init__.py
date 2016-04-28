"""
openPMD-viewer

Usage
-----
See the class OpenPMDTimeSeries to open a set of openPMD files
"""
# Make the OpenPMDTimeSeries object accessible from outside the package
from .openpmd_timeseries import OpenPMDTimeSeries, FieldMetaInformation
__all__ = ['OpenPMDTimeSeries', 'FieldMetaInformation']

# Package version number
__version__ = '0.2.0'
