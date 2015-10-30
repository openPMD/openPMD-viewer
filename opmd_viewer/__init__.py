# Make the OpenPMDTimeSeries object accessible from outside the package
from .openpmd_timeseries import OpenPMDTimeSeries, FieldMetaInformation
__all__ = ['OpenPMDTimeSeries', 'FieldMetaInformation']
