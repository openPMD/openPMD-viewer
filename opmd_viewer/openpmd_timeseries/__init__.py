# Make the OpenPMDTimeSeries accessible from outside the file main
from .main import OpenPMDTimeSeries, ParticleTracker
from .data_reader.field_metainfo import FieldMetaInformation
__all__ = ['OpenPMDTimeSeries', 'FieldMetaInformation', 'ParticleTracker']
