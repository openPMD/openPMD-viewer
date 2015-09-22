# Class that inherits from OpenPMDTimeSeries, and implements
# some standard diagnostics (emittance, etc.)
# For the moment the class is empty
from openpmd_timeseries import OpenPMDTimeSeries

class LpaDiagnostics( OpenPMDTimeSeries ):
    
    def __init__(self, path_to_dir):
        """
        Initialize an OpenPMD time series

        Parameter
        ---------
        path_to_dir : string
            The path to the directory where the openPMD files are.
            For the moment, only HDF5 files are supported. There should be
            one file per iteration, and the name of the files should end
            with the iteration number, followed by '.h5' (e.g. data0005000.h5)
        """
        OpenPMDTimeSeries.__init__( path_to_dir )
