from .data_reader import DataReader, available_backends

# Try to import mpiDataReader (only available if mpi4py is installed)
try:
    from .mpi_data_reader import mpiDataReader
    __all__ = ['DataReader', 'mpiDataReader', 'available_backends']
except ImportError:
    __all__ = ['DataReader', 'available_backends']
