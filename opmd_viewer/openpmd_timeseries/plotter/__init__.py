# This file imports the class `Plotter`, which is used
# by the class `OpenPMDTimeSeries` for plotting.
# The class `Plotter` should define the methods hist1d, hist2d, show_field_1d,
# and show_field_2d. If matplotlib is not installed, an exception will be
# raised when these methods are called.

try:
    # Check wether the plotter can be loaded
    from .plotter import Plotter
except ImportError:
    # Otherwise, create a dummy class which returns an error message
    # when the method `slider` is called.
    class Plotter(object):

        def __init__( self, t, iterations ):
            pass

        def hist1d( self, *args, **kwargs ):
            self.raise_error()

        def hist2d( self, *args, **kwargs ):
            self.raise_error()

        def show_field_1d( self, *args, **kwargs ):
            self.raise_error()

        def show_field_2d( self, *args, **kwargs ):
            self.raise_error()

        def raise_error( self ):
            raise RuntimeError(
                'Failed to load the openPMD-viewer plotter.\n'
                '(Make sure that matplotlib is installed.)\n')

__all__ = ['Plotter']
