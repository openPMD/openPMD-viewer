"""
This file is part of the openPMD-viewer.

It defines a set of methods which are useful for plotting
(and labeling the plots).

Copyright 2015-2016, openPMD-viewer contributors
Author: Remi Lehe
License: 3-Clause-BSD-LBNL
"""
import numpy as np
import math
try:
    import warnings
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib_installed = True
except ImportError:
    matplotlib_installed = False

from .numba_wrapper import numba_installed
if numba_installed:
    from .utilities import histogram_cic_1d, histogram_cic_2d

# Redefine the default matplotlib formatter for ticks
if matplotlib_installed:
    from matplotlib.ticker import ScalarFormatter

    class PowerOfThreeFormatter( ScalarFormatter ):
        """
        Formatter for matplotlib's axes ticks,
        that prints numbers as e.g. 1.5e3, 3.2e6, 0.2e-9,
        where the exponent is always a multiple of 3.

        This helps a human reader to quickly identify the closest units
        (e.g. nanometer) of the plotted quantity.

        This class derives from `ScalarFormatter`, which
        provides a nice `offset` feature.
        """
        def __init__( self, *args, **kwargs ):
            ScalarFormatter.__init__( self, *args, **kwargs )
            # Do not print the order of magnitude on the side of the axis
            self.set_scientific(False)
            # Reduce the threshold for printing an offset on side of the axis
            self._offset_threshold = 2

        def __call__(self, x, pos=None):
            """
            Function called for each tick of an axis (for matplotlib>=3.1)
            Returns the string that appears in the plot.
            """
            return self.pprint_val( x, pos )

        def pprint_val( self, x, pos=None):
            """
            Function called for each tick of an axis (for matplotlib<3.1)
            Returns the string that appears in the plot.
            """
            # Calculate the exponent (power of 3)
            xp = (x - self.offset)
            if xp != 0:
                exponent = int(3 * math.floor( math.log10(abs(xp)) / 3 ))
            else:
                exponent = 0
            # Show 3 digits at most after decimal point
            mantissa = round( xp * 10**(-exponent), 3)
            # After rounding the exponent might change (e.g. 0.999 -> 1.)
            if mantissa != 0 and math.log10(abs(mantissa)) == 3:
                exponent += 3
                mantissa /= 1000
            string = "{:.3f}".format( mantissa )
            if '.' in string:
                # Remove trailing zeros and ., for integer mantissa
                string = string.rstrip('0')
                string = string.rstrip('.')
            if exponent != 0:
                string += "e{:d}".format( exponent )
            return string

    tick_formatter = PowerOfThreeFormatter()


class Plotter(object):

    """
    Class which is used for plotting particles and fields
    (and labeling the plots)
    """

    def __init__(self, t, iterations):
        """
        Initialize the object

        Parameters
        ----------
        t: 1darray of floats (seconds)
           Time for each available iteration of the timeseries

        iterations: 1darray of ints
           Iteration number for each available iteration of the timeseries
        """
        # Default fontsize
        self.fontsize = 12

        # Register the time array and iterations array
        # (Useful when labeling the figures)
        self.t = t
        self.iterations = iterations

    def hist1d(self, q1, w, quantity1, species, current_i, nbins, hist_range,
               cmap='Blues', vmin=None, vmax=None, deposition='cic', **kw):
        """
        Plot a 1D histogram of the particle quantity q1
        Sets the proper labels

        Parameters
        ----------
        q1: 1darray of floats
            An array with one element per macroparticle, representing
            the quantity to be plotted.

        w: 1darray of floats
            An array with one element per macroparticle, representing
            the number of real particles that correspond to each macroparticle

        quantity1: string
            The name of the quantity to be plotted (for labeling purposes)

        species: string
            The name of the species from which the data is taken

        current_i: int
            The index of this iteration, within the iterations list

        nbins : int
           Number of bins for the histograms

        hist_range : list contains 2 lists of 2 floats
           Extent of the histogram along each direction

        deposition : string
            Either `ngp` (Nearest Grid Point) or `cic` (Cloud-In-Cell)
            When plotting the particle histogram, this determines how
            particles affects neighboring bins.
            `cic` (which is the default) leads to smoother results than `ngp`.

        **kw : dict, otional
           Additional options to be passed to matplotlib's bar function
        """
        # Check if matplotlib is available
        check_matplotlib()

        # Find the iteration and time
        iteration = self.iterations[current_i]
        time = self.t[current_i]

        # Check deposition method
        if deposition == 'cic' and not numba_installed:
            print_cic_unavailable()
            deposition = 'ngp'

        # Bin the particle data
        q1 = q1.astype( np.float64 )
        if deposition == 'ngp':
            binned_data, _ = np.histogram(q1, nbins, hist_range[0], weights=w)
        elif deposition == 'cic':
            binned_data = histogram_cic_1d(
                q1, w, nbins, hist_range[0][0], hist_range[0][1])
        else:
            raise ValueError('Unknown deposition method: %s' % deposition)

        # Do the plot
        bin_size = (hist_range[0][1] - hist_range[0][0]) / nbins
        bin_coords = hist_range[0][0] + bin_size * ( 0.5 + np.arange(nbins) )
        plt.bar( bin_coords, binned_data, width=bin_size, **kw )
        plt.xlim( hist_range[0] )
        plt.ylim( hist_range[1] )
        plt.xlabel(quantity1, fontsize=self.fontsize)
        plt.title("%s:   t =  %.2e s    (iteration %d)"
                  % (species, time, iteration), fontsize=self.fontsize)
        # Format the ticks
        ax = plt.gca()
        ax.get_xaxis().set_major_formatter( tick_formatter )
        ax.get_yaxis().set_major_formatter( tick_formatter )

    def hist2d(self, q1, q2, w, quantity1, quantity2, species, current_i,
                nbins, hist_range, cmap='Blues', vmin=None, vmax=None,
                deposition='cic', **kw):
        """
        Plot a 2D histogram of the particle quantity q1
        Sets the proper labels

        Parameters
        ----------
        q1: 1darray of floats
            An array with one element per macroparticle, representing
            the quantity to be plotted.

        w: 1darray of floats
            An array with one element per macroparticle, representing
            the number of real particles that correspond to each macroparticle

        quantity1, quantity2: strings
            The name of the quantity to be plotted (for labeling purposes)

        species: string
            The name of the species from which the data is taken

        current_i: int
            The index of this iteration, within the iterations list

        nbins : list of 2 ints
           Number of bins along each direction, for the histograms

        hist_range : list contains 2 lists of 2 floats
           Extent of the histogram along each direction

        deposition : string
            Either `ngp` (Nearest Grid Point) or `cic` (Cloud-In-Cell)
            When plotting the particle histogram, this determines how
            particles affects neighboring bins.
            `cic` (which is the default) leads to smoother results than `ngp`.

        **kw : dict, otional
           Additional options to be passed to matplotlib's imshow function
        """
        # Check if matplotlib is available
        check_matplotlib()

        # Find the iteration and time
        iteration = self.iterations[current_i]
        time = self.t[current_i]

        # Check deposition method
        if deposition == 'cic' and not numba_installed:
            print_cic_unavailable()
            deposition = 'ngp'

        # Bin the particle data
        q1 = q1.astype( np.float64 )
        q2 = q2.astype( np.float64 )
        if deposition == 'ngp':
            binned_data, _, _ = np.histogram2d(
                q1, q2, nbins, hist_range, weights=w)
        elif deposition == 'cic':
            binned_data = histogram_cic_2d( q1, q2, w,
                nbins[0], hist_range[0][0], hist_range[0][1],
                nbins[1], hist_range[1][0], hist_range[1][1] )
        else:
            raise ValueError('Unknown deposition method: %s' % deposition)

        # Do the plot
        plt.imshow( binned_data.T, extent=hist_range[0] + hist_range[1],
             origin='lower', interpolation='nearest', aspect='auto',
             cmap=cmap, vmin=vmin, vmax=vmax, **kw )
        plt.colorbar()
        plt.xlabel(quantity1, fontsize=self.fontsize)
        plt.ylabel(quantity2, fontsize=self.fontsize)
        plt.title("%s:   t =  %.2e s   (iteration %d)"
                  % (species, time, iteration), fontsize=self.fontsize)
        # Format the ticks
        ax = plt.gca()
        ax.get_xaxis().set_major_formatter( tick_formatter )
        ax.get_yaxis().set_major_formatter( tick_formatter )

    def show_field_1d( self, F, info, field_label, current_i, plot_range,
                            vmin=None, vmax=None, **kw ):
        """
        Plot the given field in 1D

        Parameters
        ----------
        F: 1darray of floats
            Contains the field to be plotted

        info: a FieldMetaInformation object
            Contains the information about the plotted field

        field_label: string
           The name of the field plotted (for labeling purposes)

        vmin, vmax: floats or None
           The amplitude of the field

        plot_range : list of lists
           Indicates the values between which to clip the plot,
           along the 1st axis (first list) and 2nd axis (second list)
        """
        # Check if matplotlib is available
        check_matplotlib()

        # Find the iteration and time
        iteration = self.iterations[current_i]
        time = self.t[current_i]

        # Get the x axis
        xaxis = getattr( info, info.axes[0] )
        # Plot the data
        if np.issubdtype(F.dtype, np.complexfloating):
            plot_data = abs(F) # For complex numbers, plot the absolute value
            title = "|%s|" %field_label
        else:
            plot_data = F
            title = "%s" %field_label

        # Get the title and labels
        title += " at %.2e s   (iteration %d)" % (time, iteration)
        plt.title(title, fontsize=self.fontsize)
        # Add the name of the axes
        plt.xlabel('$%s \;(m)$' % info.axes[0], fontsize=self.fontsize)

        plt.plot( xaxis, plot_data )
        # Get the limits of the plot
        # - Along the first dimension
        if (plot_range[0][0] is not None) and (plot_range[0][1] is not None):
            plt.xlim( plot_range[0][0], plot_range[0][1] )
        else:
            plt.xlim( xaxis.min(), xaxis.max() )  # Full extent of the box
        # - Along the second dimension
        if (plot_range[1][0] is not None) and (plot_range[1][1] is not None):
            plt.ylim( plot_range[1][0], plot_range[1][1] )
        # Format the ticks
        ax = plt.gca()
        ax.get_xaxis().set_major_formatter( tick_formatter )
        ax.get_yaxis().set_major_formatter( tick_formatter )

    def show_field_2d(self, F, info, slice_across, m, field_label, geometry,
                        current_i, plot_range, **kw):
        """
        Plot the given field in 2D

        Parameters
        ----------
        F: 2darray of floats
            Contains the field to be plotted

        info: a FieldMetaInformation object
            Contains the information about the plotted field

        slice_across : str, optional
           Only used for 3dcartesian geometry
           The direction across which the data is sliced

        m: int
           Only used for thetaMode geometry
           The azimuthal mode used when plotting the fields

        field_label: string
           The name of the field plotted (for labeling purposes)

        geometry: string
           Either "2dcartesian", "3dcartesian" or "thetaMode"

        plot_range : list of lists
           Indicates the values between which to clip the plot,
           along the 1st axis (first list) and 2nd axis (second list)
        """
        # Check if matplotlib is available
        check_matplotlib()

        # Find the iteration and time
        iteration = self.iterations[current_i]
        time = self.t[current_i]

        # Plot the data
        if np.issubdtype(F.dtype, np.complexfloating):
            plot_data = abs(F)
            title = "|%s|" %field_label
        else:
            plot_data = F
            title = "%s" %field_label
        plt.imshow(plot_data, extent=info.imshow_extent, origin='lower',
                   interpolation='nearest', aspect='auto', **kw)
        plt.colorbar()

        # Get the title and labels
        # Cylindrical geometry
        if geometry == "thetaMode":
            mode = str(m)
            title += " in the mode %s at %.2e s   (iteration %d)" \
                      % (mode, time, iteration)
        # 2D Cartesian geometry
        else:
            title += " at %.2e s   (iteration %d)" % (time, iteration)
        plt.title(title, fontsize=self.fontsize)

        # Add the name of the axes
        plt.xlabel('$%s \;(m)$' % info.axes[1], fontsize=self.fontsize)
        plt.ylabel('$%s \;(m)$' % info.axes[0], fontsize=self.fontsize)

        # Get the limits of the plot
        # - Along the first dimension
        if (plot_range[0][0] is not None) and (plot_range[0][1] is not None):
            plt.xlim( plot_range[0][0], plot_range[0][1] )
        # - Along the second dimension
        if (plot_range[1][0] is not None) and (plot_range[1][1] is not None):
            plt.ylim( plot_range[1][0], plot_range[1][1] )
        # Format the ticks
        ax = plt.gca()
        ax.get_xaxis().set_major_formatter( tick_formatter )
        ax.get_yaxis().set_major_formatter( tick_formatter )


def print_cic_unavailable():
    warnings.warn(
        "\nCIC particle histogramming is unavailable because \n"
        "Numba is not installed. NGP histogramming is used instead.\n"
        "Please considering installing numba (e.g. `pip install numba`)")


def check_matplotlib():
    """Raise error messages or warnings when potential issues when
    potenial issues with matplotlib are detected."""

    if not matplotlib_installed:
        raise RuntimeError( "Failed to import the openPMD-viewer plotter.\n"
            "(Make sure that matplotlib is installed.)")

    elif ('MacOSX' in matplotlib.get_backend()):
        warnings.warn("\n\nIt seems that you are using the matplotlib MacOSX "
        "backend. \n(This typically obtained when typing `%matplotlib`.)\n"
        "With recent version of Jupyter, the plots might not appear.\nIn this "
        "case, switch to `%matplotlib notebook` and restart the notebook.")
