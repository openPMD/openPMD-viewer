"""
This file is part of the openPMD-viewer.

It defines a set of methods which are useful for plotting
(and labeling the plots).

Copyright 2015-2016, openPMD-viewer contributors
Author: Remi Lehe
License: 3-Clause-BSD-LBNL
"""
import numpy as np
try:
    import warnings
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib_installed = True
except ImportError:
    matplotlib_installed = False
try:
    from .cython_function import histogram_cic_1d, histogram_cic_2d
    cython_function_available = True
except ImportError:
    cython_function_available = False


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

        hist_range : list of 2 floats
           Extent of the histogram

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
        time_fs = 1.e15 * self.t[current_i]

        # Check deposition method
        if deposition == 'cic' and not cython_function_available:
            print_cic_unavailable()
            deposition = 'ngp'

        # Bin the particle data
        q1 = q1.astype( np.float64 )
        if deposition == 'ngp':
            binned_data, _ = np.histogram(q1, nbins, hist_range, weights=w)
        elif deposition == 'cic':
            binned_data = histogram_cic_1d(
                q1, w, nbins, hist_range[0], hist_range[1])
        else:
            raise ValueError('Unknown deposition method: %s' % deposition)

        # Do the plot
        bin_size = (hist_range[1] - hist_range[0]) / nbins
        bin_coords = hist_range[0] + bin_size * ( 0.5 + np.arange(nbins) )
        plt.bar( bin_coords, binned_data, width=bin_size, **kw )
        plt.xlim( hist_range )
        plt.xlabel(quantity1, fontsize=self.fontsize)
        plt.title("%s:   t =  %.0f fs    (iteration %d)"
                  % (species, time_fs, iteration), fontsize=self.fontsize)

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
        time_fs = 1.e15 * self.t[current_i]

        # Check deposition method
        if deposition == 'cic' and not cython_function_available:
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
        plt.title("%s:   t =  %.1f fs   (iteration %d)"
                  % (species, time_fs, iteration), fontsize=self.fontsize)

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
        time_fs = 1.e15 * self.t[current_i]

        # Get the title and labels
        plt.title("%s at %.1f fs   (iteration %d)"
                % (field_label, time_fs, iteration), fontsize=self.fontsize)

        # Add the name of the axes
        plt.xlabel('$%s \;(\mu m)$' % info.axes[0], fontsize=self.fontsize)
        # Get the x axis in microns
        xaxis = 1.e6 * getattr( info, info.axes[0] )
        # Plot the data
        plt.plot( xaxis, F )
        # Get the limits of the plot
        # - Along the first dimension
        if (plot_range[0][0] is not None) and (plot_range[0][1] is not None):
            plt.xlim( plot_range[0][0], plot_range[0][1] )
        else:
            plt.xlim( xaxis.min(), xaxis.max() )  # Full extent of the box
        # - Along the second dimension
        if (plot_range[1][0] is not None) and (plot_range[1][1] is not None):
            plt.ylim( plot_range[1][0], plot_range[1][1] )

    def show_field_2d(self, F, info, slicing_dir, m, field_label, geometry,
                        current_i, plot_range, **kw):
        """
        Plot the given field in 2D

        Parameters
        ----------
        F: 2darray of floats
            Contains the field to be plotted

        info: a FieldMetaInformation object
            Contains the information about the plotted field

        slicing_dir : str, optional
           Only used for 3dcartesian geometry
           The direction along which the data is sliced

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
        time_fs = 1.e15 * self.t[current_i]

        # Get the title and labels
        # Cylindrical geometry
        if geometry == "thetaMode":
            mode = str(m)
            plt.title("%s in the mode %s at %.1f fs   (iteration %d)"
                      % (field_label, mode, time_fs, iteration),
                      fontsize=self.fontsize)
        # 2D Cartesian geometry
        elif geometry == "2dcartesian":
            plt.title("%s at %.1f fs   (iteration %d)"
                      % (field_label, time_fs, iteration),
                      fontsize=self.fontsize)
        # 3D Cartesian geometry
        elif geometry == "3dcartesian":
            slice_plane = info.axes[0] + '-' + info.axes[1]
            plt.title("%s sliced in %s at %.1f fs  (iteration %d)"
                      % (field_label, slice_plane, time_fs, iteration),
                      fontsize=self.fontsize)

        # Add the name of the axes
        plt.xlabel('$%s \;(\mu m)$' % info.axes[1], fontsize=self.fontsize)
        plt.ylabel('$%s \;(\mu m)$' % info.axes[0], fontsize=self.fontsize)

        # Plot the data
        plt.imshow(F, extent=1.e6 * info.imshow_extent, origin='lower',
                   interpolation='nearest', aspect='auto', **kw)
        plt.colorbar()

        # Get the limits of the plot
        # - Along the first dimension
        if (plot_range[0][0] is not None) and (plot_range[0][1] is not None):
            plt.xlim( plot_range[0][0], plot_range[0][1] )
        # - Along the second dimension
        if (plot_range[1][0] is not None) and (plot_range[1][1] is not None):
            plt.ylim( plot_range[1][0], plot_range[1][1] )


def print_cic_unavailable():
    warnings.warn(
        "\nCIC particle histogramming is unavailable because \n"
        "Cython is not installed. NGP histogramming is used instead.\n"
        "For CIC histogramming: \n"
        " - make sure that Cython is installed \n"
        " - then reinstall openPMD-viewer")


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
