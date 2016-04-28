"""
This file is part of the openPMD-viewer.

It defines a set of methods which are useful for plotting
(and labeling the plots).

Copyright 2015-2016, openPMD-viewer contributors
Author: Remi Lehe
License: 3-Clause-BSD-LBNL
"""

import matplotlib.pyplot as plt


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
        self.fontsize = 18

        # Register the time array and iterations array
        # (Useful when labeling the figures)
        self.t = t
        self.iterations = iterations

    def hist1d(self, q1, w, quantity1, species, current_i, nbins,
               cmap='Blues', vmin=None, vmax=None, **kw):
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

        nbins : int, optional
           Number of bins for the histograms

        **kw : dict, otional
           Additional options to be passed to matplotlib's hist
        """
        # Find the iteration and time
        iteration = self.iterations[current_i]
        time_fs = 1.e15 * self.t[current_i]

        # Do the plot
        plt.hist(q1, bins=nbins, weights=w, **kw)
        plt.xlabel(quantity1, fontsize=self.fontsize)
        plt.title("%s:   t =  %.0f fs    (iteration %d)"
                  % (species, time_fs, iteration), fontsize=self.fontsize)

    def hist2d(self, q1, q2, w, quantity1, quantity2, species, current_i,
               nbins, cmap='Blues', vmin=None, vmax=None, **kw):
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

        nbins : int, optional
           Number of bins for the histograms

        **kw : dict, otional
           Additional options to be passed to matplotlib's hist
        """
        # Find the iteration and time
        iteration = self.iterations[current_i]
        time_fs = 1.e15 * self.t[current_i]

        # Do the plot
        plt.hist2d(q1, q2, bins=nbins, cmap=cmap,
                   vmin=vmin, vmax=vmax, weights=w, **kw)
        plt.colorbar()
        plt.xlabel(quantity1, fontsize=self.fontsize)
        plt.ylabel(quantity2, fontsize=self.fontsize)
        plt.title("%s:   t =  %.1f fs   (iteration %d)"
                  % (species, time_fs, iteration), fontsize=self.fontsize)

    def show_field(self, F, info, slicing_dir, m,
                   field_label, geometry, current_i, **kw):
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
        """
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
            plt.title("%s sliced across %s at %.1f fs  (iteration %d)"
                      % (field_label, slicing_dir, time_fs, iteration),
                      fontsize=self.fontsize)

        # Add the name of the axes
        plt.xlabel('$%s \;(\mu m)$' % info.axes[1], fontsize=self.fontsize)
        plt.ylabel('$%s \;(\mu m)$' % info.axes[0], fontsize=self.fontsize)

        # Plot the data
        plt.imshow(F, extent=1.e6 * info.imshow_extent, origin='lower',
                   interpolation='nearest', aspect='auto', **kw)
        plt.colorbar()
