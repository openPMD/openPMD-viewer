"""
This file is part of the OpenPMD viewer.

It defines a set of methods which are useful for plotting.
"""

def Plotter(object):
    """
    Class which is used for plotting particles and fields
    """

    def __init__( t, iterations ):
        """
        Initialize the object

        Parameters
        ----------
        t: 1darray of floats
           Time (in seconds) for each available iteration of the timeseries

        iterations: 1darray of ints
           Iteration number for each available iteration of the timeseries        
        """
        # Default fontsize
        self.fontsize = 18

        # Register the time array and iterations array
        # (Useful when labeling the figures)
        self.t = t
        self.iterations = iterations


    def hist1d( q1, w, quantity1, current_i, nbins, cmap='Blues',
                vmin=None, vmax=None, **kw ):
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

        current_i: int
            The index of this iteration, within the iterations list

        nbins : int, optional
           Number of bins for the histograms

        **kw : dict, otional
           Additional options to be passed to matplotlib's hist
        """
        # Find the iteration and time
        iteration = self.iterations[ current_i ]
        time_fs = 1.e15*self.t[ current_i ]

        # Do the plot
        plt.hist(q1, bins=nbins, weights=w, **kw )
        plt.xlabel(quantity1, fontsize=self.fontsize)
        plt.title("t =  %.0f fs    (iteration %d)" \
                %(time_fs, iteration), fontsize=self.fontsize )


    def hist2d( q1, q2, w, quantity1, quantity2, current_i,
                nbins, cmap='Blues', vmin=None, vmax=None, **kw ):
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

        quantity1: string
            The name of the quantity to be plotted (for labeling purposes)

        current_i: int
            The index of this iteration, within the iterations list

        nbins : int, optional
           Number of bins for the histograms

        **kw : dict, otional
           Additional options to be passed to matplotlib's hist
        """
        # Find the iteration and time
        iteration = self.iterations[ current_i ]
        time_fs = 1.e15*self.t[ current_i ]

        # Do the plot
        plt.hist2d(q1, q2, bins=nbins, cmap=cmap,
                    vmin=vmin, vmax=vmax, weights=w, **kw )
        plt.colorbar()
        plt.xlabel(quantity1, fontsize=self.fontsize)
        plt.ylabel(quantity2, fontsize=self.fontsize)
        plt.title("t =  %.1f fs   (iteration %d)"  \
                %(time_fs, iteration ), fontsize=self.fontsize )
                    
