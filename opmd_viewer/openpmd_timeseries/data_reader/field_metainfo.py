"""
This file is part of the openPMD viewer.

It defines the main FieldMetaInformation class, which
is returned by `get_field` along with the array of field values,
and gathers information collected from the openPMD file.

__authors__ = "Remi Lehe"
__copyright__ = "Copyright 2015-2016, openPMD viewer contributors"
__license__ = "3-Clause-BSD-LBNL"
"""

import numpy as np


class FieldMetaInformation(object):
    """
    An object that is typically returned along with an array of field
    values, and which contains meta-information about the grid.

    Attributes
    ----------
    - axes: dict
        A dictionary of the form {0:'x', 1:'z'}, which indicates the name
        of the coordinate along each axis of the field array.
        For instance, in the case of {0:'x', 1:'y'}, the first axis of field
        array corresponds to the 'x' coordinate, and the second to 'y'.

    - xmin, xmax, zmin, zmax: double
        Scalars that indicate the position of the first grid point and
        last grid point along each axis.
        Notice that the name of these variables change according to
        the values in `axes`. For instance, if `axes` is {0: 'x', 1: 'y'},
        then these variables will be called xmin, xmax, ymin, ymax.

     - x, z: 1darrays of double
        The position of all the gridpoints, along each axis
        Notice that the name of these variables change according to
        the values in `axes`. For instance, if `axes` is {0: 'x', 1: 'y'},
        then these variables will be called x, y.

    - imshow_extent: 1darray
        An array of 4 elements that can be passed as the `extent` in
        matplotlib's imshow function.
        Because of the API of the imshow function, the coordinates are
        'swapped' inside imshow_extent. For instance, if axes is
        {0: 'x', 1: 'y'}, then imshow_extent will be [ymin, ymax, xmin, xmax].

        (NB: in the details, imshow_extent contains slightly different values
        than ymin, ymax, xmin, xmax: these values are shifted by half a cell.
        The reason for this is that imshow plots a finite-width square for each
        value of the field array.)
    """

    def __init__(self, axes, shape, grid_spacing,
                 global_offset, grid_unitSI, position, thetaMode=False):
        """
        Create a FieldMetaInformation object

        The input arguments correspond to their openPMD standard definition
        """
        # Register important initial information
        self.axes = axes
        self.imshow_extent = []

        # Create the elements
        for axis in sorted(axes.keys()):
            # Create the coordinates along this axis
            step = grid_spacing[axis] * grid_unitSI
            n_points = shape[axis]
            start = global_offset[axis] * grid_unitSI + position[axis] * step
            end = start + (n_points - 1) * step
            axis_points = np.linspace(start, end, n_points, endpoint=True)
            # Register the results in the object
            axis_name = axes[axis]
            setattr(self, axis_name, axis_points)
            setattr(self, axis_name + 'min', axis_points[0])
            setattr(self, axis_name + 'max', axis_points[-1])
            # Fill the imshow_extent in reverse order, so as to match
            # the syntax of imshow ; add a half step on each side since
            # imshow plots a square of finite width for each field value
            self.imshow_extent = \
                [start - 0.5 * step, end + 0.5 * step] + self.imshow_extent

        # Create the points below the axis if thetaMode is true
        if thetaMode:
            self.r = np.concatenate((-self.r[::-1], self.r))
            # The axis now extends from -rmax to rmax
            self.rmin = -self.rmax
            self.imshow_extent[2] = -self.imshow_extent[3]

        # Finalize imshow_extent by converting it from list to array
        self.imshow_extent = np.array(self.imshow_extent)

    def restrict_to_1Daxis(self, axis):
        """
        Suppresses the information that correspond to other axes than `axis`

        Parameters
        ----------
        axis: string
            The axis to keep
            This has to be one of the keys of the self.axes dictionary
        """
        # Check if axis is a valid key
        if (axis in self.axes.values()) is False:
            raise ValueError('`axis` is not one of the coordinates '
                             'that are present in this object.')

        # Loop through the coordinates and suppress them
        for obsolete_axis in self.axes.values():
            if obsolete_axis != axis:
                delattr(self, obsolete_axis)
                delattr(self, obsolete_axis + 'min')
                delattr(self, obsolete_axis + 'max')

        # Suppress imshow_extent and replace the dictionary
        delattr(self, 'imshow_extent')
        self.axes = {0: axis}
