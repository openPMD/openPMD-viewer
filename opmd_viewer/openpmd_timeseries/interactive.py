"""
This file is part of the openPMD viewer.

It defines an interactive interface for the viewer,
based on the IPython notebook functionalities

__authors__ = "Remi Lehe, Axel Huebl"
__copyright__ = "Copyright 2015-2016, openPMD viewer contributors"
__license__ = "3-Clause-BSD-LBNL"
"""

try:
    from ipywidgets import widgets
except ImportError:
    # If ipywidgets is not available, use the deprecated package
    # IPython.html.widgets, so that the GUI still works
    from IPython.html import widgets

from IPython.display import display, clear_output
import math
import matplotlib
import matplotlib.pyplot as plt
from functools import partial


class InteractiveViewer(object):

    def __init__(self):
        pass

    def slider(self, figsize=(8, 8),
               exclude_particle_records=['charge', 'mass'], **kw):
        """
        Navigate the simulation using a slider

        Parameters:
        -----------
        figsize: tuple
            Size of the figures

        exclude_particle_records: list of strings
            List of particle quantities that should not be displayed
            in the slider (typically because they are less interesting)

        kw: dict
            Extra arguments to pass to matplotlib's imshow
        """
        # -----------------------
        # Define useful functions
        # -----------------------

        def refresh_field(change=None, force=False):
            """
            Refresh the current field figure

            Parameters :
            ------------
            change: dictionary
                Dictionary passed by the widget to a callback functions
                whenever a change of a widget happens
                (see docstring of ipywidgets.Widget.observe)
                This is mainline a place holder ; not used in this function

            force: bool
                Whether to force the update
            """
            # Determine whether to do the refresh
            do_refresh = False
            if (self.avail_fields is not None):
                if force or fld_refresh_toggle.value:
                    do_refresh = True
            # Do the refresh
            if do_refresh:
                plt.figure(fld_figure_button.value, figsize=figsize)
                plt.clf()

                # When working in inline mode, in an ipython notebook,
                # clear the output (prevents the images from stacking
                # in the notebook)
                if 'inline' in matplotlib.get_backend():
                    clear_output()

                if fld_use_button.value:
                    i_power = fld_magnitude_button.value
                    vmin = fld_range_button.value[0] * 10 ** i_power
                    vmax = fld_range_button.value[1] * 10 ** i_power
                else:
                    vmin = None
                    vmax = None

                self.get_field(t=self.current_t, output=False, plot=True,
                    field=fieldtype_button.value, coord=coord_button.value,
                    m=convert_to_int(mode_button.value),
                    slicing=slicing_button.value, theta=theta_button.value,
                    slicing_dir=slicing_dir_button.value,
                    vmin=vmin, vmax=vmax, cmap=fld_color_button.value)

        def refresh_ptcl(change=None, force=False):
            """
            Refresh the current particle figure

            Parameters :
            ------------
            change: dictionary
                Dictionary passed by the widget to a callback functions
                whenever a change of a widget happens
                (see docstring of ipywidgets.Widget.observe)
                This is mainline a place holder ; not used in this function

            force: bool
                Whether to force the update
            """
            # Determine whether to do the refresh
            do_refresh = False
            if self.avail_species is not None:
                if force or ptcl_refresh_toggle.value:
                    do_refresh = True
            # Do the refresh
            if do_refresh:
                plt.figure(ptcl_figure_button.value, figsize=figsize)
                plt.clf()

                # When working in inline mode, in an ipython notebook,
                # clear the output (prevents the images from stacking
                # in the notebook)
                if 'inline' in matplotlib.get_backend():
                    clear_output()

                if ptcl_use_button.value:
                    i_power = ptcl_magnitude_button.value
                    vmin = ptcl_range_button.value[0] * 10 ** i_power
                    vmax = ptcl_range_button.value[1] * 10 ** i_power
                else:
                    vmin = None
                    vmax = None

                if ptcl_yaxis_button.value == 'None':
                    # 1D histogram
                    self.get_particle(t=self.current_t, output=False,
                        var_list=[ptcl_xaxis_button.value],
                        select=ptcl_select_widget.to_dict(),
                        species=ptcl_species_button.value, plot=True,
                        vmin=vmin, vmax=vmax, cmap=ptcl_color_button.value,
                        nbins=ptcl_bins_button.value)
                else:
                    # 2D histogram
                    self.get_particle(t=self.current_t, output=False,
                        var_list=[ptcl_xaxis_button.value,
                                  ptcl_yaxis_button.value],
                        select=ptcl_select_widget.to_dict(),
                        species=ptcl_species_button.value, plot=True,
                        vmin=vmin, vmax=vmax, cmap=ptcl_color_button.value,
                        nbins=ptcl_bins_button.value)

        def refresh_field_type(change):
            """
            Refresh the field type and disable the coordinates buttons
            if the field is scalar.

            Parameter
            ---------
            change: dictionary
                Dictionary passed by the widget to a callback functions
                whenever a change of a widget happens
                (see docstring of ipywidgets.Widget.observe)
            """
            if self.avail_fields[change['new']] == 'scalar':
                coord_button.disabled = True
            elif self.avail_fields[change['new']] == 'vector':
                coord_button.disabled = False
            refresh_field()

        def refresh_species(change=None):
            """
            Refresh the particle species buttons by populating them
            with the available records for the current species

            Parameter
            ---------
            change: dictionary
                Dictionary passed by the widget to a callback functions
                whenever a change of a widget happens
                (see docstring of ipywidgets.Widget.observe)
            """
            # Deactivate the particle refreshing to avoid callback
            # while modifying the widgets
            saved_refresh_value = ptcl_refresh_toggle.value
            ptcl_refresh_toggle.value = False

            # Get available records for this species
            avail_records = [q for q in
                             self.avail_record_components[
                                 ptcl_species_button.value]
                             if q not in exclude_particle_records]
            # Update the plotting buttons
            ptcl_xaxis_button.options = avail_records
            ptcl_xaxis_button.value = avail_records[0]
            ptcl_yaxis_button.options = avail_records + ['None']
            ptcl_yaxis_button.value = 'None'

            # Update the selection widgets
            for dropdown_button in ptcl_select_widget.quantity:
                dropdown_button.options = avail_records

            # Put back the previous value of the refreshing button
            ptcl_refresh_toggle.value = saved_refresh_value

        def change_t(change):
            "Plot the result at the required time"
            self.current_t = 1.e-15 * change['new']
            refresh_field()
            refresh_ptcl()

        def step_fw(b):
            "Plot the result one iteration further"
            if self.current_i < len(self.t) - 1:
                self.current_t = self.t[self.current_i + 1]
            else:
                self.current_t = self.t[self.current_i]
            slider.value = self.current_t * 1.e15

        def step_bw(b):
            "Plot the result one iteration before"
            if self.current_t > 0:
                self.current_t = self.t[self.current_i - 1]
            else:
                self.current_t = self.t[self.current_i]
            slider.value = self.current_t * 1.e15

        # ---------------
        # Define widgets
        # ---------------

        # Slider
        slider = widgets.FloatSlider(
            min=math.ceil(1.e15 * self.tmin),
            max=math.ceil(1.e15 * self.tmax),
            step=math.ceil(1.e15 * (self.tmax - self.tmin)) / 20.,
            description="t (fs)")
        slider.observe( change_t, names='value', type='change')

        # Forward button
        button_p = widgets.Button(description="+")
        button_p.on_click(step_fw)

        # Backward button
        button_m = widgets.Button(description="-")
        button_m.on_click(step_bw)

        # Display the time widgets
        container = widgets.HBox(children=[button_m, button_p, slider])
        display(container)

        # Field widgets
        # -------------
        if (self.avail_fields is not None):

            # Field type
            # ----------
            # Field button
            fieldtype_button = widgets.ToggleButtons(
                description='Field:',
                options=sorted(self.avail_fields.keys()))
            fieldtype_button.observe( refresh_field_type, 'value', 'change' )

            # Coord button
            if self.geometry == "thetaMode":
                coord_button = widgets.ToggleButtons(
                    description='Coord:', options=['x', 'y', 'z', 'r', 't'])
            elif self.geometry in ["2dcartesian", "3dcartesian"]:
                coord_button = widgets.ToggleButtons(
                    description='Coord:', options=['x', 'y', 'z'])
            coord_button.observe( refresh_field, 'value', 'change')
            # Mode and theta button (for thetaMode)
            mode_button = widgets.ToggleButtons(description='Mode:',
                                                options=self.avail_circ_modes)
            mode_button.observe( refresh_field, 'value', 'change')
            theta_button = widgets.FloatSlider(width=140, value=0.,
                description=r'Theta:', min=-math.pi / 2, max=math.pi / 2)
            theta_button.observe( refresh_field, 'value', 'change')
            # Slicing buttons (for 3D)
            slicing_dir_button = widgets.ToggleButtons(value='y',
                description='Slicing direction:', options=['x', 'y', 'z'])
            slicing_dir_button.observe( refresh_field, 'value', 'change')
            slicing_button = widgets.FloatSlider(width=150,
                description='Slicing:', min=-1., max=1., value=0.)
            slicing_button.observe( refresh_field, 'value', 'change')

            # Plotting options
            # ----------------
            # Figure number
            fld_figure_button = widgets.IntText(description='Figure ',
                                                value=0, width=50)
            # Range of values
            fld_range_button = widgets.FloatRangeSlider(
                min=-10, max=10, width=220)
            fld_range_button.observe( refresh_field, 'value', 'change')
            # Order of magnitude
            fld_magnitude_button = widgets.IntText(
                description='x 10^', value=9, width=50)
            fld_magnitude_button.observe( refresh_field, 'value', 'change')
            # Use button
            fld_use_button = widgets.Checkbox(
                description=' Use this range', value=False)
            fld_use_button.observe( refresh_field, 'value', 'change')
            # Colormap button
            fld_color_button = widgets.Select(
                options=sorted(plt.cm.datad.keys()), height=50, width=200,
                value='jet')
            fld_color_button.observe( refresh_field, 'value', 'change' )
            # Resfresh buttons
            fld_refresh_toggle = widgets.ToggleButton(
                description='Always refresh', value=True)
            fld_refresh_button = widgets.Button(
                description='Refresh now!')
            fld_refresh_button.on_click( partial(refresh_field, force=True) )

            # Containers
            # ----------
            # Field type container
            if self.geometry == "thetaMode":
                container_fields = widgets.VBox(width=260,
                    children=[
                        fieldtype_button, coord_button,
                        mode_button, theta_button])
            elif self.geometry == "2dcartesian":
                container_fields = widgets.VBox(width=260,
                    children=[fieldtype_button, coord_button])
            elif self.geometry == "3dcartesian":
                container_fields = widgets.VBox(width=260,
                    children=[
                        fieldtype_button, coord_button,
                        slicing_dir_button, slicing_button])
            # Plotting options container
            container_fld_plots = widgets.VBox(width=260,
                children=[
                    fld_figure_button, fld_range_button,
                    widgets.HBox(
                        children=[
                            fld_magnitude_button,
                            fld_use_button],
                        height=50),
                    fld_color_button])
            # Accordion for the field widgets
            accord1 = widgets.Accordion(
                children=[container_fields, container_fld_plots])
            accord1.set_title(0, 'Field type')
            accord1.set_title(1, 'Plotting options')
            # Complete field container
            container_fld = widgets.VBox(width=300,
                children=[accord1, widgets.HBox(
                    children=[fld_refresh_toggle, fld_refresh_button])])

        # Particle widgets
        # ----------------
        if (self.avail_species is not None):

            # Particle quantities
            # -------------------
            # Species selection
            ptcl_species_button = widgets.Dropdown( width=250,
                                                   options=self.avail_species )
            ptcl_species_button.observe( refresh_species, 'value', 'change')
            # Get available records for this species
            avail_records = [q for q in
                             self.avail_record_components[
                                 ptcl_species_button.value]
                             if q not in exclude_particle_records]
            # Particle quantity on the x axis
            ptcl_xaxis_button = widgets.ToggleButtons(options=avail_records)
            ptcl_xaxis_button.observe( refresh_ptcl, 'value', 'change')
            # Particle quantity on the y axis
            ptcl_yaxis_button = widgets.ToggleButtons(
                options=avail_records + ['None'], value='None')
            ptcl_yaxis_button.observe( refresh_ptcl, 'value', 'change')

            # Particle selection
            # ------------------
            # 3 selection rules at maximum
            ptcl_select_widget = ParticleSelectWidget(3,
                                 avail_records, refresh_ptcl)

            # Plotting options
            # ----------------
            # Figure number
            ptcl_figure_button = widgets.IntText(description='Figure ',
                                                 value=1, width=50)
            # Number of bins
            ptcl_bins_button = widgets.IntText(description='nbins:',
                                               value=100, width=100)
            ptcl_bins_button.observe( refresh_ptcl, 'value', 'change')
            # Colormap button
            ptcl_color_button = widgets.Select(
                options=sorted(plt.cm.datad.keys()), height=50, width=200,
                value='Blues')
            ptcl_color_button.observe( refresh_ptcl, 'value', 'change')
            # Range of values
            ptcl_range_button = widgets.FloatRangeSlider(
                min=0, max=10, width=220, value=(0, 5))
            ptcl_range_button.observe( refresh_ptcl, 'value', 'change')
            # Order of magnitude
            ptcl_magnitude_button = widgets.IntText(
                description='x 10^', value=9, width=50)
            ptcl_magnitude_button.observe( refresh_ptcl, 'value', 'change')
            # Use button
            ptcl_use_button = widgets.Checkbox(
                description=' Use this range', value=False)
            ptcl_use_button.observe( refresh_ptcl, 'value', 'change')
            # Resfresh buttons
            ptcl_refresh_toggle = widgets.ToggleButton(
                description='Always refresh', value=True)
            ptcl_refresh_button = widgets.Button(
                description='Refresh now!')
            ptcl_refresh_button.on_click( partial(refresh_ptcl, force=True) )

            # Containers
            # ----------
            # Particle quantity container
            container_ptcl_quantities = widgets.VBox(width=310,
                children=[
                    ptcl_species_button, ptcl_xaxis_button,
                    ptcl_yaxis_button])
            # Particle selection container
            container_ptcl_select = ptcl_select_widget.to_container()
            # Plotting options container
            container_ptcl_plots = widgets.VBox(width=310,
                children=[
                    ptcl_figure_button, ptcl_bins_button, ptcl_range_button,
                    widgets.HBox(
                        children=[
                            ptcl_magnitude_button,
                            ptcl_use_button],
                        height=50),
                    ptcl_color_button])
            # Accordion for the field widgets
            accord2 = widgets.Accordion(
                children=[container_ptcl_quantities, container_ptcl_select,
                          container_ptcl_plots])
            accord2.set_title(0, 'Particle quantities')
            accord2.set_title(1, 'Particle selection')
            accord2.set_title(2, 'Plotting options')
            # Complete particle container
            container_ptcl = widgets.VBox(width=370,
                children=[accord2, widgets.HBox(
                    children=[ptcl_refresh_toggle, ptcl_refresh_button])])

        # Global container
        if (self.avail_fields is not None) and \
                (self.avail_species is not None):
            global_container = widgets.HBox(
                children=[container_fld, container_ptcl])
            display(global_container)
        elif self.avail_species is None:
            display(container_fld)
        elif self.avail_fields is None:
            display(container_ptcl)


def convert_to_int(m):
    """
    Convert the string m to an int, except if m is 'all' or None
    """
    if (m == 'all') or (m is None):
        return(m)
    else:
        return(int(m))


class ParticleSelectWidget(object):

    """
    Class that groups the particle selection widgets.
    """

    def __init__(self, n_rules, avail_records, refresh_ptcl):
        """
        Initialize a set of particle selection widgets

        Parameters:
        -----------
        n_rules: int
            The number of selection rules to display

        avail_records: list of strings
            The list of available records for the current species

        refresh_ptcl: callable
            The callback function to execute when the widget is changed
        """
        self.n_rules = n_rules

        # Create widgets that determines whether the rule is used
        self.active = [widgets.Checkbox(value=False)
                       for i in range(n_rules)]
        # Create widgets that determines the quantity on which to select
        # (The Dropdown menu is empty, but is later populated by the
        # function refresh_species)
        self.quantity = [widgets.Dropdown(options=avail_records,
            description='Select ') for i in range(n_rules)]
        # Create widgets that determines the lower bound and upper bound
        self.low_bound = [widgets.FloatText(value=-1.e-1, width=90,
            description='from ') for i in range(n_rules)]
        self.up_bound = [widgets.FloatText(value=1.e-1, width=90,
            description='to ') for i in range(n_rules)]

        # Add the callback function refresh_ptcl to each widget
        for i in range(n_rules):
            self.active[i].observe( refresh_ptcl, 'value', 'change' )
            self.quantity[i].observe( refresh_ptcl, 'value', 'change' )
            self.low_bound[i].observe( refresh_ptcl, 'value', 'change' )
            self.up_bound[i].observe( refresh_ptcl, 'value', 'change' )

    def to_container(self):
        """
        Return a widget container, where all the particle selection
        widgets are placed properly, with respect to each other.
        """
        containers = []
        for i in range(self.n_rules):
            containers.append(widgets.HBox(height=40,
                children=[self.active[i], self.quantity[i]]))
            containers.append(widgets.HBox(height=50,
                children=[self.low_bound[i], self.up_bound[i]]))

        return(widgets.VBox(children=containers, width=310))

    def to_dict(self):
        """
        Return a selection dictionary of the form
        {'uz': [-0.1, 2.], 'x':[-10., 10.]}
        depending on the values of the widgets.
        """
        rule_dict = {}
        # Go through the selection rules and add the active rules
        for i in range(self.n_rules):
            if self.active[i].value is True:
                rule_dict[ self.quantity[i].value ] = \
                    [self.low_bound[i].value, self.up_bound[i].value]

        # If any rule is active, return a dictionary
        if len(rule_dict) != 0:
            return(rule_dict)
        # If no rule is active, return None
        else:
            return(None)
