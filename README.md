# openPMD viewer


## Overview

This package contains a set of tools to load and visualize the
contents of a set of [openPMD](http://www.openpmd.org/#/start) files
(typically, a timeseries).

## Installation

To install this package :

- Clone this repository using `git`

- `cd` into the directory `openPMD-viewer` and run `python setup.py install`

**Warning:** The **interactive GUI** for IPython Notebook is not
operational by default.  
This is because it requires dependencies that may be difficult to
install on some systems. If you wish to have the interactive GUI
working, install the following dependencies by hand:

- [IPython Notebook](http://ipython.org/notebook.html)  (version 4.0
or higher)  
`pip install "ipython[notebook]"` or `pip install --upgrade "ipython[notebook]"`

- [ipywidgets](https://pypi.python.org/pypi/ipywidgets/4.0.2)  
`pip install ipywidgets`

NB: For [NERSC](http://www.nersc.gov/) users, it is not necessary to
install the above two packages, as NERSC provides replacements for
them, when logging to
[https://ipython.nersc.gov](https://ipython.nersc.gov).
Therefore, NERSC users only need to install the `openPMD-viewer`
package itself.

## Usage

The routines of openPMD viewer can be used in two ways :

- Use the **Python API**, in order to write a script that loads the
  data and produces a set of pre-defined plots.

- Use the **interactive GUI inside the IPython Notebook**, in order to interactively
visualize the data.

#### Tutorials

The notebooks in the folder `tutorials/` demonstrate how to use both
the API and the interactive GUI. You can view these notebooks online
[here](https://github.com/openPMD/openPMD-viewer/tree/master/tutorials),
or, alternatively, you can run them on your local computer by typing:

`ipython notebook tutorials/`

NB: For [NERSC](http://www.nersc.gov/) users, you can run the tutorials on a
remote machine by logging in at
[https://ipython.nersc.gov](https://ipython.nersc.gov), and by
navigating to your personal copy of the directory `openPMD-viewer/tutorials`.

#### Notebook quick-starter

If you wish to use the **interactive GUI**, the installation of `openPMD-viewer` provides
a convenient executable which automatically
**creates a new pre-filled notebook** and **opens it in a
browser**. To use this executable, simply type in a regular terminal:

`openPMD_notebook`

(This executable is installed by default when running `python setup.py install`.)

## Contributing to the openPMD-viewer

We welcome contributions to the code! Please read [this page](https://github.com/openPMD/openPMD-viewer/blob/master/CONTRIBUTING.md) for
guidelines on how to contribute.
