# openPMD viewer


## Overview

This package contains a set of tools to load and visualize the contents of an
[openPMD](http://www.openpmd.org/#/start) file simulation.

## Installation

To install this package :

- Clone this repository using `git`

- `cd` into the directory `opmd_viewer` and run `python setup.py install`

**Warning:** The **interactive interface** for IPython Notebook is not
operational by default.  
This is because it requires dependencies that may be difficult to
install on some systems (e.g. Edison@NERSC).  
If you wish to have the interactive interface working, install the
following dependencies by hand:

- [IPython Notebook](http://ipython.org/notebook.html)  (version 4.0
or higher)  
`pip install "ipython[notebook]"` or `pip install --upgrade "ipython[notebook]"`

- [ipywidgets](https://pypi.python.org/pypi/ipywidgets/4.0.2)  
`pip install ipywidgets`

## Usage

The routines of openPMD viewer can be used in two ways :

- Use the **Python API**, in order to write a script that loads the
  data and produces a set of pre-defined plots.

- Use the **interactive GUI inside the IPython Notebook**, in order to interactively
visualize the data.

The notebooks in the folder `tutorials/` demonstrate how to use these
routines. You can be view these notebooks online
[here](https://github.com/openPMD/openPMD-viewer/tree/master/tutorials),
or, alternatively, you can run them on your local computer by typing:

`ipython notebook tutorials/`

