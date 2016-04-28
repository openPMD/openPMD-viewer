# openPMD-viewer

[![Build Status master](https://img.shields.io/travis/openPMD/openPMD-viewer/master.svg?label=master)](https://travis-ci.org/openPMD/openPMD-viewer/branches)
[![Build Status dev](https://img.shields.io/travis/openPMD/openPMD-viewer/dev.svg?label=dev)](https://travis-ci.org/openPMD/openPMD-viewer/branches)
[![pypi version](https://img.shields.io/pypi/v/openPMD-viewer.svg)](https://pypi.python.org/pypi/openPMD-viewer)
[![Number of PyPI Downloads](https://img.shields.io/pypi/dm/openPMD-viewer.svg)](https://pypi.python.org/pypi/openPMD-viewer)
[![License](https://img.shields.io/pypi/l/openPMD-viewer.svg)](LICENSE.txt)

## Overview

This package contains a set of tools to load and visualize the
contents of a set of [openPMD](http://www.openpmd.org/#/start) files
(typically, a timeseries).

## Installation

### Basic installation

To install this package :

- Clone this repository using `git`
```
git clone https://github.com/openPMD/openPMD-viewer.git
```

- `cd` into the directory `openPMD-viewer` and run
```
python setup.py install
```

### Installing the interactive GUI

The **interactive GUI** for IPython Notebook is not
operational by default.  
This is because it requires dependencies that may be difficult to
install on some systems. If you wish to have the interactive GUI
working, install the
[IPython Notebook](http://ipython.org/notebook.html)
(now part of the [Jupyter project](http://jupyter.org/)) by hand:  
`conda install jupyter` (for the
[Anaconda](https://www.continuum.io/downloads)
distribution) or `pip install jupyter` (for the other Python distributions)

NB: For [NERSC](http://www.nersc.gov/) users, it is not necessary to
install the above package, as NERSC provides it when logging to
[https://ipython.nersc.gov](https://ipython.nersc.gov).
Therefore, NERSC users only need to install the `openPMD-viewer`
package itself.

## Usage

The routines of openPMD-viewer can be used in two ways :

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
