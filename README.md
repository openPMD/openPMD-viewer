Opmd viewer
=======

Overview
-------

This package contains a set of tools to load and visualize the contents of an
[openPMD](http://www.openpmd.org/#/start) file simulation.

Usage
-----
The `opmd_viewer` routines can be used in two ways :

- Use the **Python API**, in order to write a script that loads the data and produces a set of
pre-defined plots.  
To do this, copy the file `examples/plotting_script.py` to your local
directory, modify it to  suit your needs, and run it with Python:  
```python plotting_script.py```

- Use the **interactive interface for IPython Notebook**, in order to interactively
visualize the data.  
To do this,  copy the file `examples/interactive_plotting.ipynb` to
your local directory, and run it with the ipython notebook:  
```ipython notebook interactive_plotting.ipynb```

All the objects and functions are internally documented, and thus their documentation can be accessed through the attribute `.__doc__` or (in an ipython shell) by using `?`.

For testing purposes, 3 example sets of openPMD files can be found [here](https://bitbucket.org/berkeleylab/opmd_viewer/downloads).

Installation
--------

To install this package :

- Clone this repository using `git`

- `cd` into the directory `opmd_viewer` and run `python setup.py install`

**Warning:** The **interactive interface** for IPython Notebook is not
installed by default.  
This is because it requires dependencies that may be difficult to
install on some systems (e.g. Edison@NERSC).  
If you wish to have the interactive interface working, install the
following dependencies by hand:

- [IPython Notebook](http://ipython.org/notebook.html)  (version 4.0
or higher)  
`pip install "ipython[notebook]"` or `pip install --upgrade "ipython[notebook]"`

- [ipywidgets](https://pypi.python.org/pypi/ipywidgets/4.0.2)  
`pip install ipywidgets`
