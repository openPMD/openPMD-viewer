# Creating a new release

This document is only relevant for maintainers of openPMD-viewer. It
explains how to create a new release. In future versions of this
packages, some of the steps below will be automatized.

## Preparing your environment for a release

Make sure that your local environment is ready for a full release on
PyPI and conda. In particular:

- you should install the package
[`twine`](https://pypi.python.org/pypi/twine).
- you should have a registered account on [PyPI](https://pypi.python.org/pypi) and [test PyPI](https://testpypi.python.org/pypi), and your `$HOME` should contain a file `.pypirc` which contains the following text:

 ```
[distutils]
index-servers=
	pypitest
	pypi

[pypitest]
repository = https://testpypi.python.org/pypi
username = <yourPypiUsername>

[pypi]
repository = https://pypi.python.org/pypi
username = <yourPypiUsername>
```

- you should have a registered account on [Anaconda.org](https://anaconda.org/)

## Creating a release on Github

- Make sure that the version number in `openpmd_viewer/__version__.py`
  correspond to the new release, and that
  the corresponding changes have been documented in `CHANGELOG.md`.

- Create a new release through the graphical interface on Github

## Uploading the package to PyPI

The code will be automatically uploaded to PyPI upon creation of a new release on Github.