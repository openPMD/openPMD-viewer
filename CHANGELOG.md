# Change Log / Release Log for openPMD-viewer

## 0.3.1

This is version `0.3.1` of openPMD-viewer. This version introduces minor changes in the way the tests are run in `setup.py`. The aim of these changes are to prepare a conda release.

## 0.3.0

This is version `0.3.0` of openPMD-viewer. This version mainly adapts the interactive GUI so that it can be used with the newer version of `ipwidgets` (`ipywidgets 5.0`), while still being compatible with previous versions of `ipwidgets`. A number of other minor changes have been introduced:

- In the method `get_particle`, the argument `species` is now optional in the case where there is only one species.
- A number of methods in the LPA addons (`LpaDiagnostics` class) now have an optional argument `plot`, which allows to directly plot the data.

## 0.2.0

This is version `0.2.0` of openPMD-viewer. A number of minor changes and fixes have been made in order to make the package more general and to prepare it for a PyPI release. Here are the main changes:

- Support for the deprecated widget package `IPython.html` has been dropped. From now on, users need to install the widget package `ipywidgets`, for the GUI to work.
- The initialization of an `OpenPMDTimeSeries` object can now be made faster by setting the optional argument `check_all_files` to `False`.
- The data reader can now support `macroWeighted` quantities. As consequence, output files from [PIConGPU](https://github.com/ComputationalRadiationPhysics/picongpu) can now be correctly read.
- The package does not assume anymore that all species contain the same particle quantities. For instance, the package will support a file that contains the positions of ions, and the positions, momenta and weighting of electrons. As part of this, the attribute `OpenPMDTimeSeries.avail_ptcl_quantities` has been replaced by a dictionary `OpenPMDTimeSeries.avail_record_components`.
- This release introduces automatic PEP8 verification as part of the automatic tests that are run on Travis CI (see CONTRIBUTING.md).
- The evaluation of the waist and duration of the laser is now based on Gaussian fit of the transverse and longtudinal profile respectively.

