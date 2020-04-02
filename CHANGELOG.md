# Change Log / Release Log for openPMD-viewer

## 1.0.1

This is a bug-fix release.

- Unreadable files are now skipped (instead of crashing the whole timeseries) ; see [#262](https://github.com/openPMD/openPMD-viewer/pull/262).

- A bug related to units (microns vs meters) was fixed in `get_emittance` and `get_current` (see [#276](https://github.com/openPMD/openPMD-viewer/pull/276))

- The quick-start notebook (`openPMD-visualization`) raised a warning saying
`Notebook validation failed` in some cases. This was fixed (see [#274](https://github.com/openPMD/openPMD-viewer/pull/274))

- When using the option `plot=True` with Python 2, openPMD-viewer crashed. This is now fixed (see [#271](https://github.com/openPMD/openPMD-viewer/pull/271)).

## 1.0.0

This version introduces major changes and breaks backward compatibility.

Here is a list of the changes:
- The import statement now uses `openpmd_viewer` instead of `opmd_viewer`, e.g.
```
from openpmd_viewer import OpenPMDTimeSeries
```
- For consistency, `ts.get_particle` now return particle positions in meters,
instead of microns. For instance, in the code below, `x`, `y`, `z` will be in
meters
```
x, y, z = ts.get_particle(['x', 'y', 'z'], iteration=1000)
```
- In `ts.get_field`, slicing can now be done in several directions, and for
1d, 2d, 3d, and circ geometries. As a consequence, this breaks backward
compatibility for 3d field:
```get_field(field=field, coord=coord, iteration=iteration)```
used to return the central slice along `y` while it now returns the full 3d field.
In addition, the name of the argument of `get_field` that triggers slicing
has been changed from `slicing_dir` to `slice_across` (and `slicing` has been
changed to `slice_relative_position`).
- `openPMD-viewer` does not rely on Cython anymore. Instead, it uses `numba`
for functions that perform a substantial amount of computation.
- A new function (`ts.iterate`) was introduced in order to quickly apply a
given function to all iterations of a time series. See the docstring of
`ts.iterate` for more information.
- The function `get_laser_envelope` does not support the argument `index` anymore
(which was effectively used in order to perform slicing). Instead, users should use
the argument `slicing_dir`. In addition, `get_laser_envelope` now supports the
argument `plot_range`.
- The function `get_laser_waist` does not support the agument `slicing_dir` anymore.

## 0.9.0

This release adds two features:
- Improved calculation of the laser envelope, using the Hilbert transform.
- Reconstruction of full 3D field from a quasi-3D dataset, when passing `theta=None`.

## 0.8.2

This is a bug-fix release. It allows the slider to work properly in JupyterLab,
by using the `%matplotlib widget` magic.

## 0.8.1

This version includes minor improvements to the viewer:
- (Experiemental) support for Windows users
- In the interactive Jupyter GUI, the user can now select the scale of the vertical axis.
- The function `get_emittance` has more options (including calculation of the slice emttance)
- The default `openPMD_notebook` now avoids warning messages about matplotlib inline, which used to occur even though `%matplotlib notebook` was used.

Many thanks to @MaxThevenet and @AngelFP for their contributions to this release!

## 0.8

This version introduces several improvements to the viewer:
- The ability to read files that contain fields in different geometries
(e.g. 3D fields and 2D slices).
- Better support for files that do not contain mesh (or do not contain
particles), including support for the openPMD 1.1.0 standard.
- Cloud-In-Cell deposition in histograms.
- Better handling of `%matplotlib notebook` for newer version of jupyter.

## 0.7.1

This version adds better support, when the local installation of matplotlib
has issues:

- The `LpaDiagnostics` can now work without matplotlib if needed.
- The `MacOSX` matplotlib backend is now avoided, since there can be issues
when using it in the latest version of Jupyter.

## 0.7.0

This version improves support for `ipywidgets` version 7, especially in
the layout of the slider.

In addition, with this version of `openPMD-viewer`, `matplotlib` is  not a
strict requirement anymore. This allows lighter installation for users that
need `openPMD-viewer` only as a data reader.

Finally, the calculation of the laser envelope in 2D has been improved
(see [PR 170](https://github.com/openPMD/openPMD-viewer/pull/170)). Note
that the function `wstd` (which is not documented in the tutorial, but
which some users might still use) has been renamed to `w_std`.

## 0.6.0

This version improves the layout of the Jupyter GUI and allows the user to
select a particular region of the plots through this GUI.

In addition, support for massless particle (e.g. photons) was added. In this
case, the momenta are returned in kg.m.s^-1, instead of using the
dimensionless momenta.

## 0.5.4

This is version `0.5.4` of openPMD-viewer.

It adds support for Python 3.4 (which erroneously dropped in the past).

## 0.5.3

This is version `0.5.3` of openPMD-viewer.

It corrects some of the issues with the size of boxes and widgets in the
interactive slider. In addition, the iteration number is now read from
the hdf5 metadata, and not the name of the file.

## 0.5.2

This is version `0.5.2` of openPMD-viewer.

It fixes some of the installation issues associated with Cython.

## 0.5.1

This is version `0.5.1` of openPMD-viewer.

It corrects a minor bug in the dependencies of the package.

## 0.5.0

This is version `0.5.0` of openPMD-viewer.

This new version includes the `ParticleTracker` object, which allows user to track individual particles across different iterations, provided that their `id` is stored in the openPMD file. Also, starting with this version, openPMD-viewer now depends on `Cython`.

For more information on how to use the `ParticleTracker`, see the tutorial notebook.

## 0.4.0

This is version `0.4.0` of openPMD-viewer.

This new version includes:
- support for 1D data
- an additional option `use_field_mesh` when plotting the particle. When set
to `True`, this option uses information from the field mesh to choose the parameters of the particle histograms (esp. the bins). This is useful in order to avoid plotting/binning artifacts (aliasing) when the particles are evenly spaced.

In addition, the module `openpmd_viewer` now has an attribute `__version__`.

## 0.3.3

This is version `0.3.3` of openPMD-viewer.

This version fixed a bug with the executable `openPMD_notebook`. More precisely, the executable was not installed, when using `pip` or `conda`. In addition, it was failing with Python 3.

## 0.3.2

This is version `0.3.2` of openPMD-viewer. The following changes were introduced:

- The conda recipe in `conda_recipe/` was simplified and streamlined.
- The documentation now explains how to install openPMD-viewer with `conda`, the instructions to release the package was put into a document `RELEASING.md`.
- A file `MANIFEST.in` was added, to avoid issues with pip and Python 3.

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
