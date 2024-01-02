# Change Log / Release Log for openPMD-viewer

## 1.10.0

* use `%matplotlib widget` by default by @BenWibking in https://github.com/openPMD/openPMD-viewer/pull/407
* BP4 Empty Skip: Inline Comment by @ax3l in https://github.com/openPMD/openPMD-viewer/pull/403
* Close #402 Fix deprecated ipython command by @RemiLehe in https://github.com/openPMD/openPMD-viewer/pull/410

**Full Changelog**: https://github.com/openPMD/openPMD-viewer/compare/1.9.0...1.10.0

## 1.9.0

* When passing `t`, choose the closest iteration by @RemiLehe in https://github.com/openPMD/openPMD-viewer/pull/347
* README: Remove Travis-CI Badges by @ax3l in https://github.com/openPMD/openPMD-viewer/pull/369
* Adding attributes to the FieldMetaInformation object by @juliettepech in https://github.com/openPMD/openPMD-viewer/pull/372
* Sphinx documentation by @RemiLehe in https://github.com/openPMD/openPMD-viewer/pull/303
* Add readthedoc configuration by @RemiLehe in https://github.com/openPMD/openPMD-viewer/pull/375
* setup.py: no upper version of openPMD-api by @ax3l in https://github.com/openPMD/openPMD-viewer/pull/378
* Docs: correct 'get_mean_gamma' description in tutorials by @IlianCS in https://github.com/openPMD/openPMD-viewer/pull/379
* Do not use pyflakes for version 3.7 by @RemiLehe in https://github.com/openPMD/openPMD-viewer/pull/385
* Impose that user always passes iteration or t by @RemiLehe in https://github.com/openPMD/openPMD-viewer/pull/383
* Docs: Improve tutorials by adding basic analysis functions by @IlianCS in https://github.com/openPMD/openPMD-viewer/pull/382
* Improve get_laser_waist by @soerenjalas in https://github.com/openPMD/openPMD-viewer/pull/359
* Python 3.8+ by @ax3l in https://github.com/openPMD/openPMD-viewer/pull/387
* Update version number and CHANGELOG.md by @RemiLehe in https://github.com/openPMD/openPMD-viewer/pull/386
* Add all record attributes to `FieldMetaInformation` by @AngelFP in https://github.com/openPMD/openPMD-viewer/pull/390
* Support complex data in `thetaMode` geometry by @AngelFP in https://github.com/openPMD/openPMD-viewer/pull/389
* fix get_data for weighting in momentum reader by @PrometheusPi in https://github.com/openPMD/openPMD-viewer/pull/393
* Simplify names of radial particle components by @AngelFP in https://github.com/openPMD/openPMD-viewer/pull/392
* Close #396 warning when reading particle data by @RemiLehe in https://github.com/openPMD/openPMD-viewer/pull/398
* Support reading of 'rt' lasy files by @RemiLehe in https://github.com/openPMD/openPMD-viewer/pull/371

**Full Changelog**: https://github.com/openPMD/openPMD-viewer/compare/1.7.0...1.9.0

## 1.8.0

- The functions `get_field` and `get_particle` now require `iteration` or `t`
  to be passed (instead of using a default iteration when none was provided).
  (See [#383](https://github.com/openPMD/openPMD-viewer/pull/383))

- The function `get_laser_waist` is more robust and does not crash when the
  laser field is 0.
  (See [#359](https://github.com/openPMD/openPMD-viewer/pull/359))

- The `FieldMEtaInformation` object has new attributes `time` and `iteration`.
  (See [#372](https://github.com/openPMD/openPMD-viewer/pull/372))

- The docstring of `get_mean_gamma` has been updated
  (See [#379](https://github.com/openPMD/openPMD-viewer/pull/379))
  and the attributes `ts.t` and `ts.iterations` are now shown in the tutorials
  (See [#382](https://github.com/openPMD/openPMD-viewer/pull/382))

## 1.7.0

This release includes a few improvements:

- The function `get_laser_waist` is more robust: it does not automatically
assume that the `z` axis is the last axis of the data. In addition, the user
can now  specify the laser propagation direction. (The default is `z`.)
(See [#345](https://github.com/openPMD/openPMD-viewer/pull/345))

- The handling of `unitSI` is now more robust. (See [#363](https://github.com/openPMD/openPMD-viewer/pull/363))

## 1.6.0

This release adds a few features:

- `openPMD-viewer` can now read complex datasets (See [#353](https://github.com/openPMD/openPMD-viewer/pull/353))

- Avoid errors in LPA diagnostics in the absence of selected particles (See [#358](https://github.com/openPMD/openPMD-viewer/pull/358))

## 1.5.0

This release fixes a few miscellaneous bugs:

- Better 3D reconstruction for `theta=None` (See [#344](https://github.com/openPMD/openPMD-viewer/pull/344))

- Better support for ADIOS data (See [#355](https://github.com/openPMD/openPMD-viewer/pull/355))

- Support for group-based encoding (See [#346](https://github.com/openPMD/openPMD-viewer/pull/346))

## 1.4.0

This new release introduces several improvements:

- The changes introduced in 1.3.0 caused a major slowdown when reading certain
types of data. This has been fixed in this new release. (See [#340](https://github.com/openPMD/openPMD-viewer/pull/340) for more details.)

- `openPMD-viewer` now supports `thetaMode` geometry with data written using
`r` as the fastest index (as written by e.g. [WarpX](https://github.com/ECP-WarpX/WarpX))
in addition to the previously supported data format which used `z` as the fastest index
(as written by e.g. [fbpic](https://github.com/fbpic/fbpic)). (See
[337](https://github.com/openPMD/openPMD-viewer/pull/337))

- `openPMD-viewer` will raise an exception if the user asks for an iteration
that is  not part of the dataset (instead of printing a message and reverting
to the first iteration, which can be confusing) (See [336](https://github.com/openPMD/openPMD-viewer/pull/336))

## 1.3.0

This new release introduces preliminary support for MR datasets
(see [#332](https://github.com/openPMD/openPMD-viewer/pull/332)).

## 1.2.0

This new release introduces several bug-fixes and miscellaneous features:

- There is a new function `get_energy_spread` that returns the energy
  spread of the beam. This is partially redundant with `get_mean_gamma`,
  which is kept for backward compatibility.
  (see [#304](https://github.com/openPMD/openPMD-viewer/pull/304)
  and [#317](https://github.com/openPMD/openPMD-viewer/pull/317))

- The 3D field reconstruction from `ThetaMode` data now has an option
  `max_resolution_3d` that limits the resolution of the final 3D array.
  This is added in order to limit the memory footprint of this array.
  (see [#307](https://github.com/openPMD/openPMD-viewer/pull/307))
  The 3D reconstruction is now also more accurate, thanks to the implementation
  of linear interpolation.
  (see [#311](https://github.com/openPMD/openPMD-viewer/pull/311))

- A bug that affected reading `ThetaMode` data with the `openpmd-api` backend
  has been fixed. (see [#313](https://github.com/openPMD/openPMD-viewer/pull/313))

- A bug that affected `get_laser_waist` has been fixed:
  (see [#320](https://github.com/openPMD/openPMD-viewer/pull/320))

## 1.1.0

This new release introduces the option to read `openPMD` files with different backends. In addition to the legacy `h5py` backend (which can read only HDF5 openPMD file), `openPMD-viewer` now has the option to use the `openpmd-api` backend (which can read both HDF5 and ADIOS openPMD files). Because the `openpmd-api` backend is thus more general, it is selected by default if available (i.e. if installed locally).
The user can override the default choice, by passing the `backend` argument when creating an `OpenPMDTimeSeries` object, and check which backend has been chosen by inspecting the `.backend` attribute of this object.

In addition, several smaller changes were introduced in this PR:
- The method `get_laser_envelope` can now take the argument `laser_propagation` in order to support lasers that do not propagates along the `z` axis.
- `openPMD-viewer` can now properly read `groupBased` openPMD  files (i.e. files that contain several iterations) [#301](https://github.com/openPMD/openPMD-viewer/pull/301).
- Users can now pass arrays of ID to the `ParticleTracker` [#283](https://github.com/openPMD/openPMD-viewer/pull/283)

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
