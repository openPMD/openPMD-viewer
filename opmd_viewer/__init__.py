"""
This is a stub that detects whether the user is attempting to use the
old import statement from openPMD-viewer 0.X (`import opmd_viewer`).

In that case, an exception is raised that prompts the user to use the new API.
"""
# Define the version number
from openpmd_viewer import __version__
__all__ = ['__version__']

raise DeprecationWarning("""
It looks like you are trying to use the API from openPMD-viewer version 0.X
but the installed version on your system is openPMD-viewer version 1.X.

* If you wish to use the new openPMD-viewer version 1.X, the import statement
should use `openpmd_viewer` instead of `opmd_viewer`, e.g.:
```
from openpmd_viewer import OpenPMDTimeSeries
```
Please have a look at the list of the changes introduced in version 1.X here:
https://github.com/openPMD/openPMD-viewer/blob/upcoming-1.0/CHANGELOG.md#10
In particular, note that `get_particle` now returns particle positions in
meters (not in microns anymore) and that the syntax for slicing fields has
changed.

* If you wish to go back to the old openPMD-viewer version 0.X, use:
```
pip install openPMD-viewer==0.9
```
""")
