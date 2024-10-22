# Contributing to the openPMD-viewer

## How to contribute

### Forking the repository

In order to contribute, please fork the [main repository](https://github.com/openPMD/openPMD-viewer):

- Click 'Fork' on the page of the main repository, in order to create a personal copy of this repository on your Github account.

- Clone this copy to your local machine:
```
git clone git@github.com:<YourUserLogin>/openPMD-viewer.git
```

### Implementing a new feature and adding it to the main repository

- Switch to the development branch
```
git checkout dev
```
and install it
```
python -m pip wheel .
python -m pip install *whl
```

- Start a new branch from the development branch, in order to
implement a new feature. (Choose a branch name that is representative of the
feature that you are implementing, e.g. `add-latex-labels` or
`fix-matplotlib-errors`)
```
git checkout -b <NewBranchName>
```

- Start coding. When your changes are ready, commit them.
```
git add <ChangedFiles>
git commit
```

- Synchronize your branch with the main repository. (It may have
  changed while you where implementing local changes.) Resolve merging
  issues if any, and commit the corresponding changes.
```
git pull git@github.com:openPMD/openPMD-viewer.git dev
```

- Test and check your code:
  - Use [pyflakes](https://pypi.python.org/pypi/pyflakes) to detect any potential bug.
  ```
  cd openPMD-viewer/
  pyflakes openpmd_viewer
  ```
  - Make sure that the tests pass (please install `wget` and `jupyter` before running the tests: `pip install wget jupyter`)
  ```
  python -m pip wheel .
  python -m pip install *whl matplotlib jupyter
  python -m pytest tests
  ```
  (Be patient: the `test_tutorials.py` can take approx. 20 seconds if
  you already downloaded the example openPMD files that are required
  in the tutorials. On the other hand, it can take several minutes if
  you have not previously downloaded these files.)

- Push the changes to your personal copy on Github
```
git push -u origin <NewBranchName>
```

- Go on your Github account and create a pull request between **your
  new feature branch** and the **dev branch of the main
  repository**. Please add some text to the pull request to describe
  what feature you just implemented and why. Please also make sure that
  the automated tests (on Github) return no error.

## Style and conventions

- Features that **modify** or **improve** the `OpenPMDTimeSeries` object
should be implemented in the
`openpmd_viewer/opempmd_timeseries` folder. Features that **build upon** the
`OpenPMDTimeSeries` object to create domain-specific analysis tools
(e.g. laser diagnostics for PIC simulations) should be implemented in
the `openpmd_viewer/addons` folder.

- Document the functions and classes that you write, by using a
  [docstring](https://www.python.org/dev/peps/pep-0257/). List the
  parameters and describe what the functions return, as in this
  example:
```python
def get_data( dset, i_slice=None, pos_slice=None ) :
    """
    Extract the data from a (possibly constant) dataset
    Slice the data according to the parameters i_slice and pos_slice

    Parameters:
    -----------
    dset: an h5py.Dataset or h5py.Group (when constant)
        The object from which the data is extracted

    i_slice: int, optional
       The index of the slice to be taken

    pos_slice: int, optional
       The position at which to slice the array
       When None, no slice is performed

    Returns:
    --------
    An np.ndarray (non-constant dataset) or a single double (constant dataset)
    """
```
Don't use documenting styles like `:param:`, `:return:`, or
`@param`, `@return`, as they are less readable.


- Lines of code should **never** have [more than 79 characters per line](https://www.python.org/dev/peps/pep-0008/#maximum-line-length).

- Names of variables, functions should be lower case (with underscore
  if needed: e.g. `get_field`). Names for classes should use the
  CapWords convention (e.g. `DataReader`). See [this page](https://www.python.org/dev/peps/pep-0008/#prescriptive-naming-conventions) for more details.
