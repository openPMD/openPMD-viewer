#!/usr/bin/env python
from setuptools import setup, find_packages

setup(name='opmd_viewer',
      version='1.0',
      author='Remi Lehe',
      author_email='remi.lehe@lbl.gov',
      description='Visualization tools for OpenPMD files',
      url='git@bitbucket.org:berkeleylab/opmd_viewer.git',
      packages = find_packages('./'),
      package_data = {'opmd_viewer':['notebook_starter/*.ipynb']},
      scripts = ['opmd_viewer/notebook_starter/openPMD_notebook'],
      install_requires=['wget', 'numpy', 'scipy', 'matplotlib', 'h5py'],
      tests_require=['pytest'],
      setup_requires=['pytest-runner']
    )
