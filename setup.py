#!/usr/bin/env python

from setuptools import setup

setup(name='opmd_viewer',
      version='1.0',
      author='Remi Lehe',
      author_email='remi.lehe@lbl.gov',
      description='Visualization tools for OpenPMD files',
      url='git@bitbucket.org:berkeleylab/opmd_viewer.git',
      install_requires=['numpy', 'scipy', 'matplotlib', 'h5py'],
      packages = ['opmd_viewer', 'opmd_viewer.openpmd_timeseries', 
                  'opmd_viewer.openpmd_timeseries.data_readers', 
                  'addons', 'addons/pic' ]
      )
