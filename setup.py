#!/usr/bin/env python
from setuptools import setup, find_packages

# Get the package requirements from the requirements.txt file
with open('requirements.txt') as f:
    install_requires = [ line.strip('\n') for line in f.readlines() ]
# Wget is installed separately, since installing it with conda fails
install_requires.append('wget')
    
# Main setup command
setup(name='opmd_viewer',
      version='1.0',
      author='Remi Lehe',
      author_email='remi.lehe@lbl.gov',
      description='Visualization tools for OpenPMD files',
      url='git@bitbucket.org:berkeleylab/opmd_viewer.git',
      packages = find_packages('./'),
      package_data = {'opmd_viewer':['notebook_starter/*.ipynb']},
      scripts = ['opmd_viewer/notebook_starter/openPMD_notebook'],
      install_requires=install_requires,
      tests_require=['pytest', 'jupyter'],
      setup_requires=['pytest-runner']
    )
