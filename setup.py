#!/usr/bin/env python
from setuptools import setup, find_packages
import opmd_viewer  # In order to extract the version number

# Get the package requirements from the requirements.txt file
with open('requirements.txt') as f:
    install_requires = [line.strip('\n') for line in f.readlines()]
# Since wget cannot be installed with conda, it is added separately here
install_requires.append('wget')

# Main setup command
setup(name='openPMD-viewer',
      version=opmd_viewer.__version__,
      maintainer='Remi Lehe',
      maintainer_email='remi.lehe@lbl.gov',
      description='Visualization tools for OpenPMD files',
      url='git@bitbucket.org:berkeleylab/opmd_viewer.git',
      packages=find_packages('./'),
      package_data={'opmd_viewer': ['notebook_starter/*.ipynb']},
      scripts=['opmd_viewer/notebook_starter/openPMD_notebook'],
      install_requires=install_requires,
      tests_require=['pytest', 'jupyter'],
      setup_requires=['pytest-runner'],
      platforms='any',
      classifiers=[
          'Programming Language :: Python',
          'Development Status :: 4 - Beta',
          'Natural Language :: English',
          'Environment :: Console',
          'Environment :: Jupyter',
          'Intended Audience :: Science/Research',
          'License :: BSD-3-Clause-LBNL',
          'Operating System :: OS Independent',
          'Topic :: Scientific/Engineering :: Physics',
          'Topic :: Scientific/Engineering :: Visualization',
          'Topic :: Database :: Front-Ends',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
      ]
      )
