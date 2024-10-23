from setuptools import setup, find_packages

# Get the long description
with open('./README.md') as f:
    long_description = f.read()
# Get the package requirements from the requirements.txt file
with open('./requirements.txt') as f:
    install_requires = [line.strip('\n') for line in f.readlines()]

# Read the version number, by executing the file openpmd_viewer/__version__.py
# This defines the variable __version__
with open('./openpmd_viewer/__version__.py') as f:
    exec( f.read() )

# Main setup command
setup(name='openPMD-viewer',
      version=__version__,
      description='Visualization tools for openPMD files',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='https://github.com/openPMD/openPMD-viewer.git',
      maintainer='Remi Lehe',
      maintainer_email='remi.lehe@lbl.gov',
      license='BSD-3-Clause-LBNL',
      packages=find_packages('.'),
      package_data={'openpmd_viewer': ['notebook_starter/*.ipynb']},
      scripts=['openpmd_viewer/notebook_starter/openPMD_notebook'],
      install_requires=install_requires,
      extras_require = {
        'all': ["ipympl", "ipywidgets", "matplotlib", "numba", "openpmd-api", "wget"],
        'GUI':  ["ipywidgets", "ipympl", "matplotlib"],
        'plot': ["matplotlib"],
        'tutorials': ["ipywidgets", "ipympl", "matplotlib", "wget"],
        'numba': ["numba"],
        'openpmd-api': ["openpmd-api"]
        },
      platforms='any',
      python_requires='>=3.8',
      classifiers=[
          'Programming Language :: Python',
          'Development Status :: 4 - Beta',
          'Natural Language :: English',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Topic :: Scientific/Engineering :: Physics',
          'Topic :: Scientific/Engineering :: Visualization',
          'Topic :: Database :: Front-Ends',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.10',
          'Programming Language :: Python :: 3.11'],
      )
