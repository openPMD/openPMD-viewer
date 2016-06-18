from setuptools import setup, find_packages
import opmd_viewer  # In order to extract the version number

# Get the long description
# If possible, use pypandoc to convert the README from Markdown
# to reStructuredText, as this is the only supported format on PyPI
try:
    import pypandoc
    long_description = pypandoc.convert( './README.md', 'rst')
except (ImportError, RuntimeError):
    long_description = open('./README.md').read()
# Get the package requirements from the requirements.txt file
with open('./requirements.txt') as f:
    install_requires = [line.strip('\n') for line in f.readlines()]

# Main setup command
setup(name='openPMD-viewer',
      version=opmd_viewer.__version__,
      description='Visualization tools for openPMD files',
      long_description=long_description,
      url='https://github.com/openPMD/openPMD-viewer.git',
      maintainer='Remi Lehe',
      maintainer_email='remi.lehe@lbl.gov',
      license='BSD-3-Clause-LBNL',
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
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Topic :: Scientific/Engineering :: Physics',
          'Topic :: Scientific/Engineering :: Visualization',
          'Topic :: Database :: Front-Ends',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5'],
      )
