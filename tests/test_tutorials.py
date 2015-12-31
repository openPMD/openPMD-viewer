"""
This test file is part of the openPMD-viewer.

It makes sure that the tutorial notebooks run without error.

Usage:
This file is meant to be run from the root directory of openPMD-viewer,
by any of the following commands
$ python tests/test_tutorials.py
$ py.test
$ python setup.py test
"""
import os, urllib, tarfile
import re

def test_tutorials():
    """Test all the tutorial notebooks"""

    # Go to the relative path where all tutorial notebooks are
    os.chdir( 'tutorials' )
    tutorial_notebooks = [ filename for filename in os.listdir('./') \
                           if filename[-6:] == '.ipynb' ]

    # Loop through the tutorials and test them
    for notebook_name in tutorial_notebooks:

        # Do a first pass where only the non-IPython features are tested.
        # (This gives better debugging information.)
        # The notebook is converted to a standard Python script and
        # run directly with `execfile`
        script_name = notebook_name[:-6] + '.py'
        os.system( 'ipython nbconvert --to=python %s' %notebook_name )
        clean_ipython_features( script_name )
        execfile( script_name )

        # Do a second pass where the full IPython features are executed.
        # The notebook is executed by the OS, through `ipython nbconvert`
        response = os.system( "ipython nbconvert --to=python "\
            + "--ExecutePreprocessor.enabled=True %s" %notebook_name)
        # Check the response from the OS
        if response != 0:
            raise OSError("The execution of the notebook %s by the OS "
                          "with `ipython nbconvert` has failed.")

        # Clean the produced script file
        os.remove( script_name )


def clean_ipython_features( script_name ):
    """
    Rewrites the Python script `script_name` by removing
    all the IPython-specific commands
    
    Parameters
    ----------
    script_name: string
        Name of the Python script
    """
    # Read the script file
    with open(script_name) as script_file:
        lines = script_file.readlines()

    # Go through the lines and replace the IPython-specific commands
    # using regular expressions
    for i in range(len(lines)):

        # Replace the lines that activate matplotlib in a notebook
        # by a line that selects the PS backend
        if re.match("[ ]*get_ipython\(\)\.magic\(u'matplotlib",
                    lines[i]) is not None:
            lines[i] = "import matplotlib; matplotlib.use('ps')\n"

        # Discard the lines that use in-notebook documentation
        if re.match("[ ]*get_ipython\(\)\.magic\(u'pinfo",
                    lines[i]) is not None:
            lines[i] = ''

        # Discard the lines that use the GUI
        if re.match("[\w]*\.slider", lines[i]) is not None:
            lines[i] = ''
                        
        # Replace the lines that call the OS by proper lines
        if re.match("[ ]*get_ipython\(\)\.system", lines[i]) is not None:
            matched = re.match("([ ]*)get_ipython\(\)\.system(.*)", lines[i])
            spaces = matched.groups()[0]
            command_line = matched.groups()[1]
            lines[i] = '%simport os; os.system%s\n' %(spaces,command_line)

    # Write the cleaned file
    with open(script_name, 'w') as script_file:
        for line in lines:
            script_file.write( line )

    

if __name__=='__main__':
    test_tutorials()
