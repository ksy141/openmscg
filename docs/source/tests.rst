Regression Tests
================

Regression tests are organized in the ``tests`` folder in the 
source folder, and built within the popular and standard 
`pytest <https://docs.pytest.org/en/latest/>`_ framework.


Run Tests
---------

Download the trajectories data from::

    https://software.rcc.uchicago.edu/mscg/downloads/data.tar.gz


Unpack the tar-ball and move it into the source folder::
    
    tar xzvf data.tar.gz
    mv data ${path_to_openmscg_package}/tests/

    
Ensure you have ``pytest`` installed in your Python environment.
Run the following command in the source folder::

    pytest

The command above is used to run through all tests in this package. To run the tests in a specific module, use the command::
    
    pytest tests/[test_file].py

Or, just run a single test function in a file::
    
    pytest tests/[test_file].py::[test_function]


    
Build New Test
--------------

1. Create a new script in the ``tests`` folder and name it as ``test_*.py``.

2. Add following heading lines to initialize the testing environment::

    import pytest
    from mscg import *

  Or, if testing new CLI module::
    
    from mscg.cli import [cli_module_name]

3. Create testing module as a function prefixed by ``test__``, for example,::

    def test_example(datafile):
        
        asset "This is a test." == ""

  ``datafile`` is a pre-defined funtion decorator to get the full path for a
  data file, for example,::

    datafile('lj_liquid.top')





