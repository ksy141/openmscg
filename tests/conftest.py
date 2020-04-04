


import pytest
import os

@pytest.fixture(scope='session')
def datafile():
    pwd = os.path.dirname(os.path.abspath(__file__))
    
    def fullpath(fn):
        return os.path.join(pwd, "data", fn)
    
    return fullpath
