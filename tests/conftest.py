


import pytest
import os

@pytest.fixture(scope='session')
def datafile():
    pwd = os.path.dirname(os.path.abspath(__file__))
    
    def fullpath(fn):
        return os.path.join(pwd, "data", fn)
    
    return fullpath


@pytest.fixture
def show():
    def __show__(name, arr):
        print(name + ": [" + ", ".join([str(_) for _ in arr]) + "]")
        
    return __show__