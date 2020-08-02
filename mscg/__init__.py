'''the mscg package provides a toolkit of data structures and API 
functions, which can be used to build up a customized or extended 
workflow of multi-scale coarse-graining modeling in Python. These
APIs are designed with object-oriented programming paradigm and
organized in sub-packages and modules. To use these APIs, simply
load the package in a Python script: ::
    
    >>> import mscg as cg
    >>> help(cg)

or, directly import sub-packages and modules into current 
namespace: ::

    >>> from mscg import *

'''

__version__ = '0.1.0'

doc_root = "https://software.rcc.uchicago.edu/mscg/docs/"

from .topology   import *
from .trajectory import *
from .pairlist   import *
from .bondlist   import *
from .tables     import *
from .matrix     import *
from .bspline    import *
from .mapper     import *

from .timer      import *
from .verbose    import *
from .cli_parser import *
from .checkpoint import *

from .tables import tables

from .model import *
from .ucg import *

from .traj_reader import *
