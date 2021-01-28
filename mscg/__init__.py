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

__version__ = '0.2.2'

doc_root = "https://software.rcc.uchicago.edu/mscg/docs/"

threads_envs = ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"]

import os

threads = 0 if 'OPENMSCG_THREADS' not in os.environ else int(os.environ['OPENMSCG_THREADS'])

if threads > 0:
    for name in threads_envs:
        os.environ[name] = str(threads)
else:
    threads = len(os.sched_getaffinity(os.getpid()))

import numpy as np

D2R = np.pi / 180.0
R2D = 180.0 / np.pi

from .verbose    import *
from .topology   import *
from .trajectory import *
from .pairlist   import *
from .bondlist   import *
from .bspline    import *
from .force      import *
from .mapper     import *
from .table      import *
from .timer      import *
from .cli_parser import *
from .checkpoint import *

from .model import *
from .ucg import *

from .traj_reader import *
