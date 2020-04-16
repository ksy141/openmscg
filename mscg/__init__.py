#

doc_root = 'https://software.rcc.uchicago.edu/mscg/docs/'

__version__ = '0.1.0'
__doc__ = doc_root

from .topology   import *
from .trajectory import *
from .pairlist   import *
from .bondlist   import *
from .tables     import *
from .matrix     import *
from .bspline    import *

from .timer      import *
from .verbose    import *
from .cli_parser import *

from .tables import tables
