from __future__ import print_function
import sys
if (sys.version_info > (3, 0)):
    from .utils import *
    # import .latex
    # import .gnuplot
else:
    from utils import *
    import latex
    import gnuplot

__all__ = ["utils", "latex", "gnuplot","progress","logs"]
