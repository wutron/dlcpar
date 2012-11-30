#
# Python module for dlcpar library
#


#
# Note, this module requires the rasmus, compbio python modules.
#

import sys, os

def load_deps(dirname="deps"):
    sys.path.append(os.path.realpath(
            os.path.join(os.path.dirname(__file__), dirname)))

# add pre-bundled dependencies to the python path,
# if they are not available already
try:
    import rasmus, compbio
except ImportError:
    load_deps()
    import rasmus, compbio


