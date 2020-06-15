"""
dlcpar
================================================
Python module for dlcpar software
"""


#
# Note, this module requires the rasmus, compbio python modules.
#

# python libraries
import os
import sys

def load_deps(dirname="deps"):
    """Add path for dependencies"""
    sys.path.append(os.path.realpath(
        os.path.join(os.path.dirname(__file__), dirname)))

# add pre-bundled dependencies to the python path,
# if they are not available already
try:
    import rasmus
    import compbio
except ImportError:
    load_deps()
    import rasmus
    import compbio


#=============================================================================
# constants

PROGRAM_NAME = u"DLCPar"
PROGRAM_VERSION_MAJOR = 2
PROGRAM_VERSION_MINOR = 0
PROGRAM_VERSION_RELEASE = 0
PROGRAM_VERSION = (PROGRAM_VERSION_MAJOR,
                   PROGRAM_VERSION_MINOR,
                   PROGRAM_VERSION_RELEASE)

if PROGRAM_VERSION_RELEASE != 0:
    PROGRAM_VERSION_TEXT = "%d.%d.%d" % (PROGRAM_VERSION_MAJOR,
                                         PROGRAM_VERSION_MINOR,
                                         PROGRAM_VERSION_RELEASE)
else:
    PROGRAM_VERSION_TEXT = "%d.%d" % (PROGRAM_VERSION_MAJOR,
                                      PROGRAM_VERSION_MINOR)
