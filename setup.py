#!/usr/bin/env python
# 
# setup for DLCpar library packages
#
# use the following to install:
#   python setup.py install
#

import os,sys
from distutils.core import setup

sys.path.insert(0, os.path.realpath(
            os.path.join(os.path.dirname(__file__), "python")))
import dlcpar
VERSION = dlcpar.PROGRAM_VERSION_TEXT

setup(
    name='dlcpar',
    version=VERSION,
    description='DLCpar',

    author='Yi-Chieh Wu',
    author_email='yjw@mit.edu',
#    url='http://compbio.mit.edu/dlcpar/',
#    download_url='http://compbio.mit.edu/dlcpar/pub/sw/dlcpar-%s.tar.gz' % VERSION,
    
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Education',
        ],
    
    package_dir = {'': 'python'},
    packages=['dlcpar',
              'dlcpar.deps.rasmus', 'dlcpar.deps.rasmus.ply',
              'dlcpar.deps.compbio',
              'dlcpar.deps.yjw', 'dlcpar.deps.yjw.bio'],
    py_modules=[],
    scripts=['bin/dlcpar', 'bin/dlcpar_search',
             'bin/dlcoal_to_dlcpar', 'bin/dlcpar_to_dlcoal'],
    ext_modules=[]
    )
