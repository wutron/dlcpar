#!/usr/bin/env python
"""
make-pkg.py
Package dlcpar for release
"""

import os
import re
import shutil
import sys
from subprocess import call, Popen, PIPE

#=============================================================================

FILES = ["bin",
         "examples",
         "python/dlcpar",

         "setup.py",
         "README.md",
         "MANUAL.md",
         "EXAMPLES.md",
         "LICENSE",
         "CHANGES"]

EXCLUDE = ["python/dlcpar/commands/view_lct.py",
           "python/dlcpar/vis/reconsvg.py",
           "examples/getexample.sh",
           ".*.pyc"]

INCLUDE = []

#=============================================================================

def main():
    """main function"""
    pkgdir = sys.argv[1]

    if os.path.exists(pkgdir):
        shutil.rmtree(pkgdir)

    exclude_expr = "|".join(EXCLUDE)

    p = Popen(["find", "-L"] + FILES, stdout=PIPE)
    all_files = [x.rstrip("\n") for x in p.stdout.readlines()]
    all_files = [x for x in all_files
                 if not re.match(exclude_expr, x)] + INCLUDE

    for fn in all_files:
        dest = os.path.join(pkgdir, fn)
        destdir = os.path.dirname(dest)
        print fn, "-->", dest

        if os.path.isfile(fn):
            # copy file
            if not os.path.exists(destdir):
                os.makedirs(destdir)
            shutil.copy(fn, dest)
        else:
            # make dir
            if not os.path.exists(dest):
                os.makedirs(dest)

    # tar
    basename = os.path.basename(pkgdir)
    print ' '.join(
        ["tar", "-C", os.path.dirname(pkgdir), "-zcvf",
         pkgdir + ".tar.gz", basename])
    call(
        ["tar", "-C", os.path.dirname(pkgdir), "-zcvf",
         pkgdir + ".tar.gz", basename])


if __name__ == "__main__":
    sys.exit(main())
