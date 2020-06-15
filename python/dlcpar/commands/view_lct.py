"""
View the labeled coalescent tree
"""

# python libraries
import argparse
import os
import sys

# dlcpar libraries
import dlcpar
from dlcpar import reconlib
from dlcpar import commands
from dlcpar.vis import reconsvg

# rasmus libraries
from rasmus import treelib
from rasmus import util

#=============================================================================

VERSION = dlcpar.PROGRAM_VERSION_TEXT

def run():
    """main program"""

    #======================
    # parser

    parser = argparse.ArgumentParser(
        usage="%(prog)s view_lct [options] <gene tree>",
        description="View the labeled coalescent tree.",
        formatter_class=commands.CustomHelpFormatter,
        add_help=False)
    parser.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)

    parser.add_argument("tree", help=argparse.SUPPRESS)

    grp_io = parser.add_argument_group("Input/Output")
    grp_io.add_argument("-s", "--stree", dest="stree",
                        metavar="<species tree>",
                        required=True,
                        help="species tree in newick format")
    grp_io.add_argument("-v", "--viewer", dest="viewer",
                        metavar="<svg viewer>",
                        default="display",
                        help="svg viewer")
    grp_io.add_argument("-o", "--output", dest="output",
                        metavar="<output>",
                        help="output file")

    grp_misc = parser.add_argument_group("Miscellaneous")
    grp_misc.add_argument("--xscale", dest="xscale",
                          metavar="<x-scaling>",
                          type=float, default=50,
                          help="x-scale factor")
    grp_misc.add_argument("--yscale", dest="yscale",
                          metavar="<y-scaling>",
                          type=float, default=50,
                          help="y-scale factor")
    grp_misc.add_argument("--names", dest="names",
                          action="store_true",
                          help="display internal node names")
    grp_misc.add_argument("--snames", dest="snames",
                          action="store_true",
                          help="display species names")

    args = parser.parse_args(sys.argv[2:])

    #=============================
    # process

    stree = treelib.read_tree(args.stree)
    inputext = ".lct.tree"
    prefix = util.replace_ext(args.tree, inputext, "")
    gtree, extra = reconlib.read_labeled_recon(prefix, stree)

    # output
    if args.output:
        stream = args.output
    else:
        stream = os.popen(args.viewer, "w")

    # label node names
    labels = {}
    if args.names:
        for node in gtree:
            if not node.is_leaf():
                labels[node.name] = "[%s]" % node.name

    # label species names
    slabels = {}
    if args.snames:
        for snode in stree:
            if not snode.is_leaf():
                slabels[snode.name] = "%s" % snode.name

    # draw tree
    reconsvg.draw_tree(gtree, stree, extra, filename=stream,
                       xscale=args.xscale, yscale=args.yscale,
                       labels=labels,
                       slabels=slabels)

    # close output
    if not args.output:
        stream.close()
