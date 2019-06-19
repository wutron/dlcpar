'''
View the labeled reconcilations.
'''
# python libraries
import os, sys
import argparse

# dlcpar libraries
import dlcpar
from dlcpar import reconlib
from dlcpar import commands
from dlcpar.vis import reconsvg

# rasmus libraries
from rasmus import treelib

#=============================================================================
VERSION = dlcpar.PROGRAM_VERSION_TEXT

def run():
    """main program"""
    #======================
    # parser
    parser = argparse.ArgumentParser(
        usage="%(prog)s view_recon [options] <prefix>",
        description="View reconcilations",
        formatter_class=commands.CustomHelpFormatter,
        add_help=False)

    parser.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)

    parser.add_argument("prefix", help=argparse.SUPPRESS)

    parser.add_argument("-s", "--stree", dest="stree",
                  metavar="<species tree>",
                  help="species tree in newick format")
    parser.add_argument("-v", "--viewer", dest="viewer",
                  metavar="<svg viewer>",
                  default="display")
    parser.add_argument("--xscale", dest="xscale",
                  metavar="<x-scaling>",
                  type=float, default=50)
    parser.add_argument("--yscale", dest="yscale",
                  metavar="<y-scaling>",
                  type=float, default=50)
    parser.add_argument("--names", dest="names",
                  action="store_true",
                  default=False,
                  help="display internal node names")
    parser.add_argument("--snames", dest="snames",
                  action="store_true",
                  default=False,
                  help="display species names")
    parser.add_argument("-o", "--output", dest="output",
                  metavar="<output>")
    args = parser.parse_args(sys.argv[2:])

    #=============================
    # process

    stree = treelib.read_tree(args.stree)
    args.prefix = args.prefix + ".lct"
    gtree, extra = reconlib.read_labeled_recon(args.prefix, stree)

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
