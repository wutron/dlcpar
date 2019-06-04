"""
Convert between reconciliation structures
"""

# python libraries
import os, sys
import argparse

# dlcpar libraries
import dlcpar
from dlcpar import common
from dlcpar import commands
from dlcpar import reconlib

# rasmus, compbio libraries
from rasmus import treelib, util
from compbio import phylo, phyloDLC

#==========================================================

VERSION = dlcpar.PROGRAM_VERSION_TEXT

def run():
    """main program"""

    #=============================
    # parser

    parser = argparse.ArgumentParser(
        usage="%(prog)s convert [options] <gene tree> ...",
        description="Convert between reconciliation structures.",
        formatter_class=commands.CustomHelpFormatter,
        add_help=False)
    parser.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)

    parser.add_argument("treefiles", nargs="+", help=argparse.SUPPRESS)

    grp_convert = parser.add_argument_group("Conversion")
    grp_convert_choices = grp_convert.add_mutually_exclusive_group()
    grp_convert_choices.add_argument("--3tree_to_lct", dest="threetree_to_lct",
                                     default=False, action="store_true",
                                     help="set to convert from 3-tree to LCT")
    grp_convert_choices.add_argument("--lct_to_3tree", dest="lct_to_threetree",
                                     default=False, action="store_true",
                                     help="set to convert from LCT to 3-Tree")

    grp_io = parser.add_argument_group("Input/Output")
    grp_io.add_argument("-s", "--stree", dest="stree",
                        metavar="<species tree>",
                        required=True,
                        help="species tree file in newick format")
    grp_io.add_argument("-S", "--smap", dest="smap",
                        metavar="<species map>",
                        required=True,
                        help="gene to species map")

    grp_misc = parser.add_argument_group("Miscellaneous (for 3tree_to_lct)")
    grp_misc.add_argument("--use-locus-recon", dest="use_locus_recon",
                          default=False, action="store_true",
                          help="if set, use locus recon file rather than MPR")
    grp_misc.add_argument("--no-delay", dest="delay",
                          default=True, action="store_false",
                          help="if set, disallow duplication between speciation and coalescence")

    args = parser.parse_args(sys.argv[2:])

    #=============================
    # check arguments

    # default extensions
    if args.threetree_to_lct:
        inputext = ".coal.tree"
        outputext = ""
    elif args.lct_to_threetree:
        inputext = ".lct.tree"
        outputext = ""

    treefiles = args.treefiles

    #=============================
    # process

    # read species tree and species map
    stree = treelib.read_tree(args.stree)
    common.check_tree(stree, args.stree)
    gene2species = phylo.read_gene2species(args.smap)

    # process genes trees
    for treefile in treefiles:
        prefix = util.replace_ext(treefile, inputext, "")
        out = util.replace_ext(treefile, inputext, outputext)

        if args.threetree_to_lct:
            # read 3-tree files
            recon = phyloDLC.Recon()
            coal_tree, extra = recon.read(prefix, stree)

            # convert
            gene_tree, labeled_recon = \
                reconlib.recon_to_labeledrecon(coal_tree, recon, stree, gene2species,
                                               locus_mpr=not args.use_locus_recon,
                                               delay=args.delay)

            # output
            labeled_recon.write(out, gene_tree)

        elif args.lct_to_threetree:
            # read lct files
            labeledrecon = reconlib.LabeledRecon()
            gene_tree, extra = labeledrecon.read(prefix, stree)

            # convert
            coal_tree, recon = reconlib.labeledrecon_to_recon(gene_tree, labeledrecon, stree)

            # output
            recon.write(out, coal_tree)