"""
Check for equality between reconciliation structures
"""

# python libraries
import argparse
import sys

# dlcpar libraries
import dlcpar
from dlcpar import commands
from dlcpar import reconlib

# rasmus, compbio libraries
from rasmus import treelib
from rasmus import util
from compbio import phylo
from compbio import phyloDLC

#==========================================================

VERSION = dlcpar.PROGRAM_VERSION_TEXT

def _equal_3tree(prefix1, prefix2, stree, gene2species,
                 use_locus_recon=True):
    """Check for equality between two 3-tree reconciliation structures"""

    #======================================
    # utilities

    def remap(tree):
        """remap internal node names based on leaves"""
        leaves = {}
        names = {}
        names_set = set()
        for node in tree.postorder():
            if node.is_leaf():
                leaves[node] = [node.name]
            else:
                leaves[node] = []
                for child in node.children:
                    leaves[node].extend(leaves[child])
                leaves[node].sort()
            names[node] = tree.unique_name(str(tuple(leaves[node])), names_set)

        for node in list(tree):
            tree.rename(node.name, names[node])

    #======================================
    # main

    recon1 = phyloDLC.Recon()
    coal_tree1, _ = recon1.read(prefix1, stree)
    if not use_locus_recon:
        recon1.locus_recon = phylo.reconcile(recon1.locus_tree, stree, gene2species)
        recon1.locus_events = phylo.label_events(recon1.locus_tree, recon1.locus_recon)
        recon1.daughters = [node for node in recon1.daughters
                            if recon1.locus_events[node.parent] == "dup"]

    recon2 = phyloDLC.Recon()
    coal_tree2, _ = recon2.read(prefix2, stree)
    if not use_locus_recon:
        recon2.locus_recon = phylo.reconcile(recon2.locus_tree, stree, gene2species)
        recon2.locus_events = phylo.label_events(recon2.locus_tree, recon2.locus_recon)
        recon2.daughters = [node for node in recon2.daughters
                            if recon2.locus_events[node.parent] == "dup"]

    # compare
    hash1 = phylo.hash_tree(coal_tree1)
    hash2 = phylo.hash_tree(coal_tree2)

    if hash1 != hash2:
        print >>sys.stderr, "coal tree mismatch"
        return False

    # process coal trees
    remap(coal_tree1)
    remap(coal_tree2)

    # compare
    return recon1 == recon2


def _equal_lct(prefix1, prefix2, stree):
    """Check for equality between two LCT reconciliation structures"""

    #======================================
    # utilities

    def remove(tree, recon):
        """remove unnecessary nodes from gene tree and recon"""
        # e.g. a node that ...
        # 1) has a single child, and either
        # 2a) is the root, or
        # 2b) is in the same species as its parent, and
        #     is in the same locus as its parent
        species_map = recon.species_map
        locus_map = recon.locus_map
        #order = recon.order

        for node in list(tree.preorder()):
            parent = node.parent
            children = node.children

            if len(children) == 1:
                if ((not parent)
                    or (parent
                        and species_map[node] == species_map[parent]
                        and locus_map[node] == locus_map[parent]
                       )
                   ):
                    # update recon
                    # these nodes are leaves and so are never in order
                    # if the node has the same locus as its parent
                    del locus_map[node]
                    del species_map[node]

                    # update gene tree
                    phylo.remove_spec_node(node, tree)

    def remap(tree):
        """remap internal node names based on leaves"""
        leaves = {}
        names = {}
        names_set = set()
        for node in tree.postorder():
            if node.is_leaf():
                leaves[node] = [node.name]
            else:
                leaves[node] = []
                for child in node.children:
                    leaves[node].extend(leaves[child])
                leaves[node].sort()
            names[node] = tree.unique_name(str(tuple(leaves[node])), names_set)

        for node in list(tree):
            tree.rename(node.name, names[node])

    #======================================
    # main

    # read files
    recon1 = reconlib.LabeledRecon()
    gene_tree1, _ = recon1.read(prefix1, stree)

    recon2 = reconlib.LabeledRecon()
    gene_tree2, _ = recon2.read(prefix2, stree)

    # compare
    tree1 = gene_tree1.copy()
    treelib.remove_single_children(tree1)

    tree2 = gene_tree2.copy()
    treelib.remove_single_children(tree2)

    hash1 = phylo.hash_tree(tree1)
    hash2 = phylo.hash_tree(tree2)

    if hash1 != hash2:
        print >>sys.stderr, "gene tree mismatch"
        return False

    # process trees and recons
    remove(gene_tree1, recon1)
    remove(gene_tree2, recon2)

    remap(gene_tree1)
    remap(gene_tree2)

    # compare
    return recon1 == recon2


def run():
    """main program"""

    #=============================
    # parser

    parser = argparse.ArgumentParser(
        usage="%(prog)s equal [options] <tree 1> <tree 2>",
        description="Check for equality of reconciliation structures.",
        formatter_class=commands.CustomHelpFormatter,
        add_help=False)
    parser.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)

    parser.add_argument("tree1", help=argparse.SUPPRESS)
    parser.add_argument("tree2", help=argparse.SUPPRESS)

    grp_format = parser.add_argument_group("Format")
    grp_format_choices = grp_format.add_mutually_exclusive_group(required=True)
    grp_format_choices.add_argument("--3t", dest="threetree",
                                    action="store_true",
                                    help="compare Three-Tree reconciliations")
    grp_format_choices.add_argument("--lct", dest="lct",
                                    action="store_true",
                                    help="compare LCT reconciliations")

    grp_io = parser.add_argument_group("Input/Output")
    grp_io.add_argument("-s", "--stree", dest="stree",
                        metavar="<species tree>",
                        required=True,
                        help="species tree file in newick format")
    grp_io.add_argument("-S", "--smap", dest="smap",
                        metavar="<species map>",
                        help="gene to species map")

    grp_ext = parser.add_argument_group("File Extensions")
    grp_ext.add_argument("-I", "--inputext", dest="inext",
                         metavar="<input file extension>",
                         help="input file extension")

    grp_misc = parser.add_argument_group("Miscellaneous [only used if --3t]")
    grp_misc.add_argument("--use-locus-lca", dest="use_locus_recon",
                          action="store_false",
                          help="set to use LCA rather than locus recon file")

    args = parser.parse_args(sys.argv[2:])

    #=============================
    # check arguments

    # default extensions
    inputext = args.inext
    if args.threetree:
        if inputext is None:
            inputext = ".coal.tree"
    elif args.lct:
        if inputext is None:
            inputext = ".lct.tree"
    else:
        parser.error("--3t or --lct required")

    if not args.use_locus_recon:
        if not args.smap:
            parser.error("-S/--smap required if --use-locus-lca set")

    #=============================
    # process

    # read files
    stree = treelib.read_tree(args.stree)
    if not args.smap:
        gene2species = None
    else:
        gene2species = phylo.read_gene2species(args.smap)

    prefix1 = util.replace_ext(args.tree1, inputext, "")
    prefix2 = util.replace_ext(args.tree2, inputext, "")

    # compare
    if args.threetree:
        eq = _equal_3tree(prefix1, prefix2, stree, gene2species,
                          use_locus_recon=args.use_locus_recon)
    else:
        eq = _equal_lct(prefix1, prefix2, stree)
    print eq
