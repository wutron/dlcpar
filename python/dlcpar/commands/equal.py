"""
Check for equality between reconciliation structures
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

def _equal_3tree(args):
    """Check for equality between two 3-tree reconciliation structures"""

    #======================================
    # utilities
    def remap(tree):
        # remap internal node names based on leaves
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
    stree = treelib.read_tree(args.stree)
    if args.smap:
        gene2species = phylo.read_gene2species(args.smap)

    recon1 = phyloDLC.Recon()
    coal_tree1, extra1 = recon1.read(args.prefix1, stree)
    if not args.use_locus_recon:
        recon1.locus_recon = phylo.reconcile(recon1.locus_tree, stree, gene2species)
        recon1.locus_events = phylo.label_events(recon1.locus_tree, recon1.locus_recon)
        recon1.daughters = filter(lambda node: recon1.locus_events[node.parent] == "dup", recon1.daughters)

    recon2 = phyloDLC.Recon()
    coal_tree2, extra2 = recon2.read(args.prefix2, stree)
    if not args.use_locus_recon:
        recon2.locus_recon = phylo.reconcile(recon2.locus_tree, stree, gene2species)
        recon2.locus_events = phylo.label_events(recon2.locus_tree, recon2.locus_recon)
        recon2.daughters = filter(lambda node: recon2.locus_events[node.parent] == "dup", recon2.daughters)

    # compare
    hash1 = phylo.hash_tree(coal_tree1)
    hash2 = phylo.hash_tree(coal_tree2)
    if hash1 != hash2:
        print >>sys.stderr, "coal tree mismatch"
        return False
    else:
        # process coal trees
        remap(coal_tree1)
        remap(coal_tree2)

        # compare
        return recon1 == recon2


def _equal_lct(args):
    """Check for equality between two LCT reconciliation structures"""

    #======================================
    # utilities

    def remove(tree, recon):
        # remove unnecessary nodes from gene tree and recon
        # e.g. a node that ...
        # 1) has a single child, and either
        # 2a) is the root, or
        # 2b) is in the same species as its parent, and
        #     is in the same locus as its parent
        species_map = recon.species_map
        locus_map = recon.locus_map
        order = recon.order

        for node in list(tree.preorder()):
            parent = node.parent
            children = node.children

            if len(children) == 1:
                child = children[0]
                if (not parent) or \
                    (parent and \
                     species_map[node] == species_map[parent] and \
                     locus_map[node] == locus_map[parent]):
                    # update recon
                    # these nodes are leaves and so are never in order
                    # if the node has the same locus as its parent
                    del locus_map[node]
                    del species_map[node]

                    # update gene tree
                    phylo.remove_spec_node(node, tree)

    def remap(tree):
        # remap internal node names based on leaves
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
    stree = treelib.read_tree(args.stree)

    recon1 = reconlib.LabeledRecon()
    args.prefix1 = args.prefix1 + ".lct"
    gene_tree1, extra1 = recon1.read(args.prefix1, stree)

    recon2 = reconlib.LabeledRecon()
    args.prefix2 = args.prefix2 + ".lct"
    gene_tree2, extra2 = recon2.read(args.prefix2, stree)

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
    else:
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
        usage="%(prog)s equal [options] <prefix 1> <prefix 2>",
        description="Check for equality of reconciliation structures.",
        formatter_class=commands.CustomHelpFormatter,
        add_help=False)
    parser.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)

    parser.add_argument("prefix1", help=argparse.SUPPRESS)
    parser.add_argument("prefix2", help=argparse.SUPPRESS)

    grp_io = parser.add_argument_group("Input/Output")
    grp_io.add_argument("--format", dest="format",
                        choices=["lct","3t"], default="lct",
                        metavar="{(lct)|3t}",
                        help="specify input format")
    grp_io.add_argument("-s", "--stree", dest="stree",
                        metavar="<species tree>",
                        required=True,
                        help="species tree file in newick format")
    grp_io.add_argument("-S", "--smap", dest="smap",
                        metavar="<species map>",
                        help="gene to species map [only used if format='3t']")
    grp_io.add_argument("--use-locus-mpr", dest="use_locus_recon",
                        default=True, action="store_false",
                        help="if set, use MPR rather than locus recon file [only used if format='3t']")

    args = parser.parse_args(sys.argv[2:])

    #=============================
    # check arguments

    if not args.use_locus_recon:
        if not args.smap:
            parser.error("-S/--smap required if --use-locus-mpr set")

    #=============================
    # process

    if args.format == "3t":
        eq = _equal_3tree(args)
    else:
        eq = _equal_lct(args)
    print eq