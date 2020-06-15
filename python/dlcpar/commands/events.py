"""
Infer event counts in a reconciliation
"""

# python libraries
import os, sys
import argparse

# dlcpar libraries
import dlcpar
from dlcpar import commands
from dlcpar import reconlib

# rasmus, compbio libraries
from rasmus import treelib, tablelib, util
from compbio import phylo, phyloDLC

#==========================================================

VERSION = dlcpar.PROGRAM_VERSION_TEXT

def _get_headers(split_coals):
    if not split_coals:
        headers = ["genes", "dup", "loss", "coal", "appear"]
    else:
        headers = ["genes", "dup", "loss", "coalspec", "coaldup", "appear"]
    return headers

def _events_3t_all(treefiles, stree, gene2species,
                   inputext, implied, split_coals, locus_lca):

    coal_trees = []
    extras = []

    for treefile in treefiles:
        prefix = util.replace_ext(treefile, inputext, "")
        coal_tree, extra = phyloDLC.read_dlcoal_recon(prefix, stree)
        coal_trees.append(coal_tree)
        extras.append(extra)

    etree = phyloDLC.count_dup_loss_coal_trees(coal_trees, extras, stree, gene2species,
                                               implied=implied,
                                               split_coals=split_coals,
                                               locus_lca=locus_lca)

    # make table
    headers = _get_headers(split_coals)
    ptable = treelib.tree2parent_table(etree, headers)

    # sort by post order
    lookup = util.list2lookup(x.name for x in stree.postorder())
    ptable.sort(key=lambda x: lookup[x[0]])

    ptable = [[str(row[0]), str(row[1]), float(row[2])] + row[3:]
              for row in ptable]

    tab = tablelib.Table(ptable,
                         headers=["nodeid", "parentid", "dist"] + headers)
    tab.write()

    return 0


def _events_3t_by_fam(treefiles, stree, gene2species,
                      inputext, implied, split_coals, locus_lca, use_famid):

    # write header
    lookup = util.list2lookup(x.name for x in stree.postorder())
    headers = _get_headers(split_coals)
    print "\t".join(["famid", "nodeid", "parentid", "dist"] + headers)

    for treefile in treefiles:
        if use_famid:
            famid = os.path.basename(os.path.dirname(treefile))
        else:
            famid = treefile

        # read files and events
        prefix = util.replace_ext(treefile, inputext, "")
        coal_tree, extra = phyloDLC.read_dlcoal_recon(prefix, stree)

        etree = phyloDLC.count_dup_loss_coal_trees([coal_tree], [extra], stree, gene2species,
                                                   implied=implied,
                                                   split_coals=split_coals,
                                                   locus_lca=locus_lca)
        ptable = treelib.tree2parent_table(etree, headers)

        # sort by post order
        ptable.sort(key=lambda x: lookup[x[0]])

        # write table
        for row in ptable:
            print "\t".join(map(str, [famid] + row))

    return 0

def _events_lct_all(treefiles, stree, gene2species,
                    inputext, implied, split_coals):

    gene_trees = []
    extras = []

    for treefile in treefiles:
        prefix = util.replace_ext(treefile, inputext, "")
        gene_tree, extra = reconlib.read_labeled_recon(prefix, stree)
        gene_trees.append(gene_tree)
        extras.append(extra)

    etree = reconlib.count_dup_loss_coal_trees(gene_trees, extras, stree, gene2species,
                                               implied=implied, split_coals=split_coals)

    # make
    headers = _get_headers(split_coals)
    ptable = treelib.tree2parent_table(etree, headers)

    # sort by post order
    lookup = util.list2lookup(x.name for x in stree.postorder())
    ptable.sort(key=lambda x: lookup[x[0]])

    ptable = [[str(row[0]), str(row[1]), float(row[2])] + row[3:]
              for row in ptable]

    tab = tablelib.Table(ptable,
                         headers=["nodeid", "parentid", "dist"] + headers)
    tab.write()

    return 0


def _events_lct_by_fam(treefiles, stree, gene2species,
                       inputext, implied, split_coals, use_famid):

    # write header
    lookup = util.list2lookup(x.name for x in stree.postorder())
    headers = _get_headers(split_coals)
    print "\t".join(["famid", "nodeid", "parentid", "dist"] + headers)

    for treefile in treefiles:
        if use_famid:
            famid = os.path.basename(os.path.dirname(treefile))
        else:
            famid = treefile

        # read files and events
        prefix = util.replace_ext(treefile, inputext, "")
        gene_tree, extra = reconlib.read_labeled_recon(prefix, stree)

        etree = reconlib.count_dup_loss_coal_trees([gene_tree], [extra], stree, gene2species,
                                                   implied=implied, split_coals=split_coals)
        ptable = treelib.tree2parent_table(etree, headers)

        # sort by post order
        ptable.sort(key=lambda x: lookup[x[0]])

        # write table
        for row in ptable:
            print "\t".join(map(str, [famid] + row))

    return 0

def run():
    """main program"""

    #=============================
    # parser

    parser = argparse.ArgumentParser(
        usage="%(prog)s events [options] <gene tree> ...",
        description="Count events in reconciliations.",
        formatter_class=commands.CustomHelpFormatter,
        add_help=False)
    parser.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)

    parser.add_argument("treefiles", nargs="+", help=argparse.SUPPRESS)

    grp_format = parser.add_argument_group("Format")
    grp_format_choices = grp_format.add_mutually_exclusive_group(required=True)
    grp_format_choices.add_argument("--3t", dest="threetree",
                                    action="store_true",
                                    help="use Three-Tree reconciliation")
    grp_format_choices.add_argument("--lct", dest="lct",
                                    action="store_true",
                                    help="use LCT reconciliation")

    grp_io = parser.add_argument_group("Input/Output")
    grp_io.add_argument("-s", "--stree", dest="stree",
                        metavar="<species tree>",
                        required=True,
                        help="species tree file in newick format")
    grp_io.add_argument("-S", "--smap", dest="smap",
                        metavar="<species map>",
                        required=True,
                        help="gene to species map")

    grp_ext = parser.add_argument_group("File Extensions")
    grp_ext.add_argument("-I","--inputext", dest="inext",
                         metavar="<input file extension>",
                         help="input file extension")

    grp_misc = parser.add_argument_group("Miscellaneous")
    grp_misc.add_argument("--split-coals", dest="split_coals",
                          action="store_true",
                          help="set to split coalescence at speciation vs duplication")
    grp_misc.add_argument("--by-fam", dest="by_fam",
                          action="store_true",
                          help="separate counts by family")
    grp_misc.add_argument("--use-famid", dest="use_famid",
                          action="store_true",
                          help="lookup famid using directory name")
    grp_misc.add_argument("--use-locus-recon", dest="use_locus_recon",
                          action="store_true", default=False,
                          help="set to use locus recon file rather than LCA [only for --3t]")
    grp_misc.add_argument("--explicit", dest="explicit",
                          action="store_true", default=False,
                          help="set to ignore extra lineages at implied speciation nodes")

    args = parser.parse_args(sys.argv[2:])

    #=============================
    # check arguments

    treefiles = args.treefiles

    #=============================
    # process

    stree = treelib.read_tree(args.stree)
    gene2species = phylo.read_gene2species(args.smap)

    if args.threetree:
        inputext = args.inext
        if inputext is None:
            inputext = ".coal.tree"

        if not args.by_fam:
            _events_3t_all(treefiles, stree, gene2species,
                           inputext,
                           implied=not args.explicit,
                           split_coals=args.split_coals,
                           locus_lca=not args.use_locus_recon)
        else:
            _events_3t_by_fam(treefiles, stree, gene2species,
                              inputext,
                              implied=not args.explicit,
                              split_coals=args.split_coals,
                              locus_lca=not args.use_locus_recon,
                              use_famid=args.use_famid)
    elif args.lct:
        inputext = args.inext
        if inputext is None:
            inputext = ".lct.tree"

        if not args.by_fam:
            _events_lct_all(treefiles, stree, gene2species,
                            inputext,
                            implied=not args.explicit,
                            split_coals=args.split_coals)

        else:
            _events_lct_by_fam(treefiles, stree, gene2species,
                               inputext,
                               implied=not args.explicit,
                               split_coals=args.split_coals,
                               use_famid=args.use_famid)

    else:
        parser.error("--3t or --lct required")
