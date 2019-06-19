"""
Infer events in a reconciliation
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

def _events_3t_all(options, treefiles, exts):

    stree = treelib.read_tree(options.stree)
    gene2species = phylo.read_gene2species(options.smap)

    coal_trees = []
    extras = []

    for treefile in treefiles:
        prefix = util.replace_ext(treefile, options.inputext, "")
        coal_tree, extra = phyloDLC.read_dlcoal_recon(prefix, stree, exts)
        coal_trees.append(coal_tree)
        extras.append(extra)

    etree = phyloDLC.count_dup_loss_coal_trees(coal_trees, extras, stree, gene2species,
                                               implied=not options.explicit,
                                               locus_mpr=not options.use_locus_recon)

    # make table
    headers = ["genes", "dup", "loss", "coal", "appear"]
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


def _events_3t_by_fam(options, treefiles, exts):

    stree = treelib.read_tree(options.stree)
    gene2species = phylo.read_gene2species(options.smap)

    # write header
    lookup = util.list2lookup(x.name for x in stree.postorder())
    headers = ["genes", "dup", "loss", "coal", "appear"]
    print "\t".join(["famid", "nodeid", "parentid", "dist"] + headers)

    for treefile in treefiles:
        if options.use_famid:
            famid = os.path.basename(os.path.dirname(treefile))
        else:
            famid = treefile

        # read files and events
        prefix = util.replace_ext(treefile, options.inputext, "")
        coal_tree, extra = phyloDLC.read_dlcoal_recon(prefix, stree, exts)

        etree = phyloDLC.count_dup_loss_coal_trees([coal_tree], [extra], stree, gene2species,
                                                   implied=not options.explicit,
                                                   locus_mpr=not options.use_locus_recon)
        ptable = treelib.tree2parent_table(etree, headers)

        # sort by post order
        ptable.sort(key=lambda x: lookup[x[0]])

        # write table
        for row in ptable:
            print "\t".join(map(str, [famid] + row))

    return 0

def _events_lct_all(options, treefiles):

    stree = treelib.read_tree(options.stree)
    gene2species = phylo.read_gene2species(options.smap)
    if options.lmap:
        gene2locus = phylo.read_gene2species(options.lmap)
    else:
        gene2locus = None

    gene_trees = []
    extras = []

    for treefile in treefiles:
        prefix = util.replace_ext(treefile, options.inputext, "")
        prefix = prefix + ".lct"
        gene_tree, extra = reconlib.read_labeled_recon(prefix, stree)
        gene_trees.append(gene_tree)
        extras.append(extra)

    etree = reconlib.count_dup_loss_coal_trees(gene_trees, extras, stree, gene2species,
                                               gene2locus, implied=not options.explicit)

    # make table
    headers = ["genes", "dup", "loss", "coal", "appear"]
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


def _events_lct_by_fam(options, treefiles):

    stree = treelib.read_tree(options.stree)
    gene2species = phylo.read_gene2species(options.smap)
    if options.lmap:
        gene2locus = phylo.read_gene2species(options.lmap)
    else:
        gene2locus = None

    # write header
    lookup = util.list2lookup(x.name for x in stree.postorder())
    headers = ["genes", "dup", "loss", "coal", "appear"]
    print "\t".join(["famid", "nodeid", "parentid", "dist"] + headers)

    for treefile in treefiles:
        if options.use_famid:
            famid = os.path.basename(os.path.dirname(treefile))
        else:
            famid = treefile

        # read files and events
        prefix = util.replace_ext(treefile, options.inputext, "")
        prefix = prefix + ".lct"
        gene_tree, extra = reconlib.read_labeled_recon(prefix, stree)

        etree = reconlib.count_dup_loss_coal_trees([gene_tree], [extra], stree, gene2species,
                                                   gene2locus, implied=not options.explicit)
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
                        required=True,
                        help="gene to species map")
    grp_io.add_argument("--lmap", dest="lmap",
                        metavar="<locus map>",
                        help="gene to locus map (species-specific) [only for lct format]")

    grp_ext = parser.add_argument_group("File Extensions")
    grp_ext.add_argument("-I","--inputext", dest="inputext",
                         metavar="<input file extension>",
                         help="input file extension")

    grp_misc = parser.add_argument_group("Miscellaneous")
    grp_misc.add_argument("--by-fam", dest="by_fam",
                          action="store_true",
                          help="")
    grp_misc.add_argument("--use-famid", dest="use_famid",
                          action="store_true",
                          help="")
    grp_misc.add_argument("--use-locus_recon", dest="use_locus_recon",
                          action="store_true", default=False,
                          help="if set, use locus recon rather than MPR [only for 3t format]")
    grp_misc.add_argument("--explicit", dest="explicit",
                          action="store_true", default=False,
                          help="set to ignore extra lineages at implied speciation nodes")

    args = parser.parse_args(sys.argv[2:])

    #=============================
    # check arguments

    treefiles = args.treefiles

    #=============================
    # process

    if args.format == "3t":
        exts = {"coal_tree": ".coal.tree",
                "coal_recon": ".coal.recon",
                "locus_tree": ".locus.tree",
                "locus_recon": ".locus.recon",
                "daughters": ".daughters"}
        if not args.by_fam:
            _events_3t_all(args, treefiles, exts)
        else:
            _events_3t_by_fam(args, treefiles, exts)
    else:
        if not args.by_fam:
            _events_lct_all(args, treefiles)

        else:
            _events_lct_by_famfam(args, treefiles)