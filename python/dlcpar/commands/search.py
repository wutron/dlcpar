"""
Solve the MPR problem using heuristic search
"""

# python libraries
import argparse
import gzip
import os
import sys
import random
import time

# dlcpar libraries
import dlcpar
from dlcpar import common
from dlcpar import commands
from dlcpar import constants
from dlcpar import reconsearch

# rasmus, compbio libraries
from rasmus import treelib
from rasmus import util
from compbio import phylo
from compbio import phyloDLC

#==========================================================

VERSION = dlcpar.PROGRAM_VERSION_TEXT

def run():
    """main program"""

    #=============================
    # parser

    parser = argparse.ArgumentParser(
        usage="%(prog)s dp [options] <gene tree> ...",
        description="Solve the MPR problem using dynamic programming.",
        formatter_class=commands.CustomHelpFormatter,
        add_help=False)
    parser.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)

    parser.add_argument("treefiles", nargs="+", help=argparse.SUPPRESS)

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
    grp_ext.add_argument("-I", "--inputext", dest="inext",
                         metavar="<input file extension>",
                         default="",
                         help="input file extension")
    grp_ext.add_argument("-O", "--outputext", dest="outext",
                         metavar="<output file extension>",
                         default=".dlcsearch",
                         help="output file extension")

    grp_costs = parser.add_argument_group("Costs")
    grp_costs.add_argument("-D", "--dupcost", dest="dupcost",
                           metavar="<dup cost>",
                           type=float, default=constants.DEFAULT_DUP_COST,
                           help="duplication cost")
    grp_costs.add_argument("-L", "--losscost", dest="losscost",
                           metavar="<loss cost>",
                           type=float, default=constants.DEFAULT_LOSS_COST,
                           help="loss cost")
    grp_costs.add_argument("-C", "--coalcost", dest="coalcost",
                           metavar="<coal cost>",
                           type=float, default=constants.DEFAULT_COAL_COST,
                           help="deep coalescence cost")

    grp_search = parser.add_argument_group("Search")
    grp_search.add_argument("-i", "--iter", dest="iter",
                            metavar="<# iterations>",
                            type=int, default=10,
                            help="number of search iterations")
    grp_search.add_argument("--nprescreen", dest="nprescreen",
                            metavar="<# prescreens>",
                            type=int, default=20,
                            help="number of prescreening iterations")
    grp_search.add_argument("--nconverge", dest="nconverge",
                            metavar="<# converge>",
                            type=int, default=None,
                            help="set to stop search after convergence -- " \
                                + "a solution has converged if it has not changed " \
                                + "for the specified number of iterations")
    grp_search.add_argument("--init-locus-tree", dest="init_locus_tree",
                            metavar="<tree file>",
                            help="initial locus tree for search")
    grp_search.add_argument("-x", "--seed", dest="seed",
                            metavar="<random seed>",
                            type=int, default=None,
                            help="random number seed")

    grp_info = parser.add_argument_group("Information")
    grp_info.add_argument("-l", "--log", dest="log",
                          action="store_true",
                          help="set to output debugging log")

    args = parser.parse_args(sys.argv[2:])

    #=============================
    # check arguments

    treefiles = args.treefiles

    # positive costs
    if args.dupcost <= 0:
        parser.error("-D/--dupcost must be positive")
    if args.losscost <= 0:
        parser.error("-L/--losscost must be positive")
    if args.coalcost <= 0:
        parser.error("-C/--coalcost must be positive")

    # search options
    if args.iter < 1:
        parser.error("-i/--iter must be a positive integer")
    if args.nprescreen < 1:
        parser.error("--nprescreen must be a positive integer")
    if args.nconverge:
        if args.nconverge < 1:
            parser.error("--nconverge must be a positive integer")
        if args.nconverge > args.iter:
            parser.error("--nconverge must be less than or equal to -i/--iter")

    #=============================
    # process

    # read species tree and species map
    stree = treelib.read_tree(args.stree)
    common.check_tree(stree, args.stree)
    gene2species = phylo.read_gene2species(args.smap)

    # random seed if not specified
    if args.seed is None:
        args.seed = int(time.time() * 100)

    # initial locus tree
    if args.init_locus_tree:
        init_locus_tree = treelib.read_tree(args.init_locus_tree)
    else:
        init_locus_tree = None

    # process genes trees
    for treefile in treefiles:
        # general output path
        out = util.replace_ext(treefile, args.inext, args.outext)

        # info file
        out_info = util.open_stream(out + ".info", 'w')

        # log file
        if args.log:
            out_log = gzip.open(out + ".log.gz", 'w')
        else:
            out_log = common.NullLog()

        # command
        cmd = "%s %s" % (os.path.basename(sys.argv[0]),
                         ' '.join(map(lambda x: x if x.find(' ') == -1 else "\"%s\"" % x,
                                      sys.argv[1:])))
        out_info.write("Version:\t%s\n" % VERSION)
        out_info.write("Command:\t%s\n\n" % cmd)
        out_info.write("Seed:\t%d\n\n" % args.seed)
        out_log.write("DLCpar version: %s\n" % VERSION)
        out_log.write("DLCpar executed with the following arguments:\n")
        out_log.write("%s\n\n" % cmd)
        out_log.write("Seed:\t%d\n\n" % args.seed)

        # read and prepare coal tree (multiple coal trees not supported)
        coal_trees = list(treelib.iter_trees(treefile))
        if len(coal_trees) > 1:
            raise Exception("unsupported: multiple coal trees per file")

        # get coal tree
        coal_tree = coal_trees[0]
        common.check_tree(coal_tree, treefile)

        # remove bootstrap and distances if they exist
        coal_tree_top = common.prepare_tree(coal_tree)

        # set random seed
        random.seed(args.seed)

        # perform reconciliation
        maxrecon, optimal_cost, runtime = reconsearch.dlc_recon(
            coal_tree_top, stree, gene2species,
            dupcost=args.dupcost, losscost=args.losscost, coalcost=args.coalcost,
            nsearch=args.iter, nprescreen=args.nprescreen, nconverge=args.nconverge,
            init_locus_tree=init_locus_tree,
            log=out_log)

        # write info
        out_info.write("Feasibility:\tfeasible\n")
        out_info.write("Runtime:\t%f sec\n" % runtime)
        out_info.write("Optimal Cost:\t%f\n" % optimal_cost)

        # write outputs
        phyloDLC.write_dlcoal_recon(out, coal_tree, maxrecon)

        # end logging
        if args.log:
            out_log.close()

        # end info
        out_info.close()
