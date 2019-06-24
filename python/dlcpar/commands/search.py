"""
Solve the MPR problem using heuristic search
"""

# python libraries
import os, sys
import argparse
import time
import random
import gzip

# dlcpar libraries
import dlcpar
from dlcpar import common
from dlcpar import commands
import dlcpar.reconsearch

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
    grp_ext.add_argument("-I","--inputext", dest="inext",
                         metavar="<input file extension>",
                         default="",
                         help="input file extension")
    grp_ext.add_argument("-O", "--outputext", dest="outext",
                         metavar="<output file extension>",
                         default=".dlcpar",
                         help="output file extension")

    grp_costs = parser.add_argument_group("Costs")
    grp_costs.add_argument("-D", "--dupcost", dest="dupcost",
                           metavar="<dup cost>",
                           type=float, default=1.0,
                           help="duplication cost")
    grp_costs.add_argument("-L", "--losscost", dest="losscost",
                           metavar="<loss cost>",
                           type=float, default=1.0,
                           help="loss cost")
    grp_costs.add_argument("-C", "--coalcost", dest="coalcost",
                           metavar="<coal cost>",
                           type=float, default=0.5,
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
                            help="set to stop search after convergence -- " +
                                 "a solution has converged if it has not changed " +
                                 "for the specified number of iterations")
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
                          help="if given, output debugging log")

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

    # initial locus tree
    if args.init_locus_tree:
        init_locus_tree = treelib.read_tree(args.init_locus_tree)
    else:
        init_locus_tree = None

    # process genes trees
    for treefile in treefiles:

        # general output path
        out = util.replace_ext(treefile, args.inext, args.outext)
        
        # start logging
        if args.log:
            log_out = gzip.open(out + ".log.gz", 'w')
        else:
            log_out = common.NullLog()

        # info file
        out_info = util.open_stream(out + ".info", 'w')

        # command
        cmd = "%s %s" % (os.path.basename(sys.argv[0]),
                         ' '.join(map(lambda x: x if x.find(' ') == -1 else "\"%s\"" % x,
                                      sys.argv[1:])))
        out_info.write("Version:\t%s\n" % VERSION)
        out_info.write("Command:\t%s\n\n" % cmd)
        log_out.write("DLCpar version: %s\n" % VERSION)
        log_out.write("DLCpar executed with the following arguments:\n")
        log_out.write("%s\n\n" % cmd)

        # set random seed
        if args.seed is None:
            args.seed = int(time.time() * 100)
        random.seed(args.seed)
        log_out.write("seed: %d\n" % args.seed)

        # read and prepare coal tree
        coal_trees = list(treelib.iter_trees(treefile))
        locus_trees = []

        for coal_tree in coal_trees:
            common.check_tree(coal_tree, treefile)

            # remove bootstrap and distances if they exist
            for node in coal_tree:
                if "boot" in node.data:
                    del node.data["boot"]
                node.dist = 0
            coal_tree.default_data.clear()

            # perform reconciliation
            maxrecon, runtime = dlcpar.reconsearch.dlc_recon(
                coal_tree, stree, gene2species,
                dupcost=args.dupcost, losscost=args.losscost, coalcost=args.coalcost, implied=True,
                nsearch=args.iter, nprescreen=args.nprescreen, nconverge=args.nconverge,
                init_locus_tree=init_locus_tree,
                log=log_out)

            # write info
            out_info.write("Feasibility:\tfeasible\n")
            out_info.write("Seed: %d\n\n" % args.seed)
            out_info.write("Runtime:\t%f sec\n" % runtime)

            # end info and log for this sample
            out_info.write("\n\n")
            log_out.write("\n\n")

        # make "consensus" reconciliation if multiple coal trees given
        if len(coal_trees) > 1:
            # make consensus locus tree
            coal_tree = phylo.consensus_majority_rule(coal_trees, rooted=True)
            phylo.ensure_binary_tree(coal_tree)
            locus_tree = phylo.consensus_majority_rule(locus_trees, rooted=True)
            phylo.ensure_binary_tree(locus_tree)
            locus_recon = phylo.reconcile(locus_tree, stree, gene2species)
            maxrecon = {
                "coal_recon": phylo.reconcile(coal_tree, locus_tree, lambda x:x),
                "locus_tree": locus_tree,
                "locus_recon": locus_recon,
                "locus_events": phylo.label_events(locus_tree, locus_recon)}
            maxrecon["daughters"] = dlcpar.simplerecon.propose_daughters(
                coal_tree, maxrecon["coal_recon"], locus_tree,
                maxrecon["locus_events"])



        # write outputs
        out = util.replace_ext(treefile, args.inext, args.outext)
        phyloDLC.write_dlcoal_recon(out, coal_tree, maxrecon)

        # end logging
        if args.log:
            log_out.close()

        # end info
        out_info.close()