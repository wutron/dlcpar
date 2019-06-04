"""
Solve the MPR problem using dynamic programming
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
from dlcpar import reconlib
import dlcpar.recondp

# rasmus, compbio libraries
from rasmus import treelib, util
from compbio import phylo

# numpy libraries
import numpy as np

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
    grp_io.add_argument("--lmap", dest="lmap",
                        metavar="<locus map>",
                        help="gene to locus map (species-specific)")
    grp_io.add_argument("-n", "--nsamples", dest="nsamples",
                        metavar="<number of reconciliations>",
                        type=int, default=1,
                        help="number of uniform random samples")

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
                           help="deep coalescence cost at speciation")
    grp_costs.add_argument("-K", "--coaldupcost", dest="coaldupcost",
                           metavar="<coal dup cost>",
                           type=float,
                           help="deep coalescence cost at duplication if different")

    grp_heur = parser.add_argument_group("Heuristics")
    grp_heur.add_argument("--no_prescreen", dest="prescreen",
                          default=True, action="store_false",
                          help="set to disable prescreen of locus maps")
    grp_heur.add_argument("--prescreen_min", dest="prescreen_min",
                          metavar="<prescreen min>",
                          type=float, default=50,
                          help="prescreen locus maps if min (forward) cost exceeds this value")
    grp_heur.add_argument("--prescreen_factor", dest="prescreen_factor",
                          metavar="<prescreen factor>",
                          type=float, default=10,
                          help="prescreen locus maps if (forward) cost exceeds this factor * min (forward) cost")
    grp_heur.add_argument("--max_loci", dest="max_loci",
                          metavar="<max # of loci>",
                          type=int, default=-1,
                          help="max # of co-existing loci (per species), " +\
                               "set to -1 for no limit")
    grp_heur.add_argument("--max_dups", dest="max_dups",
                          metavar="<max # of dups>",
                          type=int, default=4,
                          help="max # of duplications (per species), " +\
                               "set to -1 for no limit")
    grp_heur.add_argument("--max_losses", dest="max_losses",
                          metavar="<max # of losses>",
                          type=int, default=4,
                          help="max # of losses (per species), " +\
                               "set to -1 for no limit")

    grp_misc = parser.add_argument_group("Miscellaneous")
    grp_misc.add_argument("-x", "--seed", dest="seed",
                          metavar="<random seed>",
                          type=int, default=None,
                          help="random number seed")
    grp_misc.add_argument("--output_format", dest="output_format",
                          choices=["lct","3t"], default="lct",
                          metavar="{(lct)|3t}",
                          help="specify output format")

    grp_info = parser.add_argument_group("Information")
    grp_info.add_argument("-l", "--log", dest="log",
                          action="store_true",
                          help="set to output debugging log")

    args = parser.parse_args(sys.argv[2:])

    #=============================
    # check arguments

    treefiles = args.treefiles

    # required options
    if args.nsamples < 1:
        parser.error("-n/--nsamples must be at least 1")

    # positive costs
    if args.dupcost <= 0:
        parser.error("-D/--dupcost must be positive")
    if args.losscost <= 0:
        parser.error("-L/--losscost must be positive")
    if args.coalcost <= 0:
        parser.error("-C/--coalcost must be positive")
    if args.coaldupcost is None:
        args.coaldupcost = args.coalcost
    elif args.coaldupcost <= 0:
        parser.error("-K/--coaldupcost must be positive")

    # prescreen options
    if args.prescreen_min <= 0:
        parser.error("--prescreen_min must be > 0")
    if args.prescreen_factor < 1:
        parser.error("--prescreen_factor must be >= 1")

    # max loci/dups/losses options
    if args.max_loci == -1:
        args.max_loci = util.INF
    elif args.max_loci < 1:
        parser.error("--max_loci must be > 0")

    if args.max_dups == -1:
        args.max_dups = util.INF
    elif args.max_dups < 1:
        parser.error("--max_dups must be > 0")

    if args.max_losses == -1:
        args.max_losses = util.INF
    elif args.max_losses < 1:
        parser.error("--max_losses must be > 0")

    #=============================
    # process

    # read species tree and species map
    stree = treelib.read_tree(args.stree)
    common.check_tree(stree, args.stree)
    gene2species = phylo.read_gene2species(args.smap)
    if args.lmap:
        gene2locus = phylo.read_gene2species(args.lmap)
    else:
        gene2locus = None

    # random seed if not specified
    if args.seed is None:
        # note: numpy requires 32-bit unsigned integer
        args.seed = int(time.time() * 100) % (2^32-1)

    # process genes trees
    for treefile in treefiles:
        # general output path
        out = util.replace_ext(treefile, args.inext, args.outext)

        # info file
        out_info = util.open_stream(out + ".info", 'w')

        # output streams
        if args.output_format == "3t":
            out_coal_tree   = util.open_stream(out + ".coal.tree", 'w')
            out_coal_recon  = util.open_stream(out + ".coal.recon", 'w')
            out_locus_tree  = util.open_stream(out + ".locus.tree", 'w')
            out_locus_recon = util.open_stream(out + ".locus.recon", 'w')
            out_daughters   =  util.open_stream(out + ".daughters", 'w')
            filestreams = {"coal_tree"  : out_coal_tree,
                           "coal_recon" : out_coal_recon,
                           "locus_tree" : out_locus_tree,
                           "locus_recon": out_locus_recon,
                           "daughters"  : out_daughters}
        else:
            out_tree  = util.open_stream(out + ".lct.tree", 'w')
            out_recon = util.open_stream(out + ".lct.recon", 'w')
            out_order = util.open_stream(out + ".lct.order", 'w')
            filestreams = {"tree" : out_tree,
                           "recon": out_recon,
                           "order": out_order}

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
        out_log.write("DLCpar version: %s\n" % VERSION)
        out_log.write("DLCpar executed with the following arguments:\n")
        out_log.write("%s\n\n" % cmd)

        # read and prepare coal tree
        coal_trees = list(treelib.iter_trees(treefile))

        # multiple coal trees not supported
        if len(coal_trees) > 1:
            raise Exception("unsupported: multiple coal trees per file")

        # get coal tree
        coal_tree = coal_trees[0]
        common.check_tree(coal_tree, treefile)

        # remove bootstrap and distances if they exist
        coal_tree_top = coal_tree.copy()
        for node in coal_tree_top:
            if "boot" in node.data:
                del node.data["boot"]
            node.dist = 0
        coal_tree_top.default_data.clear()

        # sample reconciliations
        for i in xrange(args.nsamples):
            # start info and log for this sample
            if args.nsamples > 1:
                out_info.write("# Solution %d\n" % i)
                out_log.write("# Solution %d\n" % i)

            # set random seed (seed increases by 100 per iteration)
            seed = args.seed + 100 * i
            random.seed(seed)
            np.random.seed(seed)
            out_log.write("seed: %d\n\n" % seed)

            # perform reconciliation (use copy)
            coal_tree_tmp = coal_tree_top.copy()
            return_vals = dlcpar.recondp.dlc_recon(
                coal_tree_tmp, stree, gene2species, gene2locus,
                dupcost=args.dupcost, losscost=args.losscost, coalcost=args.coalcost, coaldupcost=args.coaldupcost,
                implied=True, delay=False,
                prescreen=args.prescreen, prescreen_min=args.prescreen_min, prescreen_factor=args.prescreen_factor,
                max_loci=args.max_loci, max_dups=args.max_dups, max_losses=args.max_losses,
                allow_both=False,
                log=out_log)

            # infeasible solution
            if return_vals[1] is None:
                out_info.write("Feasibility:\tinfeasible\n")
                out_info.write("Runtime:\tnull\n")
                out_info.write("Optimal Cost:\tnull\n")
                out_info.write("Number of Solutions:\tnull")

                # stop run
                break

            # write info
            gene_tree, labeled_recon, nsoln, optimal_cost, runtime = return_vals
            out_info.write("Feasibility:\tfeasible\n")
            out_info.write("Runtime:\t%f sec\n" % runtime)
            out_info.write("Optimal Cost:\t%f\n" % optimal_cost)
            out_info.write("Number of Solutions:\t%d" % nsoln)

            # end info and log for this sample
            if args.nsamples > 1:
                out_info.write("\n\n")
                out_log.write("\n\n")

            # write outputs
            if args.nsamples > 1:
                for stream in filestreams.itervalues():
                    stream.write("# Solution %d\n" % i)

            if args.output_format == "3t":
                coal_tree, recon = reconlib.labeledrecon_to_recon(gene_tree, labeled_recon, stree)
                recon.write(out, coal_tree, filestreams=filestreams)
            else:
                labeled_recon.write(out, gene_tree, filestreams=filestreams)

            if args.nsamples > 1:
                for stream in filestreams.itervalues():
                    stream.write("\n\n")

        # end log
        if args.log:
            out_log.close()

        # end outputs
        for stream in filestreams.itervalues():
            stream.close()

        # end info
        out_info.close()
