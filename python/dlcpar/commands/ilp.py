"""
Solve the MPR problem using integer linear programming
"""

# python libraries
import os, sys, shutil
import argparse
import time
import gzip

# dlcpar libraries
import dlcpar
from dlcpar import common
from dlcpar import commands
from dlcpar import reconlib
import dlcpar.ilprecon

# rasmus, compbio libraries
from rasmus import treelib, util
from compbio import phylo

#==========================================================

VERSION = dlcpar.PROGRAM_VERSION_TEXT

def run():
    """main program"""

    #=============================
    # parser

    parser = argparse.ArgumentParser(
        usage="%(prog)s ilp [options] <gene tree> ...",
        description="Solve the MPR problem using integer linear programming.",
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
#    grp_io.add_argument("--lmap", dest="lmap",
#                        metavar="<locus map>",
#                        help="gene to locus map (species-specific)")

    grp_ext = parser.add_argument_group("File Extensions")
    grp_ext.add_argument("-I","--inputext", dest="inext",
                         metavar="<input file extension>",
                         default="",
                         help="input file extension")
    grp_ext.add_argument("-O", "--outputext", dest="outext",
                         metavar="<output file extension>",
                         default=".dlcilp",
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

    grp_solver = parser.add_argument_group("Solver")
    grp_solver.add_argument("--solver", dest="solver",
                            choices=["CBC_CMD", "CPLEX_PY"], default="CBC_CMD",
                            help="ILP solver")
    grp_solver.add_argument("-t", "--time_limit", dest="time_limit",
                            metavar="<time limit>",
                            type=float, default=None,
                            help="ILP solver time limit in seconds")
    grp_solver.add_argument("-m", "--mem_limit", dest="mem_limit",
                            metavar="<memory limit>",
                            type=float, default=None,
                            help="ILP solver memory limit in MB")

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

    # positive time limit
    if (args.time_limit is not None) and (args.time_limit <= 0):
        parser.error("-t/--time_limit must be positive")
    if (args.mem_limit is not None) and (args.mem_limit <= 0):
        parser.error("-m/--mem_limit must be positive")

    #=============================
    # process

    # read species tree and species map
    stree = treelib.read_tree(args.stree)
    common.check_tree(stree, args.stree)
    gene2species = phylo.read_gene2species(args.smap)

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

            # create temp directory for pulp
            out_tmp = out + ".pulp"
            if os.path.exists(out_tmp):
                shutil.rmtree(out_tmp)
            os.mkdir(out_tmp)
        else:
            out_log = common.NullLog()
            out_tmp = None

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
        coal_tree_top = coal_tree.copy()
        for node in coal_tree_top:
            if "boot" in node.data:
                del node.data["boot"]
            node.dist = 0
        coal_tree_top.default_data.clear()

        # perform reconciliation
        return_vals = dlcpar.ilprecon.dlc_recon(
            coal_tree, stree, gene2species,
            dupcost=args.dupcost, losscost=args.losscost,
            coalcost=args.coalcost, coaldupcost=args.coaldupcost,
            implied=True, delay=False,
            solver=args.solver, seed=args.seed, time_limit=args.time_limit, mem_limit=args.mem_limit,
            log=out_log, info_log=out_info, tmp=out_tmp)

        # write info
        gene_tree, labeled_recon, optimal_cost, \
            runtime, runtime_setup, runtime_solve = return_vals
        out_info.write("Runtime:\t%f sec\n" % runtime)
        out_info.write("Runtime (Setup):\t%f sec\n" % runtime_setup)
        out_info.write("Runtime (Solve):\t%f sec\n" % runtime_solve)
        out_info.write("Optimal Cost:\t%f\n" % optimal_cost)

        if args.output_format == "3t":
            coal_tree, recon = reconlib.labeledrecon_to_recon(gene_tree, labeled_recon, stree)
            recon.write(out, coal_tree, filestreams=filestreams)
        else:
            labeled_recon.write(out, gene_tree, filestreams=filestreams)

        # end log
        if args.log:
            out_log.close()

        # end outputs
        for stream in filestreams.itervalues():
            stream.close()

        # end info
        out_info.close()
