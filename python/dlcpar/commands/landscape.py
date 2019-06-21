'''Find MPR landscapes across range of event counts'''
# python libraries
import os, sys
import argparse
import time
import random
import gzip

# geometry libraries
from shapely import geometry

# dlcpar libraries
import dlcpar
from dlcpar import common
from dlcpar import reconlib
from dlcpar import commands
import dlcpar.reconscape

# rasmus, compbio libraries
from rasmus import treelib, util
from compbio import phylo

# numpy libraries
import numpy as np

#==========================================================
# parser

VERSION = dlcpar.PROGRAM_VERSION_TEXT

class ParseRangeAction(argparse.Action):
        # callback to parse cost ranges
    def __call__(self,parser,namespace,values,options_string=None):
        try:
            if "-" not in values:
                a = float(values)
                b = 1/a
                if a < b:
                    low, high = a, b
                else:
                    low, high = b, a
            else:
                low, high = map(float, values.split("-"))
            assert (low > 0) and (high > 0) and (low < high)
        except:
            raise argparse.ArgumentError("invalid %s: %s" % (option, value))
        setattr(namespace, self.dest, (low, high))


def run():
    '''main program'''

    #=============================
    # parser

    parser = argparse.ArgumentParser(
        usage = "%(prog)s landscape [options] <gene tree> ...",

        description =
        "%prog is a program for finding MPR landscapes across range of event" +
        "counts.  " +
        "See http://compbio.mit.edu/dlcpar for details.",

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
    grp_io.add_argument("--events", dest="events",
                      metavar="<output events>",
                      choices=[None,"I","U"], default=None,
                      help="set to output (I)ntersection or (U)nion of events (default: None)")
    grp_io.add_argument("--draw_regions", dest="draw_regions",
                      default=False,
                      action="store_true",
                      help="set to draw regions to screen")

    grp_ext = parser.add_argument_group("File Extensions")
    grp_ext.add_argument("-I","--inputext", dest="inext",
                       metavar="<input file extension>",
                       default="",
                       help="input file extension (default: \"\")")
    grp_ext.add_argument("-O", "--outputext", dest="outext",
                       metavar="<output file extension>",
                       default=".dlcscape",
                       help="output file extension (default: \".dlcscape\")")

    grp_costs = parser.add_argument_group("Costs")
    grp_costs.add_argument("-D", "--duprange", dest="duprange",
                         metavar="<dup range>",
                         default=dlcpar.reconscape.DEFAULT_RANGE,
                         action=ParseRangeAction,
                         help="duplication range (default: %s)" % dlcpar.reconscape.DEFAULT_RANGE_STR)
    grp_costs.add_argument("-L", "--lossrange", dest="lossrange",
                         metavar="<loss range>",
                         default=dlcpar.reconscape.DEFAULT_RANGE,
                         action=ParseRangeAction,
                         help="loss range (default: %s)" % dlcpar.reconscape.DEFAULT_RANGE_STR)

    grp_heur = parser.add_argument_group("Heuristics")
    grp_heur.add_argument("--max_loci", dest="max_loci",
                        metavar="<max # of loci>",
                        default=-1, type=int,
                        help="maximum # of co-existing loci (in each ancestral species), " +\
                             "set to -1 for no limit (default: -1)")
    grp_heur.add_argument("--max_dups", dest="max_dups",
                        metavar="<max # of dups>",
                        default=4, type=int,
                        help="maximum # of duplications (in each ancestral species), " +\
                             "set to -1 for no limit (default: 4)")
    grp_heur.add_argument("--max_losses", dest="max_losses",
                        metavar="<max # of losses>",
                        default=4, type=int,
                        help="maximum # of losses (in each ancestral species), " +\
                             "set to -1 for no limit (default: 4)")

    grp_misc = parser.add_argument_group("Miscellaneous")
    grp_misc.add_argument("-x", "--seed", dest="seed",
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

    # read species tree and species map (and possibly locus map)
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
        out_info = util.open_stream(out + ".info", "w")

        # log file
        if args.log:
            out_log = gzip.open(out + ".log.gz", 'w')
        else:
            out_log = common.NullLog()

        # command
        cmd = "%s %s\n\n" % (os.path.basename(sys.argv[0]),
                             ' '.join(map(lambda x: x if x.find(' ') == -1 else "\"%s\"" % x,
                                          sys.argv[1:])))
        out_info.write("DLCscape version:\t%s\n" % VERSION)
        out_info.write("DLCscape command:\t%s\n\n" % cmd)
        out_info.write("Seed: %d\n\n" % args.seed)
        out_log.write("DLCscape version: %s\n" % VERSION)
        out_log.write("DLCscape executed with the following arguments:\n")
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

        # set random seed
        random.seed(args.seed)
        np.random.seed(args.seed)
        out_log.write("seed: %d\n\n" % args.seed)

        compute_events = (args.events is not None)

        # perform reconciliation
        return_vals = dlcpar.reconscape.dlcscape_recon(
            coal_tree_top, stree, gene2species, gene2locus,
            duprange=args.duprange, lossrange=args.lossrange,
            max_loci=args.max_loci, max_dups=args.max_dups, max_losses=args.max_losses,
            compute_events=compute_events,
            log=out_log)
        out_log.write("\n")

        # infeasible
        if return_vals[1] is None:
            out_info.write("Feasibility:\tinfeasible\n")
            out_info.close()
            if args.log:
                out_log.close()
            break

        # feasible
        cvs, runtime = return_vals

        # determine landscape
        regions = dlcpar.reconscape.get_regions(
            cvs, duprange=args.duprange, lossrange=args.lossrange,
            log=out_log)

        # write the runtime to the logging files
        out_info.write("Runtime:\t%f sec\n" % runtime)

        # write regions
        dlcpar.reconscape.write_regions(out + ".regions", regions,
                                        duprange=args.duprange, lossrange=args.lossrange,
                                        close=True)

        # write events
        if compute_events:
            dlcpar.reconscape.write_events(out + ".events", regions,
                                           intersect=(args.events=="I"),
                                           close=True)

        # draw regions if desired
        if args.draw_regions:
            dlcpar.reconscape.draw_landscape(regions, duprange=args.duprange, lossrange=args.lossrange)

        # end logging
        if args.log:
            out_log.close()

        # end info
        out_info.close()

