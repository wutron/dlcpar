"""
View landscape of equivalent regions
"""

# python libraries
import os, sys
import argparse

# dlcpar libraries
import dlcpar
from dlcpar import common
from dlcpar import commands
import dlcpar.reconscape

#=============================================================================
VERSION = dlcpar.PROGRAM_VERSION_TEXT

def run():
    """main program"""

    #==========================
    # parser

    parser = argparse.ArgumentParser(
        usage = "%(prog)s view_landscape [options] <regions>",
        description="View landscape of equivalent regions.",
        formatter_class=commands.CustomHelpFormatter,
        add_help=False)
    parser.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)

    parser.add_argument("regions", help=argparse.SUPPRESS)

    grp_display = parser.add_argument_group("Display")
    grp_display.add_argument("-o", "--output", dest="output",
                             metavar="<output>",
                             help="output file")
    grp_display.add_argument("--linear", dest="linear",
                             action="store_true",
                             help="set to use linear axes")

    args = parser.parse_args(sys.argv[2:])

    infile = args.regions

    #===========================================
    # process

    # read regions
    regions, duprange, lossrange = dlcpar.reconscape.read_regions(infile)

    # draw
    dlcpar.reconscape.draw_landscape(regions, duprange, lossrange,
                                     filename=args.output,
                                     log=not args.linear)
