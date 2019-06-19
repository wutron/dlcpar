'''Views landscape'''
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
    '''main program'''

    #==========================
    # parser

    parser = argparse.ArgumentParser(
        usage = "%(prog)s view_landscape [options] <landscape>",
        description="View landscape",
        formatter_class=commands.CustomHelpFormatter,
        add_help=False)
    
    parser.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)

    parser.add_argument("landscape", help=argparse.SUPPRESS)

    grp_display = parser.add_argument_group("Display")
    grp_display.add_argument("-o", "--output", dest="output",
                           metavar="<output>")
    grp_display.add_argument("--linear", dest="linear",
                           default=False, action="store_true")

    grp_info = parser.add_argument_group("Information")
  

    args = parser.parse_args(sys.argv[2:])

    landscape = args.landscape

    

        


    #===========================================
    # process

    # read regions
    regions, duprange, lossrange = dlcpar.reconscape.read_regions(landscape)

    # draw
    dlcpar.reconscape.draw_landscape(regions, duprange, lossrange,
                                    filename=args.output,
                                    log=not args.linear)


    