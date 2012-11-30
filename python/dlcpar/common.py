import sys
from rasmus import treelib, util

#=============================
# input

def add_common_options(parser,
                       infiles=False ,reroot=False,
                       stree=False, smap=False,
		       treeext=False, alignext=False,
		       clade=False):
    """Add common options to parser"""
    if infiles:
        parser.add_option("-i", "--input", dest="input",
                          action="append",
                          metavar="<input file>",
			  help="list of input files, one per line")
    if stree:
       parser.add_option("-s", "--stree", dest="stree",
                         metavar="<species tree>",
			 help="species tree file in newick format")
    if smap:
       parser.add_option("-S", "--smap", dest="smap",
                         metavar="<species map>",
                         help="gene to species map")

def move_option(parser, opt_str, opt_grp):
    """Move option 'opt_str' from 'parser' to 'opt_grp'"""
    if parser.has_option(opt_str):
        opt = parser.get_option(opt_str)
        parser.remove_option(opt_str)
        opt_grp.add_option(opt)

def check_req_options(parser, options,
                      species=True):
    """Check if required options are present"""
    if species and ((not options.stree) or (not options.smap)):
        parser.error("-s/--stree and -S/--smap are required")

def get_input_files(parser, options, args):
    """Determine input files from options"""
    infiles = []
    if options.input:
        for arg in options.input:
            if arg == "-":
                infiles.append(sys.stdin)
            else:
                infiles.append(open(arg))

    # determine all input lines
    files = args
    for infile in infiles:
        files.extend(map(lambda fn:fn.rstrip("\n"),infile.readlines()))
    if len(files) == 0:
        parser.error("must specify input file(s)")
        
    return files

#=============================
#  utilities

class NullLog (object):

    def __init__(self):
        pass

    def write(self, text):
        pass

    def flush(self):
        pass

def rename_nodes(tree, prefix="n"):
    """Rename nodes that all names are strings"""
    for node in list(tree.postorder()):
        if isinstance(node.name, int):
            name2 = prefix + str(node.name)
            while name2 in tree.nodes:
                name2 = prefix + str(tree.new_name())
            tree.rename(node.name, name2)
