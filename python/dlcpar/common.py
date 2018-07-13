# python libraries
import sys

# rasmus libraries
from rasmus import treelib, util

#=============================
# input

def move_option(parser, opt_str, opt_grp):
    """Move option 'opt_str' from 'parser' to 'opt_grp'"""
    if parser.has_option(opt_str):
        opt = parser.get_option(opt_str)
        parser.remove_option(opt_str)
        opt_grp.add_option(opt)


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
    """Rename nodes so that all names are strings"""
    for node in list(tree.postorder()):
        if isinstance(node.name, int):
            name2 = prefix + str(node.name)
            while name2 in tree.nodes:
                name2 = prefix + str(tree.new_name())
            tree.rename(node.name, name2)

def check_tree(tree, name=""):
    """Ensure that tree is rooted and binary"""
    if not (treelib.is_rooted(tree) and treelib.is_binary(tree)):
        raise Exception("tree must be rooted and binary: %s" % name)


def random_choice(a, p=None):
    # wrapper around numpy.random.choice
    # note: numpy cannot take list of lists so use indices
    import numpy as np
    a = list(a)
    ndx = np.random.choice(range(len(a)), p=p)
    it = a[ndx]
    return it
