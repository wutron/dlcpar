# python libraries
import sys

# rasmus libraries
from rasmus import treelib, util

#=============================
# input

def move_option(parser, opt_str, opt_grp):
    """TODO: delete"""
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


#=============================
#  tree utilities

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

def find_path(node1, node2):
    """Find the path between two nodes in a tree

    Returns node names along path from each node up to (but excluding) lca.
    Based on treelib.find_dist(...).
    """
    # keep track of input nodes for error checking
    n1, n2 = node1, node2

    # find root path for node1 [node1, ..., root]
    path1 = [node1.name]
    while node1.parent is not None:
        node1 = node1.parent
        path1.append(node1.name)

    # find root path for node2 [node2, ..., root]
    path2 = [node2.name]
    while node2.parent is not None:
        node2 = node2.parent
        path2.append(node2.name)

    # find when paths diverge (pathX[-i+1] is the lca)
    i = 1
    while i <= len(path1) and i <= len(path2) and (path1[-i] == path2[-i]):
        i += 1
    assert path1[-i+1] == path2[-i+1] == treelib.lca((n1, n2)).name, \
        (n1.name, n2.name, path1[-i+1], path2[-i+1], treelib.lca((n1, n2)).name, i)

    return (path1[-i::-1], path2[-i::-1])

def random_choice(a, p=None):
    """Return a random choice based on probabilities

    wrapper around numpy.random.choice
    note: numpy cannot take list of lists so use indices
    """
    import numpy as np
    a = list(a)
    ndx = np.random.choice(range(len(a)), p=p)
    it = a[ndx]
    return it
