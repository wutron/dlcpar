"""
common.py
Common utilities for dlcpar
"""

# numpy libraries
import numpy as np

# rasmus libraries
from rasmus import treelib

#=============================
# utilities

class NullLog(object):
    """Null log"""

    def __init__(self):
        pass

    def write(self, text):
        """Write to log"""
        pass

    def flush(self):
        """Flush log"""
        pass


def random_choice(a, p=None, normalize=False):
    """Return a random choice based on probabilities

    wrapper around numpy.random.choice
    note: numpy cannot take list of lists so use indices

    if norm, normalize probabilities
    """
    if normalize:
        tot = sum(p)
        p = [float(val) / tot for val in p]

    a = list(a)
    ndx = np.random.choice(range(len(a)), p=p)
    it = a[ndx]
    return it


#=============================
# tree utilities

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


def prepare_tree(tree):
    """Remove bootstrap and distances if they exist"""
    tree_top = tree.copy()
    for node in tree_top:
        if "boot" in node.data:
            del node.data["boot"]
        node.dist = 0
    tree_top.default_data.clear()
    return tree_top


def find_path(node1, node2):
    """Find the path between two nodes in a tree

    Returns node names along path from each node up to (but excluding) lca.
    Based on treelib.find_dist(...).
    """

    # keep track of input nodes for error checking
    lca = treelib.lca((node1, node2))
    name1, name2 = node1.name, node2.name

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
    assert path1[-i+1] == path2[-i+1] == lca.name, \
        (name1, name2, path1[-i+1], path2[-i+1], lca.name, i)

    return (path1[-i::-1], path2[-i::-1])
