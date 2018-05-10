
# integer linear programming
from pulp import *
from rasmus import treelib


def find_path_without_lca(node1, node2):
    """
    :return list[treelib.TreeNode]: the nodes on the path from node1 to node2 except for
    their LCA.
    """
    lca = treelib.lca([node1, node2])

    def path_up_to_lca(node):
        if node == lca:
            return []
        else:
            return [node] + path_up_to_lca(node.parent)

    return path_up_to_lca(node1) + path_up_to_lca(node2)


def descendants_in_species(gnode, srecon):
    """
    :param srecon: like species_map in dlcpar.reconolib.LabeledRecon
    :return list[TreeNode]: list of descendants of gnode that are in the same species as gnode. To be clear, that
    gnode is not included
    """
    children_in_species = [child for child in gnode.children if srecon[gnode] == srecon[child]]
    indirect_descendants = sum([descendants_in_species(child, srecon) for child in children_in_species], [])

    return children_in_species + indirect_descendants


def is_leaves_of_same_species(srecon, g1, g2):
    return srecon[g1] == srecon[g2] and is_species_leaf(g1, srecon) and is_species_leaf(g2, srecon)


def is_species_leaf(gnode, srecon):
    return len(gnode.children) == 0 or srecon[gnode] != srecon[gnode.children[0]]


def group_by_species(gnodes, srecon):
    gnodes_by_species = collections.defaultdict(list)
    for gnode in gnodes:
        gnodes_by_species[srecon[gnode]].append(gnode)

    return gnodes_by_species


def pairs_in_species(gnodes_grouped_by_species):
    """
    :return: all pairs of gnodes such that both gnodes in the pair is in the same species
    """
    pairs_for_each_species = [list(combination(gnodes, 2)) for gnodes in gnodes_grouped_by_species.values()]
    # combine the lists into one big list of pairs
    return sum(pairs_for_each_species, [])
