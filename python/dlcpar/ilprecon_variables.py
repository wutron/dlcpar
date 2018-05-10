
# integer linear programming
from pulp import *
from dlcpar import dlclp_utils


class IlpReconVariables(object):
    """
    contains all of the ilp variables, and other relevant structures including gtree, gnodes_by_species and more
    See the DLCLP paper for descriptiosn of the ilp variables
    """

    def __init__(self, gtree, srecon, subtrees):
        self.gtree = gtree
        self.srecon = srecon

        self.bottom_nodes = collections.defaultdict(list)
        self.top_nodes_with_child = collections.defaultdict(list)
        self.top_nodes = collections.defaultdict(list)
        self._find_tops_bottoms_and_tops_with_child(subtrees)

        self.dup_vars = LpVariable.dicts("dup", list(gtree.preorder()), 0, 1, LpInteger)
        """duplication variables"""
        self.loss_vars = LpVariable.dicts("loss", _create_loss_keys(self.top_nodes), 0, 1, LpInteger)
        """"keys are (snode, gene node at the top of snode)"""
        self.coal_vars = LpVariable.dicts("coal", _create_coal_keys(self.top_nodes_with_child), 0, 1, LpInteger)
        """keys are (snode, gene node at the top of snode)"""
        self.helper_vars = LpVariable.dicts("helper", self.loss_vars.keys(), 0, 1, LpInteger)
        """
        loss helper variables - corresponds to h_g in DLCLP paper
        keys are (snode, gene node at the top of snode)
        """
        all_pairs = list(combination(list(gtree.preorder()), 2))
        self.path_vars = LpVariable.dicts("path", all_pairs, 0, 1, LpInteger)
        """
        keys are (gene node 1, gene node 2).
        NOTE: this is combinations, not permutations, since order doesn't matter
        """
        self.gnodes_by_species = dlclp_utils.group_by_species(gtree, srecon)
        """A dictionary. keys are species nodes. For each species node the corresponding value is
        the gnodes in that species"""
        pairs_in_species = dlclp_utils.pairs_in_species(self.gnodes_by_species)
        self.delta_vars = LpVariable.dicts("coal_dup", _create_delta_keys(pairs_in_species), 0, 1, LpInteger)
        """keys are (gene node 1, gene node 2). Order matters so for each (g1, g2) key there is also a
        distinct (g2, g1) key"""
        self.orders_from_topology = _infer_orders_from_topology(gtree.root, srecon)
        """keys are pairs of gnodes that are in the same species such that one gnode is an ancestor of the other.
        For each key (g1, g2), the corresponding value is 1 if g2 more recent than g1 and 0 otherwise. 
        These are combinations not permutations"""
        self.order_vars = LpVariable.dicts("order", _create_order_keys(pairs_in_species, self.orders_from_topology),
                                           0, 1, LpInteger)
        """keys are (g1, g2) such g1 and g2 are in the same species. These are combinations not permutations"""

        self.r_vars = LpVariable.dicts("r", list(gtree.preorder()), 0, None, LpInteger)

    def get_order(self, g1, g2):
        if (g1, g2) in self.order_vars:
            return self.order_vars[(g1, g2)]
        elif (g2, g1) in self.order_vars:
            return 1 - self.order_vars[(g2, g1)]
        elif (g1, g2) in self.orders_from_topology:
            return self.orders_from_topology[(g1, g2)]
        elif (g2, g1) in self.orders_from_topology:
            return 1 - self.orders_from_topology[(g2, g1)]
        elif self.srecon[g2] in self.srecon[g1].children:
            return 1
        else:
            raise Exception('despite my best efforts I could not find order for nodes', g1, g2)

    def get_path(self, g1, g2):
        if (g1, g2) in self.path_vars:
            return self.path_vars[(g1, g2)]
        elif g1 == g2:
            return 0
        else:
            return self.path_vars[(g2, g1)]

    def _find_tops_bottoms_and_tops_with_child(self, subtrees):
        # find top and bottom nodes for each snode
        # key - snode
        # value - list of top (or bottom) nodes for that snode
        for snode, subtrees in subtrees.iteritems():
            for subtree in subtrees:
                # the root of each subtree is a top node for this snode
                self.top_nodes[snode].append(subtree[0])
                if subtree[2]:
                    self.bottom_nodes[snode].extend(subtree[2])
                if subtree[1]:
                    self.top_nodes_with_child[snode].append(subtree[0])


def _create_loss_keys(top_nodes):
    top_pairs = []
    for snode, tops in top_nodes.iteritems():
        for top in tops:
            top_pairs.append((snode, top))

    return top_pairs


def _create_coal_keys(top_nodes_with_child):
    top_pairs_with_child = []
    for snode, tops in top_nodes_with_child.iteritems():
        for top in tops:
            top_pairs_with_child.append((snode, top))

    return top_pairs_with_child


def _create_delta_keys(pairs_in_species):
    pairs_in_species_with_flipped = pairs_in_species + [(g2, g1) for (g1, g2) in pairs_in_species]

    delta_vars_keys = [(g1, g2) for g1, g2 in pairs_in_species_with_flipped
                       if g1.parent is not None and g2.parent is not None and
                       g1.parent != g2]

    return delta_vars_keys


def _create_order_keys(pairs_in_species, orders_from_topology):
    return [(g1, g2) for (g1, g2) in pairs_in_species
            if (g1, g2) not in orders_from_topology and (g2, g1) not in orders_from_topology]


def _infer_orders_from_topology(root, srecon):
    order = {}

    # infers all the orders that only involve nodes in this subtree
    def infer_orders(sub_tree_root):
        for descendant in dlclp_utils.descendants_in_species(sub_tree_root, srecon):
            order[(sub_tree_root, descendant)] = 1

        for child in sub_tree_root.children:
            infer_orders(child)

    infer_orders(root)

    return order
