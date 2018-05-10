
# integer linear programming
from pulp import *
from dlcpar import dlclp_utils
from dlcpar.ilprecon_variables import IlpReconVariables


def solve(gtree, srecon, subtrees, dupcost, losscost, coalcost, coaldupcost):
    ilp = LpProblem("dup placement", LpMinimize)

    lpvars = IlpReconVariables(gtree, srecon, subtrees)
    ilp += _create_objective_func(lpvars, dupcost, losscost, coalcost, coaldupcost)
    ilp += _create_constraints(lpvars)

    ilp.solve()

    return ilp, lpvars


def _create_objective_func(lpvars, dupcost, losscost, coalcost, coaldupcost):
    num_dups = lpSum(lpvars.dup_vars.values())
    num_losses = lpSum(lpvars.loss_vars.values())
    num_coals = lpSum(lpvars.coal_vars.values())
    num_coaldups = lpSum(lpvars.r_vars.values())

    return dupcost * num_dups + losscost * num_losses + coalcost * num_coals + coaldupcost * num_coaldups


def _create_constraints(lpvars):
    consr = LpConstraint()
    
    for snode, local_top in lpvars.coal_vars:
        # create the path constraints - if there's a dup on a given path, then that path var is 1, otherwise 0
        for (g1, g2), path_var in lpvars.path_vars.iteritems():
            path = dlclp_utils.find_path_without_lca(g1, g2)
            consr += lpSum([lpvars.dup_vars[node] for node in path]) - path_var >= 0
            consr += lpSum([lpvars.dup_vars[node] for node in path]) - len(path) * path_var <= 0

        # use the srecon to create a dictionary
        # key - snode
        # value - list of the gene leaves mapped there
        gleaves_by_species = dlclp_utils.group_by_species(lpvars.gtree.leaves(), lpvars.srecon)

        # create the dup constraints - every path between each gene leaf mapped to the same species leaf must have a dup
        for gleaves in gleaves_by_species.values():
            for g1, g2 in combination(gleaves, 2):
                path = dlclp_utils.find_path_without_lca(g1, g2)
                consr += lpSum([lpvars.dup_vars[node] for node in path]) >= 1

        # create the loss constraints -
        for snode, local_top in lpvars.loss_vars:
            local_loss = lpvars.loss_vars[(snode, local_top)]
            local_helper = lpvars.helper_vars[(snode, local_top)]

            # find the corresponding path variables to the bottom nodes
            local_bottoms = lpvars.bottom_nodes[snode]
            local_paths = [lpvars.get_path(local_top, local_bottom) for local_bottom in local_bottoms]

            # NOTE: I flipped the inequality sign of this constraint
            consr += lpSum(local_paths) - local_loss - local_helper <= len(local_bottoms) - 1

        # create the helper constraints -
        for snode, local_top in lpvars.helper_vars:
            local_helper = lpvars.helper_vars[(snode, local_top)]

            # get all top_nodes_at_snode that are lexicographically "before" this one
            top_nodes_in_snode = lpvars.top_nodes[snode]
            prev_nodes = top_nodes_in_snode[:top_nodes_in_snode.index(local_top)]

            # get all previous path vars
            prev_paths = [lpvars.get_path(local_top, g) for g in prev_nodes]

            consr += lpSum(prev_paths) + local_helper <= len(prev_nodes)
            consr += lpSum(prev_paths) + len(lpvars.top_nodes[snode]) * local_helper >= len(prev_nodes)

        # create the coal constraints -
        local_coal = lpvars.coal_vars[(snode, local_top)]

        # get all nodes that have a child (might contribute to a coal) lexicographically "before" this one
        tops_with_child_in_s = lpvars.top_nodes_with_child[snode]
        prev_nodes = tops_with_child_in_s[:tops_with_child_in_s.index(local_top)]

        # get all previous path vars
        prev_paths = [lpvars.get_path(local_top, node) for node in prev_nodes]

        consr += lpSum(prev_paths) + local_coal <= len(prev_nodes)
        consr += lpSum(prev_paths) + len(lpvars.top_nodes[snode]) * local_coal >= len(prev_nodes)

    # create the delta constraints -
    for g1, g2 in lpvars.delta_vars:
        g1_at_time_of_dup_at_g2 = lpvars.get_order(g1.parent, g2) + lpvars.get_order(g2, g1)
        g1_and_g2_at_same_locus = 1 - lpvars.get_path(g1.parent, g2.parent)

        consr += lpvars.delta_vars[(g1, g2)] >= \
                g1_at_time_of_dup_at_g2 + g1_and_g2_at_same_locus + lpvars.dup_vars[g2] - 3

    # create the transitive property constraints for order -
    for gnodes in list(combination(list(lpvars.gtree.preorder()), 3)):
        g1, g2, g3 = gnodes
        if lpvars.srecon[g1] == lpvars.srecon[g2] == lpvars.srecon[g3]:
            leaves = [g for g in gnodes if dlclp_utils.is_species_leaf(g, lpvars.srecon)]

            if len(leaves) < 2:
                consr += lpvars.get_order(g1, g3) >= lpvars.get_order(g1, g2) + lpvars.get_order(g2, g3) - 1
            elif len(leaves) == 2:
                non_leaf = [g for g in gnodes if g not in leaves][0]
                consr += lpvars.get_order(non_leaf, leaves[0]) == 1
                consr += lpvars.get_order(non_leaf, leaves[1]) == 1
            else:
                # The transitive property is not neccessary for leaves so do nothing when gnodes is all leaves
                pass

    # create the dup constraints for order of leaves -
    for (g1, g2), order in lpvars.order_vars.iteritems():
        if dlclp_utils.is_leaves_of_same_species(lpvars.srecon, g1, g2):
            consr += order >= lpvars.dup_vars[g1] - lpvars.dup_vars[g2]
            consr += order <= lpvars.dup_vars[g1] - lpvars.dup_vars[g2] + 1

    # create the r constraints -
    for g2 in lpvars.gtree.preorder():
        gnodes_in_species = lpvars.gnodes_by_species[lpvars.srecon[g2]]
        other_delta_vars_in_same_species = [lpvars.delta_vars[(g1, g2)] for g1 in gnodes_in_species
                                            if (g1, g2) in lpvars.delta_vars]
        consr += lpvars.r_vars[g2] >= lpSum(other_delta_vars_in_same_species) - 1

    return consr
