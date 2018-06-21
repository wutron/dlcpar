"""
   Code for the DLC ILP Reconciliation
   (duplications, losses, and coalescence)
"""
# python libraries
import collections
import sys

# rasmus, compbio libraries
from rasmus import treelib, util
from compbio import phylo

# dlcpar libraries
from dlcpar import common, reconlib, dlclp_utils
from dlcpar import recon
from recon import _random_choice

# integer linear programming
import pulp
from dlcpar import ilpreconlib

# The following attributes in DLCLPRecon correspond to variables described DLCLP paper
# gtree = T_g (but with implied speciation nodes)
# stree = T_s
# srecon = M
# lrecon = L

def ilp_recon(tree, stree, gene2species,
              dupcost=1, losscost=1, coalcost=1, coaldupcost=None,
              delay=True, log=sys.stdout):
    """Perform reconciliation using DLCoal model with parsimony costs"""

    reconer = DLCLPRecon(tree, stree, gene2species,
                         dupcost=dupcost, losscost=losscost, coalcost=coalcost, coaldupcost=coaldupcost,
                         delay=delay, log=log)
    return reconer.recon()


class DLCLPRecon(recon.DLCRecon):

    def __init__(self, gtree, stree, gene2species,
                 dupcost=1, losscost=1, coalcost=1, coaldupcost=None,
                 delay=True, name_internal="n", log=sys.stdout):

        # rename gene tree nodes
        common.rename_nodes(gtree, name_internal)

        self.gtree = gtree
        self.stree = stree
        self.gene2species = gene2species

        assert (dupcost >= 0) and (losscost >= 0) and (coalcost >= 0) and (coaldupcost >= 0), (dupcost, losscost, coalcost, coaldupcost)
        self.dupcost = dupcost
        self.losscost = losscost
        self.coalcost = coalcost  # actually coalspeccost, using coalcost for backwards compatibility
        self.coaldupcost = coaldupcost if coaldupcost is not None else coalcost

        if delay:
            raise Exception("delay=True not allowed")
        self.delay = delay

        self.name_internal = name_internal
        self.log = util.Timer(log)

        # these attributes are assigned when performing reconciliation using self.recon()
        #   self.srecon
        #   self.lrecon
        #   self.order
        #   self.cost


    def recon(self):
        """Perform reconciliation"""

        self.log.start("Reconciling")

        # log input gene and species trees
        # log_tree(self.gtree, self.log, func=treelib.draw_tree_names)
        # log_tree(self.stree, self.log, func=treelib.draw_tree_names)

        self.log.start("Inferring species map")
        print("hi")
        self._infer_species_map()
        self.log.stop()

        # add implied speciation nodes but first start the species tree at the right root
        substree = treelib.subtree(self.stree, self.srecon[self.gtree.root])
        subsrecon = util.mapdict(self.srecon, val=lambda snode: substree.nodes[snode.name])

        # switch internal storage with subtrees
        self.stree, substree = substree, self.stree
        self.srecon, subsrecon = subsrecon, self.srecon

        # add implied nodes (standard speciation, speciation from duplication, delay nodes)
        # then relabel events (so that factor_tree works)
        self.log.start("Adding implied nodes")
        reconlib.add_implied_nodes(self.gtree, self.stree, self.srecon, self.sevents, delay=self.delay)
        self.sevents = phylo.label_events(self.gtree, self.srecon)
        common.rename_nodes(self.gtree, self.name_internal)
        self.log.stop()

        # log gene tree (with species map)
        # log_tree(self.gtree, self.log, func=draw_tree_srecon, srecon=self.srecon)

        # find subtrees
        self.log.start("Finding subtrees")
        self.subtrees = reconlib.factor_tree(self.gtree, self.stree, self.srecon, self.sevents)
        self.log.stop()

        # create and solve ILP
        ilp, lpvars = self._solve()
        self.cost = pulp.value(ilp.objective)

        # log gene tree (with species map and locus map)
        # log_tree(self.gtree, self.log, func=draw_tree_recon, srecon=self.srecon, lrecon=self.lrecon)

        # revert to use input species tree
        self.stree = substree
        self.srecon = util.mapdict(self.srecon, val=lambda snode: self.stree.nodes[snode.name])
        
        # convert to LabeledRecon format
        self.log.start("Converting to LCT")
        labeled_recon = ilpreconlib.ilp_to_lct(self.gtree, self.srecon, lpvars)
        self.log.stop()
        recon.draw_tree_recon(self.gtree, self.srecon, labeled_recon.locus_map, minlen=15, maxlen=15)

        runtime = self.log.stop()

        return self.gtree, labeled_recon, runtime, self.cost


    def _solve(self):
        """Create and solve ILP problem"""

        self.log.start("Creating ilp variables")
        ilp = pulp.LpProblem("dup placement", pulp.LpMinimize)
        lpvars = ilpreconlib.IlpReconVariables(self.gtree, self.srecon, self.subtrees)
        self.log.stop()

        self.log.start("Building ilp constraints")
        ilp += self._create_objective_func(lpvars)
        self._add_constraints(ilp, lpvars)
        self.log.stop()

        self.log.start("Solving ilp")
        #ilp.solve()
        ilp.solve(solver=pulp.CPLEX_CMD)
        self.log.stop()

        return ilp, lpvars


    def _create_objective_func(self, lpvars):
        """ This method returns the cost function for DLC reconciliation that the ilp will minimize. """ 
        num_dups = pulp.lpSum(lpvars.dup_vars.values())
        num_losses = pulp.lpSum(lpvars._loss_vars.values())
        num_coals = pulp.lpSum(lpvars._coalspec_vars.values())
        num_coaldups = pulp.lpSum(lpvars._coaldup_vars.values())

        return self.dupcost * num_dups + self.losscost * num_losses + self.coalcost * num_coals + self.coaldupcost * num_coaldups


    def _add_constraints(self, ilp, lpvars):
        """ This function adds all the constraints as described in the dlclp paper to the ilp so that the values of the ilp
        variables will be correctly found for most parsimonious reconciliation """
        # create the path constraints - if there's a dup on a given path, then that path var is 1, otherwise 0
        for (g1, g2), path_var in lpvars._path_vars.iteritems():
            path = dlclp_utils.find_path_without_lca(g1, g2)
            ilp += pulp.lpSum([lpvars.dup_vars[node] for node in path]) - path_var >= 0
            ilp += pulp.lpSum([lpvars.dup_vars[node] for node in path]) - len(path) * path_var <= 0

        # use the srecon to create a dictionary
        # key - snode
        # value - list of the gene leaves mapped there
        gleaves_by_species = dlclp_utils.group_by_species(self.gtree.leaves(), self.srecon)

        # create the dup constraints - every path between each gene leaf mapped to the same species leaf must have a dup
        for gleaves in gleaves_by_species.values():
            for g1, g2 in pulp.combination(gleaves, 2):
                path = dlclp_utils.find_path_without_lca(g1, g2)
                ilp += pulp.lpSum([lpvars.dup_vars[node] for node in path]) >= 1

        # create the loss constraints -
        for snode, local_top in lpvars._loss_vars:
            local_loss = lpvars._loss_vars[(snode, local_top)]
            local_helper = lpvars._helper_vars[(snode, local_top)]

            # find the corresponding path variables to the bottom nodes
            local_bottoms = lpvars._bottom_nodes[snode]
            local_paths = [lpvars.get_path(local_top, local_bottom) for local_bottom in local_bottoms]

            # NOTE: I flipped the inequality sign of this constraint
            ilp += pulp.lpSum(local_paths) - local_loss - local_helper <= len(local_bottoms) - 1

        # create the helper constraints -
        for snode, local_top in lpvars._helper_vars:
            local_helper = lpvars._helper_vars[(snode, local_top)]

            # get all top_nodes_at_snode that are lexicographically "before" this one
            top_nodes_in_snode = lpvars._top_nodes[snode]
            prev_nodes = top_nodes_in_snode[:top_nodes_in_snode.index(local_top)]

            # get all previous path vars
            prev_paths = [lpvars.get_path(local_top, g) for g in prev_nodes]

            ilp += pulp.lpSum(prev_paths) + local_helper <= len(prev_nodes)
            ilp += pulp.lpSum(prev_paths) + len(lpvars._top_nodes[snode]) * local_helper >= len(prev_nodes)

        # create the coal constraints -
        for snode, local_top in lpvars._coalspec_vars:
            local_coal = lpvars._coalspec_vars[(snode, local_top)]

            # get all nodes that have a child (might contribute to a coal) lexicographically "before" this one
            tops_with_child_in_s = lpvars._top_nodes_with_child[snode]
            prev_nodes = tops_with_child_in_s[:tops_with_child_in_s.index(local_top)]

            # get all previous path vars
            prev_paths = [lpvars.get_path(local_top, node) for node in prev_nodes]

            ilp += pulp.lpSum(prev_paths) + local_coal <= len(prev_nodes)
            ilp += pulp.lpSum(prev_paths) + len(lpvars._top_nodes[snode]) * local_coal >= len(prev_nodes)

        # create the delta constraints -
        for g1, g2 in lpvars._delta_vars:
            g1_at_time_of_dup_at_g2 = lpvars.get_order(g1.parent, g2) + lpvars.get_order(g2, g1)
            g1_and_g2_at_same_locus = 1 - lpvars.get_path(g1.parent, g2.parent)

            ilp += lpvars._delta_vars[(g1, g2)] >= \
                    g1_at_time_of_dup_at_g2 + g1_and_g2_at_same_locus + lpvars.dup_vars[g2] - 3

        # create the transitive property constraints for order -
        for gnodes in list(pulp.combination(list(self.gtree.preorder()), 3)):
            g1, g2, g3 = gnodes
            if self.srecon[g1] == self.srecon[g2] == self.srecon[g3]:
                leaves = [g for g in gnodes if dlclp_utils.is_species_leaf(g, self.srecon)]

                if len(leaves) < 2:
                    ilp += lpvars.get_order(g1, g3) >= lpvars.get_order(g1, g2) + lpvars.get_order(g2, g3) - 1
                elif len(leaves) == 2:
                    non_leaf = [g for g in gnodes if g not in leaves][0]
                    ilp += lpvars.get_order(non_leaf, leaves[0]) == 1
                    ilp += lpvars.get_order(non_leaf, leaves[1]) == 1
                else:
                    # The transitive property is not neccessary for leaves so do nothing when gnodes is all leaves
                    pass

        # create the dup constraints for order of leaves -
        for (g1, g2), order in lpvars.order_vars.iteritems():
            if dlclp_utils.is_leaves_of_same_species(self.srecon, g1, g2):
                ilp += order >= lpvars.dup_vars[g1] - lpvars.dup_vars[g2]
                ilp += order <= lpvars.dup_vars[g1] - lpvars.dup_vars[g2] + 1

        # create the r constraints -
        for g2 in self.gtree.preorder():
            gnodes_in_species = lpvars._gnodes_by_species[self.srecon[g2]]
            other_delta_vars_in_same_species = [lpvars._delta_vars[(g1, g2)] for g1 in gnodes_in_species
                                                if (g1, g2) in lpvars._delta_vars]
            ilp += lpvars._coaldup_vars[g2] >= pulp.lpSum(other_delta_vars_in_same_species) - 1
