"""
   Code for the DLC ILP Reconciliation
   (duplications, losses, and coalescence)
"""

# python libraries
import sys
import collections

# rasmus, compbio libraries
from rasmus import treelib, util
from compbio import phylo

# dlcpar libraries
from dlcpar import common, reconlib

# integer linear programming
import pulp
from dlcpar import ilpreconlib
try:            # default to use CPLEX_PY if available
    import cplex
    solver_name = "CPLEX_PY"
except:         # otherwise use pulp's default solver
    solver_name = "CBC_CMD"

# The following attributes in DLCLPRecon correspond to variables described DLCLP paper
# gtree = T_g (with implied speciation nodes)
# stree = T_s
# srecon = M
# lrecon = L

#==========================================================

def ilp_recon(tree, stree, gene2species,
              dupcost=1, losscost=1, coalcost=1, coaldupcost=None,
              timeLimit=None, delay=True, optimize=True, log=sys.stdout):
    """Perform reconciliation using DLCoal model with parsimony costs"""

    reconer = DLCLPRecon(tree, stree, gene2species,
                         dupcost=dupcost, losscost=losscost, coalcost=coalcost, coaldupcost=coaldupcost,
                         timeLimit=timeLimit, delay=delay, optimize=optimize, log=log)
    return reconer.recon()


class DLCLPRecon(object):

    def __init__(self, gtree, stree, gene2species,
                 dupcost=1, losscost=1, coalcost=1, coaldupcost=None, timeLimit=None,
                 delay=True, optimize=True,
                 name_internal="n", log=sys.stdout):

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

        self.timeLimit = timeLimit
        self.optimize = optimize
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

    #=============================
    # main methods

    def recon(self):
        """Perform reconciliation"""

        self.log.start("Reconciling")

        # log input gene and species trees
        self.log.log("gene tree\n")
        reconlib.log_tree(self.gtree, self.log, func=treelib.draw_tree_names)
        self.log.log("species tree\n")
        reconlib.log_tree(self.stree, self.log, func=treelib.draw_tree_names)

        # infer species map
        self.log.start("Inferring species map")
        self.srecon = phylo.reconcile(self.gtree, self.stree, self.gene2species)
        self.log.stop()
        self.log.log("\n\n")

        # add implied speciation nodes but first start the species tree at the right root
        substree = treelib.subtree(self.stree, self.srecon[self.gtree.root])
        subsrecon = util.mapdict(self.srecon, val=lambda snode: substree.nodes[snode.name])

        # switch internal storage with subtrees
        self.stree, substree = substree, self.stree
        self.srecon, subsrecon = subsrecon, self.srecon

        # add implied nodes (standard speciation, speciation from duplication, delay nodes)
        # then relabel events (so that factor_tree works)
        self.log.start("Adding implied nodes")
        sevents = phylo.label_events(self.gtree, self.srecon)
        reconlib.add_implied_nodes(self.gtree, self.stree, self.srecon, sevents, delay=self.delay)
        common.rename_nodes(self.gtree, self.name_internal)
        self.log.stop()

        # log gene tree (with species map)
        reconlib.log_tree(self.gtree, self.log, func=reconlib.draw_tree_recon, srecon=self.srecon)

        # infer locus map
        self.log.start("Inferring locus map")
        ilp, lpvars, setup_runtime, solve_runtime = self._infer_locus_map()
        self.log.stop()
        self.log.log("\n\n")
        
        # revert to use input species tree
        self.stree = substree
        self.srecon = util.mapdict(self.srecon, val=lambda snode: self.stree.nodes[snode.name])

        # convert to LabeledRecon data structure
        self.log.start("Converting to LCT")
        labeled_recon = ilpreconlib.ilp_to_lct(self.gtree, lpvars)
        self.log.stop()

        # log gene tree (with species map and locus map)
        self.log.log("gene tree (with species and locus map)\n")
        reconlib.log_tree(self.gtree, self.log, func=reconlib.draw_tree_recon,
                          srecon=labeled_recon.species_map, lrecon=labeled_recon.locus_map)
    
        # calculate runtime
        runtime = self.log.stop()

        return self.gtree, labeled_recon, runtime, setup_runtime, solve_runtime, self.cost


    def _infer_locus_map(self):
        """Infer (and assign) locus map"""

        # create and solve ILP
        self.log.start("Creating ilp variables")
        ilp = pulp.LpProblem("dlcilp", pulp.LpMinimize)
        lpvars = ilpreconlib.IlpReconVariables(self.gtree, self.stree, self.srecon)
        self.log.stop()

        self.log.start("Building ilp constraints")
        ilp += self._create_objective_func(lpvars)
        self._add_constraints(ilp, lpvars)
        setup_runtime = self.log.stop()

        if solver_name == "CPLEX_PY":
            ilpsolver = pulp.solvers.CPLEX_PY(timeLimit=self.timeLimit)
        else:
            ilpsolver = pulp.solvers.PULP_CBC_CMD(maxSeconds=self.timeLimit)
        
        self.log.start("Solving ilp")
        ilp.solve(solver=ilpsolver)
        solve_runtime = self.log.stop()
        self.log.log("Solver: " + solver_name)
        self.cost = pulp.value(ilp.objective)

        return ilp, lpvars, setup_runtime, solve_runtime


    def _create_objective_func(self, lpvars):
        """Return cost function for ILP formulation."""
        num_dups = pulp.lpSum(lpvars.dup_vars.values())
        num_losses = pulp.lpSum(lpvars._loss_vars.values())
        num_coals = pulp.lpSum(lpvars._coalspec_vars.values())
        num_coaldups = pulp.lpSum(lpvars._coaldup_vars.values())

        return self.dupcost * num_dups + self.losscost * num_losses + \
            self.coalcost * num_coals + self.coaldupcost * num_coaldups


    def _add_constraints(self, ilp, lpvars):
        """Add constraints for ILP formulation."""

        all_gnodes = list(self.gtree.preorder())
        
        # create dup constraints
        leaves_by_species = collections.defaultdict(list)
        for leaf in self.gtree.leaves():
            leaves_by_species[self.srecon[leaf]].append(leaf)
        for leaves in leaves_by_species.itervalues():
            for g1, g2 in pulp.combination(leaves, 2):
                path1, path2 = common.find_path(g1, g2)
                path = path1 + path2
                nodes = [self.gtree[name] for name in path]
                ilp += pulp.lpSum([lpvars.dup_vars[node] for node in nodes]) >= 1

        # create duplication optimization constraints (from a gnode g with 2 children g' and g'', only 1 will have a duplication)
        #this is an optimization constraint which is not crucial --MORGAN
        if self.optimize:
            for gnode in all_gnodes:
                if not gnode.is_leaf():
                    assert(len(gnode.children) == 2, "ilprecon only takes binary gene trees")
                    ilp += pulp.lpSum([lpvars.dup_vars[node] for node in gnode.children]) <= 1

                
        # create the path constraints - if there is dup on given path, then that path var is 1, otherwise 0
        for (g1, g2), path_var in lpvars._path_vars.iteritems():
            path1, path2 = common.find_path(g1, g2)
            path = path1 + path2
            nodes = [self.gtree[name] for name in path]
            ilp += path_var <= pulp.lpSum([lpvars.dup_vars[node] for node in nodes])
            ilp += pulp.lpSum([lpvars.dup_vars[node] for node in nodes]) <= len(path) * path_var

        # create loss constraints
        for (snode, gnode), local_loss in lpvars._loss_vars.iteritems():
            local_helper = lpvars._helper_vars[(snode, gnode)]
            local_bottoms = lpvars._bottom_nodes[snode]
            local_paths = [lpvars.get_path(gnode, local_bottom) for local_bottom in local_bottoms]
            ilp += local_loss >= pulp.lpSum(local_paths) - len(local_bottoms) + 1 - local_helper

        # create helper constraints
        for (snode, gnode), local_helper in lpvars._helper_vars.iteritems():
            # get all top_nodes_at_snode that are lexicographically "before" this one
            top_nodes_in_snode = lpvars._top_nodes[snode]
            prev_nodes = top_nodes_in_snode[:top_nodes_in_snode.index(gnode)]

            # get all previous path vars
            prev_paths = [lpvars.get_path(gnode, other) for other in prev_nodes]

            ilp += local_helper <= len(prev_nodes) - pulp.lpSum(prev_paths)
            ilp +=  len(prev_nodes) - pulp.lpSum(prev_paths) <= len(lpvars._top_nodes[snode]) * local_helper

        # create coal constraints
        for (snode, gnode), local_coal in lpvars._coalspec_vars.iteritems():
            # get all nodes that have a child (might contribute to a coal) lexicographically "before" this one
            tops_with_child_in_s = lpvars._top_nodes_with_child[snode]
            prev_nodes = tops_with_child_in_s[:tops_with_child_in_s.index(gnode)]

            # get all previous path vars
            prev_paths = [lpvars.get_path(gnode, other) for other in prev_nodes]

            ilp += local_coal <= len(prev_nodes) - pulp.lpSum(prev_paths)
            ilp += len(prev_nodes) - pulp.lpSum(prev_paths) <=  len(lpvars._top_nodes[snode]) * local_coal

        # create delta constraints
        for g1, g2 in lpvars._delta_vars:
            # check if g1.parent < g2 and g2 < g1; set to 0 if neither, 1 if either, 2 if both
            g1_at_time_of_dup_at_g2 = lpvars.get_order(g1.parent, g2) + lpvars.get_order(g2, g1)
            g1_and_g2_at_same_locus = 1 - lpvars.get_path(g1.parent, g2.parent)

            ilp += lpvars._delta_vars[(g1, g2)] >= \
                    g1_at_time_of_dup_at_g2 + g1_and_g2_at_same_locus + lpvars.dup_vars[g2] - 3

        # create order constraints (transitive property)
        for gnodes in pulp.combination(all_gnodes, 3):
            g1, g2, g3 = gnodes
            if self.srecon[g1] == self.srecon[g2] == self.srecon[g3]:
                snode = self.srecon[g1]
                leaves = [g for g in gnodes if g in lpvars._bottom_nodes[snode]]
                non_leaves = [g for g in gnodes if g not in leaves]
                nleaves = len(leaves)

                if nleaves == 0:
                    ilp += lpvars.get_order(g1, g3) >= lpvars.get_order(g1, g2) + lpvars.get_order(g2, g3) - 1
                elif nleaves == 1:                    
                    ilp += lpvars.get_order(non_leaves[0], leaves[0]) == 1
                    ilp += lpvars.get_order(non_leaves[1], leaves[0]) == 1
                if nleaves == 2:
                    ilp += lpvars.get_order(non_leaves[0], leaves[0]) == 1
                    ilp += lpvars.get_order(non_leaves[0], leaves[1]) == 1
                else:
                    # transitive property is not neccessary if all gnodes are leaves
                    pass

        # create order constraints (duplication property, e.g. duplications as early as possible)
        #this is an optimization constraint which is not crucial --MORGAN
        if self.optimize:
            for (g1, g2), local_order in lpvars.order_vars.iteritems():
                if self.srecon[g1] == self.srecon[g2]:
                    snode = self.srecon[g1]
                    if (g1 in lpvars._bottom_nodes[snode]) and (g2 in lpvars._bottom_nodes[snode]):
                        ilp += local_order >= lpvars.dup_vars[g1] - lpvars.dup_vars[g2]
                        ilp += local_order <= lpvars.dup_vars[g1] - lpvars.dup_vars[g2] + 1

        # create r constraints
        for g2 in all_gnodes:
            gnodes_in_species = lpvars._gnodes_by_species[self.srecon[g2]]
            other_delta_vars_in_same_species = [lpvars._delta_vars[(g1, g2)] for g1 in gnodes_in_species
                                                if (g1, g2) in lpvars._delta_vars]
            ilp += lpvars._coaldup_vars[g2] >= pulp.lpSum(other_delta_vars_in_same_species) - 1
