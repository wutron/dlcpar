"""
   Code for the DLC ILP Reconciliation
   (duplications, losses, and coalescence)
"""

import gzip

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

def ilp_recon(tree, stree, gene2species, seed,
              dupcost=1, losscost=1, coalcost=0.5, coaldupcost=None,
              time_limit=None, delay=True, log=sys.stdout): 
    """Perform reconciliation using DLCoal model with parsimony costs"""

    reconer = DLCLPRecon(tree, stree, gene2species, seed,
                         dupcost=dupcost, losscost=losscost, coalcost=coalcost, coaldupcost=coaldupcost,
                         time_limit=time_limit, delay=delay, log=log)
    return reconer.recon()


class DLCLPRecon(object):

    def __init__(self, gtree, stree, gene2species, seed,
                 dupcost=1, losscost=1, coalcost=1, coaldupcost=None, time_limit=None,  
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

        self.time_limit = time_limit
        if delay:
            raise Exception("delay=True not allowed")
        self.delay = delay

        self.name_internal = name_internal
        self.log = util.Timer(log)

        self.seed = seed

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

            ilpsolver = pulp.CPLEX_PY()        
            ilpsolver.buildSolverModel(ilp)

            #set time_limit and seed
            if self.time_limit:
                ilpsolver.solverModel.parameters.timelimit.set(int(self.time_limit))
            ilpsolver.solverModel.parameters.randomseed.set(self.seed)

            self.log.start("Solving ilp")
            ilpsolver.callSolver(ilp)
            ilpsolver.findSolutionValues(ilp)
        else:  
            options = ['randomSeed ' + str(self.seed), 'randomCbcSeed ' + str(self.seed)]
            self.log.start("Solving ilp")
            ilp.solve(pulp.solvers.PULP_CBC_CMD(maxSeconds=self.time_limit, options=options))
            
        solve_runtime = self.log.stop()
        self.log.log("Solver: " + solver_name)
        self.cost = pulp.value(ilp.objective)

        #print all variables to log file
        self.log.log("\n")
        self.log.log("Variables after solving")

        self.log.log("\nDuplication Variables (value = 1 if duplication on edge between gnode and gnode's parent, 0 otherwise)")
        self._log_var(lpvars.dup_vars, False)

        self.log.log("\nPath Variables (value = 1 if path between gene nodes has at least one duplication, 0 otherwise)")
        all_gnodes = list(self.gtree.preorder())
        all_pairs = list(pulp.combination(all_gnodes, 2)) #key list
        for gtuple in all_pairs:          
            _path_var = lpvars._path_vars[gtuple]
            self.log.log( "\t", gtuple, ": ", _path_var.varValue) 

        self.log.log("\nOrder Variables (value = 1 if second is more recent than first, 0 otherwise)")
        self._log_var(lpvars.order_vars, False)    

        self.log.log("\nLoss Variables (value = 1 if gnode creates a loss in snode, 0 otherwise)")
        self._log_var(lpvars._loss_vars, True)

        self.log.log("\nCoalescence at Speciation Variables (value = 1 if gnode on same locus as another gnode at top of snode with child, 0 otherwise)")
        self._log_var(lpvars._coalspec_vars, True)

        self.log.log("\nCoalescence at Duplication Delta Variables (delta_{g1,g2} variables, see paper for description)")
        self._log_var(lpvars._delta_vars, False)

        self.log.log("\nCoalescence at Duplication r Variables (number of coalescenses due to dup at g)")
        self._log_var(lpvars._coaldup_vars, False)           

        self.log.log("\nHelper Variables (value = 1 if g on same locus as gnode2 mapped to top of snode and gnode2 < gnode, 0 otherwise)")
        self._log_var(lpvars._helper_vars, True)

        self.log.log("\nOrders from Tree")
        for gtuple, order_var in lpvars._orders_from_tree.iteritems():
            self.log.log( "\t", gtuple, ": ", order_var) 

        #print optimal cost from ilp
        self.log.log("\nOptimal Cost:\t%f" % self.cost)
        self.log.log("Dups:\t%f" % sum(var.varValue for var in lpvars.dup_vars.values()))
        self.log.log("Losses:\t%f" % sum(var.varValue for var in lpvars._loss_vars.values()))
        self.log.log("CoalSpecs:\t%f" % sum(var.varValue for var in lpvars._coalspec_vars.values()))
        self.log.log("CoalDups:\t%f\n" % sum(var.varValue for var in lpvars._coaldup_vars.values()))

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
        
        #=============================
        # create dup constraints
        leaves = collections.defaultdict(list)
        for node, snode in self.srecon.iteritems():
            if node.is_leaf():
                leaves[snode].append(node)
            leaves[snode].sort(key=lambda node: node.name)

        for snode in self.stree.preorder():
            for g1, g2 in pulp.combination(leaves[snode], 2):
                path1, path2 = common.find_path(g1, g2)
                path = path1 + path2
                nodes = [self.gtree[name] for name in path] 
                ilp += pulp.lpSum([lpvars.dup_vars[node] for node in nodes]) >= 1
        
        #=============================
        # create the path constraints - if there is dup on given path, then that path var is 1, otherwise 0
        all_pairs = list(pulp.combination(all_gnodes, 2)) #key list
        for (g1, g2) in all_pairs:
            path_var = lpvars._path_vars[(g1, g2)]
            path1, path2 = common.find_path(g1, g2)
            path = path1 + path2
            nodes = [self.gtree[name] for name in path]
            ilp += path_var <= pulp.lpSum([lpvars.dup_vars[node] for node in nodes])
            ilp += pulp.lpSum([lpvars.dup_vars[node] for node in nodes]) <= len(path) * path_var
            
        # create loss constraints

        for (snode, gnode) in sorted(lpvars._loss_vars.keys(), key = lambda node: str(node)):
            local_loss = lpvars._loss_vars[(snode, gnode)]

            local_helper = lpvars._helper_vars[(snode, gnode)]
            local_bottoms = lpvars._bottom_nodes[snode]
            local_paths = [lpvars.get_path(gnode, local_bottom) for local_bottom in local_bottoms]
            ilp += local_loss >= pulp.lpSum(local_paths) - len(local_bottoms) + 1 - local_helper 
            
            #ilp += len(local_bottoms) + local_helper >= pulp.lpSum(local_paths) + 1 - len(local_bottoms) * (local_loss)

        # create helper constraints
       
        for (snode, gnode) in sorted(lpvars._helper_vars.keys(), key = lambda node: str(node)):
            # get all top_nodes_at_snode that are lexicographically "before" this one
            local_helper = lpvars._helper_vars[(snode, gnode)]
            
            top_nodes_in_snode = lpvars._top_nodes[snode]
            prev_nodes = top_nodes_in_snode[:top_nodes_in_snode.index(gnode)]

            # get all previous path vars
            prev_paths = [lpvars.get_path(gnode, other) for other in prev_nodes]

            ilp += local_helper <= len(prev_nodes) - pulp.lpSum(prev_paths)
            ilp += len(prev_nodes) - pulp.lpSum(prev_paths) <= len(lpvars._top_nodes[snode]) * local_helper

        # create coal constraints

        for (snode, gnode) in sorted(lpvars._coalspec_vars.keys(), key = lambda node: str(node)):
            # get all nodes that have a child (might contribute to a coal) lexicographically "before" this one
            local_coal = lpvars._coalspec_vars[(snode, gnode)]
            
            tops_with_child_in_s = lpvars._top_nodes_with_child[snode]
            prev_nodes = tops_with_child_in_s[:tops_with_child_in_s.index(gnode)]

            # get all previous path vars
            prev_paths = [lpvars.get_path(gnode, other) for other in prev_nodes]
            ilp += local_coal <= len(prev_nodes) - pulp.lpSum(prev_paths)
            ilp += len(prev_nodes) - pulp.lpSum(prev_paths) <=  len(lpvars._top_nodes[snode]) * local_coal

        # create delta constraints

        for g1, g2 in sorted(lpvars._delta_vars.keys(), key = lambda node: str(node)):
            
            g1_at_time_of_dup_at_g2 = lpvars.get_order(g1.parent, g2, False) + lpvars.get_order(g2, g1, False)
            g1_and_g2_at_same_locus = 1 - lpvars.get_path(g1.parent, g2.parent)

            ilp += lpvars._delta_vars[g1,g2] >= \
                    g1_at_time_of_dup_at_g2 + g1_and_g2_at_same_locus + lpvars.dup_vars[g2] - 3
            ilp += lpvars._delta_vars[g1,g2] <= lpvars.get_order(g1.parent, g2, False)
            ilp += lpvars._delta_vars[g1,g2] <= lpvars.get_order(g2, g1, False)
            ilp += lpvars._delta_vars[g1,g2] <= lpvars.dup_vars[g2]
            ilp += lpvars._delta_vars[g1,g2] <= g1_and_g2_at_same_locus

        # create order constraints (transitive property)

        for gnodes in pulp.combination(all_gnodes, 3):
            g1, g2, g3 = gnodes
            if not (self.srecon[g1] == self.srecon[g2] == self.srecon[g3]):
                continue
            
            snode = self.srecon[g1]
            leaves = [g for g in gnodes if g in lpvars._bottom_nodes[snode]]
            non_leaves = [g for g in gnodes if g not in leaves]
            nleaves = len(leaves)
            assert (nleaves >=0 and nleaves < 4), ("Impossible Number of Leaves", nleaves)

            #transitivity of order
            ilp += lpvars.get_order(g1, g3, False) >= lpvars.get_order(g1, g2, False) + lpvars.get_order(g2, g3, False) - 1
            if nleaves == 1:                    
                ilp += lpvars.get_order(non_leaves[0], leaves[0], False) == 1
                ilp += lpvars.get_order(non_leaves[1], leaves[0], False) == 1
            if nleaves == 2:
                ilp += lpvars.get_order(non_leaves[0], leaves[0], False) == 1
                ilp += lpvars.get_order(non_leaves[0], leaves[1], False) == 1

        # create zeta constraints
                
        for g1, g2 in sorted(lpvars._zeta_vars.keys(), key = lambda node: str(node)):
            if (g1, g2) not in lpvars._orders_from_tree and (g2, g1) not in lpvars._orders_from_tree:
                zeta_val = lpvars._zeta_vars[(g1, g2)]
                ilp += zeta_val <= 1 - lpvars.dup_vars[g1]
                ilp += zeta_val <= lpvars.dup_vars[g2]
                ilp += zeta_val >=  lpvars.dup_vars[g2] + (1 - lpvars.dup_vars[g1]) - 1
                ilp += lpvars.get_order(g2, g1, False) >= zeta_val

        # create r constraints
        for g2 in all_gnodes:
            species_nodes = lpvars._gnodes[self.srecon[g2]] + lpvars._cnodes[self.srecon[g2]]
            other_delta_vars_in_same_species = [lpvars._delta_vars[(g1, g2)] for g1 in species_nodes
                                                if (g1, g2) in lpvars._delta_vars]
            ilp += lpvars._coaldup_vars[g2] >= pulp.lpSum(other_delta_vars_in_same_species) - 1
            # TODO: less than or equal constraint for r variable


    def _log_var(self, lpvar_dict, key_type):
        """ logs dictionaries that are sorted in the same order as when they are used in the constraints
            key_type is a boolean to determine the type of key: 
            True is key = (snode, gnode), False otherwise """
        
        sorted_key_list = sorted(lpvar_dict.keys(), key = lambda node: str(node))
        
        if key_type:
            for (snode, gnode) in sorted_key_list:
                var = lpvar_dict[(snode, gnode)]
                self.log.log( "\t", gnode, "in", snode, ": ", var.varValue) 
        else:
            for key in sorted_key_list:
                var = lpvar_dict[key]
                self.log.log( "\t", key, ": ", var.varValue)