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
              dupcost=1, losscost=1, coalcost=1, coaldupcost=None,
              time_limit=None, delay=True, mpr_constraints=True, log=sys.stdout): 
    """Perform reconciliation using DLCoal model with parsimony costs"""

    reconer = DLCLPRecon(tree, stree, gene2species, seed,
                         dupcost=dupcost, losscost=losscost, coalcost=coalcost, coaldupcost=coaldupcost,
                         time_limit=time_limit, delay=delay, mpr_constraints=mpr_constraints, log=log)
    return reconer.recon()


class DLCLPRecon(object):

    def __init__(self, gtree, stree, gene2species, seed,
                 dupcost=1, losscost=1, coalcost=1, coaldupcost=None, time_limit=None,  
                 delay=True, mpr_constraints=True,
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

        self.time_limit = time_limit
        self.mpr_constraints = mpr_constraints
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

            #set time_limit to CPLEX default
            if self.time_limit is None:
                self.time_limit = 1e+75

            ilpsolver = pulp.CPLEX_PY()        
            ilpsolver.buildSolverModel(ilp)

            #log results, used to check if seed has changed
            cplex_out_log = "/home/julia/Desktop/CS/projects/dlcpar/examples/results.dlcilp.cplex.log"
            ilpsolver.solverModel.set_results_stream(cplex_out_log)
            ilpsolver.solverModel.set_warning_stream(cplex_out_log)
            ilpsolver.solverModel.set_error_stream(cplex_out_log)
            ilpsolver.solverModel.set_log_stream(cplex_out_log)

            #set time_limit, seed, and output
            ilpsolver.solverModel.parameters.timelimit.set(int(self.time_limit))
            ilpsolver.solverModel.parameters.randomseed.set(self.seed)
            ilpsolver.solverModel.parameters.mip.display.set(5)

            self.log.start("Solving ilp")
            #ilp.solve(solver=ilpsolver) - this actually creates a new solver
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

        self.log.log("\nTopology Order Variables")
        for gtuple, order_var in lpvars._orders_from_tree.iteritems():
            self.log.log( "\t", gtuple, ": ", order_var) 

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

        sorted_leaves_by_species_keys = self._sort_keys(leaves_by_species.keys())

        #original: for leaves in leaves_by_species.itervalues():
        for key in sorted_leaves_by_species_keys:
            leaves = leaves_by_species[key]
            for g1, g2 in pulp.combination(leaves, 2):
                path1, path2 = common.find_path(g1, g2)
                path = path1 + path2
                nodes = [self.gtree[name] for name in path]
                ilp += pulp.lpSum([lpvars.dup_vars[node] for node in nodes]) >= 1
            
        # create duplication optimization constraints (from a gnode g with 2 children g' and g'', only 1 will have a duplication)
        #this is an optimization constraint which is not crucial --MORGAN
        if self.mpr_constraints:
            for gnode in all_gnodes:
                if not gnode.is_leaf():
                    assert len(gnode.children) == 2 #"ilprecon only takes binary gene trees"
                    ilp += pulp.lpSum([lpvars.dup_vars[node] for node in gnode.children]) <= 1
        # create the path constraints - if there is dup on given path, then that path var is 1, otherwise 0
        all_pairs = list(pulp.combination(all_gnodes, 2)) #key list
        for key in all_pairs:
            g1, g2 = key
            path_var = lpvars._path_vars[key]
            path1, path2 = common.find_path(g1, g2)
            path = path1 + path2
            nodes = [self.gtree[name] for name in path]
            ilp += path_var <= pulp.lpSum([lpvars.dup_vars[node] for node in nodes])
            ilp += pulp.lpSum([lpvars.dup_vars[node] for node in nodes]) <= len(path) * path_var
            
        # create loss constraints

        sorted_loss_keys = self._sort_keys(lpvars._loss_vars.keys())

        for (snode, gnode) in sorted_loss_keys:
            local_loss = lpvars._loss_vars[(snode, gnode)]

            local_helper = lpvars._helper_vars[(snode, gnode)]
            local_bottoms = lpvars._bottom_nodes[snode]
            local_paths = [lpvars.get_path(gnode, local_bottom) for local_bottom in local_bottoms]
            ilp += local_loss >= pulp.lpSum(local_paths) - len(local_bottoms) + 1 - local_helper

        # create helper constraints
       
        sorted_helper_keys = self._sort_keys(lpvars._helper_vars.keys())
       
        for (snode, gnode) in sorted_helper_keys:
            # get all top_nodes_at_snode that are lexicographically "before" this one
            local_helper = lpvars._helper_vars[(snode, gnode)]
            
            top_nodes_in_snode = lpvars._top_nodes[snode]
            prev_nodes = top_nodes_in_snode[:top_nodes_in_snode.index(gnode)]

            # get all previous path vars
            prev_paths = [lpvars.get_path(gnode, other) for other in prev_nodes]

            ilp += local_helper <= len(prev_nodes) - pulp.lpSum(prev_paths)
            ilp += len(prev_nodes) - pulp.lpSum(prev_paths) <= len(lpvars._top_nodes[snode]) * local_helper

        # create coal constraints

        sorted_coalspec_keys = self._sort_keys(lpvars._coalspec_vars.keys())

        for (snode, gnode) in sorted_coalspec_keys:
            # get all nodes that have a child (might contribute to a coal) lexicographically "before" this one
            local_coal = lpvars._coalspec_vars[(snode, gnode)]
            
            tops_with_child_in_s = lpvars._top_nodes_with_child[snode]
            prev_nodes = tops_with_child_in_s[:tops_with_child_in_s.index(gnode)]

            # get all previous path vars
            prev_paths = [lpvars.get_path(gnode, other) for other in prev_nodes]
            ilp += local_coal <= len(prev_nodes) - pulp.lpSum(prev_paths)
            ilp += len(prev_nodes) - pulp.lpSum(prev_paths) <=  len(lpvars._top_nodes[snode]) * local_coal

        # create delta constraints

        sorted_delta_keys = self._sort_keys(lpvars._delta_vars.keys())

        for g1, g2 in sorted_delta_keys:
        # check if g1.parent < g2 and g2 < g1; set to 0 if neither, 1 if either, 2 if both
            
            g1_at_time_of_dup_at_g2 = lpvars.get_order(g1.parent, g2) + lpvars.get_order(g2, g1)
            g1_and_g2_at_same_locus = 1 - lpvars.get_path(g1.parent, g2.parent)

            ilp += lpvars._delta_vars[g1,g2] >= \
                    g1_at_time_of_dup_at_g2 + g1_and_g2_at_same_locus + lpvars.dup_vars[g2] - 3
            ilp += lpvars._delta_vars[g1,g2] <= lpvars.get_order(g1.parent, g2)
            ilp += lpvars._delta_vars[g1,g2] <= lpvars.get_order(g2, g1)
            ilp += lpvars._delta_vars[g1,g2] <= lpvars.dup_vars[g2]
            ilp += lpvars._delta_vars[g1,g2] <= g1_and_g2_at_same_locus

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
        # this is an optimization constraint which is not crucial --MORGAN
        
        if self.mpr_constraints:
            sorted_order_keys = self._sort_keys(lpvars.order_vars.keys())
            for g1, g2 in sorted_order_keys:
                local_order = lpvars.order_vars[(g1, g2)]
                if self.srecon[g1] == self.srecon[g2]:
                    snode = self.srecon[g1]
                    if (g1 in lpvars._bottom_nodes[snode]) and (g2 in lpvars._bottom_nodes[snode]):
                        ilp += local_order >= lpvars.dup_vars[g1] - lpvars.dup_vars[g2]
                        ilp += local_order <= lpvars.dup_vars[g1] - lpvars.dup_vars[g2] + 1

        # # create zeta constraints
        
        sorted_zeta_keys = self._sort_keys(lpvars._zeta_vars.keys())
        
        for g1, g2 in sorted_zeta_keys:
            if (g1, g2) not in lpvars._orders_from_tree and (g2, g1) not in lpvars._orders_from_tree:
                zeta_val = lpvars._zeta_vars[(g1, g2)]
                ilp += zeta_val <= 1 - lpvars.dup_vars[g1]
                ilp += zeta_val <= lpvars.dup_vars[g2]
                ilp += zeta_val >=  lpvars.dup_vars[g2] + (1 - lpvars.dup_vars[g1]) - 1
                ilp += lpvars.get_order(g2, g1) >= zeta_val

        # create r constraints
        for g2 in all_gnodes:
            gnodes_in_species = lpvars._gnodes_by_species[self.srecon[g2]]
            other_delta_vars_in_same_species = [lpvars._delta_vars[(g1, g2)] for g1 in gnodes_in_species
                                                if (g1, g2) in lpvars._delta_vars]
            ilp += lpvars._coaldup_vars[g2] >= pulp.lpSum(other_delta_vars_in_same_species) - 1


    def _sort_keys(self, keys):      
        """sorts a list of keys"""
        for k in keys:
            # "bubble sort" the list
            for i in range(len(keys)):
                for j in range(0,len(keys)-i-1):
                    if str(keys[j]) > str(keys[j+1]):
                        keys[j], keys[j+1] = keys[j+1], keys[j]
        return keys

    def _log_var(self, lpvar_dict, key_type):
        """ logs dictionaries that are sorted in the same order as when they are used in the constraints
            key_type is a boolean to determine the type of key: 
            True is key = (snode, gnode), False otherwise """
        sorted_key_list = self._sort_keys(lpvar_dict.keys())
        if key_type:
            for (snode, gnode) in sorted_key_list:
                var = lpvar_dict[(snode, gnode)]
                self.log.log( "\t", gnode, "in", snode, ": ", var.varValue) 
        else:
            for key in sorted_key_list:
                var = lpvar_dict[key]
                self.log.log( "\t", key, ": ", var.varValue)

    def _find_sum(self, lpvar_dict):
        sum = 0
        for var in lpvar_dict.values():
            sum += var.varValue
        return sum