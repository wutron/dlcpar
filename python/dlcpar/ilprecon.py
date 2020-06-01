"""

   Code to solve the DLC MPR Problem using ILP

"""

# python libraries
import os, sys, shutil
import gzip
import collections

# rasmus, compbio libraries
from rasmus import treelib, util
from compbio import phylo

# dlcpar libraries
from dlcpar import common
from dlcpar import reconlib

# pulp libraries
import pulp
import pulp.apis
try:
    import cplex
except:
    pass

from dlcpar import ilpreconlib

#==========================================================
# globals

ILPNAME = "dlcilp"

#==========================================================

def dlc_recon(tree, stree, gene2species,
              dupcost=1, losscost=1, coalcost=1, coaldupcost=None,
              implied=True, delay=True,
              solver="CBC_CMD", seed=None,
              time_limit=None, mem_limit=None, num_threads=None,
              log=sys.stdout, info_log=sys.stdout, tmp=None):
    """Perform reconciliation using DLCoal model with parsimony costs"""

    reconer = DLCRecon(tree, stree, gene2species,
                       dupcost=dupcost, losscost=losscost, coalcost=coalcost, coaldupcost=coaldupcost,
                       implied=implied, delay=delay,
                       solver=solver, seed=seed,
                       time_limit=time_limit, mem_limit=mem_limit, num_threads=num_threads,
                       log=log, info_log=info_log, tmp=tmp)
    return reconer.recon()



class DLCRecon(object):

    def __init__(self, gtree, stree, gene2species,
                 dupcost=1, losscost=1, coalcost=0.5, coaldupcost=None,
                 implied=True, delay=True,
                 solver="CBC_CMD", seed=None,
                 time_limit=None, mem_limit=None, num_threads=None,
                 name_internal="n",
                 log=sys.stdout, info_log=sys.stdout, tmp=None):

        # rename gene tree nodes
        common.rename_nodes(gtree, name_internal)

        self.gtree = gtree
        self.stree = stree
        self.gene2species = gene2species

        assert (dupcost > 0) and (losscost >= 0) and\
               (coalcost > 0) and (coaldupcost is None or coaldupcost > 0)
        self.dupcost = dupcost
        self.losscost = losscost
        self.coalcost = coalcost  # actually coalspeccost, using coalcost for backwards compatibility
        self.coaldupcost = coaldupcost if coaldupcost is not None else coalcost

        self.solver = solver
        self.seed = seed
        self.time_limit = time_limit
        self.mem_limit = mem_limit
        self.num_threads = num_threads

        if not implied:
            raise Exception("implied=False not allowed")
        if delay:
            raise Exception("delay=True not allowed")
        self.delay = delay

        self.name_internal = name_internal
        self.log = util.Timer(log)
        self.info_log = util.Timer(info_log)
        self.tmp = tmp

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
        self._infer_species_map()
        self.log.log("\n\n")

        # add implied speciation nodes but first start the species tree at the right root
        substree = treelib.subtree(self.stree, self.srecon[self.gtree.root])
        subsrecon = util.mapdict(self.srecon, val=lambda snode: substree.nodes[snode.name])

        # switch internal storage with subtrees
        self.stree, substree = substree, self.stree
        self.srecon, subsrecon = subsrecon, self.srecon

        # add implied nodes (standard speciation, speciation from duplication, delay nodes)
        # then relabel events (so that factor_tree works)
        reconlib.add_implied_nodes(self.gtree, self.stree, self.srecon, self.sevents, delay=self.delay)
        self.sevents = phylo.label_events(self.gtree, self.srecon)
        common.rename_nodes(self.gtree, self.name_internal)

        # log gene tree (with species map)
        self.log.log("gene tree (with species map)\n")
        reconlib.log_tree(self.gtree, self.log, func=reconlib.draw_tree_recon, srecon=self.srecon)

        # infer locus map
        ilp, lpvars, runtime_setup, runtime_solve = self._infer_locus_map()
        self.log.log("\n\n")

        # convert to LabeledRecon data structure
        self.log.start("Converting to LCT")
        _, labeled_recon = ilpreconlib.ilp_to_lct(self.gtree, lpvars)
        self.log.stop()

        # log gene tree (with species map and locus map)
        self.log.log("gene tree (with species and locus map)\n")
        self.lrecon = labeled_recon.locus_map
        self.order = labeled_recon.order
        reconlib.log_tree(self.gtree, self.log, func=reconlib.draw_tree_recon, srecon=self.srecon, lrecon=self.lrecon)

        # revert to use input species tree
        self.stree = substree
        self.srecon = util.mapdict(self.srecon, val=lambda snode: self.stree.nodes[snode.name])
        self.order = util.mapdict(self.order, key=lambda snode: self.stree.nodes[snode.name])

        # calculate runtime
        runtime = self.log.stop()

        return self.gtree, labeled_recon, self.cost, runtime, runtime_setup, runtime_solve


    def _infer_species_map(self):
        """Infer (and assign) species map"""

        self.log.start("Inferring species map")
        self.srecon = phylo.reconcile(self.gtree, self.stree, self.gene2species)
        self.sevents = phylo.label_events(self.gtree, self.srecon)
        self.log.stop()


    def _infer_locus_map(self):
        """Infer (and assign) locus map"""

        self.log.start("Inferring locus map")

        # create ILP
        self.log.start("Creating ilp variables")
        ilp = pulp.LpProblem(ILPNAME, pulp.LpMinimize)
        lpvars = ilpreconlib.IlpReconVariables(self.gtree, self.stree, self.srecon)
        runtime_create_vars = self.log.stop()

        self.log.start("Building ilp constraints")
        ilp += self._create_objective_func(lpvars)
        self._add_constraints(ilp, lpvars)
        runtime_create_constraints = self.log.stop()

        runtime_setup = runtime_create_vars + runtime_create_constraints

        # solve ILP
        if self.solver == "CPLEX_PY":
            ilpsolver = pulp.CPLEX_PY()
            ilpsolver.buildSolverModel(ilp)

            # log stdout
            if self.tmp:
                cplex_out_log = os.path.join(self.tmp, ILPNAME + ".mip.display")
                ilpsolver.solverModel.set_results_stream(cplex_out_log)
                ilpsolver.solverModel.set_warning_stream(cplex_out_log)
                ilpsolver.solverModel.set_error_stream(cplex_out_log)
                ilpsolver.solverModel.set_log_stream(cplex_out_log)
                ilpsolver.solverModel.parameters.mip.display.set(5)
            else:
                ilpsolver.solverModel.parameters.mip.display.set(0)

            # set parameters
            if self.seed:
                ilpsolver.solverModel.parameters.randomseed.set(self.seed)
            if self.time_limit:
                ilpsolver.solverModel.parameters.timelimit.set(int(self.time_limit))
            if self.mem_limit:
                ilpsolver.solverModel.parameters.workmem.set(int(self.mem_limit))
            if self.num_threads:
                ilpsolver.solverModel.parameters.threads.set(self.num_threads)
            if self.tmp:
                ilp.writeLP(os.path.join(self.tmp, ILPNAME + ".lp"))

            # solve ilp
            self.log.start("Solving ilp")
            ilpsolver.callSolver(ilp)
            runtime_solve = self.log.stop()

            if self.tmp:
                # write solution html file to temp dir
                ilpsolver.solverModel.solution.write(
                    os.path.join(self.tmp, ILPNAME + "-pulp.sol"))

            ilpsolver.findSolutionValues(ilp)

        else:
            # set options
            options = []
            if self.seed:
                options.append('randomCbcSeed ' + str(self.seed))
            if self.mem_limit:
                raise Exception("cbc memory limit not implemented")

            # log mps and sol files if log option is selected
            if not self.tmp:
                keepFiles=False
            else:
                keepFiles=True
                ilp.writeLP(os.path.join(self.tmp, ILPNAME + ".lp"))

            # solve
            self.log.start("Solving ilp")
            ilp.solve(pulp.apis.PULP_CBC_CMD(maxSeconds=self.time_limit, threads=self.num_threads,
                                             options=options, keepFiles=keepFiles))
            runtime_solve = self.log.stop()

            # move logs
            if self.tmp:
                shutil.move(ILPNAME + "-pulp.mps", self.tmp)
                shutil.move(ILPNAME + "-pulp.sol", self.tmp)

        # round variable values to ensure binary integrality
        self._round_variables(lpvars)

        # get minimum cost
        self.cost = pulp.value(ilp.objective)

        # log solver and status
        self.log.log("Solver: %s" % self.solver)
        self.info_log.log("Solver: %s" % self.solver)
        self.log.log("Solver status:\t%s" % pulp.LpStatus[ilp.status])
        self.info_log.log("Solver status:\t%s" % pulp.LpStatus[ilp.status])
        self.log.log("\n\n")

        # print all variables to log file
        self.log.log("Solved Variables")
        self._log_vars(lpvars)

        # log optimal cost from ilp
        self.log.log()
        self.log.log("Optimal Cost:\t%f" % self.cost)
        self.log.log("Dups:\t%d" % sum(var.varValue for var in lpvars.dup_vars.values()))
        self.log.log("Losses:\t%d" % sum(var.varValue for var in lpvars._loss_vars.values()))
        self.log.log("CoalSpecs:\t%d" % sum(var.varValue for var in lpvars._coalspec_vars.values()))
        self.log.log("CoalDups:\t%d" % sum(var.varValue for var in lpvars._coaldup_vars.values()))
        self.log.log()

        self.log.stop()

        return ilp, lpvars, runtime_setup, runtime_solve


    def _create_objective_func(self, lpvars):
        """Return cost function for ILP formulation"""

        num_dups = pulp.lpSum(lpvars.dup_vars.values())
        num_losses = pulp.lpSum(lpvars._loss_vars.values())
        num_coals = pulp.lpSum(lpvars._coalspec_vars.values())
        num_coaldups = pulp.lpSum(lpvars._coaldup_vars.values())

        return self.dupcost * num_dups + self.losscost * num_losses \
            + self.coalcost * num_coals + self.coaldupcost * num_coaldups


    def _add_constraints(self, ilp, lpvars):
        """Add constraints for ILP formulation"""

        gtree = self.gtree
        stree = self.stree
        srecon = self.srecon

        all_gnodes = list(self.gtree.preorder())


        #=============================
        # dup constraints

        leaves = collections.defaultdict(list)
        for node, snode in srecon.iteritems():
            if node.is_leaf():
                leaves[snode].append(node)
            leaves[snode].sort(key=lambda node: node.name)

        # leaves in same species are paralogous
        for snode in stree:
            for g1, g2 in pulp.combination(leaves[snode], 2):
                path1, path2 = common.find_path(g1, g2)
                path = path1 + path2
                ilp += pulp.lpSum([lpvars.dup_vars[gname] for gname in path]) >= 1

                # self.log.log("dup constraint s:", snode.name, "\tpulp_sum:", pulp.lpSum([lpvars.dup_vars[gname] for gname in path]), ">=", 1)


        #=============================
        # order constraints

        # transitivity of order
        for snode in stree:
            for (gname1, gname2, gname3) in pulp.permutation(lpvars._names[snode.name], 3):
                if ((gname1, gname3) not in lpvars.order_vars.keys()) or \
                   ((gname1, gname2) not in lpvars.order_vars.keys()) or \
                   ((gname2, gname3) not in lpvars.order_vars.keys()):
                    continue

                ilp += lpvars.order_vars[gname1, gname3] >= lpvars.order_vars[gname1, gname2] + lpvars.order_vars[gname2, gname3] - 1

                # self.log.log("gnode triple: ", (gname1, gname2, gname3))
                # self.log.log("transitivity constraint: ", lpvars.order_vars[gname1, gname3], " >= ", lpvars.order_vars[gname1, gname2], " + ", lpvars.order_vars[gname2, gname3], " - 1")

        # opposite orders
        for gname1, gname2 in sorted(lpvars.order_vars.keys(), key=lambda x: str(x)):
            ilp += lpvars.order_vars[gname1, gname2] == 1 - lpvars.order_vars[gname2, gname1]

            # self.log.log("opposite orders: ", lpvars.order_vars[gname1, gname2], " == 1 - ", lpvars.order_vars[gname2, gname1])

        # orders from topology
        for gname1, gname2 in lpvars._orders_from_tree:
            ilp += lpvars.order_vars[gname1, gname2] == 1.0
            ilp += lpvars.order_vars[gname2, gname1] == 0.0

            # self.log.log("orders from tree:", lpvars.order_vars[gname1, gname2], "== 1")
            # self.log.log("orders from tree:", lpvars.order_vars[gname2, gname1], "== 0")


        #=============================
        # loss constraints

        for (sname, gname) in sorted(lpvars._loss_vars.keys(), key=lambda x: str(x)):
            local_loss = lpvars._loss_vars[sname, gname]
            local_lambda = lpvars._lambda_vars[sname, gname]
            local_bottoms = lpvars._bottom_names[sname]
            local_paths = [lpvars._path_vars[gname, local_bottom] for local_bottom in local_bottoms]

            ilp += local_loss >= 1 - local_lambda + pulp.lpSum(local_paths) - len(local_bottoms)

            # self.log.log("loss local_loss: ", local_loss, "\ts: ", sname, "\tg: ", gname)
            # self.log.log("loss constraint: ", local_loss, ">=", "1 -", local_lambda ", +, " pulp.lpSum(local_paths), " - ", len(local_bottoms)


        #=============================
        # path constraints

        for g1, g2 in pulp.combination(all_gnodes, 2):
            gname1 = g1.name
            gname2 = g2.name
            local_path = lpvars._path_vars[gname1, gname2]

            # find path
            path1, path2 = common.find_path(g1, g2)
            path = path1 + path2

            # consistent with duplication variables
            ilp += local_path <= pulp.lpSum([lpvars.dup_vars[gname] for gname in path])
            ilp += pulp.lpSum([lpvars.dup_vars[gname] for gname in path]) <= len(path) * local_path

            # path equality constraints
            # (not in paper - paper uses permutations)
            ilp += local_path == lpvars._path_vars[g2.name, g1.name]

            # self.log.log("path constraint gtuple: ", (g1.name, g2.name))
            # self.log.log("path constraint one: ", path_var, " <= ", pulp.lpSum([lpvars.dup_vars[gname] for gname in path]))
            # self.log.log("path constraint two: ", pulp.lpSum([lpvars.dup_vars[gname] for gname in path]), " <= ", len(path), "*", path_var)

        # create path constraint that no path from a gnode to itself can have a duplication
        # (not in paper - required in pulp for lpSum over empty list)
        for g in all_gnodes:
            ilp += lpvars._path_vars[g.name, g.name] == 0


        #=============================
        # lambda constraints

        for (sname, gname) in sorted(lpvars._lambda_vars.keys(), key=lambda x: str(x)):
            local_lambda = lpvars._lambda_vars[sname, gname]

            # all top nodes of snode that precede this one
            tops = lpvars._top_names[sname]
            prev = tops[:tops.index(gname)]
            prev_paths = [lpvars._path_vars[gname, other] for other in prev]

            ilp += local_lambda <= len(prev) - pulp.lpSum(prev_paths)
            ilp += len(prev) - pulp.lpSum(prev_paths) <= len(tops) * local_lambda

            # self.log.log("lambda local_lambda: ", local_lambda, "\ts: ", sname, "\tg: ", gname)
            # self.log.log("lambda top_nodes: ", tops)
            # self.log.log("lambda prev_nodes: ", prev)
            # self.log.log("lambda constraint one: ", local_lambda, " <= ", len(prev_names), " - ", pulp.lpSum(prev_paths))
            # self.log.log("lambda constraint two: ", len(prev_names), " - ", pulp.lpSum(prev_paths), " <= ", len(top_names), "*", local_lambda)


        #=============================
        # coalspec constraints

        for (sname, gname) in sorted(lpvars._coalspec_vars.keys(), key=lambda x: str(x)):
            local_coal = lpvars._coalspec_vars[sname, gname]

            # all survived nodes ofsnode that precede this one
            survived = lpvars._survived_names[sname]
            prev = survived[:survived.index(gname)]
            prev_paths = [lpvars._path_vars[gname, other] for other in prev]

            ilp += len(prev) - pulp.lpSum(prev_paths) <=  len(survived) * local_coal

            # self.log.log("coalspec local_coalspec: ", local_coal, "\ts: ", sname, "\tg: ", gname)
            # self.log.log("coalspec survived_nodes: ", survived)
            # self.log.log("coalspec prev_nodes: ", prev)
            # self.log.log("coalspec constraint one: ", local_coal, " <= ", len(prev), " - ", pulp.lpSum(prev_paths))
            # self.log.log("coalspec constraint two: ", len(prev), "-", pulp.lpSum(prev_paths), " <= ", len(survived), "*", local_coal)


        #=============================
        # omega constraints

        for gname1, gname2 in sorted(lpvars._omega_vars.keys(), key=lambda x: str(x)):
            local_omega = lpvars._omega_vars[gname1, gname2]

            ilp += local_omega <= 1 - lpvars.dup_vars[gname1]
            ilp += local_omega <= lpvars.dup_vars[gname2]
            ilp += local_omega >=  lpvars.dup_vars[gname2] + (1 - lpvars.dup_vars[gname1]) - 1
            ilp += lpvars.order_vars[gname2, gname1] >= local_omega

            # self.log.log("omega gtuple: ", (gname1, gname2))
            # self.log.log("omega constraint omega_val: ", omega_val, " lpvars.dup_vars[g2]: ", lpvars.dup_vars[gname2], \
            #                                                         " lpvars.dup_vars[g1]: ", lpvars.dup_vars[gname1], \
            #                                                         " lpvars.order_vars[g2, g1] ", lpvars.order_vars[gname2, gname1])

        #=============================
        # kappa constraints

        for gname1, gname2 in sorted(lpvars._kappa_vars.keys(), key=lambda x: str(x)):
            pname1 = gtree.nodes[gname1].parent.name
            pname2 = gtree.nodes[gname2].parent.name

            g1_and_g2_at_same_locus = 1 - lpvars._path_vars[pname1, pname2]
            g1_at_time_of_dup_at_g2 = lpvars.order_vars[pname1, gname2] + lpvars.order_vars[gname2, gname1]

            ilp += lpvars._kappa_vars[gname1, gname2] <= lpvars.dup_vars[gname2]
            ilp += lpvars._kappa_vars[gname1, gname2] <= g1_and_g2_at_same_locus
            ilp += lpvars._kappa_vars[gname1, gname2] <= lpvars.order_vars[pname1, gname2]
            ilp += lpvars._kappa_vars[gname1, gname2] <= lpvars.order_vars[gname2, gname1]
            ilp += lpvars._kappa_vars[gname1, gname2] >= \
                lpvars.dup_vars[gname2] + g1_and_g2_at_same_locus + g1_at_time_of_dup_at_g2 - 3

            # self.log.log("kappa gtuple: ", (gname1, gname2))
            # self.log.log("kappa constraint: ", lpvars._kappa_vars[gname1, gname2], ">=", \
            #     lpvars.dup_vars[gname2], "+", g1_and_g2_at_same_locus, "+", g1_at_time_of_dup_at_g2, "- 3")

        #=============================
        # k_g constraints

        for g2 in all_gnodes:
            gname2 = g2.name
            sname = srecon[g2].name
            other = lpvars._names[sname] + lpvars._bottom_children_names[sname]
            other_kappa_vars = [lpvars._kappa_vars[gname1, gname2] for gname1 in other if (gname1, gname2) in lpvars._kappa_vars]

            ilp += lpvars._coaldup_vars[gname2] >= pulp.lpSum(other_kappa_vars) - 1

            # self.log.log("coaldup g2: ", g2, " constraint: ", lpvars._coaldup_vars[gname2], " >= ", pulp.lpSum(other_kappa_vars), " - 1")


    def _round_variables(self, lpvars):
        # ensure every variable is binary or an int (for _coaldup_vars)

        lpvar_dicts = lpvars._get_binary_dicts()
        for lpvar_dict in lpvar_dicts:
            for key in lpvar_dict.keys():
                if lpvar_dict[key].varValue > 0.5:
                    lpvar_dict[key].varValue = 1.0
                else:
                    lpvar_dict[key].varValue = 0.0

        for key in lpvars._coaldup_vars.keys():
            lpvars._coaldup_vars[key].varValue = round(lpvars._coaldup_vars[key].varValue)


    def _log_vars_helper(self, lpvar_dict, key_type):
        """Logs dictionary that are sorted in the same order as when they are used in the constraints
           key_type is an integer to determine the type of key:
               0 : key = gname
               1 : key = (gname1, gname2)
               2 : key = (sname, gname)"""

        sorted_key_list = sorted(lpvar_dict.keys(), key=lambda x: str(x))

        if key_type == 0:
            for key in sorted_key_list:
                var = lpvar_dict[key]
                self.log.log( "\t", key, ": ", var.varValue)
                assert var.varValue.is_integer(), (var.varValue, " must be an integer value")

        elif key_type == 1:
            for (gname1, gname2) in sorted_key_list:
                var = lpvar_dict[gname1, gname2]
                self.log.log( "\t", "(%s, %s)" % (gname1, gname2), ": ", var.varValue)
                assert var.varValue.is_integer(), (var.varValue, " must be an integer value")

        elif key_type == 2:
            for (sname, gname) in sorted_key_list:
                var = lpvar_dict[sname, gname]
                self.log.log( "\t", gname, "in", sname, ": ", var.varValue)
                assert var.varValue.is_integer(), (var.varValue, " must be an integer value")

        else:
            raise Exception("invalid key_type")



    def _log_vars(self, lpvars):
        """Logs dictionaries"""

        gtree = self.gtree

        self.log.log("Duplication Variables")
        self._log_vars_helper(lpvars.dup_vars, 0)
        self.log.log("\n")

        self.log.log("Order Variables")
        self._log_vars_helper(lpvars.order_vars, 1)
        self.log.log("\n")

        self.log.log("Omega Variables")
        self._log_vars_helper(lpvars._omega_vars, 1)
        self.log.log("\n")

        self.log.log("Loss Variables")
        self._log_vars_helper(lpvars._loss_vars, 2)
        self.log.log("\n")

        self.log.log("Path Variables")
        self._log_vars_helper(lpvars._path_vars, 1)
        self.log.log("\n")

        self.log.log("Lambda Variables")
        self._log_vars_helper(lpvars._lambda_vars, 2)
        self.log.log("\n")

        self.log.log("Coalescence at Speciation Variables")
        self._log_vars_helper(lpvars._coalspec_vars, 2)
        self.log.log("\n")

        self.log.log("Kappa Variables")
        self._log_vars_helper(lpvars._kappa_vars, 1)
        self.log.log("\n")

        self.log.log("Coalescence at Duplication Variables")
        self._log_vars_helper(lpvars._coaldup_vars, 0)
        self.log.log("\n")

        self.log.log("Orders from Tree")
        for (gname1, gname2) in sorted(lpvars._orders_from_tree, key=lambda gtuple: str(gtuple)):
            self.log.log( "\t", "(%s, %s)" % (gname1, gname2), ": ", lpvars.order_vars[gname1, gname2].varValue)
        self.log.log("\n")

        self.log.log("\nOrders not from Tree")
        for (gname1, gname2) in sorted(lpvars.order_vars.keys(), key=lambda x: str(x)):
            if ((gname1, gname2) not in lpvars._orders_from_tree) and ((gname2, gname1) not in lpvars._orders_from_tree):
                self.log.log( "\t", "(%s, %s)" % (gname1, gname2), ": ", lpvars.order_vars[gname1, gname2].varValue)
        self.log.log("\n")
