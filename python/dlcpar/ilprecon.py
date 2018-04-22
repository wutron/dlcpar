"""

   Code for the DLC ILP Reconciliation
   (duplications, losses, and coalescence)

"""

# python libraries
import sys
import random
import collections
import StringIO
import itertools
import recon

# rasmus libraries
from rasmus import treelib, util
from compbio import phylo

# dlcpar libraries
from dlcpar import common
from dlcpar import reconlib

# numpy
import numpy as np

# integer linear programming
from pulp import *

from dlcpar.recon import *
from recon import _random_choice

# gtree = G (with implied speciation nodes)
# stree = S
# srecon = M
# lrecon = L

#==========================================================
# globals

INF = util.INF
INIT_LOCUS = 1

import operator
def prod(iterable):
    return reduce(operator.mul, iterable, 1)

#==========================================================$
# utilities

def _random_choice(a, p=None):
    # wrapper around numpy.random.choice
    # note: numpy cannot take list of lists so use indices
    a = list(a)
    ndx = np.random.choice(range(len(a)), p=p)
    it = a[ndx]
    return it

#==========================================================

def ilp_recon(tree, stree, gene2species, gene2locus=None,
              dupcost=1, losscost=1, coalcost=1, coaldupcost=None,
              implied=True, delay=True,
              prescreen=False, prescreen_min=INF, prescreen_factor=INF,
              max_loci=INF, max_dups=INF, max_losses=INF, allow_both=False,
              log=sys.stdout):
    """Perform reconciliation using DLCoal model with parsimony costs"""

    reconer = DLCLPRecon(tree, stree, gene2species, gene2locus,
                       dupcost=dupcost, losscost=losscost, coalcost=coalcost, coaldupcost=coaldupcost,
                       implied=implied, delay=delay,
                       prescreen=prescreen, prescreen_min=prescreen_min, prescreen_factor=prescreen_factor,
                       max_loci=max_loci, max_dups=max_dups, max_losses=max_losses, allow_both=allow_both,
                       log=log)
    return reconer.recon()



class DLCLPRecon(DLCRecon):

    def __init__(self, gtree, stree, gene2species, gene2locus=None,
                 dupcost=1, losscost=1, coalcost=1, coaldupcost=None,
                 implied=True, delay=True,
                 prescreen=False, prescreen_min=INF, prescreen_factor=INF,
                 max_loci=INF, max_dups=INF, max_losses=INF, allow_both=False,
                 name_internal="n", log=sys.stdout):

        # rename gene tree nodes
        common.rename_nodes(gtree, name_internal)

        self.gtree = gtree
        self.stree = stree
        self.gene2species = gene2species
        self.gene2locus = gene2locus

        assert (dupcost >= 0) and (losscost >= 0) and (coalcost >= 0) and (coaldupcost >= 0), (dupcost, losscost, coalcost, coaldupcost)
        self.dupcost = dupcost
        self.losscost = losscost
        self.coalcost = coalcost  # actually coalspeccost, using coalcost for backwards compatibility
        self.coaldupcost = coaldupcost if coaldupcost is not None else coalcost

        self.implied = implied
        self.delay = delay

        assert (prescreen_min > 0) and (prescreen_factor >= 1)
        self.prescreen = prescreen
        self.prescreen_min = prescreen_min
        self.prescreen_factor = prescreen_factor

        assert (max_loci > 0) and (max_dups > 0) and (max_losses > 0)
        self.max_loci = max_loci
        self.max_dups = max_dups
        self.max_losses = max_losses
        self.allow_both = allow_both

        self.name_internal = name_internal
        self.log = util.Timer(log)

        # these attributes are assigned when performing reconciliation using self.recon()
        #   self.ilp
        #   self.srecon
        #   self.lrecon
        #   self.order
        #   self.nsoln
        #   self.cost

    #=============================
    # main methods

    def get_g_node(self, name_part):
        return [g for g in list(self.gtree.preorder()) if name_part in str(g)][0]


    def recon(self):
        """Perform reconciliation"""

        """
        self.log.start("Checking feasibility")
        feasible = self._are_constraints_consistent()
        self.log.stop()
        if not feasible:
            self.log.log("tree not feasible")
            return self.gtree, None
        """

        self.log.start("Reconciling")

        # log input gene and species trees
        #self.log.log("gene tree\n")
        #log_tree(self.gtree, self.log, func=treelib.draw_tree_names)
        #self.log.log("species tree\n")
        #log_tree(self.stree, self.log, func=treelib.draw_tree_names)

        # infer species map
        self._infer_species_map()
        #self.log.log("\n\n")

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
        #self.log.log("gene tree (with species map)\n")
        #log_tree(self.gtree, self.log, func=draw_tree_srecon, srecon=self.srecon)

        #treelib.draw_tree(self.gtree)

        # find subtrees
        self.subtrees = reconlib.factor_tree(self.gtree, self.stree, self.srecon, self.sevents)
        # compute the leaves of those subtrees
        self.sorted_leaves = self._find_sorted_leaves(self.subtrees)

        # create the ILP problem
        self.ilp = LpProblem("dup placement", LpMinimize)

        dup_vars, r_vars, delta_vars, order_vars, path_vars = self._ilpize()

        #self.dup_approx_constraint(dup_vars)

        # run the ilp!
        self.ilp.solve()


        dup_placement = [node for node in dup_vars if dup_vars[node].varValue == 1.0]
        deltas = [var.varValue for var in delta_vars.values() if var.varValue is not None]
        print('delta vars sum', len(deltas), sum(deltas))
        delta_var_1s = [pair for pair in delta_vars if delta_vars[pair].varValue == 1]
        print('delta vars 1', delta_var_1s)
        r_var_1s = [g for g in r_vars if r_vars[g].varValue == 1]
        print('r vars sum', r_var_1s)
        dup_var_1s = [g for g in dup_vars if dup_vars[g].varValue == 1]
        print('dup var 1s', dup_var_1s)

        def get_path(g1, g2):
            if (g1, g2) in path_vars:
                return path_vars[(g1, g2)].varValue
            elif g1 == g2:
                return 0
            else:
                return path_vars[(g2, g1)].varValue
        # for g1, g2 in delta_var_1s:
        #     print('g2 has dup', g1, g2, dup_vars[g2].varValue, get_path(g1.parent, g2.parent),
        #           'order')


        #for variable in self.ilp.variables():
        #    print variable.name, "=", variable.varValue

        #print "Total Cost (D, L, C): ", value(self.ilp.objective)


        self.cost = value(self.ilp.objective)
        print('the cost is', self.cost)

        #print dup_placement

        # infer locus map
        self.lrecon = self._infer_locus_map(dup_placement)
        #self.log.log("\n\n")

        # postprocessing - given the dups that were parsimonius for DLC, find an ordering
        # that minimizes coal_dup
        ncoal_dup = 0
        self.order = {}
        #if self.coaldupcost == 0.0:
        #    self.order = {}
        #else:
        #    self.order, ncoal_dup = self._infer_opt_order(dup_placement)
        #self.order, ncoal_dup = self._infer_opt_order(dup_placement)

        self.cost += self.coaldupcost * ncoal_dup
        # print('ncoal_dup, ', ncoal_dup, sum([r_vars[node].varValue for node in r_vars]))
        # log gene tree (with species map and locus map)
        #self.log.log("gene tree (with species and locus map)\n")
        #log_tree(self.gtree, self.log, func=draw_tree_recon, srecon=self.srecon, lrecon=self.lrecon)

        #print "TOTAL COST: ", self.cost 

        # revert to use input species tree
        self.stree = substree
        self.srecon = util.mapdict(self.srecon, val=lambda snode: self.stree.nodes[snode.name])
        self.order = util.mapdict(self.order, key=lambda snode: self.stree.nodes[snode.name])

        # convert to LabeledRecon data structure
        #labeled_recon = reconlib.LabeledRecon(self.srecon, self.lrecon, self.order)
        labeled_recon = None
        # print('lrecon', self.lrecon, len(set(self.lrecon.values())))

        recon.draw_tree_recon(self.gtree, self.srecon, self.lrecon, minlen= 10, maxlen = 10)

        #print self.srecon
        #print self.lrecon
        #print self.order

        # calculate runtime
        runtime = self.log.stop()

        # this was the original code. I changed it
        # return self.gtree, labeled_recon, runtime, self.cost
        return self.gtree, self.lrecon, runtime, self.cost

    def _infer_species_map(self):
        """Infer (and assign) species map"""

        self.log.start("Inferring species map")
        self.srecon = phylo.reconcile(self.gtree, self.stree, self.gene2species)
        self.sevents = phylo.label_events(self.gtree, self.srecon)
        self.log.stop()


    def _ilpize(self):
        """Convert the input to ILP parameters"""
        # find top and bottom nodes for each snode
        # key - snode
        # value - list of top (or bottom) nodes for that snode
        top_nodes = collections.defaultdict(list)
        top_nodes_with_child = collections.defaultdict(list)
        bottom_nodes = collections.defaultdict(list)
        for snode, subtrees in self.subtrees.iteritems():
            for subtree in subtrees:
                # the root of each subtree is a top node for this snode
                top_nodes[snode].append(subtree[0])
                if subtree[2]:
                    bottom_nodes[snode].extend(subtree[2])
                if subtree[1]:
                    top_nodes_with_child[snode].append(subtree[0])

        top_pairs = []
        for snode, tops in top_nodes.iteritems():
            for top in tops:
                top_pairs.append((snode, top))

        top_pairs_with_child = []
        for snode, tops in top_nodes_with_child.iteritems():
            for top in tops:
                top_pairs_with_child.append((snode, top))

        # duplication variables - for each node, whether or not there is a dup above it
        # key - gene node
        # value - the corresponding dup LpVariable
        dup_vars = LpVariable.dicts("dup", list(self.gtree.preorder()), 0, 1, LpInteger)

        # loss variables - for each top node, whether or not it is a lost top locus
        # key - (snode, gene node at the top of snode)
        # value - the corresponding loss LpVariable
        loss_vars = LpVariable.dicts("loss", top_pairs, 0, 1, LpInteger)

        # coal variables - for each top node, whether or not it is an extra lineage
        # key - (snode, gene node at the top of snode)
        # value - the corresponding coal LpVariable
        coal_vars = LpVariable.dicts("coal", top_pairs_with_child, 0, 1, LpInteger)

        # loss helper variables - for each top node, whether or not it is the first to be lost
        # key - (snode, gene node at the top of snode)
        # value - the corresponding coal LpVariable
        helper_vars = LpVariable.dicts("helper", top_pairs, 0, 1, LpInteger)

        # path variables - for each pair of nodes - is there a dup along the path between them?
        # key - (gene node 1, gene node 2)
        # value - the corresponding LpVariable
        # NOTE: this is combinations, not permutations, since order doesn't matter
        # that means you might need to check both possible orders for the right key
        all_pairs = list(combination(list(self.gtree.preorder()), 2))
        path_vars = LpVariable.dicts("path", all_pairs, 0, 1, LpInteger)

        gnodes_by_species = collections.defaultdict(list)
        for gnode in self.gtree.preorder():
            gnodes_by_species[self.srecon[gnode]].append(gnode)

        pairs_for_each_species = [list(combination(gnodes, 2)) for gnodes in gnodes_by_species.values()]

        # combine the lists into one big list of pairs
        pairs_in_species = sum(pairs_for_each_species, [])

        pairs_in_species_with_flipped = pairs_in_species + [(g2, g1) for (g1, g2) in pairs_in_species]

        delta_vars_keys = [(g1, g2) for g1, g2 in pairs_in_species_with_flipped
                           if g1.parent is not None and self.srecon[g1] == self.srecon[g1.parent]]
        delta_vars = LpVariable.dicts("coal_dup", delta_vars_keys, 0, 1, LpInteger)


        orders_from_topology = infer_orders_from_topology(self.gtree.root, self.srecon)
        # order_keys has all the pairs in a species that require an order variable
        order_keys = [(g1, g2) for (g1, g2) in pairs_in_species
                      if (g1, g2) not in orders_from_topology and (g2, g1) not in orders_from_topology
                      and not is_leaves_of_same_species(self.srecon, g1, g2)]

        order_vars = LpVariable.dicts("order", order_keys, 0, 1, LpInteger)

        r_vars = LpVariable.dicts("r", list(self.gtree.preorder()), 0, None, LpInteger)

        # objective function - dup cost * number of dups + loss cost * number of losses + coal cost * number of coals
        self.ilp += lpSum([self.dupcost * dup_vars[i] for i in dup_vars])\
                    + lpSum([self.losscost * loss_vars[i] for i in loss_vars])\
                    + lpSum([self.coalcost * coal_vars[i] for i in coal_vars])\
                    + lpSum([self.coaldupcost * r for r in r_vars.values()])


        def get_order(g1, g2):
            if (g1, g2) in order_vars:
                return order_vars[(g1, g2)]
            elif (g2, g1) in order_vars:
                return 1 - order_vars[(g2, g1)]
            elif (g1, g2) in orders_from_topology:
                return orders_from_topology[(g1, g2)]
            else:
                return 1 - orders_from_topology[(g2, g1)]

        def get_path(g1, g2):
            if (g1, g2) in path_vars:
                return path_vars[(g1, g2)]
            elif g1 == g2:
                return 0
            else:
                return path_vars[(g2, g1)]

        # create the path constraints - if there's a dup on a given path, then that path var is 1, otherwise 0
        for genes, path_var in path_vars.iteritems():
            path = self._find_path(genes[0], genes[1])
            self.ilp += lpSum([dup_vars[node] for node in path]) - path_var >= 0
            self.ilp += lpSum([dup_vars[node] for node in path]) - len(path) * path_var <= 0

        # use the smap to create a dictionary
        # key - snode
        # value - list of the gene leaves mapped there
        mapping_dict = collections.defaultdict(list)
        for gleaf in self.gtree.leaves():
            mapping_dict[self.gene2species(gleaf.name)].append(gleaf)

        # create the dup constraints - every path between each gene leaf mapped to the same species leaf must have a dup
        for gleaves in mapping_dict.itervalues():
            for gene1, gene2 in combination(gleaves, 2):
                path = self._find_path(gene1, gene2)
                self.ilp += lpSum([dup_vars[node] for node in path]) >= 1

        # create the loss constraints - 
        for snode, local_top in top_pairs:
            local_loss = loss_vars[(snode, local_top)]
            local_helper = helper_vars[(snode, local_top)]
            
            # find the corresponding path variables to the bottom nodes
            local_paths = []
            local_bottoms = bottom_nodes[snode]
            for local_bottom in local_bottoms:
                # try both orders of keys
                if (local_top, local_bottom) in path_vars:
                    local_paths.append(path_vars[(local_top, local_bottom)])
                else:
                    local_paths.append(path_vars[(local_bottom, local_top)])

            # NOTE: I flipped the inequality sign of this constraint
            self.ilp += lpSum(local_paths) - local_loss - local_helper <= len(local_bottoms) - 1

        # create the helper constraints -
        for snode, local_top in top_pairs:
            local_helper = helper_vars[(snode, local_top)]

            # get all nodes lexicographically "before" this one
            prev_nodes = []
            for node in top_nodes[snode]:
                if node == local_top:
                    break
                else:
                    prev_nodes.append(node)

            # get all previous path vars
            prev_paths = []
            for node in prev_nodes:
                if (local_top, node) in path_vars:
                    prev_paths.append(path_vars[(local_top, node)])
                else:
                    prev_paths.append(path_vars[(node, local_top)])

            self.ilp += lpSum(prev_paths) + local_helper <= len(prev_nodes)
            self.ilp += lpSum(prev_paths) + len(top_nodes[snode]) * local_helper >= len(prev_nodes)

        # create the coal constraints -
        for snode, local_top in top_pairs_with_child:
            local_coal = coal_vars[(snode, local_top)]

            # get all nodes that have a child (might contribute to a coal) lexicographically "before" this one
            prev_nodes = []
            for node in top_nodes_with_child[snode]:
                if node == local_top:
                    break
                else:
                    prev_nodes.append(node)

            # get all previous path vars
            prev_paths = []
            for node in prev_nodes:
                if (local_top, node) in path_vars:
                    prev_paths.append(path_vars[(local_top, node)])
                else:
                    prev_paths.append(path_vars[(node, local_top)])

            self.ilp += lpSum(prev_paths) + local_coal <= len(prev_nodes)
            self.ilp += lpSum(prev_paths) + len(top_nodes[snode]) * local_coal >= len(prev_nodes)

        for g1, g2 in delta_vars_keys:
            if g1.parent is not None and g2.parent is not None and \
                            g1.parent != g2:

                # orderg1parent_g2 should be 1 if g2 is more recent than g2.parent
                if self.srecon[g2] in self.srecon[g1.parent].children:
                    orderg1parent_g2 = 1
                else:
                    orderg1parent_g2 = get_order(g1.parent, g2)

                # orderg2_g1 should be 1 if g1 is more recent than g2 or g1 and g2 must be
                # at approximately the same time because they are leaves of the same species
                if is_leaves_of_same_species(self.srecon, g1, g2):
                    # they are at approximately the same time so order constraint is satisfied
                    orderg2_g1 = 1
                else:
                    orderg2_g1 = get_order(g2, g1)

                self.ilp += delta_vars[(g1, g2)] >= orderg1parent_g2 + orderg2_g1 + \
                                                    (1 - get_path(g1.parent, g2.parent)) + dup_vars[g2] - 3

        for gnodes in list(combination(list(self.gtree.preorder()), 3)):
            g1, g2, g3 = gnodes
            if self.srecon[g1] == self.srecon[g2] == self.srecon[g3]:
                leaves = [g for g in gnodes if len(g.children) == 0 or self.srecon[g] != self.srecon[g.children[0]]]

                if len(leaves) < 2:
                    self.ilp += get_order(g1, g3) >= get_order(g1, g2) + get_order(g2, g3) - 1
                elif len(leaves) == 2:
                    non_leaf = [g for g in gnodes if g not in leaves][0]
                    self.ilp += get_order(non_leaf, leaves[0]) == 1
                    self.ilp += get_order(non_leaf, leaves[1]) == 1

                    # if all the nodes are leaves so we don't order them

        for g2 in self.gtree.preorder():
            gnodes_in_species = gnodes_by_species[self.srecon[g2]]
            relevant_delta_vars = [delta_vars[(g1, g2)] for g1 in gnodes_in_species
                                   if g1 != g2 and (g1, g2) in delta_vars]
            self.ilp += r_vars[g2] >= lpSum(relevant_delta_vars) - 1


        return dup_vars, r_vars, delta_vars, order_vars, path_vars

    #TODO: don't need this one-liner...
    def _infer_locus_map(self, dups):
        """Infer (and assign) locus map"""
        return self._generate_locus_map(self.gtree, self.gtree.root, self.gtree.leaves(), dups)

    def _infer_opt_order(self, dups):
        """minimize coal_dup based on dup placement"""

        # key - snode
        # value - list of the dup nodes as an order of those nodes
        order = {}
        coal_dups = 0

        for snode in self.stree.preorder():
            # restrict to local lrecon
            # local nodes includes all nodes mapped to this snode, as well as the top loci of this snode
            local_nodes = []
            for root, rootchild, leaves in self.subtrees[snode]:
                if rootchild:
                    for node in self.gtree.preorder(root, is_leaf=lambda x: x in leaves):
                        local_nodes.append(node)
                else:
                    local_nodes.append(root)
                    
            local_lrecon = util.subdict(self.lrecon, [node for node in self.gtree if node in local_nodes])
            # find the dups that occurred in this snode - these are the ones we must order
            local_dups = [node for node in dups if self.srecon[node] == snode] 

            nodefunc = lambda node: node
            start, min_orders, nsoln = self._find_min_coal_dup(local_lrecon, self.subtrees[snode], nodefunc=nodefunc,
                                                        dup_nodes=local_dups, all_leaves=self.sorted_leaves[snode])
            min_order = {}
            for locus, orderings in min_orders.iteritems():
                min_order[locus] = _random_choice(orderings)
            ncoal_dup = self._count_coal_dup(local_lrecon, min_order, start, nodefunc=nodefunc)
            order[snode] = min_order
 
            #print order, ncoal_dup
            coal_dups += ncoal_dup

        # print("TOTAL COAL_DUPS: ", coal_dups)

        return order

    #=============================
    # utilities

    def _infer_min_dup_loss(self):
        """Infer cost if only duplications and loss are allowed"""
        recon = phylo.reconcile(self.gtree, self.stree, self.gene2species)
        events = phylo.label_events(self.gtree, recon)

        ndup = phylo.count_dup(self.gtree, events)
        nloss = phylo.count_loss(self.gtree, self.stree, recon)

        return ndup, nloss

    def _infer_trivial_locus_map(self):
        """For debugging only"""
        self.lrecon = dict([(node, 1) for node in self.gtree])
        self.order = dict()

    def _find_all_leaves(self, subtrees):
        all_leaves = []
        for (root, rootchild, leaves) in subtrees:
            if leaves is not None:
                all_leaves.extend(leaves)
        return all_leaves

    def _find_path(self, node1, node2):
        """Returns a list of nodes, the nodes on the path from node1 to node2.
           Does not include the LCA - this is meant to include nodes which might
           have a possible duplication along that path. Since the LCA's incident
           edge is not part of the path, it is not included."""
        path = []
        lca = treelib.lca([node1, node2])

        # walk up from node1 to the LCA, adding those nodes to path
        nodewalk = node1
        while nodewalk != lca:
            path.append(nodewalk)
            nodewalk = nodewalk.parent
        
        # walk up from node2
        nodewalk = node2
        while nodewalk != lca:
            path.append(nodewalk)
            nodewalk = nodewalk.parent

        return path

    def _generate_locus_map(self, tree, root, leaves, dups):
        """returns a dictionary which maps the nodes of tree to loci based on the dup placement"""
        # dup_placement is a list of nodes whose parents are mapped to different loci
        # key - node
        # value - locus (starts at 0)
        locus = 1
        lrecon = {}

        for node in tree.preorder(root, is_leaf = lambda x: x in leaves):
            # if it's the root, it has the first locus
            if node == root:
                lrecon[node] = locus
            # if it has a dup, it has a different locus
            elif node in dups:
                locus += 1
                lrecon[node] = locus
            # if there's no dup, it's the same as the parent
            else:
                lrecon[node] = lrecon[node.parent]

        return lrecon

    def _find_sorted_leaves(self, subtrees):
        # get node ids (for sorting)
        ids = dict((node.name, gid) for (gid, node) in enumerate(self.gtree.preorder()))
        sorted_leaves = {}
        for snode in self.stree.preorder():
            subtrees_snode = subtrees[snode]

            if len(subtrees_snode) == 0:
                # handle root separately
                if snode is self.stree.root:
                    sorted_leaves[snode] = [self.gtree.root]
                continue

            subtrees_snode.sort(key=lambda (root, rootchild, leaves): ids[root.name])
            leaves_snode = []
            for (root, rootchild, leaves) in subtrees_snode:
                if leaves is not None:
                    leaves_snode.extend(leaves)
            leaves_snode.sort(key=lambda node: ids[node.name])
            sorted_leaves[snode] = leaves_snode
        return sorted_leaves

    # THIS CONSTRAINT IS WRONG :(
    def dup_approx_constraint(self, dup_vars):
        """uses the dup-loss model to compute an upper bound on the number of dups"""
        recon = phylo.reconcile(self.gtree, self.stree, self.gene2species)
        events = phylo.label_events(self.gtree, recon)
        max_dups = phylo.count_dup(self.gtree, events)
        self.ilp += lpSum(dup_vars) <= max_dups

a = None

def infer_orders_from_topology(root, srecon):
    order = {}

    # infers all the orders that only involve nodes in this subtree
    def infer_orders(sub_tree_root):
        for descendant in descendants_in_species(sub_tree_root, srecon):
            order[(sub_tree_root, descendant)] = 1

        for child in sub_tree_root.children:
            infer_orders(child)

    infer_orders(root)

    return order


def descendants_in_species(gnode, srecon):
    children_in_species = filter(lambda child: srecon[gnode] == srecon[child], gnode.children)

    return children_in_species + \
           sum([descendants_in_species(child, srecon) for child in children_in_species], [])


def is_leaves_of_same_species(srecon, g1, g2):
    return srecon[g1] == srecon[g2] and \
           (len(g1.children) == 0 or srecon[g1] != srecon[g1.children[0]]) and \
           (len(g2.children) == 0 or srecon[g2] != srecon[g2.children[0]])
