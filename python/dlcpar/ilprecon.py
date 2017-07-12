"""

   Code for the DLC Parsimony Reconciliation
   (duplications, losses, and coalescence)

"""

# python libraries
import sys
import random
import collections
import StringIO
import itertools

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
        self.coaldupcost = coaldupcost if coaldupcost else coalcost

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
        #self.log.log("gene tree (with species map)\n")
        #log_tree(self.gtree, self.log, func=draw_tree_srecon, srecon=self.srecon)

        #treelib.draw_tree(self.gtree)

        # find subtrees
        self.subtrees = reconlib.factor_tree(self.gtree, self.stree, self.srecon, self.sevents)
        # compute the leaves of those subtrees
        self.sorted_leaves = self._find_sorted_leaves(self.subtrees)

        # create the ILP problem
        self.ilp = LpProblem("dup placement", LpMinimize)

        dup_vars = self._ilpize()

        # run the ilp!
        self.ilp.solve()

        dup_placement = [node for node in dup_vars if dup_vars[node].varValue == 1.0]

        for variable in self.ilp.variables():
            print variable.name, "=", variable.varValue

        print "Total Cost (D, L, C): ", value(self.ilp.objective)

        self.cost = value(self.ilp.objective)

        #print dup_placement

        # infer locus map
        self.lrecon = self._infer_locus_map(dup_placement)
        #self.log.log("\n\n")

        # postprocessing - given the dups that were parsimonius for DLC, find an ordering
        # that minimizes coal_dup
        self.order, ncoal_dup = self._infer_opt_order(dup_placement)

        self.cost += self.coaldupcost * ncoal_dup

        # log gene tree (with species map and locus map)
        #self.log.log("gene tree (with species and locus map)\n")
        #log_tree(self.gtree, self.log, func=draw_tree_recon, srecon=self.srecon, lrecon=self.lrecon)

        print "TOTAL COST: ", self.cost 

        # revert to use input species tree
        self.stree = substree
        self.srecon = util.mapdict(self.srecon, val=lambda snode: self.stree.nodes[snode.name])
        self.order = util.mapdict(self.order, key=lambda snode: self.stree.nodes[snode.name])

        # convert to LabeledRecon data structure
        labeled_recon = reconlib.LabeledRecon(self.srecon, self.lrecon, self.order)

        # calculate runtime
        runtime = self.log.stop()

        return self.gtree, labeled_recon, runtime, self.cost

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

        # objective function - dup cost * number of dups + loss cost * number of losses + coal cost * number of coals
        self.ilp += lpSum([self.dupcost * dup_vars[i] for i in dup_vars]) + lpSum([self.losscost * loss_vars[i] for i in loss_vars]) + lpSum([self.coalcost * coal_vars[i] for i in coal_vars])

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

        return dup_vars

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
            # find the dups that occurred in this snode - these are the ones we must order
            local_dups = [node for node in dups if self.srecon[node] == snode] 

            nodefunc = lambda node: node
            start, min_orders, nsoln = self._find_min_coal_dup(self.lrecon, self.subtrees[snode], nodefunc=nodefunc,
                                                        dup_nodes=local_dups, all_leaves=self.sorted_leaves[snode])
            min_order = {}
            for locus, orderings in min_orders.iteritems():
                min_order[locus] = _random_choice(orderings)
            ncoal_dup = self._count_coal_dup(self.lrecon, min_order, start, nodefunc=nodefunc)
            order[snode] = min_order
 
            #print order, ncoal_dup
            coal_dups += ncoal_dup

        #print "TOTAL COAL_DUPS: ", coal_dups

        return order, coal_dups

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
        locus = 0
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

