#!/usr/bin/env python

"""

   Code for the DLC Parsimony Reconciliation
   (duplications, losses, and coalescence)

"""

import sys
import random
import collections
import StringIO
import math

from rasmus import treelib, util
from compbio import phylo

from yjw import combinatorics

from dlcpar import common
from dlcpar import reconlib

# gtree = G (with implied speciation nodes)
# stree = S
# srecon = M
# lrecon = L

#==========================================================
# utilities

import operator
def prod(iterable):
    return reduce(operator.mul, iterable, 1)

#==========================================================

def dlc_recon(tree, stree, gene2species,
              dupcost=1, losscost=1, coalcost=1,
              implied=True, delay=True,
              prescreen=False, prescreen_min=5, prescreen_factor=2,
              max_loci=util.INF, max_dups=util.INF, max_losses=util.INF,
              log=sys.stdout):
    """Perform reconciliation using DLCoal model with parsimony costs"""

    reconer = DLCRecon(tree, stree, gene2species,
                       dupcost=dupcost, losscost=losscost, coalcost=coalcost,
                       implied=implied, delay=delay,
                       prescreen=prescreen, prescreen_min=prescreen_min, prescreen_factor=prescreen_factor,
                       max_loci=max_loci, max_dups=max_dups, max_losses=max_losses,
                       log=log)
    return reconer.recon()

    

class DLCRecon(object):

    def __init__(self, gtree, stree, gene2species,
                 dupcost=1, losscost=1, coalcost=1,
                 implied=True, delay=True,
                 prescreen=False, prescreen_min=5, prescreen_factor=2,
                 max_loci=util.INF, max_dups=util.INF, max_losses=util.INF,
                 name_internal="n", log=sys.stdout):

        # rename gene tree nodes
        common.rename_nodes(gtree, name_internal)
        
        self.gtree = gtree
        self.stree = stree
        self.gene2species = gene2species
        
        assert (dupcost >= 0) and (losscost >= 0) and (coalcost >= 0)
        self.dupcost = dupcost
        self.losscost = losscost
        self.coalcost = coalcost
        
        self.implied = implied
        self.delay = delay
        assert (prescreen_min > 0) and (prescreen_factor > 0)
        self.prescreen = prescreen
        self.prescreen_min = prescreen_min
        self.prescreen_factor = prescreen_factor
        assert (max_loci > 0) and (max_dups > 0) and (max_losses > 0)
        self.max_loci = max_loci
        self.max_dups = max_dups
        self.max_losses = max_losses


        self.name_internal = name_internal
        self.log = util.Timer(log)

    #=============================
    # main methods

    def recon(self):
        """Perform reconciliation"""

        self.log.start("Reconciling")

        # log input gene and species trees
        self.log.log("gene tree\n")
        log_tree(self.gtree, self.log, func=treelib.draw_tree_names)
        self.log.log("species tree\n")
        log_tree(self.stree, self.log, func=treelib.draw_tree_names)
               
        # infer species map
        self._infer_species_map()
        self.log.log("\n\n")

        # add implied speciation nodes but first start the species tree at the right root
        substree = treelib.subtree(self.stree, self.srecon[self.gtree.root])
        subrecon = util.mapdict(self.srecon, val=lambda snode: substree.nodes[snode.name])
        
        # switch internal storage with subtrees
        self.stree, subtree = substree, self.stree
        self.srecon, subrecon = subrecon, self.srecon

        # add implied nodes (standard speciation, speciation from duplication, delay nodes)
        # then relabel events (so that factor_tree works)
        reconlib.add_implied_nodes(self.gtree, self.stree, self.srecon, self.sevents, delay=self.delay)
        self.sevents = phylo.label_events(self.gtree, self.srecon)
        common.rename_nodes(self.gtree, self.name_internal)

        # log gene tree (with species map)
        self.log.log("gene tree (with species map)\n")
        log_tree(self.gtree, self.log, func=draw_tree_srecon, srecon=self.srecon)

        # infer locus map
        self._infer_locus_map()
        self.log.log("\n\n")        

        # log gene tree (with species map and locus map)
        self.log.log("gene tree (with species and locus map)\n")
        log_tree(self.gtree, self.log, func=draw_tree_recon, srecon=self.srecon, lrecon=self.lrecon)

        # revert to use input species tree
        self.stree = subtree
        self.srecon = util.mapdict(self.srecon, val=lambda snode: self.stree.nodes[snode.name])
        self.order = util.mapdict(self.order, key=lambda snode: self.stree.nodes[snode.name])

        # convert to LabeledRecon data structure
        labeled_recon = reconlib.LabeledRecon(self.srecon, self.lrecon, self.order)

        self.log.stop()
        
        return self.gtree, labeled_recon


    def _infer_species_map(self):
        """Infer (and assign) species map"""
        
        self.log.start("Inferring species map")
        self.srecon = phylo.reconcile(self.gtree, self.stree, self.gene2species)
        self.sevents = phylo.label_events(self.gtree, self.srecon)
        self.log.stop()


    def _infer_locus_map(self):
        """Infer (and assign) locus map"""

        gtree = self.gtree
        stree = self.stree

        # get node ids (for sorting)
        ids = dict((node.name, gid) for (gid, node) in enumerate(gtree.preorder()))

        # find speciation subtrees and sort by id of root
        # also find sorted leaves
        subtrees = reconlib.factor_tree(gtree, stree, self.srecon, self.sevents)
        sorted_leaves = {}
        for snode in self.stree.preorder(stree.root):
            subtrees_snode = subtrees[snode]

            if len(subtrees_snode) == 0:
                # handle root separately
                sorted_leaves[snode] = [gtree.root]
                continue
            
            subtrees_snode.sort(key=lambda (root, rootchild, leaves): ids[root.name])
            leaves_snode = []
            for (root, rootchild, leaves) in subtrees_snode:
                if leaves is not None:
                    leaves_snode.extend(leaves)
            leaves_snode.sort(key=lambda node: ids[node.name])
            sorted_leaves[snode] = leaves_snode
                    
        # enumerate valid locus maps, then infer optimal
        self.log.start("Inferring locus map")
        locus_maps = self._enumerate_locus_maps(subtrees, sorted_leaves)
        self.log.log("\n\n")
        self._infer_opt_locus_map(locus_maps, subtrees)
##        self._infer_trivial_locus_map()
        self.log.stop()
        

    #=============================
    # event/cost methods -- these operate at the species branch level

    def _find_all_leaves(self, subtrees):
        all_leaves = []
        for (root, rootchild, leaves) in subtrees:
            if leaves is not None:
                all_leaves.extend(leaves)
        return all_leaves


    def _count_events(self, lrecon, subtrees, nodefunc=lambda node: node.name,
                      min_cost=util.INF, return_cost=True,
                      all_leaves=None,
		      max_dups=util.INF, max_losses=util.INF):
        """
        Count number of dup, loss, coal events
        
        Returns
        - number of duplications
        - number of losses
        - minimum number of extra lineages over all internal node orderings
        - the optimal internal node ordering
        """

        extra = {"species_map" : self.srecon, "locus_map" : lrecon}

        # defaults
        ndup, nloss, ncoal_spec, ncoal_dup, order, cost = util.INF, util.INF, util.INF, util.INF, {}, util.INF

        # duplications
##        ndup = reconlib.count_dup_snode(self.gtree, self.stree, extra, snode=None,
##                                        subtrees_snode=subtrees,
##                                        nodefunc=nodefunc)
        dup_nodes = reconlib.find_dup_snode(self.gtree, self.stree, extra, snode=None,
                                            subtrees_snode=subtrees,
                                            nodefunc=nodefunc)
        ndup = len(dup_nodes)
        # check if dups exceeds mincost (allow equality to select from random)
        if (ndup > max_dups) or (self._compute_cost(ndup, 0, 0) > min_cost):
            if return_cost:
                return ndup, nloss, ncoal_spec, ncoal_dup, order, cost
            else:
                return ndup, nloss, ncoal_spec, ncoal_dup, order

        # losses
        nloss = reconlib.count_loss_snode(self.gtree, self.stree, extra, snode=None,
                                          subtrees_snode=subtrees,
                                          nodefunc=nodefunc)
        # check if dups + losses exceeds mincost (allow equality to select from random)
        if (nloss > max_losses) or (self._compute_cost(ndup, nloss, 0) > min_cost):
            if return_cost:
                return ndup, nloss, ncoal_spec, ncoal_dup, order, cost
            else:
                return ndup, nloss, ncoal_spec, ncoal_dup, order
        
        # extra lineages at speciations
        ncoal_spec = reconlib.count_coal_snode_spec(self.gtree, self.stree, extra, snode=None,
                                                    subtrees_snode=subtrees,
                                                    nodefunc=nodefunc,
                                                    implied=self.implied)
        # check if dups + losses + coal (spec) exceeds mincost (allow equality to select from random)
        if self._compute_cost(ndup, nloss, ncoal_spec) > min_cost:
            if return_cost:
                return ndup, nloss, ncoal_spec, ncoal_dup, order, cost
            else:
                return ndup, nloss, ncoal_spec, ncoal_dup, order
        
        # extra lineages at duplications
        ncoal_dup, order = self._count_min_coal_dup(lrecon, subtrees, nodefunc=nodefunc, dup_nodes=dup_nodes)

        if return_cost:
            cost = self._compute_cost(ndup, nloss, ncoal_spec + ncoal_dup)
            return ndup, nloss, ncoal_spec, ncoal_dup, order, cost
        else:
            return ndup, nloss, ncoal_spec, ncoal_dup, order


    def _count_min_coal_dup_nosearch(self, lrecon, subtrees, nodefunc=lambda node: node.name,
                                     dup_nodes=None, all_leaves=None):
        """
        Over all internal node orderings, find minimum number of inferred extra lineages due to duplications.
        Finds the optimum node ordering WITHOUT enumerating over all internal node orderings.
        """
        
        extra = {"species_map" : self.srecon, "locus_map" : lrecon}

        # extra lineages at duplications
        start = self._find_locus_order_start(lrecon, subtrees, nodefunc=nodefunc,
                                             dup_nodes=dup_nodes, all_leaves=all_leaves)
        min_order = self._find_locus_order(lrecon, subtrees, start, nodefunc=nodefunc,
                                           dup_nodes=dup_nodes, all_leaves=all_leaves)
        min_ncoal_dup = self._count_coal_dup(lrecon, min_order, start, nodefunc=nodefunc)
               
        return min_ncoal_dup, min_order
    

    def _count_min_coal_dup_search(self, lrecon, subtrees, nodefunc=lambda node: node.name,
                                   dup_nodes=None, all_leaves=None):
        """
        Over all internal node orderings, find minimum number of inferred extra lineages due to duplications.
        Finds the optimum node ordering by enumerating over all internal node orderings.
        """
        
        extra = {"species_map" : self.srecon, "locus_map" : lrecon}

        # extra lineages at duplications - try all orderings of internal nodes
        min_ncoal_dup = util.INF
        min_order_lst = None
        start = self._find_locus_order_start(lrecon, subtrees, nodefunc=nodefunc,
                                             dup_nodes=dup_nodes, all_leaves=all_leaves)
        for order in self._enumerate_locus_order(lrecon, subtrees, start, nodefunc=nodefunc,
                                                 all_leaves=all_leaves):
            ncoal_dup = self._count_coal_dup(lrecon, order, start, nodefunc=nodefunc)

            if ncoal_dup < min_ncoal_dup:
                min_ncoal_dup = ncoal_dup
                min_order_lst = [order]
            elif ncoal_dup == min_ncoal_dup:
                min_order_lst.append(order)
                
        # handle 1) no orderings (== no duplications) and 2) multiple optimal orderings (select one randomly)
        if min_order_lst is None:
            min_ncoal_dup = 0
            min_order = None
        else:
            min_order = random.choice(min_order_lst)
        
        return min_ncoal_dup, min_order


    _count_min_coal_dup = _count_min_coal_dup_nosearch


    def _find_locus_order_start(self, lrecon, subtrees, nodefunc=lambda node: node.name,
                                dup_nodes=None, all_leaves=None):
        """Helper function for _count_coal"""

        extra = {"species_map" : self.srecon, "locus_map" : lrecon}

        if dup_nodes is None:
            dup_nodes = reconlib.find_dup_snode(self.gtree, self.stree, extra, snode=None,
                                                subtrees=subtrees, nodefunc=nodefunc)

        if all_leaves is None:
            all_leaves = self._find_all_leaves(subtrees)
        
        # rest of code is adapted from first routine in reconlib.count_coal_snode_dup        
        start = collections.defaultdict(list)
        parent_loci = set()
        for node in dup_nodes:
            if node.parent:
                parent_loci.add(lrecon[nodefunc(node.parent)])
        # for each locus found, if this locus is a "parent locus",
        # add the children if the dup node is not a leaf
        # (leaves never incur extra lineages in this species branch)
        for node in dup_nodes:
            locus = lrecon[nodefunc(node)]
            if (node in all_leaves) or (locus not in parent_loci):
                continue
            for child in node.children:
                start[locus].append(child)
        # for each locus that exists at the top of the species branch,
        # if this locus is a "parent locus",
        # add the child if it exists (assume immediate loss if child does not exist)
        for (root, rootchild, leaves) in subtrees:
            locus = lrecon[nodefunc(root)]
            if locus not in parent_loci:
                continue

            # handle root separately
            if not root.parent:
                for child in root.children:
                    if nodefunc(child) in lrecon:   # DIFFERENT from reconlib: ensure node in sbranch
                        start[locus].append(child)
            else:
                if rootchild:
                    start[locus].append(rootchild)

        return start


    def _count_coal_dup(self, lrecon, order, start, nodefunc=lambda node: node.name):
        """Helper function for _count_coal"""

        assert set(start) == set(order), (dict(start), order)
        
        # code is adapted from second routine in reconlib.count_coal_snode_dup
        ncoal = 0
        count = {}
        for plocus, nodes in order.iteritems():
            current = start[plocus][:]   # DIFFERENT from reconlib: use copy!!
            num_lineages = len(current)
            for next in nodes:
                assert num_lineages == len(current), (num_lineages, nodes)
                assert next in current, (next, current)

                # locus of next node
                next_locus = lrecon[nodefunc(next)]

                # keep if leaf and locus does not change : leaves (extant genes) exist to present time
                if (next.is_leaf()) and (plocus == next_locus):
                    pass
                else:              
                    current.remove(next)
                    num_lineages -= 1

                # update lineage count and list of nodes
                if plocus == next_locus:
                    # deep coalescence - keep even if next in leaves to allow for delay btwn coalescence and speciation
                    for child in next.children:
                        current.append(child)
                        num_lineages += 1
                else:
                    # duplication
                    if num_lineages > 1:
                        ncoal += num_lineages - 1
        return ncoal


    def _find_locus_order(self, lrecon, subtrees, start=None, nodefunc=lambda node: node.name,
                          dup_nodes=None, all_leaves=None):
        """
        Find internal node ordering that minimizes the number of inferred extra lineages due to duplications.

        TODO: merge with _count_coal_dup
        """
        
        gtree = self.gtree
        extra = {"species_map" : self.srecon, "locus_map" : lrecon}

        if dup_nodes is None:
            dup_nodes = reconlib.find_dup_snode(self.gtree, self.stree, extra, snode=None,
                                                subtrees=subtrees, nodefunc=nodefunc)
            
        if all_leaves is None:
            all_leaves = self._find_all_leaves(subtrees) 

        def get_local_order(start, locus):           
            # find nodes with parent locus = locus
            # for each node, also find path from "start" node (non-inclusive of end-points)
            paths = {}
            for node in start:
                # recur over subtree
                for node2 in gtree.preorder(node, is_leaf=lambda x: x in all_leaves):
                    if lrecon[nodefunc(node2.parent)] == locus:
                        if node2 is node:
                            paths[node2] = collections.deque()
                        else:
                            paths[node2] = collections.deque(paths[node2.parent]); paths[node2].append(node2.parent)

            # keep track of all nodes, current nodes (equal to start), and duplicated nodes that have this parent locus
            all_nodes = set(paths.keys())
            current = start[:]
            dup = filter(lambda node: lrecon[nodefunc(node.parent)] == locus, dup_nodes)

            # retain non-trivial paths
            for node, nodes in paths.items():
                if len(nodes) == 0:
                    del paths[node]

            # get local order
            # 1) pick duplicated node with shortest path
            #    a) if this node can be chosen next, choose it
            #    b) otherwise, choose all nodes in the path
            # 2) update paths
            # 3) recur until no duplicated node is left
            # 4) choose rest of nodes (at this point, number of children no longer matters)
            local_order = []
            while len(dup) > 0:
                # find next duplicated node
                items = filter(lambda node: node in current, dup)
                if len(items) == 0:
                    items, min_dist = util.minall(dup, minfunc=lambda node: len(paths[node]))
                next_dup = random.choice(items)
                dup.remove(next_dup)

                # add nodes in path first
                if next_dup in paths:
                    path_nodes = list(paths[next_dup])
                    del paths[next_dup]
                    for next in path_nodes:
                        assert next in current, (next, current)
                        assert lrecon[nodefunc(next)] == locus, (next.name, lrecon[nodefunc(next)], locus)
                        local_order.append(next)
                        current.remove(next)
                        current.extend(next.children)

                        # update paths
                        for node, nodes in paths.items():
                            if nodes[0] is next:
                                nodes.popleft()
                            assert next not in nodes
                            if len(nodes) == 0:
                                del paths[node]
                            else:
                                paths[node] = nodes

                # now add the duplicated node
                next = next_dup
                assert next in current, (next, current)
                local_order.append(next)
                current.remove(next)
            all_nodes.difference_update(local_order)
            while len(all_nodes) > 0:
                tochoose = filter(lambda node: node in current, all_nodes)
                next = random.choice(tochoose)
                all_nodes.remove(next)
                
                assert next in current, (next, current)
                local_order.append(next)
                current.remove(next)
                current.extend(next.children)
            return local_order

        order = {}
        for locus, nodes in start.iteritems():
            order[locus] = get_local_order(nodes, locus)

        return order


    def _enumerate_locus_order(self, lrecon, subtrees, start=None, nodefunc=lambda node: node.name,
                               all_leaves=None):
        """
        Find all internal node orderings.
        
        TODO: merge with _count_coal_dup
        """

        extra = {"species_map" : self.srecon, "locus_map" : lrecon}

        if all_leaves is None:
            all_leaves = self._find_all_leaves(subtrees)
            
        if start is None:
            start = self._find_locus_order_start(lrecon, subtrees, nodefunc=nodefunc,
                                                 all_leaves=all_leaves)

        def helper(lst, tochoose, locus):
            if len(tochoose) == 0:
                yield lst
            else:
                for child in tochoose:
                    next_lst = lst[:]
                    next_lst.append(child)
                    
                    next_tochoose = tochoose[:]
                    next_tochoose.remove(child)
                    if child not in all_leaves and lrecon[nodefunc(child)] == locus:
                        next_tochoose.extend(child.children)
                        
                    for local_order in helper(next_lst, next_tochoose, locus):
                        yield local_order

        order = collections.defaultdict(list)
        for locus, nodes in start.iteritems():
            for local_order in helper([], nodes, locus):
                order[locus].append(local_order)

        choices = [order[locus] for locus in start]
        for choice in combinatorics.product(*choices):
            order = {}
            for ndx, locus in enumerate(start):
                order[locus] = choice[ndx]
            yield order


    def _compute_cost(self, ndup, nloss, ncoal):
        """Find reconciliation cost"""
        return ndup*self.dupcost + nloss*self.losscost + ncoal*self.coalcost



    #=============================
    # locus (and locus map) methods

    def _find_unique_loci(self, lrecon, leaves, start=1):
        """Returns unique leaf loci"""

        # unique mapping
        x = map(lambda node: lrecon[node.name], leaves)
        y = []
        m = {}
        next = start
        for locus in x:
            if locus not in m:
                m[locus] = next
                next += 1
            y.append(m[locus])
        y = tuple(y)
            
        return y, m


    def _evolve_subtree(self, root, leaves, state, start=1, next=None):
        """Given state changes, find the locus at the nodes"""

        gtree = self.gtree

        if next is None:
            next = start + 1        
        assert next > start, (start, next)
        
        lrecon = {}
        for node in gtree.preorder(root, is_leaf=lambda x: x in leaves):
            if node is root:
                lrecon = {node.name : start}
            else:
                if state[node.name]:
                    lrecon[node.name] = next
                    next += 1
                else:
                    lrecon[node.name] = lrecon[node.parent.name]
        return lrecon


    def _find_locus_states_sbranch(self, subtrees, is_leaf,
                                   max_loci=util.INF, max_dups=util.INF):
        """Finds the valid locus states for each subtree in the species branch"""
        
        gtree = self.gtree
        
        # storage for states for all subtrees
        # each node of the tree is labeled with T/F corresponding to whether a dup occurred along the branch
        all_states = []

        # iterate through subtrees in sbranch to find leaf states
        # start at rootchild since duplication along (root, rootchild) will be handled when combining subtrees
        for (root, rootchild, leaves) in subtrees:
            if not rootchild:
                # loss
                all_states.append([None])
                continue

            # find maximum number of dup in subtree
            # K can be less than Ksub due to speciation nodes 
            max_dups_subtree = len(leaves)
            if max_dups < max_dups_subtree:
                max_dups_subtree = max_dups

            # initialize states storage and rootchild node
            # TODO: more efficient to use matrix to store states
            states_down = collections.defaultdict(list)
            state = {rootchild.name : False}
            states_down[rootchild].append(state)

            # recur down subtree to find leaf states 
            for node in gtree.preorder(rootchild, is_leaf=lambda x: x in leaves):
                if node in leaves:
                    continue

                for state in states_down[node]:
                    ndup = util.counteq(True, state.values())
                    children = node.children
                    nchildren = len(children)

                    if nchildren == 0:
                        raise Exception("invalid number of children: %s" % node.name)
                    elif nchildren == 1:
                        child = children[0]
                        s1 = state.copy();  s1[child.name] = False
                        states_down[child].append(s1)
                        if ndup < max_dups_subtree:
                            s2 = state.copy();  s2[child.name] = True
                            states_down[child].append(s2)
                    elif nchildren == 2:
                        left, right = children
                        s1 = state.copy();      s1[left.name] = False;  s1[right.name] = False
                        states_down[left].append(s1)
                        if ndup < max_dups_subtree:
                            # daughter lineage matters
                            s2 = state.copy();  s2[left.name] = False;  s2[right.name] = True
                            states_down[left].append(s2)
                            s3 = state.copy();  s3[left.name] = True;   s3[right.name] = False
                            states_down[left].append(s3)

                            # duplication along both child lineages incurs a deep coalescence so
                            # it is always less parsimonious than dup from parent followed by dup in one child lineage
                            # (dup in parent and both children incurs a loss so again is less parsimonious)
                            ##if ndup + 1 < max_dups_subtree:
                            ##    s4 = h.copy(); s4[left.name] = True; s4[right.name] = True
                            ##    states_down[left].append(s4)
                        states_down[right] = states_down[left]
                    else:
                        raise Exception("invalid number of children: %s" % node.name)

            # recur up subtree to combine leaf states
            # TODO: check if copies are needed
            states_up = collections.defaultdict(list)
            for node in gtree.postorder(rootchild, is_leaf=lambda x: x in leaves):
                # base case (leaf)
                if node in leaves:
                    for s in states_down[node]:
                        states_up[node].append(s)
                    continue

                children = node.children
                nchildren = len(children)
            
                if nchildren == 0:
                    raise Exception("invalid number of children: %s" % node.name)
                elif nchildren == 1:
                    child = children[0]
                    for s in states_up[child]:
                        states_up[node].append(s)
                elif nchildren == 2:
                    left, right = children
                    
                    # another base case (both left and right children are leaves)
                    if (left in leaves) and (right in leaves):
                        assert states_down[left] == states_down[right], (states_down[left], states_down[right])
                        for s in states_down[left]:
                            states_up[node].append(s)
                        continue

                    # combine children
                    for sleft in states_up[left]:
                        sleftset = set(sleft)
                        for sright in states_up[right]:
                            intersect = sleftset.intersection(sright)
                            if all([sleft[name] == sright[name] for name in intersect]):
                                s = sleft.copy()
                                s.update(sright)
                                states_up[node].append(s)
                else:
                    raise Exception("invalid number of children: %s" % node.name)

            # ensure valid states if sbranch is leaf branch
            # also, for all species branches, allow duplication along root branch
            states = []       # state for one subtree

            # iterate over states
            for i, state in enumerate(states_up[rootchild][:]):
                valid = True

                if is_leaf:
                    # find locus at leaves (in this subtree), then check if valid
                    lrecon = self._evolve_subtree(rootchild, leaves, state)                    
                    leaf_loci = [lrecon[node.name] for node in leaves]
                    valid = len(set(leaf_loci)) == len(leaf_loci)
                else:
                    if max_loci != util.INF:
                        # find locus at leaves (in this subtree), then check if valid
                        lrecon = self._evolve_subtree(rootchild, leaves, state)                    
                        leaf_loci = [lrecon[node.name] for node in leaves]
                        valid = len(set(leaf_loci)) <= max_loci
                        
                if valid:
                    # store original (without duplication along root branch)
                    states.append(state)

                    # allow duplication along root branch
                    if root is not rootchild:
                        ndup = util.counteq(True, state.values())
                        if ndup < max_dups_subtree:
                            s1 = state.copy(); s1[rootchild.name] = True
                            states.append(s1)
                        
            # store
            all_states.append(states)
            
        return all_states


    def _enumerate_locus_maps(self, subtrees, sorted_leaves, init_locus=1):
        """
        Enumerates the valid locus maps by recurring down the species tree.
        For each combination of locus assignments at the top and bottom of an sbranch,
        also find the optimal locus map and order.
        (This saves memory compared to storing all valid locus maps.)
        """

        self.log.start("Enumerating locus maps")

        gtree = self.gtree
        stree = self.stree
        gene2species = self.gene2species
        
        # find maximum number of dup - does not assume that species map is MPR
        recon = phylo.reconcile(gtree, stree, gene2species)
        events = phylo.label_events(gtree, recon)
        max_dups = phylo.count_dup(gtree, events)
        max_loci_sbranch = self.max_loci
        max_dups_sbranch = min(max_dups, self.max_dups)
        max_losses_sbranch = self.max_losses
        self.log.log("Max # loci: %s" % str(max_loci_sbranch))
        self.log.log("Max # dup: %d" % max_dups)
        self.log.log("Max # dup per sbranch: %s" % str(max_dups_sbranch))
        self.log.log("Max # loss per sbranch: %s" % str(max_losses_sbranch))
        self.log.log()

        # partitions at each sbranch
        # key1 = snode, key2 = bottom_loci, key3 = top_loci, value = (lrecon, order, cost)
        PS = {}

        if self.prescreen:
            # key1 = snode, key2 = bottom_loci, value = min cost-to-go (from root) required to assign bottom_loci to sbranch
            GS = collections.defaultdict(dict)

	# recur down species tree
        sroot = stree.root  # stree.root == srecon[gtree.root] due to assignment of substree
        for snode in stree.preorder(sroot):
            self.log.start("Working on snode %s" % snode.name)
            is_leaf = snode.is_leaf()

            # get subtrees in this sbranch
            subtrees_snode = subtrees[snode]

            # nothing happened along the branch if there are no subtrees
            # still have to initialize loci if at root though
            if len(subtrees_snode) == 0:
                # handle root separately
                if snode is stree.root:
                    top_loci = bottom_loci = (init_locus,)
                    lrecon = {gtree.root.name: init_locus}
                    order = {}
                    cost = 0
                    PS[snode] = collections.defaultdict(dict)
                    PS[snode][bottom_loci][top_loci] = (lrecon, order, cost)
                    if self.prescreen:
                        GS[snode][bottom_loci] = 0
                self.log.log("Empty sbranch")
                self.log.stop()
                continue
            
            # initialize storage for this sbranch
            PS[snode] = collections.defaultdict(dict)
            FS = collections.defaultdict(dict)  # key = (bottom_loci, top_loci), val = (mincost, mindup)
            
            # hash subtrees using root
            subtrees_hash = {}
            for (root, rootchild, leaves) in subtrees_snode:
                subtrees_hash[root] = (rootchild, leaves)

            # get leaves for this sbranch
            leaves_snode = sorted_leaves[snode]

            # get states (changed/unchanged) for each branch of each subtree in sbranch
            states = self._find_locus_states_sbranch(subtrees_snode, is_leaf,
                                                     max_loci=max_loci_sbranch,
                                                     max_dups=util.INF if is_leaf else max_dups_sbranch)

            # top of this sbranch is the bottom of the parent sbranch
            if snode.parent in PS:
                top_loci_lst = PS[snode.parent]
                top_leaves = sorted_leaves[snode.parent]
            else:
                top_loci_lst = {(init_locus,): None}
                top_leaves = [gtree.root]

            # TODO: incorporate max_dup from sbranch to sbranch
            self.log.log("top nodes: %s" % ','.join(map(lambda node: node.name, top_leaves)))
            self.log.log("bottom nodes: %s" % ','.join(map(lambda node: node.name, leaves_snode)))
##            self.log.log("top_loci: %s" % ';'.join(map(str, top_loci_lst.keys())))
            self.log.log("number of assignments at top nodes: %d" % len(top_loci_lst))
            states_len = map(len, states)
            self.log.log("number of states: %s" % ','.join(map(str, states_len)))
            self.log.log("number of state combinations: %d" % prod(states_len))
            self.log.log("")

            # combine subtrees
            for ndx1, (top_loci, _) in enumerate(top_loci_lst.iteritems()):
                # start with loci at top of sbranch
                assert len(top_loci) == len(states), (len(top_loci), len(states))

                # initialize next loci with top of sbranch
                init_next = max(top_loci) + 1
                init_lrecon = {}
                for i, start in enumerate(top_loci):
                    init_lrecon[top_leaves[i].name] = start
                
                for ndx2, state in enumerate(combinatorics.product(*states)):
                    next = init_next
                    lrecon = init_lrecon.copy()

                    #================================================================================
                    # find next lrecon using states
                    for i, start in enumerate(top_loci):
                        s = state[i]
                        if s is not None:
                            root = top_leaves[i]
                            rootchild, leaves = subtrees_hash[root]

                            # evolve from root to rootchild
                            if s[rootchild.name]:
                                rootchild_loci = next
                                next += 1
                            else:
                                rootchild_loci = start

                            # evolve from rootchild to leaves
                            l = self._evolve_subtree(rootchild, leaves, state=s, start=rootchild_loci, next=next)
                            assert all([name not in lrecon for name in l if name != root.name]), (l, lrecon)
                            lrecon.update(l)
                            next = max(init_next, max(lrecon.values())) + 1

                    #================================================================================
                    # (unique) loci at bottom of sbranch
                    bottom_loci, mapping = self._find_unique_loci(lrecon, leaves_snode)

                    #================================================================================
                    # check validity
                    leaf_loci = [lrecon[node.name] for node in leaves_snode]
                    if is_leaf:
                        # if leaf, see if all leaves belong to distinct loci
                        # checks for validity across all leaves, not just leaves in a subtree
                        if len(set(leaf_loci)) != len(leaf_loci):
##                            self.log.log("\t%s -> %s : [%d, %d] (skipped - leaf loci)" % (str(top_loci), str(bottom_loci), ndx1, ndx2))
##                            self.log.log("\t\tlrecon: %s" % str(lrecon))
##                            self.log.log()
                            continue
                    else:
                        # skip if does not satisfy ancestral lineage count
                        if len(set(leaf_loci)) > max_loci_sbranch:
                            self.log.log("\t%s -> %s : [%d, %d] (skipped - locus count)" % (str(top_loci), str(bottom_loci), ndx1, ndx2))
                            self.log.log("\t\tlrecon: %s" % str(lrecon))
                            self.log.log()
                            continue

                    #================================================================================                        
                    # find optimal cost (and order) for this lrecon
                    if (bottom_loci, top_loci) in FS:
                        mincost, mindup = FS[(bottom_loci,top_loci)]
                    else:
                        mincost, mindup = util.INF, util.INF
                    ndup, nloss, ncoal_spec, ncoal_dup, order, cost = self._count_events(lrecon, subtrees_snode,
                                                                                         min_cost=mincost,
                                                                                         return_cost=True,
                                                                                         all_leaves=leaves_snode,
											 max_dups=util.INF if is_leaf else max_dups_sbranch,
											 max_losses=util.INF if is_leaf else max_losses_sbranch)
                    ncoal = ncoal_spec + ncoal_dup

                    # skip if exceeds max # of dup/loss for (ancestral) sbranch
                    if not is_leaf:
                        if ndup > max_dups_sbranch:
                            self.log.log("\t%s -> %s : [%d, %d] (skipped - dup)" % (str(top_loci), str(bottom_loci), ndx1, ndx2))
                            self.log.log("\t\tlrecon: %s" % str(lrecon))
                            self.log.log()
                            continue
                        if nloss > max_losses_sbranch:
                            self.log.log("\t%s -> %s : [%d, %d] (skipped - loss)" % (str(top_loci), str(bottom_loci), ndx1, ndx2))
                            self.log.log("\t\tlrecon: %s" % str(lrecon))
                            self.log.log()
                            continue

                    # log
                    self.log.log("\t%s -> %s : [%d, %d]" % (str(top_loci), str(bottom_loci), ndx1, ndx2))
                    self.log.log("\t\tlrecon: %s" % str(lrecon))
                    self.log.log("\t\torder: %s" % str(order))
                    self.log.log("\t\tmapping: %s" % str(mapping))
                    self.log.log("\t\tndup: %s, nloss: %s, ncoal: %s (spec: %s, dup: %s), cost: %g" % \
                                 (str(ndup), str(nloss), str(ncoal), str(ncoal_spec), str(ncoal_dup), cost))
                    self.log.log()

                    # update storage
                    if top_loci not in PS[snode][bottom_loci]:
                        PS[snode][bottom_loci][top_loci] = []

                    # update optimum if better
                    # a solution is better if 1) it has lower cost, or 2) it has equal cost and lower ndups
                    item = (lrecon, order)                    
                    if cost < mincost:
                        FS[(bottom_loci, top_loci)] = (cost, ndup)
                        PS[snode][bottom_loci][top_loci] = [item]
                    if cost == mincost:
                        if ndup < mindup:
                            FS[(bottom_loci, top_loci)] = (cost, ndup)
                            PS[snode][bottom_loci][top_loci] = [item]
                        elif ndup == mindup:
                            FS[(bottom_loci, top_loci)] = (cost, ndup)
                            PS[snode][bottom_loci][top_loci].append(item)

                    #================================================================================

            # for each locus assignment at top and bottom of sbranch,
            # determine single optimal (lrecon, order), also keep track of cost
            self.log.log("optimal costs")            
            for bottom_loci, d in PS[snode].iteritems():
                for top_loci, lst in d.iteritems():
                    # choose a single optimum
                    if len(lst) == 1:
                        item = lst[0]
                    else:
                        item = random.choice(lst)
                    cost, ndup = FS[(bottom_loci, top_loci)]
                    PS[snode][bottom_loci][top_loci] = item + (cost,)

                    # determine cost-to-go (from root)
                    if self.prescreen:
                        if not snode.parent:
                            GS[snode][bottom_loci] = cost
                        else:
                            GS[snode][bottom_loci] = GS[snode.parent][top_loci] + cost

                    # log
                    lrecon, order = item
                    self.log.log("\t%s -> %s" % (str(top_loci), str(bottom_loci)))
                    self.log.log("\t\tlrecon: %s" % str(lrecon))
                    self.log.log("\t\torder: %s" % str(order))
                    self.log.log("\t\tcost: %g" % FS[(bottom_loci, top_loci)][0])
                    self.log.log()

            if self.prescreen:
                self.log.log("prescreen")

                # determine min cost across all bottom_loci
                mincost = min(GS[snode].values())
                self.log.log("\tmin/max cost: %d/%d" % (mincost, max(GS[snode].values())))

                if mincost <= self.prescreen_min:
                    self.log.log("\tskipping prescreen")
                else:
                    ntotal = 0
                    npruned = 0
                    thr = self.prescreen_factor * mincost
                    for bottom_loci, cost in GS[snode].iteritems():
                        ntotal += 1
                        if cost > thr:    # prescreen bottom_loci from this sbranch
                            npruned += 1
                            del PS[snode][bottom_loci]
                            self.log.log("\tpruned %s : %d" % (str(bottom_loci), cost))
                    self.log.log("\tpruned %d/%d assignments" % (npruned, ntotal))
                self.log.log()

            self.log.stop()
        self.log.stop()
            
        return PS


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
    

    def _infer_opt_locus_map(self, locus_maps, subtrees):
        """Find optimal locus map by recurring up the species tree (through dynamic programming)"""

        # locus_maps is a multi-dimensional dict with the structure
        # key1 = snode, key2 = bottom_loci, key3 = top_loci, value = (lrecon, order, cost)

        self.log.start("Inferring optimal locus map")

        gtree = self.gtree
        stree = self.stree
       
        # dynamic programming storage
        F = {}      # cost-to-go
                    # key1 = snode, key2 = top_loci, val = (bottom_loci, cost-to-go)
                    # keep track of bottom_loci for traceback

        # recur up species tree to determine (optimal) cost-to-go
        self.log.log("DP table format: key = top_loci, val = (bottom_loci, cost-to-go)")
        null_value = ((), 0)
        for snode in stree.postorder():
            self.log.start("Working on snode %s" % snode.name)
            F[snode] = {}

            if snode not in locus_maps:
                # nothing in this sbranch
                F[snode][()] = null_value
                self.log.log("Empty sbranch")
                self.log.stop()
                continue

            # get stored values for the sbranch
            locus_maps_snode = locus_maps[snode]
            subtrees_snode = subtrees[snode]

            if snode.is_leaf():
                # leaf base case
                for bottom_loci, d in locus_maps_snode.iteritems():
                    for top_loci, (lrecon, order, cost) in d.iteritems():                        
                        F[snode][top_loci] = (bottom_loci, cost)
            else:
                if len(snode.children) != 2:
                    raise Exception("non-binary species tree")

                # find cost-to-go of assigning top_loci to top of sbranch
                # = cost of top_loci to bottom_loci along sbranch
                #   + cost of bottom_loci at top of left child
                #   + cost of bottom_loci at top of right child
                sleft, sright = snode.children
                costs = collections.defaultdict(list)   # separate list for each assignment of top_loci for this sbranch
                for bottom_loci, d in locus_maps_snode.iteritems():
                    # find cost-to-go in children
                    if not self.prescreen:
                        cost_left = F[sleft][bottom_loci][1]
                        cost_right = F[sright][bottom_loci][1]
                    else:
                        # locus assignment may have been removed due to prescreening
                        if bottom_loci in F[sleft]:
                            cost_left = F[sleft][bottom_loci][1]
                        else:
                            cost_left = util.INF
                        if bottom_loci in F[sright]:
                            cost_right = F[sright][bottom_loci][1]
                        else:
                            cost_right = util.INF
                    children_cost = cost_left + cost_right

                    # add cost in this sbranch
                    for top_loci, (lrecon, order, cost) in d.iteritems():
                        cost_to_go = cost + children_cost
                        costs[top_loci].append((bottom_loci, cost_to_go))
                # select optimum cost-to-go for assigning top_loci to top of sbranch
                for top_loci, lst in costs.iteritems():
                    items, mincost = util.minall(lst,
                                                 minfunc=lambda (bottom_loci, cost_to_go): cost_to_go)
                    F[snode][top_loci] = random.choice(items)

            self.log.log("DP table")
            for top_loci, (bottom_loci, cost_to_go) in F[snode].iteritems():
                self.log.log("%s -> %s : %g" % (str(top_loci), str(bottom_loci), cost_to_go))
            self.log.stop()

        # termination
        # not necessary since cost along root sbranch already determined,
        # and by design, F[sroot] is always assigned locus = init_locus
        assert len(F[stree.root]) == 1, F[stree.root]
        top_loci, (bottom_loci, cost_to_go) = F[stree.root].items()[0]
        self.log.log("")
        self.log.log("Optimal cost: %g" % cost_to_go)
        self.log.log("")

        # recur down species tree to assign loci
        next_locus = 1
        lrecon = {}
        order = collections.defaultdict(dict)
        G = {}      # used in traceback, key = snode, val = bottom_loci of sbranch

        # initialize root
        lrecon[gtree.root] = next_locus
        next_locus += 1

        for snode in stree.preorder():
            self.log.start("Working on snode %s" % snode.name)

            # determine top_loci and bottom_loci
            if snode is stree.root:
                # root base case
                top_loci, (bottom_loci, cost_to_go) = F[snode].items()[0]
            else:
                top_loci = G[snode.parent]
                (bottom_loci, cost_to_go) = F[snode][top_loci]

            # update traceback
            G[snode] = bottom_loci
               
            # find loci and order
            if top_loci == ():
                pass
            else:
                # get stored values for the sbranch
                locus_maps_snode = locus_maps[snode]
                subtrees_snode = subtrees[snode]

                # get optimum for sbranch from DP table
                self.log.log("%s -> %s : %g" % (str(top_loci), str(bottom_loci), F[snode][top_loci][1]))
                (local_lrecon, local_order, cost) = locus_maps_snode[bottom_loci][top_loci]
                self.log.log("lrecon: %s" % str(local_lrecon))
                self.log.log("order: %s" % local_order)

                # update lrecon
                for (root, rootchild, leaves) in subtrees_snode:
                    if not rootchild:
                        continue
                    
                    for node in gtree.preorder(rootchild, is_leaf=lambda x: x in leaves):
                        if node.parent:
                            if local_lrecon[node.name] != local_lrecon[node.parent.name]:
                                lrecon[node] = next_locus
                                next_locus += 1
                            else:
                                lrecon[node] = lrecon[node.parent]

                # update order
                for plocus, lst in local_order.iteritems():
                    new_plocus = lrecon[lst[0].parent]
                    order[snode][new_plocus] = lst

            self.log.stop()

        self.lrecon = lrecon
        self.order = dict(order)
            
        self.log.stop()

#==========================================================
# tree logging

def log_tree(gtree, log, func=None, *args, **kargs):
    """print tree to log"""
           
    treeout = StringIO.StringIO()
    if not func:
        gtree.write(treeout, oneline=True, *args, **kargs)
    else:
        func(gtree, out=treeout, minlen=20, maxlen=20, *args, **kargs)
    log.log("\n%s\n" % treeout.getvalue())
    treeout.close()
    

def draw_tree_srecon(tree, srecon, *args, **kargs):
    labels = {}
    for node in tree.nodes.values():
        if not node.is_leaf():
            labels[node.name] = "%s [%s]" % (node.name, srecon[node].name)
    
    treelib.draw_tree(tree, labels, *args, **kargs)


def draw_tree_lrecon(tree, lrecon, *args, **kargs):
    labels = {}
    for node in tree.nodes.values():
        if not node.is_leaf():
            labels[node.name] = "%s" % node.name
        else:
            labels[node.name] = ""
        labels[node.name] += " (%s)" % lrecon[node]
    
    treelib.draw_tree(tree, labels, *args, **kargs)


def draw_tree_recon(tree, srecon, lrecon, *args, **kargs):
    labels = {}
    for node in tree.nodes.values():
        if not node.is_leaf():
            labels[node.name] = "%s [%s]" % (node.name, srecon[node].name)
        else:
            labels[node.name] = ""
        labels[node.name] += " (%s)" % lrecon[node]
    
    treelib.draw_tree(tree, labels, *args, **kargs)
