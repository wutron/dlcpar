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

#==========================================================

def dlc_recon(tree, stree, gene2species, gene2locus=None,
              dupcost=1, losscost=1, coalcost=1, coaldupcost=None,
              implied=True, delay=False,
              prescreen=False, prescreen_min=INF, prescreen_factor=INF,
              max_loci=INF, max_dups=INF, max_losses=INF, allow_both=False,
              log=sys.stdout):
    """Perform reconciliation using DLCoal model with parsimony costs"""

    reconer = DLCRecon(tree, stree, gene2species, gene2locus,
                       dupcost=dupcost, losscost=losscost, coalcost=coalcost, coaldupcost=coaldupcost,
                       implied=implied, delay=delay,
                       prescreen=prescreen, prescreen_min=prescreen_min, prescreen_factor=prescreen_factor,
                       max_loci=max_loci, max_dups=max_dups, max_losses=max_losses, allow_both=allow_both,
                       log=log)
    return reconer.recon()



class DLCRecon(object):

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

        assert (dupcost >= 0) and (losscost >= 0) and (coalcost >= 0) and (coaldupcost >= 0)
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
        #   self.srecon
        #   self.lrecon
        #   self.order
        #   self.nsoln
        #   self.cost

    #=============================
    # main methods

    def recon(self):
        """Perform reconciliation"""

        self.log.start("Checking feasibility")
        feasible = self._are_constraints_consistent()
        self.log.stop()
        if not feasible:
            self.log.log("tree not feasible")
            return self.gtree, None

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
        log_tree(self.gtree, self.log, func=draw_tree_srecon, srecon=self.srecon)

        # infer locus map
        self._infer_locus_map()
        self.log.log("\n\n")

        # log gene tree (with species map and locus map)
        self.log.log("gene tree (with species and locus map)\n")
        log_tree(self.gtree, self.log, func=draw_tree_recon, srecon=self.srecon, lrecon=self.lrecon)

        # revert to use input species tree
        self.stree = substree
        self.srecon = util.mapdict(self.srecon, val=lambda snode: self.stree.nodes[snode.name])
        self.order = util.mapdict(self.order, key=lambda snode: self.stree.nodes[snode.name])

        # convert to LabeledRecon data structure
        labeled_recon = reconlib.LabeledRecon(self.srecon, self.lrecon, self.order)

        # calculate runtime
        runtime = self.log.stop()

        return self.gtree, labeled_recon, self.nsoln, self.cost, runtime


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

        # find speciation subtrees and sort by id of top, also find sorted bottoms
        subtrees = reconlib.factor_tree(gtree, stree, self.srecon, self.sevents)
        sorted_bottoms = {}
        for snode in stree.preorder():
            subtrees_snode = subtrees[snode]

            if len(subtrees_snode) == 0:
                # handle root separately
                if snode is stree.root:
                    sorted_bottoms[snode] = [gtree.root]
                continue

            subtrees_snode.sort(key=lambda (root, rootchild, leaves): ids[root.name])
            bottoms_snode = []
            for (root, rootchild, leaves) in subtrees_snode:
                if leaves is not None:
                    bottoms_snode.extend(leaves)
            bottoms_snode.sort(key=lambda node: ids[node.name])
            sorted_bottoms[snode] = bottoms_snode

        # enumerate valid locus maps, then infer optimal
        self.log.start("Inferring locus map")
        locus_maps = self._enumerate_locus_maps(subtrees, sorted_bottoms)
        self.log.log("\n\n")
        self._infer_opt_locus_map(locus_maps, subtrees)
##        self._infer_trivial_locus_map()
        self.log.stop()

        # return optimal tiles for each species (needed for reconscape)
        return locus_maps


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


    def _find_bottoms(self, subtrees):
        bottoms = []
        for (root, rootchild, leaves) in subtrees:
            if leaves is not None:
                bottoms.extend(leaves)
        return bottoms


    #=============================
    # event/cost methods -- these operate at the species branch level

    def _compute_cost(self, ndup, nloss, ncoalspec, ncoaldup):
        """Find reconciliation cost"""
        return ndup*self.dupcost + nloss*self.losscost + ncoalspec*self.coalcost + ncoaldup*self.coaldupcost


    def _count_events(self, lrecon, subtrees, nodefunc=lambda node: node.name,
                      bottoms=None,
                      max_dups=INF, max_losses=INF,
                      min_cost=INF):
        """
        Count number of dup, loss, coal events and infer optimal partial order

        Return
        - number of duplications
        - number of losses
        - number of extra lineages at speciations
        - minimum number of extra lineages at duplications over all internal node orderings
        - an optimal internal node ordering
        - the number of equally parsimonious partial orderings
        - events for this tile (always None, required for reconscape)

        Wrap this information as a single element in a list.
        This is required because reconscape returns a list of the above tuples, one tuple per optimal order.
        """

        extra = {"species_map" : self.srecon, "locus_map" : lrecon}

        # defaults
        ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, cost, events = INF, INF, INF, INF, {}, INF, INF, None

        # duplications
##        ndup = reconlib.count_dup_snode(self.gtree, self.stree, extra, snode=None,
##                                        subtrees_snode=subtrees,
##                                        nodefunc=nodefunc)
        dup_nodes = reconlib.find_dup_snode(self.gtree, self.stree, extra, snode=None,
                                            subtrees_snode=subtrees,
                                            nodefunc=nodefunc)
        ndup = len(dup_nodes)
        # check if dups exceeds mincost (allow equality to select from random)
        if (ndup > max_dups) or (self._compute_cost(ndup, 0, 0, 0) > min_cost):
            return [(ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events)]

        # losses
        nloss = reconlib.count_loss_snode(self.gtree, self.stree, extra, snode=None,
                                          subtrees_snode=subtrees,
                                          nodefunc=nodefunc)
        # check if dups + losses exceeds mincost (allow equality to select from random)
        if (nloss > max_losses) or (self._compute_cost(ndup, nloss, 0, 0) > min_cost):
            return [(ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events)]

        # extra lineages at speciations
        ncoal_spec = reconlib.count_coal_spec_snode(self.gtree, self.stree, extra, snode=None,
                                                    subtrees_snode=subtrees,
                                                    nodefunc=nodefunc,
                                                    implied=self.implied)
        # check if dups + losses + coal (spec) exceeds mincost (allow equality to select from random)
        if self._compute_cost(ndup, nloss, ncoal_spec, 0) > min_cost:
            return [(ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events)]

        # extra lineages at duplications
        start = self._find_locus_orders_start(lrecon, subtrees, nodefunc=nodefunc,
                                              dup_nodes=dup_nodes, bottoms=bottoms)
        orders, nsoln = self._find_locus_orders(lrecon, subtrees, start, nodefunc=nodefunc,
                                                dup_nodes=dup_nodes, bottoms=bottoms)
        # sample optimal partial order
        order = {}
        for locus, orderings in orders.iteritems():
            order[locus] = common.random_choice(orderings)
        ncoal_dup = self._count_coal_dup(lrecon, order, start, nodefunc=nodefunc)

        return [(ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events)]


    def _count_coal_dup(self, lrecon, order, start, nodefunc=lambda node: node.name):
        """Count number of coalescences at duplications"""

        # code is adapted from reconlib.count_coal_snode_dup
        ncoal = 0
        coals = self._find_coal_dup(lrecon, order, start, nodefunc=nodefunc)
        for coal in coals:
            assert len(coal) > 1
            ncoal += len(coal) - 1
        return ncoal


    def _find_coal_dup(self, lrecon, order, start, nodefunc=lambda node: node.name):
        """Find extra lineages at duplications"""

        assert set(start) == set(order), (dict(start), order)

        # code is adapted from second routine in reconlib.find_coal_snode_dup
        coals = []
        for plocus, nodes in order.iteritems():
            current = start[plocus][:]   # DIFFERENT from reconlib: use copy!!
            num_lineages = len(current)
            for next_node in nodes:
                assert num_lineages == len(current), (num_lineages, nodes)
                assert next_node in current, (next_node, current)

                # locus of next node
                next_locus = lrecon[nodefunc(next_node)]

                # keep if leaf and locus does not change : leaves (extant genes) exist to present time
                # note: this special case may not be necessary since leaf nodes no longer in order
                if (next_node.is_leaf()) and (plocus == next_locus):
                    pass
                else:
                    current.remove(next_node)
                    num_lineages -= 1

                # update lineage count and list of nodes
                if plocus == next_locus:
                    # deep coalescence
                    # note: keep even if next_node in bottoms to allow for delay btwn coalescence and speciation
                    #       this special case may not be necessary since leaf nodes no longer in order
                    for child in next_node.children:
                        current.append(child)
                        num_lineages += 1
                else:
                    # duplication
                    if num_lineages > 1:
                        assert len(current) == num_lineages # sanity check
                        coals.append(current[:])            # use copy because current will be updated

        return coals


    def _find_locus_orders_start(self, lrecon, subtrees, nodefunc=lambda node: node.name,
                                 dup_nodes=None, bottoms=None):
        """Helper function to find starting lineages for each locus"""

        extra = {"species_map" : self.srecon, "locus_map" : lrecon}

        if dup_nodes is None:
            dup_nodes = reconlib.find_dup_snode(self.gtree, self.stree, extra, snode=None,
                                                subtrees=subtrees, nodefunc=nodefunc)

        if bottoms is None:
            bottoms = self._find_bottoms(subtrees)

        # rest of code is adapted from first routine in reconlib.count_coal_snode_dup
        start = collections.defaultdict(list)
        parent_loci = set()
        for node in dup_nodes:
            pnode = node.parent
            if pnode:
                parent_loci.add(lrecon[nodefunc(pnode)])
        # for each locus found, if this locus is a "parent locus",
        # add the children if the dup node is not a leaf
        # (leaves never incur extra lineages in this species branch)
        for node in dup_nodes:
            locus = lrecon[nodefunc(node)]
            if (node in bottoms) or (locus not in parent_loci):
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
                for child in root.children:         # add all children of root
                    if nodefunc(child) in lrecon:   # DIFFERENT from reconlib: ensure node in sbranch
                        start[locus].append(child)
            else:
                if rootchild:
                    start[locus].append(rootchild)

        return start


    def _find_locus_orders(self, lrecon, subtrees, start=None, nodefunc=lambda node: node.name,
                           dup_nodes=None, bottoms=None):
        """
        Find internal node orderings that minimize number of inferred extra lineages due to duplications.

        TODO: merge with _count_coal_dup
        """

        gtree = self.gtree
        extra = {"species_map" : self.srecon, "locus_map" : lrecon}

        if dup_nodes is None:
            dup_nodes = reconlib.find_dup_snode(self.gtree, self.stree, extra, snode=None,
                                                subtrees=subtrees, nodefunc=nodefunc)
        if bottoms is None:
            bottoms = self._find_bottoms(subtrees)

        #=============================
        # utilities

        def get_local_order(start, locus):

            def get_cost(order):
                # calculate the cost -- locus map is fixed so can ignore count changes due to ...
                # a) multiple branches having same locus at top of species branch
                # b) count decreases at duplications

                n = len(dups)
                if n == 0:
                    return 0

                cost = 0
                ndup = 0          # keep track of how many dups seen so far
                last_dup_ndx = -1 # keep track of last seen dup

                for i, node in enumerate(order):
                    if node in dups:
                        extra_lineages = i - last_dup_ndx - 1  # number of extra lineages at this dup ...
                        cost += extra_lineages * (n - ndup)    # multiplied by number of times this cost is incurred
                        ndup += 1
                        last_dup_ndx = i
                assert ndup == n, (ndup, n, order, dups)

                return cost


            def get_order_helper(curr_order):
                # recursively get optimal orderings and calculate costs for each
                # efficient alternative to enumerating all permutations of dups

                placed_dups = dups.intersection(curr_order)
                placed_nondups = set(filter(lambda node: node not in dups, curr_order))

                if len(placed_dups) == len(dups):
                    # base case: no duplications left to add
                    # add rest of nodes according to canonical order
                    # that is, sort tops of subtrees in alphnumeric order, then preorder within each subtree

                    # find roots of subtrees that are left (i.e. nodes that can be chosen immediately)
                    # a) children of placed non-dup nodes that have not been added and are not dups
                    # b) starting nodes / lineages that are not dups
                    roots = []
                    for node in placed_nondups:
                        roots.extend([child for child in node.children
                                      if child not in placed_nondups and child not in dups])
                    for node in start:
                        if node not in placed_nondups and node not in dups:
                            roots.append(node)

                    # follow canonical form specified
                    assert len(roots) == len(set(roots)), roots
                    roots.sort(key=lambda node: node.name)
                    for root in roots:
                        for node in gtree.preorder(root, is_leaf=lambda x: x in bottoms):
                            assert lrecon[nodefunc(node)] == locus, (node.name, lrecon[nodefunc(node)], locus)

                            # skip if leaf node
                            if node.is_leaf() or len(node.children) == 1:
                                continue
                            curr_order.append(node)

                    # yield order
                    yield curr_order

                else:
                    # recursive case: duplication left to add
                    # add duplications with minimum constraints, then recur

                    # update paths
                    updated_paths = util.mapdict(dup_paths, val=lambda lst: filter(lambda node: node not in placed_nondups, lst))

                    # list of dups with minimum constraints considering already placed non-dups
                    unplaced_dups = dups.difference(placed_dups)
                    best_dup_nodes, _ = util.minall(unplaced_dups, minfunc=lambda dup_node: len(updated_paths[dup_node]))
                    for dup_node in best_dup_nodes:
                        new_curr_order = curr_order[:]  # new list to not affect curr_order in other iterations
                        nodes = updated_paths[dup_node] # nodes to add above duplication (to satisfy temporal constraints)
                        new_curr_order.extend(nodes)    # add non-duplication nodes
                        new_curr_order.append(dup_node) # add duplication

                        # recur to add remaining nodes
                        for order in get_order_helper(new_curr_order):
                            yield order


            #=============================
            # main function of get_local_order

            # find nodes with parent locus = locus
            # for each node node, also find path from "start" node (non-inclusive of end-points)
            paths = {}
            for node in start:
                # recur over subtree
                for node2 in gtree.preorder(node, is_leaf=lambda x: x in bottoms):
                    pnode2 = node2.parent
                    if lrecon[nodefunc(pnode2)] == locus:
                        if node2 is node:
                            paths[node2] = collections.deque()
                        else:
                            paths[node2] = collections.deque(paths[pnode2])   # path ending before parent
                            paths[node2].append(pnode2)                       # path ending at parent

            # note: faster to consider only most parsimonious ordering above last duplication node
            #       this implementation adds nodes after last duplication even if not parsimonious

            dups = set(filter(lambda node: lrecon[nodefunc(node.parent)] == locus, dup_nodes))
            dup_paths = util.subdict(paths, dups)

            # get optimal orders (and associated costs) by recurring down optimal duplication paths
            # still need to compute costs because duplications may share ancestors
            order_count = collections.defaultdict(set)  # key = cost, val = set of duporders
            for order in get_order_helper([]):
                cost = get_cost(order)
                order_count[cost].add(tuple(order))
            optimal_dup_orders = order_count[min(order_count.keys())]

            # count optimal orders
            local_orders = [list(order) for order in optimal_dup_orders]
            nsoln = len(optimal_dup_orders)

            return local_orders, nsoln

        #=============================
        # find optimal partial orders for each locus
        orders = {}
        nsoln = 1
        for locus, nodes in start.iteritems():
            lorders, lnsoln = get_local_order(nodes, locus)
            orders[locus] = lorders
            nsoln *= lnsoln
        return orders, nsoln


    #=============================
    # constraint methods (for species-specific locus maps)

    def __find_path(self, node1, node2):
        """Find the path between two nodes in a tree

        Returns node names along path from each node up to (but excluding) lca.
        Based on treelib.find_dist(...).
        """
        # keep track of input nodes for error checking
        n1, n2 = node1, node2

        # find root path for node1 [node1, ..., root]
        path1 = [node1.name]
        while node1.parent is not None:
            node1 = node1.parent
            path1.append(node1.name)

        # find root path for node2 [node2, ..., root]
        path2 = [node2.name]
        while node2.parent is not None:
            node2 = node2.parent
            path2.append(node2.name)

        # find when paths diverge (pathX[-i+1] is the lca)
        i = 1
        while i <= len(path1) and i <= len(path2) and (path1[-i] == path2[-i]):
            i += 1
        assert path1[-i+1] == path2[-i+1] == treelib.lca((n1, n2)).name, \
            (n1.name, n2.name, path1[-i+1], path2[-i+1], treelib.lca((n1, n2)).name, i)

        return (path1[-i::-1], path2[-i::-1])


    def _find_constraints_nodups(self, leaves):
        """Determine invalid branches for duplications based on species-specific loci"""

        stree = self.stree
        gene2locus = self.gene2locus

        constraints = set()
        for snode in stree.leaves():
            # skip if no leaves
            if snode not in leaves:
                continue

            # find (species-specific) loci for this species
            # and for each locus, the set of associated genes
            loci = collections.defaultdict(set)
            for leaf in leaves[snode]:
                loci[gene2locus(leaf.name)].add(leaf)

            # for each locus, check gene pairs
            for locus, genes in loci.iteritems():
                for gene1, gene2 in itertools.combinations(genes, 2):
                    path1, path2 = self.__find_path(gene1, gene2)
                    constraints.update(path1 + path2)
        return constraints


    def _find_constraints_dups(self, leaves):
        """Determine necessary paths for duplications based on species-specific loci"""

        stree = self.stree
        gene2locus = self.gene2locus

        paths = []
        for snode in stree.leaves():
            # skip if no leaves
            if snode not in leaves:
                continue

            # find (species-specific) loci for this species
            # and for each locus, the set of associated genes
            loci = collections.defaultdict(set)
            for leaf in leaves[snode]:
                loci[gene2locus(leaf.name)].add(leaf)

            # for each pair of loci, check gene pairs
            for locus1, locus2 in itertools.combinations(loci, 2):
                genes1 = loci[locus1]
                genes2 = loci[locus2]

                for gene1, gene2 in itertools.product(genes1, genes2):
                    path1, path2 = self.__find_path(gene1, gene2)
                    paths.append((path1, path2))
        return paths


    def _are_constraints_consistent(self):
        """Return whether constraints implied by species-specific loci are consistent"""

        if self.gene2locus is None:
            return True

        gtree = self.gtree
        stree = self.stree
        gene2species = self.gene2species
        recon = phylo.reconcile(gtree, stree, gene2species)

        leaves = collections.defaultdict(list)

        for leaf in gtree.leaves():
            leaves[recon[leaf]].append(leaf)

        constraints_nodups = self._find_constraints_nodups(leaves)
        constraints_dups = self._find_constraints_dups(leaves)

        # iterate through paths that require a duplication
        for path1, path2 in constraints_dups:
            # check if no duplication is allowed along entire path
            if all([name in constraints_nodups for name in path1 + path2]):
                return False

        return True

    #=============================
    # locus (and locus map) methods

    def _find_unique_loci(self, lrecon, nodes, start=1):
        """Return unique loci"""

        # unique mapping
        x = [lrecon[node.name] for node in nodes]
        y = []
        m = {}
        next_locus = start
        for locus in x:
            if locus not in m:
                m[locus] = next_locus
                next_locus += 1
            y.append(m[locus])
        y = tuple(y)

        return y, m


    def _evolve_subtree(self, root, leaves, state, start_locus=INIT_LOCUS, next_locus=None):
        """Given state changes, find the locus at the nodes"""

        gtree = self.gtree

        if next_locus is None:
            next_locus = start_locus + 1
        assert next_locus > start_locus, (start_locus, next_locus)

        lrecon = {}
        for node in gtree.preorder(root, is_leaf=lambda x: x in leaves):
            if node is root:
                lrecon = {node.name : start_locus}
            else:
                if state[node.name]:
                    lrecon[node.name] = next_locus
                    next_locus += 1
                else:
                    lrecon[node.name] = lrecon[node.parent.name]
        return lrecon


    def _find_locus_states_sbranch(self, subtrees, is_leaf,
                                   max_loci=INF, max_dups=INF,
                                   constraints=None):
        """Find the valid locus states for each subtree in the species branch"""

        gtree = self.gtree

        # no constraints on duplications
        if constraints is None:
            constraints = set()

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
                            if child.name not in constraints:
                                s2 = state.copy();  s2[child.name] = True
                                states_down[child].append(s2)
                    elif nchildren == 2:
                        left, right = children
                        s1 = state.copy();      s1[left.name] = False;  s1[right.name] = False
                        states_down[left].append(s1)
                        if ndup < max_dups_subtree:
                            # daughter lineage matters
                            if right.name not in constraints:
                                s2 = state.copy();  s2[left.name] = False;  s2[right.name] = True
                                states_down[left].append(s2)
                            if left.name not in constraints:
                                s3 = state.copy();  s3[left.name] = True;   s3[right.name] = False
                                states_down[left].append(s3)

                            # dup along both child lineages (with or without dup from parent)
                            # is always equal or less parsimonious than dup from parent followed by dup in one child lineage
                            if self.allow_both and ndup + 1 < max_dups_subtree:
                                if (left.name not in constraints) and (right.name not in constraints):
                                    s4 = state.copy(); s4[left.name] = True; s4[right.name] = True
                                    states_down[left].append(s4)
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

            # ensure valid states (for single subtree) if sbranch is leaf branch
            # also, for all species branches, allow duplication along root branch
            states = []       # state for one subtree

            if is_leaf:
                gene2locus = self.gene2locus
                if gene2locus is None:
                    gene2locus = phylo.gene2species
                lrecon = dict([(leaf.name, gene2locus(leaf.name)) for leaf in leaves])
                bottom_loci_from_locus_map, _ = self._find_unique_loci(lrecon, leaves)

            # iterate over states
            for i, state in enumerate(states_up[rootchild][:]):
                if is_leaf:
                    # find locus at leaves (in this subtree), then check if valid
                    lrecon = self._evolve_subtree(rootchild, leaves, state)
                    bottom_loci, _ = self._find_unique_loci(lrecon, leaves)
                    if bottom_loci != bottom_loci_from_locus_map:
                        continue
                else:
                    if max_loci != INF:
                        # find locus at bottoms (in this subtree), then check if valid
                        lrecon = self._evolve_subtree(rootchild, leaves, state)
                        bottom_loci = [lrecon[node.name] for node in leaves]
                        if len(set(bottom_loci)) > max_loci:
                            continue

                # store original (without duplication along root branch)
                states.append(state)

                # allow duplication along root branch
                if root is not rootchild:
                    ndup = util.counteq(True, state.values())
                    if ndup < max_dups_subtree and rootchild.name not in constraints:
                        s1 = state.copy(); s1[rootchild.name] = True
                        states.append(s1)

            # store
            all_states.append(states)

        return all_states


    def _enumerate_locus_maps(self, subtrees, sorted_bottoms):
        """
        Enumerate the valid locus maps by recurring down the species tree.
        For each combination of locus assignments at the top and bottom of an sbranch,
        also find the optimal locus map and order.
        (This saves memory compared to storing all valid locus maps.)
        """

        self.log.start("Enumerating locus maps")

        gtree = self.gtree
        stree = self.stree
        gene2species = self.gene2species
        gene2locus = self.gene2locus

        # locus constraints
        if gene2locus is not None:
            constraints = self._find_constraints_nodups(sorted_bottoms)
        else:
            constraints = None

        # find maximum number of dup - does not assume that species map is MPR
        recon = phylo.reconcile(gtree, stree, gene2species)
        events = phylo.label_events(gtree, recon)
        max_loci_sbranch = self.max_loci
        max_dups_sbranch = self.max_dups
        max_losses_sbranch = self.max_losses
        self.log.log("Max # loci per sbranch: %s" % max_loci_sbranch)
        self.log.log("Max # dup per sbranch: %s" % max_dups_sbranch)
        self.log.log("Max # loss per sbranch: %s" % max_losses_sbranch)
        self.log.log()

        # partitions at each sbranch
        # Originally:
        #     key1 = snode, key2 = bottom_loci, key3 = top_loci
        #     val = list of items (lrecon, order, ndup, nloss, ncoalspec, ncoaldup, cost, nsoln)
        #           single item   (lrecon, order, cost, nsoln)
        #     [val is list generally and single item after filtering partitions]
        # Changing to:
        #     key1 = snode, key2 = (bottom_loci, top_loci)
        #     val = list of items (lrecon, order, ndup, nloss, ncoalspec, ncoaldup, cost, nsoln)
        #           single item   (lrecon, order, cost, nsoln)
        #     [val is list generally and single item after filtering partitions]

        PS = {}

        if self.prescreen:
            GS = collections.defaultdict(dict)

        # recur down species tree
        for snode in stree.preorder():
            self.log.start("Working on snode %s" % snode.name)
            parent_snode = snode.parent
            is_leaf = snode.is_leaf()

            # get subtrees in this sbranch
            subtrees_snode = subtrees[snode]

            # TODO: need to update if-statement to match current structure
            # nothing happened along the branch if there are no subtrees
            if len(subtrees_snode) == 0:
                # handle root separately (initialize loci)
                if snode is stree.root:
                    PS[snode] = {}
                    top_loci = bottom_loci = (INIT_LOCUS,)
                    self._initialize_partitions_sroot(PS[snode], bottom_loci, top_loci)
                    if self.prescreen:
                        self._initialize_prescreen_sroot(GS[snode], bottom_loci, top_loci)
                self.log.log("Empty sbranch")
                self.log.stop()
                continue

            # Changed to match new PS format
            # initialize storage for this sbranch
            PS[snode] = {}

            # hash subtrees using top node
            subtrees_hash = {}
            for (root, rootchild, leaves) in subtrees_snode:
                subtrees_hash[root] = (rootchild, leaves)

            # get bottom nodes for this sbranch
            bottoms = sorted_bottoms[snode]

            # get bottom loci for this sbranch
            if is_leaf:
                if gene2locus is not None:
                    lrecon = dict([(bottom.name, gene2locus(bottom.name)) for bottom in bottoms])
                    bottom_loci_from_locus_map, _ = self._find_unique_loci(lrecon, bottoms)

            # get states (changed/unchanged) for each branch of each subtree in sbranch
            states = self._find_locus_states_sbranch(subtrees_snode, is_leaf,
                                                     max_loci=max_loci_sbranch,
                                                     max_dups=INF if is_leaf else max_dups_sbranch,
                                                     constraints=constraints)

            # top of this sbranch is the bottom of the parent sbranch
            # Changed to match new PS format
            # Changed top_loci_lst from a dictionary to a set 
            if parent_snode in PS:
                top_loci_lst = set()
                for (top_loci, _) in PS[parent_snode]:
                    top_loci_lst.add(top_loci)
                tops = sorted_bottoms[parent_snode]
                # assert len(top_loci_lst) == len(set(top_loci_lst)), top_loci_lst
                assert top_loci_lst
            else:
                top_loci_lst = set([(INIT_LOCUS,)])
                tops = [gtree.root]

            # TODO: incorporate max_dup from sbranch to sbranch
            self.log.log("top nodes: %s" % ','.join(map(lambda node: node.name, tops)))
            self.log.log("bottom nodes: %s" % ','.join(map(lambda node: node.name, bottoms)))
##            self.log.log("top_loci: %s" % ';'.join(map(str, top_loci_lst.keys())))
            self.log.log("number of assignments at top nodes: %d" % len(top_loci_lst)) # TODO: make sure this is still valid for new implementation of top_loci_lst
            states_len = map(len, states)
            self.log.log("number of states: %s" % ','.join(map(str, states_len)))
            self.log.log("number of state combinations: %d" % prod(states_len))
            self.log.log("")

            # combine subtrees
            # for loop statement changed to reflect new format of top_loci_lst
            for ndx1, top_loci in enumerate(top_loci_lst): 
                # start with loci at top of sbranch
                assert len(top_loci) == len(states), (len(top_loci), len(states))

                # initialize next locus with top of sbranch
                init_next_locus = max(top_loci) + 1
                init_lrecon = {}
                for i, start in enumerate(top_loci):
                    init_lrecon[tops[i].name] = start

                for ndx2, state in enumerate(itertools.product(*states)):
                    next_locus = init_next_locus
                    lrecon = init_lrecon.copy()

                    #=============================
                    # find next lrecon using states

                    for i, start in enumerate(top_loci):
                        s = state[i]
                        if s is not None:
                            root = tops[i]
                            rootchild, leaves = subtrees_hash[root]

                            # evolve from root to rootchild
                            if s[rootchild.name]:
                                rootchild_loci = next_locus
                                next_locus += 1
                            else:
                                rootchild_loci = start

                            # evolve from rootchild to leaves
                            l = self._evolve_subtree(rootchild, leaves, state=s,
                                                     start_locus=rootchild_loci, next_locus=next_locus)
                            assert all([name not in lrecon for name in l if name != root.name]), (l, lrecon)
                            lrecon.update(l)
                            next_locus = max(init_next_locus, max(lrecon.values())) + 1

                    #=============================
                    # (unique) loci at bottom of sbranch

                    bottom_loci, mapping = self._find_unique_loci(lrecon, bottoms)

                    #=============================
                    # check validity

                    if is_leaf:
                        # checks for validity across all leaves, not just leaves in a subtree
                        if gene2locus is None:
                            if len(set(bottom_loci)) != len(bottom_loci):
                                continue
                        else:
                            if bottom_loci != bottom_loci_from_locus_map:
                                continue
                    else:
                        if len(set(bottom_loci)) > max_loci_sbranch:  # exceeds max # lineages for (ancestral) sbranch
                            self.log.log("\t%s -> %s : [%d, %d] (skipped - locus count)" % (top_loci, bottom_loci, ndx1, ndx2))
                            self.log.log("\t\tlrecon: %s" % lrecon)
                            self.log.log()
                            continue

                    #=============================
                    # find optimal cost (and order) for this lrecon
                    solns = self._find_optimal_cost(PS[snode], bottom_loci, top_loci,
                                                    lrecon, subtrees=subtrees_snode, bottoms=bottoms,
                                                    max_dups=INF if is_leaf else max_dups_sbranch,
                                                    max_losses=INF if is_leaf else max_losses_sbranch,
                                                    snode=snode)
                    # all solns should have same cost, use first as proxy for logging info
                    ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events = solns[0]

                    # skip?
                    skip = ""
                    if not is_leaf:
                        if ndup > max_dups_sbranch:         # exceeds max # of dup for (ancestral) sbranch
                            skip = "dup"
                        elif nloss > max_losses_sbranch:    # exceeds max # of loss for (ancestral) sbranch
                            skip = "loss"
                    if (skip == "") and ((ndup == INF) or (nloss == INF) or (ncoal_spec == INF) or (ncoal_dup == INF)):    # non-optimal
                        skip = "non-optimal"

                    if skip != "":
                        self.log.log("\t%s -> %s : [%d, %d] (skipped - %s)" % (top_loci, bottom_loci, ndx1, ndx2, skip))
                        self.log.log("\t\tlrecon: %s" % lrecon)
                        self.log.log()
                        continue

                    # log
                    self.log.log("\t%s -> %s : [%d, %d]" % (top_loci, bottom_loci, ndx1, ndx2))
                    self.log.log("\t\tlrecon: %s" % lrecon)
                    self.log.log("\t\torder: %s" % order)
                    self.log.log("\t\tmapping: %s" % mapping)
                    self.log.log("\t\tnsoln: %g" % nsoln)
                    self.log.log("\t\tndup: %s, nloss: %s, ncoal: %s (spec: %s, dup: %s)" % \
                                 (ndup, nloss, ncoal_spec + ncoal_dup, ncoal_spec, ncoal_dup))
                    self.log.log()

                    # update storage
                    for solution in solns:
                        ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events = solution

                        self._update_partitions(PS[snode], bottom_loci, top_loci,
                                                lrecon, order,
                                                ndup, nloss, ncoal_spec, ncoal_dup, nsoln, events)

                    #=============================

            # do any valid locus maps exist?
            if not PS[snode]:
                raise Exception("no valid locus maps exist")

            # for each locus assignment at top and bottom of sbranch,
            # filter storage
            self._filter_partitions(PS[snode])

            # prescreen
            if self.prescreen:
                self._prescreen(PS[snode], GS[snode], GS.get(parent_snode, None))

            self.log.stop()
        self.log.stop()

        return PS


    def _infer_opt_locus_map(self, locus_maps, subtrees):
        """
        Find optimal locus map by first recurring up the species tree (through dynamic programming) to determine optimal cost
        then tracing back down the species tree to determine optimal locus map
        """

        self.log.start("Inferring optimal locus map")

        # compute DP table
        F = self._dp_table(locus_maps, subtrees)

        # terminate
        self._dp_terminate(F)

        # traceback
        self._dp_traceback(locus_maps, subtrees, F)

        self.log.stop()

    #=============================
    # locus partition methods (used by _enumerate_locus_maps)
    # partition structure:
    #     key1 = bottom_loci, key2 = top_loci
    #     value = (lrecon, order, ndup, nloss, ncoalspec, ncoaldup, cost, nsoln) => (lrecon, order, cost, nsoln)
    # (note: partition is not a good descriptor)

    def _initialize_partitions_sroot(self, partitions, bottom_loci, top_loci):
        # Changed to match new PS format
        
        lrecon = {self.gtree.root.name: INIT_LOCUS}
        order = {}
        cost = 0
        nsoln = 1
        partitions[(bottom_loci,top_loci)] = (lrecon, order, cost, nsoln)


    def _find_optimal_cost(self, partitions, bottom_loci, top_loci,
                           lrecon, subtrees, bottoms=None,
                           max_dups=INF, max_losses=INF, snode=None):

        # function changed to reflect new structure of PS
        mincost = INF
        if (bottom_loci,top_loci) in partitions:
            # lst contains items (lrecon, order, ndup, nloss, ncoalspec, ncoaldup, cost, nsoln)
            mincost = min([item[6] for item in partitions[(bottom_loci,top_loci)]])

        solns = self._count_events(lrecon, subtrees, bottoms=bottoms,
                                   max_dups=max_dups, max_losses=max_losses,
                                   min_cost=mincost)

        return solns

    def _update_partitions(self, partitions, bottom_loci, top_loci,
                           lrecon, order,
                           ndup, nloss, ncoalspec, ncoaldup, nsoln, events):
        # a solution is better if 1) it has lower cost, or 2) it has equal cost and lower ndups

        # Changed to match new PS format
        
        mincost, mindup = INF, INF
        if (bottom_loci,top_loci) in partitions:
            # lst contains items (lrecon, order, ndup, nloss, ncoalspec, ncoaldup, cost, nsoln)
            mincost, mindup = min([(item[6], item[2]) for item in partitions[(bottom_loci,top_loci)]])

        cost = self._compute_cost(ndup, nloss, ncoalspec, ncoaldup)
        item = (lrecon, order, ndup, nloss, ncoalspec, ncoaldup, cost, nsoln)

        if (bottom_loci,top_loci) not in partitions:
            partitions[(bottom_loci,top_loci)] = []

        if cost < mincost:
            partitions[(bottom_loci,top_loci)] = [item]
        elif cost == mincost:
            if ndup < mindup:
                partitions[(bottom_loci,top_loci)] = [item]
            elif ndup == mindup:
                assert item not in partitions[(bottom_loci,top_loci)]
                partitions[(bottom_loci,top_loci)].append(item)


    def _filter_partitions(self, partitions):
        """
        JQ - updates the value of locus_maps to be val = (lrecon, order, cost, total_nsoln)
        """
        # partitions is PS[snode].
        # JQ-partitions: key1 = (bottom_loci, top_loci), value = (lrecon, order, ndup, nloss, ncoalspec, ncoaldup, cost, nsoln)
        # JQ-partitions: key1 = (bottom_loci, top_loci), value = (lrecon, order, cost, total_nsoln)

        self.log.log("optimal costs")

        #JQ - need to change this loop:
        #MT-commented out old version
        """
        for bottom_loci, d in partitions.iteritems():
            for top_loci, lst in d.iteritems():
        """ 
        for (bottom_loci, top_loci) in partitions.iteritems():     
            # lst contains items (lrecon, order, ndup, nloss, ncoalspec, ncoaldup, cost, nsoln)
            # where nsoln is the number of partial orderings that are optima for the lrecon

            # if multiple optima exist (len(lst) > 1):
            # a) add number of solutions across multiple optima
            # b) choose single optimum (lrecon, order) based on weights
            #    weigh the choices according to number of optimal partial orderings per locus map
            # note: numpy cannot take list of lists so use indices
            nsoln = [item[7] for item in lst]
            total_nsoln = sum(nsoln)
            weights = map(lambda val: float(val) / total_nsoln, nsoln)
            item = common.random_choice(lst, p=weights)

            # set the chosen optimum back into partitions
            # update number of solutions to be total across all optimal locus maps
            lrecon, order, ndup, nloss, ncoalspec, ncoaldup, cost, nsoln = item
            item = (lrecon, order, cost, total_nsoln)

            #JQ - change this line
            #MT-commented out old version

            partitions[(bottom_loci, top_loci)] = item
            """
            partitions[bottom_loci][top_loci] = item
            """

            # log
            self.log.log("\t%s -> %s" % (top_loci, bottom_loci))
            self.log.log("\t\tlrecon: %s" % lrecon)
            self.log.log("\t\torder: %s" % order)
            self.log.log("\t\tcost: %g" % cost)
            self.log.log("\t\tnsoln: %g" % total_nsoln)
            self.log.log()

    #=============================
    # prescreen methods (used by _enumerate_locus_maps)
    # key1 = snode, key2 = bottom_loci, value = min cost-to-go (from root) required to assign bottom_loci to sbranch

    def _initialize_prescreen_sroot(self, prescreen, bottom_loci, top_loci):
        prescreen[bottom_loci] = 0


    def _prescreen(self, partitions, prescreen, prescreen_parent):
        self.log.log("prescreen")

        # update min cost-to-go (from root)
        # TODO: need to collapse these for loops into one that iterates over (bottom_loci,top_loci)
        # GS structure: GS[snode][bottom_loci] = cost
        """
        # Original version
        for bottom_loci, d in partitions.iteritems():
            cost_lst = []
            for top_loci, item in d.iteritems():
                lrecon, order, cost, nsoln = item
                if prescreen_parent is not None:
                    parent_cost = prescreen_parent[top_loci]
                else:
                    parent_cost = 0
                cost_lst.append(parent_cost + cost)
            prescreen[bottom_loci] = min(cost_lst)
        """

        for (bottom_loci,top_loci), items in partitions.iteritems():
            lrecon, order, cost, nsoln = item
            if prescreen_parent is not None:
                parent_cost = prescreen_parent[top_loci]
            else:
                parent_cost = 0
            total_cost = parent_cost + cost
            if bottom_loci in prescreen:
                if prescreen[bottom_loci] > total_cost:
                    prescreen[bottom_loci] = total_cost
            else:
                prescreen[bottom_loci] = total_cost


        # determine min cost across all bottom_loci
        mincost = min(prescreen.values())
        self.log.log("\tmin cost: %d" % mincost)

        if mincost <= self.prescreen_min:
            self.log.log("\tskipping prescreen")
        else:
            ntotal = 0
            npruned = 0
            thr = self.prescreen_factor * mincost
            for bottom_loci, cost in prescreen.iteritems():
                ntotal += 1
                if cost > thr:    # prescreen bottom_loci from this sbranch
                    npruned += 1
                    del partitions[bottom_loci]
                    self.log.log("\tpruned %s : %d" % (bottom_loci, cost))
            self.log.log("\tpruned %d/%d assignments" % (npruned, ntotal))
        self.log.log()


    #=============================
    # DP table methods (used by _infer_opt_locus_map)

    def _dp_table(self, locus_maps, subtrees):
        # locus_maps is a multi-dimensional dict with the structure
        # JQ-old: key1 = snode, key2 = bottom_loci, key3 = top_loci, value = (lrecon, order, cost, nsoln)
        
        # JQ-refactored: key1 = snode, key2 = (bottom_loci, top_loci), value = (lrecon, order, cost, nsoln)
        # JQ-F is unchanged

        stree = self.stree
        gene2locus = self.gene2locus

        # dynamic programming storage
        F = {}      # key1 = snode, key2 = top_loci
                    # val = (bottom_loci, cost-to-go, nsoln)

        for snode in stree.postorder():
            self.log.start("Working on snode %s" % snode.name)
            F[snode] = {}

            if snode not in locus_maps: #JQ-should be good, still has snode as first key
                # nothing in this sbranch
                F[snode][()] = ((), 0, 1)
                self.log.log("Empty sbranch")
                self.log.stop()
                continue

            # get stored values for the sbranch
            locus_maps_snode = locus_maps[snode] 
            #JQ-locus_maps_snode: currently is dict with key1 = bottom_loci, key2 = top_loci, value = (lrecon, order, cost, nsoln)
            #JQ-now is dict with key1 = (bottom_loci, top_loci), value = (lrecon, order, cost, nsoln)
            
            subtrees_snode = subtrees[snode]

            if snode.is_leaf():
                # leaf base case

                #JQ-this for loop has to be changed
                #JQ-possible refactoring:
                #MT-commented out old version
                
                for (bottom_loci, top_loci),(lrecon, order, cost, nsoln) in locus_maps_snode.iteritems():
                    assert top_loci not in F[snode]
                    F[snode][top_loci] = (bottom_loci, cost, nsoln)
                """
                for bottom_loci, d in locus_maps_snode.iteritems():
                    for top_loci, (lrecon, order, cost, nsoln) in d.iteritems():
                        assert top_loci not in F[snode]
                        F[snode][top_loci] = (bottom_loci, cost, nsoln)
                """
            else:
                if len(snode.children) != 2:
                    raise Exception("non-binary species tree")

                # find
                # a) cost-to-go of assigning top_loci to top of sbranch
                #    = cost of top_loci to bottom_loci along sbranch
                #      + cost of bottom_loci at top of left child
                #      + cost of bottom_loci at top of right child
                # b) nsoln-to-go of assigning top_loci to top of sbranch
                #    = nsoln of top_loci to bottom_loci along sbranch
                #      * nsoln of bottom_loci at top of left child
                #      * nsoln of bottom_loci at top of right child

                sleft, sright = snode.children
                costs = collections.defaultdict(list)   # separate list for each assignment of top_loci for this sbranch
                
                
                #JQ-have to change for loop. currently looping over bottom_loci and then top_loci
                #proposed refactoring:
                #MT-commented out old version
                

                for (bottom_loci, top_loci),(lrecon, order, cost, nsoln) in locus_maps_snode.iteritems():
                    # find cost-to-go and nsoln-to-go in children
                    # locus assignment may have been removed due to search heuristics
                    _, cost_left, nsoln_left = F[sleft].get(bottom_loci, (None, INF, 0))
                    _, cost_right, nsoln_right = F[sright].get(bottom_loci, (None, INF, 0))

                    # update cost-to-go and nsoln-to-go based on children
                    children_cost = cost_left + cost_right
                    children_nsoln = nsoln_left * nsoln_right

                    cost_to_go = cost + children_cost
                    nsoln_to_go = nsoln * children_nsoln
                    item = (bottom_loci, cost_to_go, nsoln_to_go)
                    assert item not in costs[top_loci], (snode, bottom_loci, top_loci, item)
                    costs[top_loci].append(item)
                """
                for bottom_loci, d in locus_maps_snode.iteritems():
                    # find cost-to-go and nsoln-to-go in children
                    # locus assignment may have been removed due to search heuristics
                    _, cost_left, nsoln_left = F[sleft].get(bottom_loci, (None, INF, 0))
                    _, cost_right, nsoln_right = F[sright].get(bottom_loci, (None, INF, 0))

                    # update cost-to-go and nsoln-to-go based on children
                    children_cost = cost_left + cost_right
                    children_nsoln = nsoln_left * nsoln_right

                    # add cost and multiply nsoln in this sbranch
                    for top_loci, (lrecon, order, cost, nsoln) in d.iteritems():
                        cost_to_go = cost + children_cost
                        nsoln_to_go = nsoln * children_nsoln
                        item = (bottom_loci, cost_to_go, nsoln_to_go)
                        assert item not in costs[top_loci], (snode, bottom_loci, top_loci, item)
                        costs[top_loci].append(item)
                """
                
                # select optimum cost-to-go for assigning top_loci to top of sbranch
                for top_loci, lst in costs.iteritems():
                    # lst contains items (bottom_loci, cost_to_go, nsoln_to_go)
                    items, mincost = util.minall(lst, minfunc=lambda item: item[1])

                    if gene2locus is None:
                        # single sample per species
                        # cannot have all solutions have infinite cost
                        # (e.g. all viable solutions removed due to heuristics)
                        assert mincost != INF, (top_loci, lst)
                    else:
                        # multiple samples per species
                        # locus maps per branch are only checked for validity once all locus maps are put together
                        # this top_loci is invalid because it never satisfied duplication constraints
                        if mincost == INF:
                            continue

                    # item = (bottom_loci, cost_to_go, nsoln_to_go)
                    nsoln = [item[2] for item in items]
                    total_nsoln = sum(nsoln)
                    weights = map(lambda val: float(val) / total_nsoln, nsoln)

                    # choose uniform sample by weighting choices based on number of solutions
                    item = common.random_choice(items, p=weights)

                    # set the chosen optimum
                    # update number of solutions to be total across all optimal choices
                    bottom_loci, cost_to_go, nsoln_to_go = item
                    item = (bottom_loci, cost_to_go, total_nsoln)
                    F[snode][top_loci] = item

            self.log.log("DP table")
            for top_loci, (bottom_loci, cost_to_go, nsoln_to_go) in F[snode].iteritems():
                self.log.log("%s -> %s : %g [%g]" % (top_loci, bottom_loci, cost_to_go, nsoln_to_go))
            self.log.stop()

        return F



    def _dp_terminate(self, F):
        stree = self.stree

        # not necessary since cost along root sbranch already determined,
        # and by design, F[sroot] is always assigned locus = INIT_LOCUS
        assert len(F[stree.root]) == 1, F[stree.root]
        bottom_loci, cost_to_go, nsoln_to_go = F[stree.root].values()[0]
        self.log.log("")
        self.log.log("Optimal cost: %g" % cost_to_go)
        self.log.log("Number of solutions: %g" % nsoln_to_go)
        self.log.log("")



    def _dp_traceback(self, locus_maps, subtrees, F):

        gtree = self.gtree
        stree = self.stree

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
                top_loci, (bottom_loci, cost_to_go, nsoln_to_go) = F[snode].items()[0]
            else:
                top_loci = G[snode.parent]
                (bottom_loci, cost_to_go, nsoln_to_go) = F[snode][top_loci]

            # update traceback
            G[snode] = bottom_loci

            # find loci and order
            if top_loci == ():
                pass
            else:
                # get stored values for the sbranch
                locus_maps_snode = locus_maps[snode] 
                #JQ-locus_maps_snode: currently is dict with key1 = bottom_loci, key2 = top_loci, value = (lrecon, order, cost, nsoln)
                #JQ-now is dict with key1 = (bottom_loci, top_loci), value = (lrecon, order, cost, nsoln) 
               
                subtrees_snode = subtrees[snode]

                # get optimum for sbranch from DP table
                self.log.log("%s -> %s : %g [%g]" % \
                             (top_loci, bottom_loci, F[snode][top_loci][1], F[snode][top_loci][2]))
                

                #MT-commented out old version
                """
                #the following line needs to be changed
                local_lrecon, local_order, cost, nsoln = locus_maps_snode[bottom_loci][top_loci]
                #proposed refactoring:
                """
                
                local_lrecon, local_order, cost, nsoln = locus_maps_snode[(bottom_loci,top_loci)]
                
                
                self.log.log("lrecon: %s" % local_lrecon)
                self.log.log("order: %s" % local_order)

                # update lrecon
                for (root, rootchild, leaves) in subtrees_snode:
                    if not rootchild:
                        continue

                    for node in gtree.preorder(rootchild, is_leaf=lambda x: x in leaves):
                        pnode = node.parent
                        if pnode:
                            if local_lrecon[node.name] != local_lrecon[pnode.name]:
                                lrecon[node] = next_locus
                                next_locus += 1
                            else:
                                lrecon[node] = lrecon[pnode]

                # update order
                for plocus, lst in local_order.iteritems():
                    new_plocus = lrecon[lst[0].parent]
                    order[snode][new_plocus] = lst

            self.log.stop()

        self.lrecon = lrecon
        self.order = dict(order)
        self.cost = F[stree.root].values()[0][1]
        self.nsoln = F[stree.root].values()[0][2]


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
