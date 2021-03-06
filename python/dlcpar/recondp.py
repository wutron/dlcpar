"""
recondp.py
Library to solve DLC MPR Problem using DP
"""


# python libraries
import collections
import itertools
import operator
import sys

# rasmus, compbio libraries
from rasmus import timer
from rasmus import treelib
from rasmus import util
from compbio import phylo

# dlcpar libraries
from dlcpar import common
from dlcpar import constants
from dlcpar import reconlib

# gtree = G (with implied speciation nodes)
# stree = S
# srecon = M
# lrecon = L

#==========================================================
# globals

INF = util.INF
INIT_LOCUS = 1

def prod(iterable):
    """Return the product of an iterable"""
    return reduce(operator.mul, iterable, 1)

#==========================================================

def dlc_recon(tree, stree, gene2species, gene2locus=None,
              dupcost=constants.DEFAULT_DUP_COST,
              losscost=constants.DEFAULT_LOSS_COST,
              coalcost=constants.DEFAULT_COAL_COST,
              coaldupcost=None,
              implied=True, delay=False,
              prescreen=False, prescreen_min=INF, prescreen_factor=INF,
              max_loci=INF, max_dups=INF, max_losses=INF,
              allow_both=False,
              log=sys.stdout):
    """Perform reconciliation using dynamic programming"""

    reconer = DLCRecon(tree, stree, gene2species, gene2locus,
                       dupcost=dupcost, losscost=losscost,
                       coalcost=coalcost, coaldupcost=coaldupcost,
                       implied=implied, delay=delay,
                       prescreen=prescreen, prescreen_min=prescreen_min, prescreen_factor=prescreen_factor,
                       max_loci=max_loci, max_dups=max_dups, max_losses=max_losses,
                       allow_both=allow_both,
                       log=log)
    return reconer.recon()



class DLCRecon(object):
    """Reconciliation class using dynamic programming"""

    def __init__(self, gtree, stree, gene2species, gene2locus=None,
                 dupcost=constants.DEFAULT_DUP_COST,
                 losscost=constants.DEFAULT_LOSS_COST,
                 coalcost=constants.DEFAULT_COAL_COST,
                 coaldupcost=None,
                 implied=True, delay=False,
                 prescreen=False, prescreen_min=INF, prescreen_factor=INF,
                 max_loci=INF, max_dups=INF, max_losses=INF,
                 allow_both=False,
                 name_internal="n", log=sys.stdout):

        # rename gene tree nodes
        common.rename_nodes(gtree, name_internal)

        self.gtree = gtree
        self.stree = stree
        self.gene2species = gene2species
        self.gene2locus = gene2locus

        assert (dupcost > 0) and (losscost > 0)
        assert (coalcost > 0) and (coaldupcost is None or coaldupcost > 0)
        self.dupcost = dupcost
        self.losscost = losscost
        self.coalcost = coalcost  # coalspeccost, using coalcost for backwards compatibility
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
        self.log = timer.Timer(log)

        # these attributes are assigned when performing reconciliation using self.recon()
        self.srecon = None
        self.lrecon = None
        self.order = None
        self.nsoln = None
        self.cost = None

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
        sevents = phylo.label_events(self.gtree, self.srecon)
        reconlib.add_implied_nodes(self.gtree, self.stree, self.srecon, sevents, delay=self.delay)
        common.rename_nodes(self.gtree, self.name_internal)

        # log gene tree (with species map)
        self.log.log("gene tree (with species map)\n")
        reconlib.log_tree(self.gtree, self.log, func=reconlib.draw_tree_recon,
                          srecon=self.srecon)

        # infer locus map
        self._infer_locus_map()
        self.log.log("\n\n")

        # log gene tree (with species map and locus map)
        self.log.log("gene tree (with species and locus map)\n")
        reconlib.log_tree(self.gtree, self.log, func=reconlib.draw_tree_recon,
                          srecon=self.srecon, lrecon=self.lrecon)

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
        self.log.stop()


    def _infer_locus_map(self):
        """Infer (and assign) locus map"""

        gtree = self.gtree
        stree = self.stree
        srecon = self.srecon

        # find speciation subtrees and sort by id of top, also find sorted bottoms
        sevents = phylo.label_events(gtree, srecon)
        subtrees = reconlib.factor_tree(gtree, stree, srecon, sevents)
        sorted_bottoms = {}
        for snode in stree.preorder():
            subtrees_snode = subtrees[snode]

            if not subtrees_snode: # len(subtrees_snode) == 0
                # handle root separately
                if snode is stree.root:
                    sorted_bottoms[snode] = [gtree.root]
                continue

            subtrees_snode.sort(key=lambda (top, topchild, bottoms): top.name)
            bottoms_snode = []
            for (_, _, bottoms) in subtrees_snode:
                if bottoms is not None:
                    bottoms_snode.extend(bottoms)
            bottoms_snode.sort(key=lambda node: node.name)
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


    def _find_all_bottoms(self, subtrees):
        all_bottoms = []
        for (_, _, bottoms) in subtrees:
            if bottoms is not None:
                all_bottoms.extend(bottoms)
        return all_bottoms


    #=============================
    # event/cost methods -- these operate at the species branch level

    def _compute_cost(self, ndup, nloss, ncoalspec, ncoaldup):
        """Find reconciliation cost"""
        return ndup*self.dupcost + nloss*self.losscost \
            + ncoalspec*self.coalcost + ncoaldup*self.coaldupcost


    def _count_events(self, lrecon, subtrees, nodefunc=lambda node: node.name,
                      all_bottoms=None,
                      max_dups=INF, max_losses=INF,
                      min_cost=INF):
        """Count number of dup, loss, coal events and infer optimal partial order

        Return
        - number of duplications
        - number of losses
        - number of extra lineages at speciations
        - minimum number of extra lineages at duplications over all internal node orderings
        - an optimal internal node ordering
        - the number of equally parsimonious partial orderings
        - events for this tile (always None, required for reconscape)

        Wrap this information as a single element in a list.
        This is required because reconscape returns a list of the above tuples,
        one tuple per optimal order.
        """

        extra = {"species_map" : self.srecon, "locus_map" : lrecon}

        # defaults
        ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln = INF, INF, INF, INF, {}, INF
        events = None

        # duplications
##        ndup = reconlib.count_dup_snode(self.gtree, self.stree, extra, snode=None,
##                                        subtrees_snode=subtrees,
##                                        nodefunc=nodefunc)
        dup_nodes = reconlib.find_dup_snode(self.gtree, self.stree, extra, snode=None,
                                            subtrees_snode=subtrees,
                                            nodefunc=nodefunc)
        ndup = len(dup_nodes)
        # check if dups exceeds mincost (allow eq to select from random)
        if (ndup > max_dups) or (self._compute_cost(ndup, 0, 0, 0) > min_cost):
            return [(ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events)]

        # losses
        nloss = reconlib.count_loss_snode(self.gtree, self.stree, extra, snode=None,
                                          subtrees_snode=subtrees,
                                          nodefunc=nodefunc)
        # check if dups + losses exceeds mincost (allow eq to select from random)
        if (nloss > max_losses) or (self._compute_cost(ndup, nloss, 0, 0) > min_cost):
            return [(ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events)]

        # extra lineages at speciations
        ncoal_spec = reconlib.count_coal_spec_snode(self.gtree, self.stree, extra, snode=None,
                                                    subtrees_snode=subtrees,
                                                    nodefunc=nodefunc,
                                                    implied=self.implied)
        # check if dups + losses + coal (spec) exceeds mincost (allow eq to select from random)
        if self._compute_cost(ndup, nloss, ncoal_spec, 0) > min_cost:
            return [(ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events)]

        # extra lineages at duplications
        start = self._find_locus_orders_start(lrecon, subtrees, nodefunc=nodefunc,
                                              dup_nodes=dup_nodes, all_bottoms=all_bottoms)
        orders, nsoln = self._find_locus_orders(lrecon, subtrees, start, nodefunc=nodefunc,
                                                dup_nodes=dup_nodes, all_bottoms=all_bottoms)

        # sample optimal partial order uniformly at random
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
        for lineages in coals:
            assert len(lineages) > 1
            ncoal += len(lineages) - 1
        return ncoal


    def _find_coal_dup(self, lrecon, order, start, nodefunc=lambda node: node.name):
        """Find contemporary lineages at duplications"""

        coals = self._find_contemporary_lineages(lrecon, order, start, nodefunc=nodefunc)
        new_coals = []
        for lineages in coals.itervalues():
            if len(lineages) > 1:
                new_coals.append(lineages)
        return new_coals


    def _find_contemporary_lineages(self, lrecon, order, start, nodefunc=lambda node: node.name):
        """Find contemporary lineages at duplications"""

        assert set(start) == set(order), (dict(start), order)

        # code is adapted from second routine in reconlib.find_coal_snode_dup
        coals = {}  # key = dup node, val = list of contermporary lineages
        for plocus, nodes in order.iteritems():
            current = start[plocus][:]   # DIFFERENT from reconlib: use copy!!
            num_lineages = len(current)
            for next_node in nodes:
                assert num_lineages == len(current), (num_lineages, nodes)
                assert next_node in current, (next_node, current)

                # locus of next node
                next_locus = lrecon[nodefunc(next_node)]

                # keep if leaf and locus does not change
                # leaves (extant genes) exist to present time
                # TODO: this case may not be necessary since leaf nodes no longer in order
                if (next_node.is_leaf()) and (plocus == next_locus):
                    pass
                else:
                    current.remove(next_node)
                    num_lineages -= 1

                # update lineage count and list of nodes
                if plocus == next_locus:
                    # deep coalescence
                    # note: keep even if next_node in bottoms to allow for delay btwn coal and spec
                    # TODO: this may not be necessary since leaf nodes no longer in order
                    for child in next_node.children:
                        current.append(child)
                        num_lineages += 1
                else:
                    # duplication
                    assert len(current) == num_lineages # sanity check
                    coals[next_node] = current[:]       # use copy because current will be updated

        return coals


    def _find_locus_orders_start(self, lrecon, subtrees, nodefunc=lambda node: node.name,
                                 dup_nodes=None, all_bottoms=None):
        """Helper function to find starting lineages for each locus"""

        extra = {"species_map" : self.srecon, "locus_map" : lrecon}

        if dup_nodes is None:
            dup_nodes = reconlib.find_dup_snode(self.gtree, self.stree, extra, snode=None,
                                                subtrees=subtrees, nodefunc=nodefunc)

        if all_bottoms is None:
            all_bottoms = self._find_all_bottoms(subtrees)

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
            if (node in all_bottoms) or (locus not in parent_loci):
                continue
            for child in node.children:
                start[locus].append(child)
        # for each locus that exists at the top of the species branch,
        # if this locus is a "parent locus",
        # add the child if it exists (assume immediate loss if child does not exist)
        for (top, topchild, _) in subtrees:
            locus = lrecon[nodefunc(top)]
            if locus not in parent_loci:
                continue

            # handle root separately
            if not top.parent:
                for child in top.children:          # add all children of top
                    if nodefunc(child) in lrecon:   # DIFFERENT from reconlib: ensure node in sbranch
                        start[locus].append(child)
            else:
                if topchild:
                    start[locus].append(topchild)

        return start


    def _find_locus_orders(self, lrecon, subtrees, start=None, nodefunc=lambda node: node.name,
                           dup_nodes=None, all_bottoms=None):
        """Find partial order that minimizes number of inferred extra lineages due to duplications.

        TODO: merge with _count_coal_dup
        """

        gtree = self.gtree
        extra = {"species_map" : self.srecon, "locus_map" : lrecon}

        if dup_nodes is None:
            dup_nodes = reconlib.find_dup_snode(self.gtree, self.stree, extra, snode=None,
                                                subtrees=subtrees, nodefunc=nodefunc)
        if all_bottoms is None:
            all_bottoms = self._find_all_bottoms(subtrees)

        #=============================
        # utilities

        def get_local_order(start, locus):
            """Find optimal partial order for single parent locus"""

            def get_cost(order):
                """calculate cost -- locus map is fixed so can ignore count changes due to ...
                a) multiple branches having same locus at top of species branch
                b) count decreases at duplications
                """

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
                """recursively get optimal orderings and calculate costs for each
                efficient alternative to enumerating all permutations of dups"""

                used_dups = dups.intersection(curr_order)
                unused_dups = dups.difference(used_dups)
                used_nondups = set([node for node in curr_order if node not in dups])

                if not unused_dups: # len(unused_dups) == 0
                    # base case: no duplications left to add
                    # add rest of nodes in canonical order
                    # (sort tops of subtrees, then preorder within each subtree)

                    # find tops of subtrees that are left (these nodes can be chosen immediately)
                    # a) children of used non-dup nodes that have not been added and are not dups
                    # b) starting nodes / lineages that are not dups
                    tops = []
                    for node in used_nondups:
                        tops.extend([child for child in node.children
                                     if child not in used_nondups and child not in dups])
                    for node in start:
                        if node not in used_nondups and node not in dups:
                            tops.append(node)

                    # follow canonical form specified
                    assert len(tops) == len(set(tops)), tops
                    tops.sort(key=lambda node: node.name) # tops in alphanumeric order
                    for top in tops:                      # preorder within each subtree
                        for node in gtree.preorder(top, is_leaf=lambda x: x in all_bottoms):
                            assert lrecon[nodefunc(node)] == locus, \
                                (node.name, lrecon[nodefunc(node)], locus)

                            # skip if leaf node
                            if node.is_leaf() or (len(node.children) == 1 and node in all_bottoms):
                                continue
                            curr_order.append(node)

                    # yield order
                    yield curr_order

                else:
                    # recursive case: duplication left to add
                    # add duplications with minimum constraints, then recur

                    # update paths
                    updated_paths = {}
                    for dup_node, lst in dup_paths.iteritems():
                        updated_paths[dup_node] = [node for node in lst if node not in used_nondups]

                    # get dups with minimum constraints considering already placed non-dups
                    best_dup_nodes, _ = \
                        util.minall(unused_dups,
                                    minfunc=lambda dup_node: len(updated_paths[dup_node]))

                    # add each dup to order
                    for dup_node in best_dup_nodes:
                        new_curr_order = curr_order[:]  # new list to not affect curr_order in other iterations
                        nodes = updated_paths[dup_node] # nodes to add above duplication (to satisfy temporal constraints)
                        new_curr_order.extend(nodes)    # add non-duplication nodes
                        new_curr_order.append(dup_node) # add duplication

                        # recur to add more duplications
                        for order in get_order_helper(new_curr_order):
                            yield order


            #=============================
            # main function of get_local_order

            # find nodes with parent locus = locus
            # for each node node, also find path from "start" node (non-inclusive of end-points)
            paths = {}
            for node in start:
                # recur over subtree
                for node2 in gtree.preorder(node, is_leaf=lambda x: x in all_bottoms):
                    pnode2 = node2.parent
                    if lrecon[nodefunc(pnode2)] == locus:
                        if node2 is node:
                            paths[node2] = collections.deque()
                        else:
                            paths[node2] = collections.deque(paths[pnode2]) # path ending before parent
                            paths[node2].append(pnode2)                     # path ending at parent

            # note: faster to consider only most parsimonious ordering above last duplication node
            #       this implementation adds nodes after last duplication even if not parsimonious

            dups = set([node for node in dup_nodes if lrecon[nodefunc(node.parent)] == locus])
            dup_paths = util.subdict(paths, dups)

            # get optimal orders by recurring down optimal duplication paths
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
            for genes in loci.itervalues():
                for gene1, gene2 in itertools.combinations(genes, 2):
                    path1, path2 = common.find_path(gene1, gene2)
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
                    path1, path2 = common.find_path(gene1, gene2)
                    paths.append(path1 + path2)
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
        for path in constraints_dups:
            # check if no duplication is allowed along entire path
            if all([name in constraints_nodups for name in path]):
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
        # each node of tree is labeled T/F corresponding to whether dup occurred along branch
        all_states = []

        # iterate through subtrees in sbranch to find leaf states
        # start at topchild since dup along (top, topchild) will be handled when combining subtrees
        for (top, topchild, bottoms) in subtrees:
            if not topchild:
                # loss
                all_states.append([None])
                continue

            # find maximum number of dup in subtree
            # K can be less than Ksub due to speciation nodes
            max_dups_subtree = len(bottoms)
            if max_dups < max_dups_subtree:
                max_dups_subtree = max_dups

            # initialize states storage and topchild node
            # TODO: more efficient to use matrix to store states
            states_down = collections.defaultdict(list)
            state = {topchild.name : False}
            states_down[topchild].append(state)

            # recur down subtree to find leaf states
            for node in gtree.preorder(topchild, is_leaf=lambda x: x in bottoms):
                if node in bottoms:
                    continue

                for state in states_down[node]:
                    ndup = util.counteq(True, state.values())
                    children = node.children
                    nchildren = len(children)

                    if nchildren == 0:
                        raise Exception("invalid number of children: %s" % node.name)

                    elif nchildren == 1:
                        child = children[0]
                        s1 = state.copy()
                        s1[child.name] = False
                        states_down[child].append(s1)

                        if ndup < max_dups_subtree:
                            if child.name not in constraints:
                                s2 = state.copy()
                                s2[child.name] = True
                                states_down[child].append(s2)

                    elif nchildren == 2:
                        left, right = children
                        s1 = state.copy()
                        s1[left.name] = False
                        s1[right.name] = False
                        states_down[left].append(s1)

                        if ndup < max_dups_subtree:
                            # daughter lineage matters
                            if right.name not in constraints:
                                s2 = state.copy()
                                s2[left.name] = False
                                s2[right.name] = True
                                states_down[left].append(s2)
                            if left.name not in constraints:
                                s3 = state.copy()
                                s3[left.name] = True
                                s3[right.name] = False
                                states_down[left].append(s3)

                            # dup along both child lineages (with or w/o dup from parent)
                            # is always equal or worse than dup from parent
                            # followed by dup in one child lineage
                            if (self.allow_both \
                                and ndup + 1 < max_dups_subtree \
                                and left.name not in constraints \
                                and right.name not in constraints):
                                s4 = state.copy()
                                s4[left.name] = True
                                s4[right.name] = True
                                states_down[left].append(s4)
                        states_down[right] = states_down[left]

                    else:
                        raise Exception("invalid number of children: %s" % node.name)

            # recur up subtree to combine leaf states
            # TODO: check if copies are needed
            states_up = collections.defaultdict(list)
            for node in gtree.postorder(topchild, is_leaf=lambda x: x in bottoms):
                # base case (leaf)
                if node in bottoms:
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

                    # another base case (both left and right children are bottoms)
                    if (left in bottoms) and (right in bottoms):
                        assert states_down[left] == states_down[right], \
                            (states_down[left], states_down[right])
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
                lrecon = dict([(bottom.name, gene2locus(bottom.name)) for bottom in bottoms])
                bottom_loci_from_locus_map, _ = self._find_unique_loci(lrecon, bottoms)

            # iterate over states
            for state in states_up[topchild][:]:
                if is_leaf:
                    # find locus at bottoms (in this subtree), then check if valid
                    lrecon = self._evolve_subtree(topchild, bottoms, state)
                    bottom_loci, _ = self._find_unique_loci(lrecon, bottoms)
                    if bottom_loci != bottom_loci_from_locus_map:
                        continue
                else:
                    if max_loci != INF:
                        # find locus at bottoms (in this subtree), then check if valid
                        lrecon = self._evolve_subtree(topchild, bottoms, state)
                        bottom_loci = [lrecon[node.name] for node in bottoms]
                        if len(set(bottom_loci)) > max_loci:
                            continue

                # store original (without duplication along root branch)
                states.append(state)

                # allow duplication along root branch
                if top is not topchild:
                    ndup = util.counteq(True, state.values())
                    if ndup < max_dups_subtree and topchild.name not in constraints:
                        s1 = state.copy()
                        s1[topchild.name] = True
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
        gene2locus = self.gene2locus

        # locus constraints
        if gene2locus is not None:
            constraints = self._find_constraints_nodups(sorted_bottoms)
        else:
            constraints = None

        # log heuristics
        self.log.log("Max # loci per sbranch: %f" % self.max_loci)
        self.log.log("Max # dup per sbranch: %f" % self.max_dups)
        self.log.log("Max # loss per sbranch: %f" % self.max_losses)
        self.log.log()

        # tiles at each sbranch
        #     key1 = snode, key2 = (bottom_loci, top_loci)
        #     val = list of items (lrecon, order, ndup, nloss, ncoalspec, ncoaldup, cost, nsoln)
        #           single item   (lrecon, order, cost, nsoln)
        #     (val is list generally and a single item after choosing tiles)
        tiles = {}

        if self.prescreen:
            tiles_ps = collections.defaultdict(dict)

        # recur down species tree
        for snode in stree.preorder():
            self.log.start("Working on snode %s" % snode.name)
            parent_snode = snode.parent
            is_leaf = snode.is_leaf()

            # get subtrees in this sbranch
            subtrees_snode = subtrees[snode]

            # TODO: need to update if-statement to match current structure
            # nothing happened along the branch if there are no subtrees
            if not subtrees_snode: # len(subtrees_snode) == 0
                # handle root separately (initialize loci)
                if snode is stree.root:
                    tiles[snode] = {}
                    self._initialize_tiles_sroot(tiles[snode])
                    if self.prescreen:
                        self._initialize_prescreen_sroot(tiles_ps[snode])
                self.log.log("Empty sbranch")
                self.log.stop()
                continue

            # get bottom nodes for this sbranch
            bottoms = sorted_bottoms[snode]

            # get states (changed/unchanged) for each branch of each subtree in sbranch
            states = self._find_locus_states_sbranch(subtrees_snode, is_leaf,
                                                     max_loci=self.max_loci,
                                                     max_dups=INF if is_leaf else self.max_dups,
                                                     constraints=constraints)

            # top of this sbranch is the bottom of the parent sbranch
            if parent_snode in tiles:
                top_loci_lst = set()
                for (top_loci, _) in tiles[parent_snode]:
                    top_loci_lst.add(top_loci)
                tops = sorted_bottoms[parent_snode]
                assert top_loci_lst
            else:
                top_loci_lst = set([(INIT_LOCUS,)])
                tops = [gtree.root]

            # TODO: incorporate max_dup from sbranch to sbranch
            self.log.log("top nodes: %s" % ','.join([node.name for node in tops]))
            self.log.log("bottom nodes: %s" % ','.join([node.name for node in bottoms]))
##            self.log.log("top_loci: %s" % ';'.join(map(str, top_loci_lst)))
            self.log.log("number of assignments at top nodes: %d" % len(top_loci_lst))
            states_len = map(len, states)
            self.log.log("number of states: %s" % ','.join(map(str, states_len)))
            self.log.log("number of state combinations: %d" % prod(states_len))
            self.log.log("")

            # combine top loci and states in sbranch
            tiles[snode] = self._enumerate_locus_maps_sbranch_helper(snode, subtrees_snode, is_leaf,
                                                                     tops, bottoms,
                                                                     top_loci_lst, states)

            # do any valid locus maps exist?
            if not tiles[snode]:
                raise Exception("no valid locus maps exist")

            # for each locus assignment at top and bottom of sbranch, filter storage
            self._choose_tile(tiles[snode])

            # prescreen
            if self.prescreen:
                self._prescreen(tiles[snode], tiles_ps[snode], tiles_ps.get(parent_snode, None))

            self.log.stop()
        self.log.stop()

        return tiles


    def _enumerate_locus_maps_sbranch_helper(self, snode, subtrees, is_leaf,
                                             tops, bottoms,
                                             top_loci_lst, states):
        """Enumerate the valid locus maps from set of top_loci and states."""

        tiles = {}

        # hash subtrees using top node
        subtrees_hash = {}
        for (stop, stopchild, sbottoms) in subtrees:
            subtrees_hash[stop] = (stopchild, sbottoms)

        # get bottom loci for this sbranch
        if is_leaf:
            gene2locus = self.gene2locus
            if gene2locus is not None:
                lrecon = dict([(bottom.name, gene2locus(bottom.name)) for bottom in bottoms])
                bottom_loci_from_locus_map, _ = self._find_unique_loci(lrecon, bottoms)

        for ndx1, top_loci in enumerate(top_loci_lst):
            # start with loci at top of sbranch
            assert len(top_loci) == len(states), (len(top_loci), len(states))

            # initialize next locus with top of sbranch
            init_next_locus = max(top_loci) + 1
            init_lrecon = {}
            for i, start in enumerate(top_loci):
                init_lrecon[tops[i].name] = start

            # combine states across subtrees
            for ndx2, state in enumerate(itertools.product(*states)):
                next_locus = init_next_locus
                lrecon = init_lrecon.copy()

                #=============================
                # find next lrecon using states

                for i, start in enumerate(top_loci):
                    s = state[i]
                    if s is not None:
                        stop = tops[i]
                        stopchild, sbottoms = subtrees_hash[stop]

                        # evolve from top to topchild
                        if s[stopchild.name]:
                            stopchild_loci = next_locus
                            next_locus += 1
                        else:
                            stopchild_loci = start

                        # evolve from topchild to bottoms
                        local_lrecon = self._evolve_subtree(stopchild, sbottoms, state=s,
                                                            start_locus=stopchild_loci,
                                                            next_locus=next_locus)
                        assert all([name not in lrecon for name in local_lrecon
                                    if name != stop.name]), (local_lrecon, lrecon)
                        lrecon.update(local_lrecon)
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
                    if len(set(bottom_loci)) > self.max_loci:
                        # exceeds max # lineages for (ancestral) sbranch
                        self.log.log("\t%s -> %s : [%d, %d] (skipped - locus count)" % \
                                        (top_loci, bottom_loci, ndx1, ndx2))
                        self.log.log("\t\tlrecon: %s" % lrecon)
                        self.log.log()
                        continue

                #=============================
                # find optimal cost (and order) for this lrecon
                solns = self._find_optimal_solns(tiles, bottom_loci, top_loci,
                                                 lrecon, subtrees=subtrees, bottoms=bottoms,
                                                 max_dups=INF if is_leaf else self.max_dups,
                                                 max_losses=INF if is_leaf else self.max_losses,
                                                 snode=snode)
                # all solns should have same cost, use first as proxy for logging info
                ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events = solns[0]

                # skip?
                skip = ""
                if not is_leaf:
                    if ndup > self.max_dups:       # exceeds max # of dup for (ancestral) sbranch
                        skip = "dup"
                    elif nloss > self.max_losses:  # exceeds max # of loss for (ancestral) sbranch
                        skip = "loss"
                if (skip == "") and ((ndup == INF) or (nloss == INF) \
                    or (ncoal_spec == INF) or (ncoal_dup == INF)):  # non-optimal
                    skip = "non-optimal"

                if skip != "":
                    self.log.log("\t%s -> %s : [%d, %d] (skipped - %s)" % \
                                    (top_loci, bottom_loci, ndx1, ndx2, skip))
                    self.log.log("\t\tlrecon: %s" % lrecon)
                    self.log.log()
                    continue

                # log
                self.log.log("\t%s -> %s : [%d, %d]" % \
                                (top_loci, bottom_loci, ndx1, ndx2))
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
                    self._update_tiles(tiles, bottom_loci, top_loci,
                                       lrecon, order,
                                       ndup, nloss, ncoal_spec, ncoal_dup, nsoln, events)

        return tiles


    def _infer_opt_locus_map(self, locus_maps, subtrees):
        """Find optimal locus map

        First recur up the species tree (through dynamic programming) to determine optimal cost,
        then trace back down the species tree to determine optimal locus map.
        """

        self.log.start("Inferring optimal locus map")

        # compute DP table
        dp_table = self._dp_compute(locus_maps)

        # terminate
        self._dp_terminate(dp_table)

        # traceback
        self._dp_traceback(locus_maps, subtrees, dp_table)

        self.log.stop()


    #=============================
    # locus tiles methods (used by _enumerate_locus_maps)
    # tiles structure:
    #     key  = (bottom_loci, top_loci)
    #     value = (lrecon, order, ndup, nloss, ncoalspec, ncoaldup, cost, nsoln) => (lrecon, order, cost, nsoln)

    def _initialize_tiles_sroot(self, tiles):
        """Initialize tiles for species root"""
        lrecon = {self.gtree.root.name: INIT_LOCUS}
        order = {}
        cost = 0
        nsoln = 1

        top_loci = bottom_loci = (INIT_LOCUS,)
        key = (bottom_loci, top_loci)
        tiles[key] = (lrecon, order, cost, nsoln)


    def _find_optimal_solns(self, tiles, bottom_loci, top_loci,
                            lrecon, subtrees, bottoms=None,
                            max_dups=INF, max_losses=INF, snode=None):
        """Find solutions for tiles"""

        mincost = INF
        key = (bottom_loci, top_loci)
        if key in tiles:
            # lst contains items (lrecon, order, ndup, nloss, ncoalspec, ncoaldup, cost, nsoln)
            mincost = min([item[6] for item in tiles[key]])

        solns = self._count_events(lrecon, subtrees, all_bottoms=bottoms,
                                   max_dups=max_dups, max_losses=max_losses,
                                   min_cost=mincost)

        return solns


    def _update_tiles(self, tiles, bottom_loci, top_loci,
                      lrecon, order,
                      ndup, nloss, ncoalspec, ncoaldup, nsoln, events):
        """Update tiles to keep with minimum reconciliation cost"""

        # a solution is better if 1) it has lower cost, or 2) it has equal cost and lower ndups

        mincost, mindup = INF, INF
        key = (bottom_loci, top_loci)

        if key in tiles:
            # lst contains items (lrecon, order, ndup, nloss, ncoalspec, ncoaldup, cost, nsoln)
            mincost, mindup = min([(item[6], item[2]) for item in tiles[key]])

        cost = self._compute_cost(ndup, nloss, ncoalspec, ncoaldup)
        item = (lrecon, order, ndup, nloss, ncoalspec, ncoaldup, cost, nsoln)

        if key not in tiles:
            tiles[key] = []

        if cost < mincost:
            tiles[key] = [item]
        elif cost == mincost:
            if ndup < mindup:
                tiles[key] = [item]
            elif ndup == mindup:
                assert item not in tiles[key]
                tiles[key].append(item)


    def _choose_tile(self, tiles):
        """Choose single tile with minimum reconciliation cost"""
        # Updates value of dictionary to be val = (lrecon, order, cost, nsoln)

        self.log.log("optimal costs")

        for (bottom_loci, top_loci), lst in tiles.iteritems():
            # lst contains items (lrecon, order, ndup, nloss, ncoalspec, ncoaldup, cost, nsoln)
            # where nsoln is the number of partial orderings that are optima for the lrecon

            # sample tile with weights proportional to number of optimal partial orderings
            nsoln = [item[7] for item in lst]
            total_nsoln = sum(nsoln)
            item = common.random_choice(lst, p=nsoln, normalize=True)

            # set the chosen optimum back into tiles
            # update number of solutions to be total across all optimal locus maps
            lrecon, order = item[:2]
            cost = item[-2]
            item = (lrecon, order, cost, total_nsoln)

            tiles[bottom_loci, top_loci] = item

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

    def _initialize_prescreen_sroot(self, prescreen):
        """Initialize prescreen dictionary"""
        bottom_loci = (INIT_LOCUS,)
        prescreen[bottom_loci] = 0


    def _prescreen(self, tiles, tiles_ps, tiles_ps_parent):
        """Prescreen"""
        self.log.log("prescreen")

        # update min cost-to-go (from root)
        for (bottom_loci, top_loci), item in tiles.iteritems():
            cost = item[2] # (lrecon, order, cost, nsoln)

            if tiles_ps_parent is not None:
                parent_cost = tiles_ps_parent[top_loci]
            else:
                parent_cost = 0
            total_cost = parent_cost + cost

            if bottom_loci in tiles_ps:
                if tiles_ps[bottom_loci] > total_cost:
                    tiles_ps[bottom_loci] = total_cost
            else:
                tiles_ps[bottom_loci] = total_cost

        # determine min cost across all bottom_loci
        mincost = min(tiles_ps.values())
        self.log.log("\tmin cost: %d" % mincost)

        if mincost <= self.prescreen_min:
            self.log.log("\tskipping prescreen")
        else:
            ntotal = 0
            npruned = 0
            thr = self.prescreen_factor * mincost
            for bottom_loci, cost in tiles_ps.iteritems():
                ntotal += 1
                if cost > thr:    # prescreen bottom_loci from this sbranch
                    npruned += 1
                    del tiles[bottom_loci]
                    self.log.log("\tpruned %s : %d" % (bottom_loci, cost))
            self.log.log("\tpruned %d/%d assignments" % (npruned, ntotal))
        self.log.log()


    #=============================
    # DP table methods (used by _infer_opt_locus_map)

    def _dp_compute(self, locus_maps):
        """Compute DP table"""
        # locus_maps is a multi-dimensional dict with the structure
        # key1 = snode, key2 = (bottom_loci, top_loci), value = (lrecon, order, cost, nsoln)

        stree = self.stree
        gene2locus = self.gene2locus

        # dynamic programming storage
        dp_table = {}   # key1 = snode, key2 = top_loci
                        # val = (bottom_loci, cost_to_go, nsoln_to_go)

        for snode in stree.postorder():
            self.log.start("Working on snode %s" % snode.name)
            dp_table[snode] = {}

            if snode not in locus_maps:
                # nothing in this sbranch
                dp_table[snode][()] = ((), 0, 1)
                self.log.log("Empty sbranch")
                self.log.stop()
                continue

            # get stored values for the sbranch
            locus_maps_snode = locus_maps[snode]

            if snode.is_leaf():
                # leaf base case

                for (bottom_loci, top_loci), (_, _, cost, nsoln) in locus_maps_snode.iteritems():
                    assert top_loci not in dp_table[snode]
                    dp_table[snode][top_loci] = (bottom_loci, cost, nsoln)
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
                costs = collections.defaultdict(list) # key = top_loci, val = (bottom_loci, cost, nsoln)

                for (bottom_loci, top_loci), (_, _, cost, nsoln) in locus_maps_snode.iteritems():
                    # find cost-to-go and nsoln-to-go in children
                    # locus assignment may have been removed due to search heuristics
                    _, cost_left, nsoln_left = dp_table[sleft].get(bottom_loci, (None, INF, 0))
                    _, cost_right, nsoln_right = dp_table[sright].get(bottom_loci, (None, INF, 0))

                    # update cost-to-go and nsoln-to-go based on children
                    children_cost = cost_left + cost_right
                    children_nsoln = nsoln_left * nsoln_right

                    # add cost and multiple nsoln in this sbranch
                    cost_to_go = cost + children_cost
                    nsoln_to_go = nsoln * children_nsoln
                    item = (bottom_loci, cost_to_go, nsoln_to_go)
                    assert item not in costs[top_loci], (snode, bottom_loci, top_loci, item)
                    costs[top_loci].append(item)

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
                        # multiple samples per species - check validity
                        # locus maps per branch only checked once all locus maps are put together
                        # this top_loci is invalid because it never satisfied dup constraints
                        if mincost == INF:
                            continue

                    # item = (bottom_loci, cost_to_go, nsoln_to_go)
                    # sample relative locus maps with weights proportional to number of solutions
                    nsoln = [item[2] for item in items]
                    total_nsoln = sum(nsoln)
                    item = common.random_choice(items, p=nsoln, normalize=True)

                    # set the chosen optimum
                    # update number of solutions to be total across all optimal choices
                    bottom_loci, cost_to_go, nsoln_to_go = item
                    item = (bottom_loci, cost_to_go, total_nsoln)
                    dp_table[snode][top_loci] = item

            self.log.log("DP table")
            for top_loci, (bottom_loci, cost_to_go, nsoln_to_go) in dp_table[snode].iteritems():
                self.log.log("%s -> %s : %g [%g]" % (top_loci, bottom_loci, cost_to_go, nsoln_to_go))
            self.log.stop()

        return dp_table


    def _dp_terminate(self, dp_table):
        """Terminate DP to find optimal cost"""
        stree = self.stree

        # not necessary since cost along root sbranch already determined,
        # and by design, dp_table[sroot] is always assigned locus = INIT_LOCUS
        assert len(dp_table[stree.root]) == 1, dp_table[stree.root]
        _, cost_to_go, nsoln_to_go = dp_table[stree.root].values()[0]
        self.log.log("")
        self.log.log("Optimal cost: %g" % cost_to_go)
        self.log.log("Number of solutions: %g" % nsoln_to_go)
        self.log.log("")


    def _dp_traceback(self, locus_maps, subtrees, dp_table):
        """Traceback DP to find optimal reconciliation"""

        gtree = self.gtree
        stree = self.stree

        # recur down species tree to assign loci
        next_locus = 1
        lrecon = {}
        order = collections.defaultdict(dict)
        tb = {}    # used in traceback, key = snode, val = bottom_loci of sbranch

        # initialize root
        lrecon[gtree.root] = next_locus
        next_locus += 1

        for snode in stree.preorder():
            self.log.start("Working on snode %s" % snode.name)

            # determine top_loci and bottom_loci
            if snode is stree.root:
                # root base case
                top_loci = dp_table[snode].keys()[0]
            else:
                top_loci = tb[snode.parent]
            bottom_loci = dp_table[snode][top_loci][0]

            # update traceback
            tb[snode] = bottom_loci

            # find loci and order
            if top_loci == ():
                pass
            else:
                # get stored values for the sbranch
                locus_maps_snode = locus_maps[snode] # (local_lrecon, local_order, cost, nsoln)
                subtrees_snode = subtrees[snode]

                # get optimum for sbranch from DP table
                self.log.log("%s -> %s : %g [%g]" % \
                             (top_loci, bottom_loci, \
                              dp_table[snode][top_loci][1], dp_table[snode][top_loci][2]))
                local_lrecon, local_order, _, _ = locus_maps_snode[bottom_loci, top_loci]

                self.log.log("lrecon: %s" % local_lrecon)
                self.log.log("order: %s" % local_order)

                # update lrecon
                for (_, topchild, bottoms) in subtrees_snode:
                    if not topchild:
                        continue

                    for node in gtree.preorder(topchild, is_leaf=lambda x: x in bottoms):
                        pnode = node.parent
                        if pnode:
                            if local_lrecon[node.name] != local_lrecon[pnode.name]:
                                lrecon[node] = next_locus
                                next_locus += 1
                            else:
                                lrecon[node] = lrecon[pnode]

                # update order
                for _, lst in local_order.iteritems():
                    new_plocus = lrecon[lst[0].parent]
                    order[snode][new_plocus] = lst

            self.log.stop()

        self.lrecon = lrecon
        self.order = dict(order)
        self.cost = dp_table[stree.root].values()[0][1]
        self.nsoln = dp_table[stree.root].values()[0][2]
