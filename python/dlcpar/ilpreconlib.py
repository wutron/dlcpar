# integer linear programming library

# python libraries
import collections

# pulp libraries
import pulp

# rasmus, compbio libraries
from rasmus import util
from compbio import phylo

# dlcpar libraries
from dlcpar import reconlib

#==========================================================
# globals

INIT_LOCUS = 1

#=============================================================================
# reconciliation data structures

class IlpReconVariables(object):
    """
    The reconciliation data structure for the DLC ILP formulation

    input variables :
        gtree (Tree) :
            gene tree
        stree (Tree) :
            species tree
        srecon (dict) :
            species map

    LCT variables (dicts - all keys are node names):
        dup_vars :
            key = g
            value = 1 if duplication on edge to g
                    0 otherwise
        order_vars :
            key = (g1, g2)
            value = 1 if g1 precedes g2 in time
                    0 otherwise

    solver variables (dicts - all keys are node names):
        _omega_vars :
            key = (g1 at bottom of snode, g2 of s)
            value = 1 if no duplication on edge to g1 and duplication on edge to g2
                    0 otherwise
        _loss_vars :
            key = (s, g at top of s)
            value = 1 if g creates a loss in s
                    0 otherwise
        _path_vars :
            key = (g1, g2)
            value = 1 if there is at least one dup on path between g1 and g2
                    0 otherwise
        _lambda_vars :
            key = (s, g at top of snode)
            value = 1 if g on same locus as some other gprime that precedes g
                    0 otherwise
        _coalspec_vars :
            key = (s, g at top of snode with child in s)
            value = 1 if child of g is extra lineage
                    0 otherwise
        _kappa_vars :
            key = (g1 in nodes or bottom children of s, g2 of s that is incomparable)
            value = 1 if duplication on edge to g2 and g1 is contemporary lineage
                    0 otherwise
        _coaldup_vars :
            key = g
            value = number of coalesences due to a duplication at g

    structure variables:
        _gnodes (dict) :
            key = sname, value = list of gnodes mapped to snode
        _top_nodes (dict) :
            key = sname, value = list of gnodes at top of snode
        _bottom_nodes (dict) :
            key = sname, value = list of gnodes at bottom of snode
        _survived_nodes (dict) :
            key = sname, value = list of gnodes at top of snode with child also mapped to snode
        _bottom_children_nodes (dict) :
            key = sname, value = list of childrens of bottoms of the snode
        _comparable_nodes (list) :
            elements are tuples (gname1, gname2) where g1 is ancestor of g2
        _orders_from_tree (list) :
            elements are tuples (gname1, gname2), where g1 is older than g2 by (gene tree or species tree) topology
    """


    def __init__(self, gtree, stree, srecon, all_vars=True):
        self.gtree = gtree
        self.stree = stree
        self.srecon = srecon

        self._create_structure_vars(all_vars)
        self._create_recon_vars(all_vars)
        if all_vars:
            self._create_solver_vars()


    def _create_structure_vars(self, all_vars):
        """Create variables for capturing structure in trees"""

        gtree = self.gtree
        stree = self.stree
        srecon = self.srecon

        #========================================
        # node sets

        self._gnodes = collections.defaultdict(list)
        for gnode in gtree:
            sname = srecon[gnode].name
            self._gnodes[sname].append(gnode)

        if all_vars:
            sevents = phylo.label_events(gtree, srecon)
            subtrees = reconlib.factor_tree(gtree, stree, srecon, sevents)

            self._top_nodes = collections.defaultdict(list)
            self._bottom_nodes = collections.defaultdict(list)
            self._survived_nodes = collections.defaultdict(list)
            self._bottom_children_nodes = collections.defaultdict(list)

            for snode, subtrees_snode in subtrees.iteritems():
                sname = snode.name
                for (top, topchild, bottoms) in subtrees_snode:
                    self._top_nodes[sname].append(top)
                    if topchild:
                        self._survived_nodes[sname].append(top)
                    if bottoms:
                        self._bottom_nodes[sname].extend(bottoms)
                        for gnode in bottoms:
                            self._bottom_children_nodes[sname].extend(gnode.children)

        #========================================
        # helper functions

        def descendants_in_species(gnode):
            """Return list of descendants of gnode in same species as gnode"""

            children_in_species = [child for child in gnode.children if srecon[gnode] == srecon[child]]
            bottom_childs = []
            for child in gnode.children:
                if srecon[child] != srecon[gnode]:
                    bottom_childs.append(child)
            indirect_descendants = sum([descendants_in_species(child) for child in children_in_species], [])
            return children_in_species + indirect_descendants + bottom_childs

        def descendants_not_in_species(snode):
            """Return list of descendants of a particular gnode from its species tree"""

            children = []
            for schild in snode.children:
                for g2 in self._gnodes[schild.name]:
                    children.append(g2)
            return children

        #========================================
        # orderings from topology

        self._comparable_nodes = []
        self._orders_from_tree = []

        for gnode in gtree.preorder():
            for descendant in descendants_in_species(gnode):
                self._comparable_nodes.append((gnode.name, descendant.name))

        self._orders_from_tree = self._comparable_nodes[:] # copy
        for gnode in gtree.preorder():
            snode = srecon[gnode]
            for descendant in descendants_not_in_species(snode):
                key = (gnode.name, descendant.name)
                if key not in self._orders_from_tree:
                    self._orders_from_tree.append(key)


    def _create_recon_vars(self, all_vars):
        """Creates variables necessary for converting to LCT and inferring evolutionary events"""

        gtree = self.gtree
        stree = self.stree
        srecon = self.srecon

        #========================================
        # duplication variables

        # d_g
        dup_keys = [node.name for node in gtree]
        self.dup_vars = pulp.LpVariable.dicts("dup", dup_keys, 0, 1, pulp.LpBinary)

        #========================================
        # order variables

        # o_{g1,g2}
        # create variables not already ordered by tree
        order_keys = []
        for (gname1, gname2) in self._orders_from_tree:
            order_keys.append((gname1, gname2))
            order_keys.append((gname2, gname1))

        for snode in stree:
            sname = snode.name
            for g1 in self._gnodes[sname]:
                for g2 in self._gnodes[sname]:
                    if (g1 != g2) and ((g1.name, g2.name) not in order_keys):
                        order_keys.append((g1.name, g2.name)) # combinations not permutations

        self.order_vars = pulp.LpVariable.dicts("order", order_keys, 0, 1, pulp.LpBinary)


    def _create_solver_vars(self):
        """Create variables for solving ILP"""

        gtree = self.gtree
        stree = self.stree
        all_gnodes = list(gtree.preorder())

        #========================================
        # order variables
        # technically required to LCT constraints but only needed when solving ILP

        # omega_{g1,g2}
        omega_keys = []
        for sname in stree.leaf_names():
            for g1 in self._bottom_nodes[sname]:
                for g2 in self._gnodes[sname]:
                    if g1 != g2:
                        omega_keys.append((g1.name, g2.name))
        self._omega_vars = pulp.LpVariable.dicts("omega", omega_keys, 0, 1, pulp.LpBinary)

        #========================================
        # loss variables

        # l_sg
        loss_keys = []
        for sname, tops in self._top_nodes.iteritems():
            for top in tops:
                loss_keys.append((sname, top.name))
        self._loss_vars = pulp.LpVariable.dicts("loss", loss_keys, 0, 1, pulp.LpBinary)

        # p_{g,g'}
        path_keys = []
        for g1, g2 in pulp.permutation(all_gnodes, 2):
            path_keys.append((g1.name, g2.name))
        for g in all_gnodes:
            path_keys.append((g.name,g.name)) # trivial paths
        self._path_vars = pulp.LpVariable.dicts("path", path_keys, 0, 1, pulp.LpBinary)

        # lambda_sg
        self._lambda_vars = pulp.LpVariable.dicts("lambda", loss_keys, 0, 1, pulp.LpBinary)


        #========================================
        # coalspec variables

        # c_sg
        coalspec_keys = []
        for sname, tops in self._survived_nodes.iteritems():
            for top in tops:
                coalspec_keys.append((sname, top.name))
        self._coalspec_vars = pulp.LpVariable.dicts("coal_spec", coalspec_keys, 0, 1, pulp.LpBinary)


        #========================================
        # coaldup variables

        # kappa_{g1,g2}
        kappa_vars_keys = []
        for snode in stree:
            sname = snode.name
            for g1 in self._gnodes[sname] + self._bottom_children_nodes[sname]:
                for g2 in self._gnodes[sname]:
                    gname1 = g1.name
                    gname2 = g2.name
                    if (gname1 != gname2) and \
                       ((gname1, gname2) not in self._comparable_nodes) and \
                       ((gname2, gname1) not in self._comparable_nodes):
                        if g1.parent and g2.parent:
                            kappa_vars_keys.append((gname1, gname2))
        self._kappa_vars = pulp.LpVariable.dicts("coal_dup_helper", kappa_vars_keys, 0, 1, pulp.LpBinary)

        # k_g
        self._coaldup_vars = pulp.LpVariable.dicts("coal_dup", list(node.name for node in all_gnodes), 0, None, pulp.LpInteger)


    #=============================================================================
    # rounding utilities for CPLEX_PY

    def _get_binary_dicts(self):
        # every dictionary except _coal_dup_vars, as it takes on integer values
        return [self.dup_vars, self.order_vars, self._omega_vars,
                self._loss_vars, self._path_vars, self._lambda_vars, \
                self._coalspec_vars, self._kappa_vars]


#=============================================================================
# conversion utilities

def ilp_to_lct(gtree, lpvars):
    """Convert from ILPReconVariables to LabeledRecon"""

    #========================================
    # find relevant ILP variables

    srecon = lpvars.srecon
    dup_vars = lpvars.dup_vars
    order_vars = lpvars.order_vars

    #========================================
    # find locus map
    # adapted from recondp._evolve_subtree

    dup_nodes = [gtree.nodes[gname] for gname, dup_var in dup_vars.iteritems() if dup_var.varValue == 1.0]

    locus = INIT_LOCUS
    lrecon = {}

    for gnode in gtree.preorder():
        if gnode == gtree.root:  # if root, it has the first locus
            lrecon[gnode] = locus
        elif gnode in dup_nodes: # if dup, it has a new locus
            locus += 1
            lrecon[gnode] = locus
        else:                    # if no dup, it has same locus as parent
            lrecon[gnode] = lrecon[gnode.parent]

    #========================================
    # find order

    order = {}

    # create (snode, plocus) pairs for which corresponding gnodes will be ordered
    parent_loci = set()
    for gnode in gtree:
        pnode = gnode.parent
        if pnode:
            locus = lrecon[gnode]
            plocus = lrecon[pnode]

            if locus != plocus:
                snode = srecon[gnode]
                parent_loci.add((snode, plocus))

    # find nodes for each species and locus
    for gnode in gtree.preorder():
        if gnode.parent:
            snode = srecon[gnode]
            locus = lrecon[gnode]
            plocus = lrecon[gnode.parent]

            if (snode, plocus) in parent_loci:
                # skip if same locus as parent and leaf node or special bottom node
                if locus == plocus and (gnode.is_leaf() or \
                    (len(gnode.children) == 1 and all([snode != srecon[child] for child in gnode.children]))):
                    continue

                order.setdefault(snode, {})
                order[snode].setdefault(plocus, [])
                order[snode][plocus].append(gnode)

    # find duplication order and partitions
    for snode, d in order.iteritems():
        for plocus, lst in d.iteritems():
            # "insertion sort" the list
            for i in xrange(1, len(lst)):
                g1 = lst[i]
                j = i-1
                while j >= 0 and order_vars[g1.name, lst[j].name].varValue == 1:
                    lst[j+1] = lst[j]
                    j -= 1
                lst[j+1] = g1

            # sanity check that all the order variables are satisfied by the order in lst (after the insertion sort)
            for g1, g2 in list(pulp.combination(lst, 2)):
                assert order_vars[g1.name, g2.name].varValue == 1.0, \
                    (((g1.name, g2.name), lst), "is not in the correct order \
                    with order value: ", order_vars[g1.name, g2.name].varValue)

    #========================================
    # put everything together

    labeled_recon = reconlib.LabeledRecon(srecon, lrecon, order)

    #========================================
    # check conversion

    stree = lpvars.stree.copy()

    loss_vars = lpvars._loss_vars
    coalspec_vars = lpvars._coalspec_vars
    coaldup_vars = lpvars._coaldup_vars

    loss_tuples = [(stree.nodes[sname], gtree.nodes[gname]) \
                   for (sname, gname), loss_var in loss_vars.iteritems() if loss_var.varValue == 1.0]
    coalspec_tuples = [(stree.nodes[sname], gtree.nodes[gname]) \
                       for (sname, gname), coalspec_var in coalspec_vars.iteritems() if coalspec_var.varValue == 1.0]
    coaldup_tuples = [(gtree.nodes[gname], coaldup_var) \
                      for gname, coaldup_var in lpvars._coaldup_vars.iteritems() if coaldup_var.varValue > 0]

    reconlib.init_dup_loss_coal_tree(stree)
    LCT_dups, LCT_losses, LCT_coalspecs, LCT_coaldups = _get_event_vars_from_LCT(gtree, stree, labeled_recon)

    assert set(dup_nodes) == set(LCT_dups) \
        and set(loss_tuples) == set(LCT_losses) \
        and set(coalspec_tuples) == set(LCT_coalspecs) \
        and set(coaldup_tuples) == set(LCT_coaldups)

    #========================================
    # finished

    return gtree, labeled_recon


def lct_to_ilp(gtree, stree, labeledrecon):
    """Convert from LabeledRecon to ILPReconVariables"""

    #========================================
    # find dup_vars

    lrecon = labeledrecon.locus_map
    dup_vars = {}
    for gnode in lrecon:
        if gnode.parent and lrecon[gnode] != lrecon[gnode.parent]:
            dup_vars[gnode] = 1.0
        else:
            dup_vars[gnode] = 0.0

    #========================================
    # find order_vars

    order = labeledrecon.order
    order_vars = {}
    for snode in order:
        for plocus in order[snode]:
            gene_lst = order[snode][plocus]
            for i, g1 in enumerate(gene_lst[:-1]):
                for g2 in gene_lst[i+1:]:
                    order_vars[g1, g2] = 1.0 # g1 is older than g2

    #========================================
    # put everything together

    lpvars = IlpReconVariables(gtree, stree, labeledrecon.species_map, all_vars=False)
    for gnode, var_value in dup_vars.iteritems():
        lpvars.dup_vars[gnode.name].varValue = var_value
    for gnodes, var_value in order_vars.iteritems():
        if gnodes in lpvars.order_vars:
            lpvars.order_vars[gnodes.name].varValue = var_value
    return lpvars


def _get_event_vars_from_LCT(gtree, stree, labeled_recon):
    """Returns event vars from LCT"""

    # adapted from reconlib.count_dup_loss_coal_tree

    # use stree to modify internal species map and order
    extra = labeled_recon.get_dict()
    srecon = util.mapdict(extra["species_map"], val=lambda snode: stree.nodes[snode.name])
    lrecon = extra["locus_map"]
    order = util.mapdict(extra["order"], key=lambda snode: stree.nodes[snode.name])

    extra = {"species_map": srecon,
             "locus_map":   lrecon,
             "order":       order}

    # factor gene tree
    events = phylo.label_events(gtree, srecon)
    subtrees = reconlib.factor_tree(gtree, stree, srecon, events)

    # get events along each species branch
    dups = []
    losses = []
    coalspecs = []
    coaldups = []

    for snode in stree:
        subtrees_snode = subtrees[snode]
        if len(subtrees_snode) == 0:
            continue

        # find dups
        dup_nodes = reconlib.find_dup_snode(gtree, stree, extra, snode,
                                     subtrees, subtrees_snode)
        dups.extend(dup_nodes)

        # find losses
        losses_snode = reconlib.find_loss_snode(gtree, stree, extra, snode,
                                               subtrees, subtrees_snode)
        if len(losses_snode) > 0:
            for loss in losses_snode:
                losses.append((snode, loss[0])) # use first node

        # find extra lineages due to coalescence at speciations
        coalspec_lineages = reconlib.find_coal_spec_snode(gtree, stree, extra, snode,
                                                          subtrees, subtrees_snode,
                                                          implied=True)
        for lineages in coalspec_lineages:
            assert len(lineages) > 1
            for i in range(1, len(lineages)):                 # skip first node
                coalspecs.append((snode, lineages[i].parent)) # use parents

        # find extra lineages due to coalescence at duplications
        coaldup_lineages = reconlib.find_coal_dup_snode(gtree, stree, extra, snode,
                                                        subtrees, subtrees_snode,
                                                        return_dups=True)
        for dup, lineages in coaldup_lineages.iteritems():
            assert len(lineages) > 1
            coaldups.append((dup, len(lineages)-1))

    return dups, losses, coalspecs, coaldups
