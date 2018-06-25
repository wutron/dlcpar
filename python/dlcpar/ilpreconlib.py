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

INIT_LOCUS = 1

#=============================================================================
# reconciliation data structures

class IlpReconVariables(object):
    """
    The reconciliation data structure for the DLC ILP formulation

    input variables:
        gtree
        stree
        srecon
    recon variables:
        dup_vars           :    key = gnode, value = 1 if duplication on edge between gnode and gnode's parent, 0 otherwise
        order_vars
    solver variables:
    structure variables:
        _gnodes_by_species :    TBD
        _orders_from_topology : TBD
        _bottom_nodes :         key = snode, value = list of gnodes mapped to bottom of snode
        _top_nodes_with_child : key = snode, value = list of gnodes mapped to top of snode with child also mapped to snode
        _top_nodes :            key = snode, value = list of gnodes mapped to top of snode
    """

    def __init__(self, gtree, stree, srecon, all_vars=True):
        self.gtree = gtree
        self.stree = stree
        self.srecon = srecon

        self._create_recon_vars()
        if all_vars:
            self._create_solver_vars()


    def _create_recon_vars(self):
        """Creates variables necessary for converting to LCT and inferring evolutionary events"""

        gtree = self.gtree
        srecon = self.srecon

        #========================================
        # duplication variables

        # d_e, key = gnode, value = 1 if duplication on edge between gnode and gnode's parent, 0 otherwise
        self.dup_vars = pulp.LpVariable.dicts("dup", list(gtree.preorder()), 0, 1, pulp.LpInteger)

        #========================================
        # order variables

        # infer required orders from gene tree topology
        # key = (g1,g2) such that g1,g2 in same species and one node is ancestor of other
        #       note that only one of (g1,g2) or (g2,g1) is in dictionary
        # value = 1 if g2 more recent than g1 and 0 otherwise
        def descendants_in_species(gnode):
            """Return list of descendants of gnode in same species as gnode"""
            children_in_species = [child for child in gnode.children if srecon[gnode] == srecon[child]]
            indirect_descendants = sum([descendants_in_species(child) for child in children_in_species], [])
            return children_in_species + indirect_descendants

        self._orders_from_topology = {}
        for gnode in gtree:
            for descendant in descendants_in_species(gnode):
                self._orders_from_topology[gnode, descendant] = 1

        # o_{g1,g2}, key = (g1,g2), value = 1 if g2 more recent than g1, 0 otherwise
        self._gnodes_by_species = collections.defaultdict(list)
        for gnode in gtree:
            self._gnodes_by_species[srecon[gnode]].append(gnode)
        pairs_in_species = []        
        for gnodes in self._gnodes_by_species.itervalues():
            pairs = list(pulp.combination(gnodes, 2))
            pairs_in_species.extend(pairs)
        
        order_keys = [(g1, g2) for (g1, g2) in pairs_in_species
                      if (g1, g2) not in self._orders_from_topology and (g2, g1) not in self._orders_from_topology]
        self.order_vars = pulp.LpVariable.dicts("order", order_keys, 0, 1, pulp.LpInteger)


    def _create_solver_vars(self):
        """Create variables for solving ILP"""

        all_gnodes = list(self.gtree.preorder())
        sevents = phylo.label_events(self.gtree, self.srecon)
        subtrees = reconlib.factor_tree(self.gtree, self.stree, self.srecon, sevents)

        #========================================
        # structures containing groups of nodes

        self._bottom_nodes = collections.defaultdict(list)
        self._top_nodes_with_child = collections.defaultdict(list)
        self._top_nodes = collections.defaultdict(list)

        # find top and bottom nodes for each snode
        for snode, subtrees_snode in subtrees.iteritems():
            for (root, rootchild, leaves) in subtrees_snode:
                self._top_nodes[snode].append(root)
                if rootchild:
                    self._top_nodes_with_child[snode].append(root)
                if leaves:
                    self._bottom_nodes[snode].extend(leaves)                

        #========================================
        # loss variables

        loss_keys = []
        for snode, tops in self._top_nodes.iteritems(): # contains every gene mapped to top of snode 
            for top in tops:
                loss_keys.append((snode, top))

        # l_sg
        # key = (snode, gnode at top of snode)
        # value = 1 if gnode creates a loss in snode, 0 otherwise
        self._loss_vars = pulp.LpVariable.dicts("loss", loss_keys, 0, 1, pulp.LpInteger) 

        #========================================
        # coalspec variables
    
        coalspec_keys = []
        for snode, tops in self._top_nodes_with_child.iteritems():
            for top in tops:
                coalspec_keys.append((snode, top))

        # c_sg variables
        # key = (snode, gene node at top of snode with child)
        # value = 1 if gnode on same locus as another gnode at top of snode with child, 0 otherwise
        self._coalspec_vars = pulp.LpVariable.dicts("coal_spec", coalspec_keys, 0, 1, pulp.LpInteger) 

        #========================================
        # constraint variables

        # h_g
        # key = (snode, gnode at top of snode)
        # value = 1 if g on same locus as gnode2 mapped to top of snode and gnode2 < gnode, 0 otherwise
        self._helper_vars = pulp.LpVariable.dicts("helper", self._loss_vars.keys(), 0, 1, pulp.LpInteger) 

        # p_{g,g'}
        # key = pair of gene nodes
        # value = 1 if path between gene nodes has at least one duplication, 0 otherwise
        all_pairs = list(pulp.combination(all_gnodes, 2))
        self._path_vars = pulp.LpVariable.dicts("path", all_pairs, 0, 1, pulp.LpInteger) 

        # delta_{g1,g2} variables, see paper for description
        delta_vars_keys = []
        for gnodes in self._gnodes_by_species.itervalues():

            for g1, g2 in pulp.combination(gnodes, 2):
                if ((g1,g2) not in self._orders_from_topology) and ((g2,g1) not in self._orders_from_topology):
                    if g1.parent and g2.parent:
                        delta_vars_keys.append((g1,g2))
                        delta_vars_keys.append((g2,g1))
        self._delta_vars = pulp.LpVariable.dicts("coal_dup_helper", delta_vars_keys, 0, 1, pulp.LpInteger) 

        # r_g
        # key = gnode
        # value = number of coalescense due to dup at g
        self._coaldup_vars = pulp.LpVariable.dicts("coal_dup", all_gnodes, 0, None, pulp.LpInteger)


    def get_order(self, g1, g2):
        """Return 1 if g2 more recent than g1, 0 otherwise.

        Checks order_vars, then orders_from_topology, then ???"""
        if (g1, g2) in self.order_vars:
            return self.order_vars[(g1, g2)]
        elif (g2, g1) in self.order_vars:
            return 1 - self.order_vars[(g2, g1)]

        elif (g1, g2) in self._orders_from_topology:
            return self._orders_from_topology[(g1, g2)]
        elif (g2, g1) in self._orders_from_topology:
            return 1 - self._orders_from_topology[(g2, g1)]

        elif self.srecon[g2] in self.srecon[g1].children:
            # g2 maps to a species node that is a child of the species node g1 maps to
            # needed for delta_vars checking g1.parent against g2
            return 1

        else:
            raise Exception("Could not find order for nodes (%s,%s)" % (g1,g2))


    def get_path(self, g1, g2):
        """Return 1 if there is a duplication on path between g1 and g2, 0 otherwise."""

        if g1 == g2:
            assert (g1,g2) not in self._path_vars
            return 0
        elif (g1, g2) in self._path_vars:
            return self._path_vars[(g1, g2)]
        elif (g2,g1) in self._path_vars:
            return self._path_vars[(g2, g1)]
        else:
            raise Exception("Could not find path variable for nodes (%s,%s)" % (g1,g2))


#=============================================================================
# conversion utilities

def ilp_to_lct(gtree, lpvars):
    """Converts ILPReconVariables to a LabeledRecon"""

    #========================================
    # find species map
    srecon = lpvars.srecon

    #========================================
    # find locus map
    # adapted from recon._evolve_subtree

    dup_nodes = [g for g, dup_var in lpvars.dup_vars.iteritems() if dup_var.varValue == 1.0]

    locus = INIT_LOCUS
    lrecon = {}
    
    for gnode in gtree.preorder():
        if gnode == gtree.root:  # if root, it has the first locus
            lrecon[gnode] = locus
        elif gnode in dup_nodes: # if it has a dup, it has a different locus
            locus += 1
            lrecon[gnode] = locus
        else:                    # if there's no dup, it's the same as the parent
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

    for gnode in gtree:
        # skip gene tree root
        if gnode.parent is None:
            continue

        snode = srecon[gnode]
        locus = lrecon[gnode]
        plocus = lrecon[gnode.parent]

        if (snode, plocus) in parent_loci:
            # skip if same locus as parent and leaf node
            if locus == plocus and (gnode.is_leaf() or \
                (len(gnode.children) == 1 and all([snode != srecon[child] for child in gnode.children]))):
                continue

            order.setdefault(snode, {})
            order[snode].setdefault(plocus, [])
            order[snode][plocus].append(gnode)

  
    for snode, d in order.iteritems():
        for plocus, lst in d.iteritems():
            # "bubble sort" the list
            for i in range(len(lst)-1):
                for j in range(len(lst)-1):
                    g1, g2 = lst[j], lst[j+1]
                    if lpvars.get_order(g2, g1) ==1:
                        # swap consecutive genes
                        lst[j], lst[j+1] = lst[j+1], lst[j]

    #========================================
    # put everything together
    return reconlib.LabeledRecon(srecon, lrecon, order)


def lct_to_ilp(gtree, stree, labeledrecon):
    """ Converts LabeledRecon into ILPReconVariables"""

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
                    order_vars[g1, g2] = 1 # g1 is older than g2

    #========================================
    # put everything together
    lpvars = IlpReconVariables(gtree, stree, labeledrecon.species_map, all_vars=False)
    for gnode, var_value in dup_vars.iteritems():
        lpvars.dup_vars[gnode].varValue = var_value
    for gnodes, var_value in order_vars.iteritems():
        if gnodes in lpvars.order_vars:
            lpvars.order_vars[gnodes].varValue = var_value
    return lpvars

