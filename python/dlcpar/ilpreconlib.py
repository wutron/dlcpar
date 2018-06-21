
# integer linear programming

# python libraries
import collections

# pulp libraries
import pulp 

# rasmus libraries
from rasmus import util

# dlcpar libraries
from dlcpar import reconlib

#=============================================================================
# reconciliation data structures

INIT_LOCUS = 1

class IlpReconVariables(object):
    """
    The reconciliation data structure for the DLC ILP formulation
    """

    def __init__(self, gtree, srecon, subtrees, all_vars=True):
        self.srecon = srecon
        self.all_vars = all_vars
        self._create_recon_vars(gtree, srecon)
        self._create_solver_vars(gtree, srecon, subtrees)

    def _create_recon_vars(self, gtree, srecon):
        """Creates variables necessary for inferring evolutionary events and making LCT"""
        #========================================
        # d_e variables in paper, key = gnode value = 1 if duplication on edge between gnode and gnode's parent, 0 otherwise
        self.dup_vars = pulp.LpVariable.dicts("dup", list(gtree.preorder()), 0, 1, pulp.LpInteger)

        #========================================
        # create order variables
        self._gnodes_by_species = collections.defaultdict(list)
        for gnode in gtree:
            self._gnodes_by_species[srecon[gnode]].append(gnode)
        pairs_for_each_species = [list(pulp.combination(gnodes, 2)) for gnodes in self._gnodes_by_species.values()]
        pairs_in_species= sum(pairs_for_each_species, [])

        def descendants_in_species(gnode, srecon):
            """returns list of descendants of gnode in same species as gnode"""
            children_in_species = [child for child in gnode.children if srecon[gnode] == srecon[child]]
            indirect_descendants = sum([descendants_in_species(child, srecon) for child in children_in_species], [])

            return children_in_species + indirect_descendants
        
        # infer required orders from gene tree topology
        # key = pair of gnodes in same species such that one gnode is an ancestor of other
        #       note that only one of (g1,g2) or (g2,g1) is in dictionary
        # value = 1 if g2 more recent than g1 and 0 otherwise
        self._orders_from_topology = {}
        def infer_orders(subtree_root):
            for descendant in descendants_in_species(subtree_root, srecon):
                self._orders_from_topology[subtree_root, descendant] = 1
            for child in subtree_root.children:
                infer_orders(child)
        infer_orders(gtree.root)   
        
        # o_{g1,g2} variables in paper, key = snode value = 1 if g2 more recent than g1, 0 otherwise
        order_keys = [(g1, g2) for (g1, g2) in pairs_in_species
                if (g1, g2) not in self._orders_from_topology and (g2, g1) not in self._orders_from_topology]

        self.order_vars = pulp.LpVariable.dicts("order", order_keys,
                                           0, 1, pulp.LpInteger)



    
    def _create_solver_vars(self, gtree, srecon, subtrees):
        """Create variables for solving ilp"""
            
        if self.all_vars:
            #========================================
            # structures containing groups of nodes
            self._bottom_nodes = collections.defaultdict(list)         # key = snode value = list of gnodes mapped to bottom of snode
            self._top_nodes_with_child = collections.defaultdict(list) # key = snode value = list of gnodes mapped to top of snode with child also mapped to snode
            self._top_nodes = collections.defaultdict(list)            # key = snode value = list of gnodes mapped to top of snode
            
            # find top and bottom nodes for each snode
            for snode, subtrees in subtrees.iteritems():
                for subtree in subtrees:
                    # the root of each subtree is a top node for this snode
                    self._top_nodes[snode].append(subtree[0])
                    if subtree[2]:
                        self._bottom_nodes[snode].extend(subtree[2])
                    if subtree[1]:
                        self._top_nodes_with_child[snode].append(subtree[0])         

            #========================================
            # creates loss variables
            loss_keys = []
            for snode, tops in self._top_nodes.iteritems(): # contains every gene mapped to top of snode 
                for top in tops:
                    loss_keys.append((snode, top))
            # l_sg variables in paper, key = (snode, gnode at top of snode), value = 1 if gnode creates a loss in snode, 0 otherwise
            self._loss_vars = pulp.LpVariable.dicts("loss", loss_keys, 0, 1, pulp.LpInteger) 

            #========================================
            # creates coalspec variables
            coalspec_keys = []
            for snode, tops in self._top_nodes_with_child.iteritems():
                for top in tops:
                    coalspec_keys.append((snode, top))
            # c_sg variables in paper, key = (snode, gene node at top of snode with child) value = 1 if gnode on same locus as another gnode at top of snode with child, 0 otherwise
            self._coalspec_vars = pulp.LpVariable.dicts("coal_spec", coalspec_keys, 0, 1, pulp.LpInteger) 

            #========================================
            # h_g variables in paper, key = (snode, gnode at top of snode) value = 1 if g on same locus as gnode2 mapped to top of snode and gnode2<gnode, 0 otherwise
            self._helper_vars = pulp.LpVariable.dicts("helper", self._loss_vars.keys(), 0, 1, pulp.LpInteger) 

            #========================================
            all_pairs = list(pulp.combination(list(gtree.preorder()), 2))
            # p_gg' variables in paper, key = pair of gene nodes, value = 1 if path between gene nodes has at least one duplication, 0 otherwise
            self._path_vars = pulp.LpVariable.dicts("path", all_pairs, 0, 1, pulp.LpInteger) 
            #========================================
            pairs_for_each_species = [list(pulp.combination(gnodes, 2)) for gnodes in self._gnodes_by_species.values()]
            pairs_in_species= sum(pairs_for_each_species, [])
            pairs_in_species_with_flipped = pairs_in_species + [(g2, g1) for (g1, g2) in pairs_in_species]
            delta_vars_keys = [(g1, g2) for g1, g2 in pairs_in_species_with_flipped
                        if g1.parent is not None and g2.parent is not None and
                        g1.parent != g2]
            # delta_g1g2 variables in paper, see paper for description
            self._delta_vars = pulp.LpVariable.dicts("coal_dup_helper", delta_vars_keys, 0, 1, pulp.LpInteger) 

            #========================================
            # r_g variables in paper, key = gnode, value = number of coalescense due to dup at g
            self._coaldup_vars = pulp.LpVariable.dicts("coal_dup", list(gtree.preorder()), 0, None, pulp.LpInteger) 

    def get_order(self, g1, g2):
        """Returns 1 if g2 more recent than g1, 0 otherwise.
        Checks order_vars, then orders_from_topology, then ???"""
        if (g1, g2) in self.order_vars:
            return self.order_vars[(g1, g2)]
        elif (g2, g1) in self.order_vars:
            return 1 - self.order_vars[(g2, g1)]
        elif (g1, g2) in self._orders_from_topology:
            return self._orders_from_topology[(g1, g2)]
        elif (g2, g1) in self._orders_from_topology:
            return 1 - self._orders_from_topology[(g2, g1)]
        elif self.srecon[g2] in self.srecon[g1].children: # g2 maps to a species node that is a child of the species node g1 maps to
            print("srecon", g1, g2)
            return 1
        else:
            raise Exception("Could not find order for nodes (%s,%s)" % (g1,g2))

    def get_path(self, g1, g2):
        """Returns 1 if there is a duplication on path between g1 and g2, 0 otherwise."""
        if g1 == g2:
            assert (g1,g2) not in self._path_vars
            return 0
        elif (g1, g2) in self._path_vars:
            return self._path_vars[(g1, g2)]
        elif (g2,g1) in self._path_vars:
            return self._path_vars[(g2, g1)]
        else:
            raise Exception("Could not find path variable for nodes (%s,%s)" % (g1,g2))



def ilp_to_lct(gtree, srecon, lpvars):
    """Converts ILPReconVariables to a LabeledRecon"""
    dup_nodes = [g for g, dup_var in lpvars.dup_vars.iteritems() if dup_var.varValue == 1.0]
    lrecon = _infer_locus_map(gtree, dup_nodes)
    order = make_order(lpvars.order_vars, srecon, lrecon)
    lct = reconlib.LabeledRecon(srecon, lrecon, order)
    return lct

def lct_to_ilp(gtree, labeledrecon):
    """ Converts LabeledRecon into ILPReconVariables """
    dup_vars = lrecon_to_dupvars(labeledrecon.locus_map)
    order_vars = order_to_ordervars(labeledrecon.order)
    lpvars = IlpReconVariables(gtree, labeledrecon.species_map, None, all_vars=False)
    for gnode, varValue in dup_vars.iteritems():
        lpvars.dup_vars[gnode] = varValue
    for gnodes, varValue in order_vars.iteritems():
        lpvars.order_vars[gnodes] = varValue 
    return lpvars

def _infer_locus_map(gtree, dup_nodes):
    """Infer locus map
    Adapted from recon._evolve_subtree"""
    locus = INIT_LOCUS
    lrecon = {}
    
    for gnode in gtree.preorder():
        # if root, it has the first locus
        if gnode == gtree.root:
            lrecon[gnode] = locus
        # if it has a dup, it has a different locus
        elif gnode in dup_nodes:
            locus += 1
            lrecon[gnode] = locus
        # if there's no dup, it's the same as the parent
        else:
            lrecon[gnode] = lrecon[gnode.parent]

    return lrecon
    

def make_order(order_vars, srecon, lrecon):
    """Returns order dictionary (see LabeledRecon)"""
    order = {}
    
    # create (snode, plocus) pairs for which corresponding gnodes will be ordered
    addSL = []
    for gnode in lrecon:
        if gnode.parent and lrecon[gnode] != lrecon[gnode.parent]:
            slPair = (srecon[gnode], lrecon[gnode.parent])
            addSL.append(slPair)

    for gnode in srecon:
        # only adds gene nodes if the gene node meets the requirement of being an internal node, or having more than 1 child or having a duplication leading to it
        # skip gene tree root
        if gnode.parent is None:
            continue
        
        # only adds gnodes if they are not a species leaf unless there is a dup from the gnode's parent to it, or the gnode has more than 1 child
        if lrecon[gnode]!=lrecon[gnode.parent] or not (len(gnode.children) == 0 or srecon[gnode] != srecon[gnode.children[0]]) or len(gnode.children)>1:
            snode = srecon[gnode]
            plocus = lrecon[gnode.parent]

            # next, the node will only be added if the species and plocus it corresponds to contain a duplication at some location
            if (snode, plocus) in addSL:                
                if snode in order:
                    if plocus in order[snode]:
                        order[snode][plocus].append(gnode)
                    else:
                        order[snode][plocus] = [gnode]
                else:
                    order[snode] = {plocus: [gnode]}

    # orders each gene list using order variables
    for snode in order:
        for plocus in order[snode]:
            order[snode][plocus] = _gnode_sort(order[snode][plocus], order_vars)
    return order

        
def _gnode_sort(gene_lst, order_vars):
    """This function uses a order variables to order a list of genes from the same species 
       and locus from least recent to most recent"""
    for i in range(len(gene_lst)-1):
        for j in range(len(gene_lst)-1):
            if (gene_lst[j+1], gene_lst[j]) in order_vars:
                    
                # if two consecutive genes are out of order, they are switched
                if order_vars[(gene_lst[j+1], gene_lst[j])] == 1:
                    temp = gene_lst[j]
                    gene_lst[j] = gene_lst[j+1]
                    gene_lst[j+1] = temp
    return gene_lst

def lrecon_to_dupvars(lrecon):
    """Makes dictionary of gene nodes and duplication variables indicating if there is a duplication lead to gnode """
    dup_vars = {}
    for gnode in lrecon:

        # if the gene node differs in locus from its parent, there must be a duplication leading to it
        if gnode.parent and lrecon[gnode] != lrecon[gnode.parent]:
            dup_vars[gnode] = 1.0

        # if the gene node has the same locus as its parent, or is the root node and has no parent, then there was no duplication
        else:
            dup_vars[gnode] = 0.0
    return dup_vars

def order_to_ordervars(order):
    """takes in order from LabeledRecon and returns a dictionary of order variables"""
    order_vars = {}
    for snode in order:
        for plocus in order[snode]:
            gene_lst = order[snode][plocus]
            for i, g1 in enumerate(len(gene_lst)-1):
                for g2 in gene_lst[i+1:]:
                    order_vars[g1, g2] = 1 # g1 is older than g2
    return order_vars