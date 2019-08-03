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
        order_vars         :    key = (gnode1, gnode2), value = 1 if gnode2 is more recent than gnode1, 0 otherwise
    solver variables:
        _loss_vars         :    key = (snode, gnode mapped to top of snode), value = 1 if gnode creates a loss in snode, 0 otherwise
        _path_vars         :    key = (gnode1, gnode2), value = 1 if there is at least one dup on path between gnode1 and gnode2, 0 otherwise
        _lambda_vars       :    key = (snode, gnode mapped to top of snode), value = 1 if g on same locus as gnode2 mapped to top of snode and gnode2 < gnode, 0 otherwise
        _coalspec_vars     :    key = (snode, gnode at top of snode with child), value = 1 if gnode on same locus as another gnode at top of snode with child, 0 otherwise
        _coaldup_vars      :    key = gnode, value = number of coalesences due to a duplication at gnode 
        _kappa_vars        :    key = (gnode1, gnode2), see paper for value description
        _omega_vars         :    key = (gnode1 mapped to bottom of a snode, gnode2 of a snode), value = 1 if branch to gnode2 duplicates and branch to gnode1 doesn't, 0 otherwise
    structure variables:
        _gnodes :    key = snode, value = list of gnodes mapped to snode
        _orders_from_tree : key = (gnode1, gnode2), value = 1 if gnode2 more recent than gnode1, 0 otherwise
        _bottom_nodes :         key = snode, value = list of gnodes mapped to bottom of snode
        _top_nodes_with_child : key = snode, value = list of gnodes mapped to top of snode with child also mapped to snode
        _top_nodes :            key = snode, value = list of gnodes mapped to top of snode
        _cnodes :    key = snode, value = list of childrens of bottoms of the snode
        _incomparable_nodes : list of keys (g1, g2) where g2 is a descendant of g1
    """

    def __init__(self, gtree, stree, srecon, all_vars=True):
        self.gtree = gtree
        self.stree = stree
        self.srecon = srecon

        self._create_recon_vars(all_vars)
        if all_vars:
            self._create_solver_vars()


    def _create_recon_vars(self,all_vars):
        """Creates variables necessary for converting to LCT and inferring evolutionary events"""

        gtree = self.gtree
        srecon = self.srecon

        #========================================
        # duplication variables

        # d_e, key = gnode, value = 1 if duplication on edge between gnode and gnode's parent, 0 otherwise
        self.dup_vars = pulp.LpVariable.dicts("dup", list(node.name for node in gtree.preorder()), 0, 1, pulp.LpBinary)
        #========================================
        # order variables
        # o_{g1,g2}, key = (g1,g2), value = 1 if g2 more recent than g1, 0 otherwise
        # order variables in paper include both orders from tree and order vars set by the ilp

        # infer required orders from gene tree topology
        # key = (g1,g2) such that g1,g2 in same species and one node is ancestor of other
        #       note that only one of (g1,g2) or (g2,g1) is in dictionary
        # value = 1 if g2 more recent than g1 and 0 otherwise

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
            # for node in self._bottom_nodes[snode]:
            #     children.extend(node.children)

            for s in snode.children:
                for g2 in self._gnodes[s]:
                    children.append(g2)
            return children 

        if all_vars:
            sevents = phylo.label_events(self.gtree, self.srecon)
            subtrees = reconlib.factor_tree(self.gtree, self.stree, self.srecon, sevents)

            #========================================
            # structures containing groups of nodes
            self._bottom_nodes = collections.defaultdict(list)
            self._top_nodes_with_child = collections.defaultdict(list)
            self._top_nodes = collections.defaultdict(list)
            for snode, subtrees_snode in subtrees.iteritems():
                for (root, rootchild, leaves) in subtrees_snode:
                    self._top_nodes[snode.name].append(root)
                    if rootchild:
                        self._top_nodes_with_child[snode.name].append(root)
                    if leaves:
                        self._bottom_nodes[snode.name].extend(leaves) 

        # create order_(g1, g2) variables not already ordered by tree
        self._gnodes = collections.defaultdict(list)
        self._cnodes = collections.defaultdict(list)
        for gnode in gtree:
            snode = srecon[gnode]
            self._gnodes[snode].append(gnode)
            if gnode in self._bottom_nodes[snode.name]:
                for child in gnode.children:
                    assert child not in self._cnodes[snode]
                    self._cnodes[snode].append(child)

        # create orders_from_tree, contains order_(gnodes, descendants) determined through the tree
        # create incomparable_nodes, used to determine keys for kappa dictionary 
        self._orders_from_tree = []
        self._incomparable_nodes = []
        for gnode in gtree:
            snode = srecon[gnode]
            for descendant in descendants_in_species(gnode):
                self._orders_from_tree.append((gnode.name, descendant.name))
                # self._orders_from_tree[gnode, descendant] = 1
                self._incomparable_nodes.append((gnode.name, descendant.name))
            for descendant in descendants_not_in_species(snode):
                # self._orders_from_tree[gnode, descendant] = 1
                if (gnode, descendant) not in self._orders_from_tree:
                    self._orders_from_tree.append((gnode.name, descendant.name))

        order_keys = []
        for (g1, g2) in self._orders_from_tree:
            order_keys.append((g1, g2))
            order_keys.append((g2, g1))
        
        for snode in self.stree:
            for g1 in self._gnodes[snode]:                
                for g2 in self._gnodes[snode]:
                    if g1 != g2 and (g1, g2) not in order_keys:
                        order_keys.append((g1.name,g2.name))

        self.order_vars = pulp.LpVariable.dicts("order", order_keys, 0, 1, pulp.LpBinary)
        
        # out_order = util.open_stream("order_keys.txt", 'w')
        # out_order.write("order keys:\t" + str(order_keys))
        # out_order.write("\n\ngnodes:\t" + str(self._gnodes))
        # out_order.write("\n\ncnodes:\t" + str(self._cnodes))
        # out_order.write("\n\norders_from_tree:\t" + str(self._orders_from_tree))


    def _create_solver_vars(self):
        """Create variables for solving ILP"""

        all_gnodes = list(self.gtree.preorder())
        #========================================
        # loss variables

        loss_keys = []
        for snode, tops in self._top_nodes.iteritems(): # contains every gene mapped to top of snode 
            for top in tops:
                loss_keys.append((snode, top))

        # l_sg
        # key = (snode, gnode at top of snode)
        # value = 1 if gnode creates a loss in snode, 0 otherwise
        self._loss_vars = pulp.LpVariable.dicts("loss", list((s,g.name) for (s,g) in loss_keys), 0, 1, pulp.LpBinary) 

        #========================================
        # coalspec variables
    
        coalspec_keys = []
        for snode, tops in self._top_nodes_with_child.iteritems():
            for top in tops:
                coalspec_keys.append((snode, top))

        # c_sg variables
        # key = (snode, gene node at top of snode with child)
        # value = 1 if gnode on same locus as another gnode at top of snode with child, 0 otherwise
        self._coalspec_vars = pulp.LpVariable.dicts("coal_spec", list((s,g.name) for (s,g) in coalspec_keys), 0, 1, pulp.LpBinary) 
        #========================================
        # constraint variables

        # lambda_g
        # key = (snode, gnode at top of snode)
        # value = 1 if gnode on same locus as gnode2 mapped to top of snode and gnode2 < gnode, 0 otherwise
        self._lambda_vars = pulp.LpVariable.dicts("lambda", self._loss_vars.keys(), 0, 1, pulp.LpBinary) 
        
        # p_{g,g'}
        # key = pair of gene nodes
        # value = 1 if path between gene nodes has at least one duplication, 0 otherwise

        # add path variable for a gnode to itself (used in kappa constraints)
        all_pairs = list(pulp.permutation(all_gnodes, 2))
        for g in all_gnodes:
            all_pairs.append((g,g))
        
        self._path_vars = pulp.LpVariable.dicts("path", list((n1.name,n2.name) for (n1,n2) in all_pairs), 0, 1, pulp.LpBinary) 

        #========================================
        # coaldup variables

        # kappa_{g1,g2} variables
        # key = pair of gene nodes
        # value = see paper for description
        kappa_vars_keys = []
        for snode in self.stree:
            for g1 in self._gnodes[snode] + self._cnodes[snode]:
                for g2 in self._gnodes[snode]:
                    if g1 != g2 and \
                       (g1.name, g2.name) not in self._incomparable_nodes and \
                       (g2.name, g1.name) not in self._incomparable_nodes:
                        if g1.parent and g2.parent:
                            kappa_vars_keys.append((g1,g2))

        self._kappa_vars = pulp.LpVariable.dicts("coal_dup_helper", list((n1.name,n2.name) for (n1,n2) in kappa_vars_keys), 0, 1, pulp.LpBinary) 
        
        # k_g
        # key = gnode
        # value = number of coalescenses due to dup at g
        self._coaldup_vars = pulp.LpVariable.dicts("coal_dup", list(node.name for node in all_gnodes), 0, None, pulp.LpInteger)

        #========================================
        # omega variables

        omega_keys = []
        leaf_species = [snode for snode in self.stree.leaves()]
        for snode in leaf_species:
            for g1 in self._bottom_nodes[snode.name]:
                for g2 in self._gnodes[snode]:
                    if g1 != g2:
                        omega_keys.append((g1, g2))

        # omega variables
        # key = (gnode1 mapped to bottom of a snode, gnode2 mapped to the bottom of that snode)
        # value = 1 if branch to gnode2 duplicates and branch to gnode1 doesn't and gnode1 is a leaf, 0 otherwise
        self._omega_vars = pulp.LpVariable.dicts("omega", list((n1.name,n2.name) for (n1,n2) in omega_keys), 0, 1, pulp.LpBinary) 

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

    # dup_nodes = [g for g, dup_var in lpvars.dup_vars.iteritems() if dup_var.varValue == 1.0]
    dup_nodes = [gtree.nodes[g] for g, dup_var in lpvars.dup_vars.iteritems() if dup_var.varValue == 1.0]

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
            # "insertion sort" the list
            # print("list\t" + str(lst))
            for i in xrange(1, len(lst)):
                g1 = lst[i]
                j = i-1
                
                # print lpvars.order_vars[g1.name, lst[j].name].varValue
                while j >= 0 and lpvars.order_vars[g1.name, lst[j].name].varValue == 1:
                    # print g1, lst[j]
                    lst[j+1] = lst[j]
                    j -= 1
                lst[j+1] = g1
            # print("list done\t" + str(lst))
            # sanity check that all the order variables are satisfied by the order in lst (after the insertion sort)
            for (g1, g2) in list(pulp.combination(lst, 2)):
                assert lpvars.order_vars[g1.name, g2.name].varValue==1, (((g1.name, g2.name),lst), "is not in the correct order")
                
            
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
        lpvars.dup_vars[gnode.name].varValue = var_value
    for gnodes, var_value in order_vars.iteritems():
        if gnodes in lpvars.order_vars:
            lpvars.order_vars[gnodes.name].varValue = var_value
    return lpvars

