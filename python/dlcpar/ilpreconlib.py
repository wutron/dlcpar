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
    LCT variables:
        dup_vars               :    key = gnode, value = 1 if duplication on edge between gnode and gnode's parent, 0 otherwise
        order_vars             :    key = (gnode1, gnode2), value = 1 if gnode2 is more recent than gnode1, 0 otherwise
    solver variables:
        _loss_vars             :    key = (snode, gnode mapped to top of snode), value = 1 if gnode creates a loss in snode, 0 otherwise
        _path_vars             :    key = (gnode1, gnode2), value = 1 if there is at least one dup on path between gnode1 and gnode2, 0 otherwise
        _lambda_vars           :    key = (snode, gnode mapped to top of snode), value = 1 if g on same locus as gnode2 mapped to top of snode and gnode2 < gnode, 0 otherwise
        _coalspec_vars         :    key = (snode, gnode at top of snode with child), value = 1 if gnode on same locus as another gnode at top of snode with child, 0 otherwise
        _omega_vars            :    key = (gnode1 mapped to bottom of a snode, gnode2 of a snode), value = 1 if branch to gnode2 duplicates and branch to gnode1 doesn't, 0 otherwise
        _kappa_vars            :    key = (gnode1, gnode2), see paper for value description
        _coaldup_vars          :    key = gnode, value = number of coalesences due to a duplication at gnode 
    structure variables:
        _gnodes                :    key = snode, value = list of gnodes mapped to snode
        _orders_from_tree      :    list of keys (g1, g2), where g2 is more recent than g1 by topology
        _bottom_nodes          :    key = snode, value = list of gnodes mapped to bottom of snode
        _top_nodes_with_child  :    key = snode, value = list of gnodes mapped to top of snode with child also mapped to snode
        _top_nodes             :    key = snode, value = list of gnodes mapped to top of snode
        _cnodes                :    key = snode, value = list of childrens of bottoms of the snode
        _incomparable_nodes    :    list of keys (g1, g2) where g2 is a descendant of g1
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

        # d_g
        # key = gnode
        # value = 1 if duplication on edge e(g) leading to g, 0 otherwise
        self.dup_vars = pulp.LpVariable.dicts("dup", list(node.name for node in gtree.preorder()), 0, 1, pulp.LpBinary)
        
        #========================================
        # order variables

        # o_{g1,g2}
        # key = (g1,g2)
        # value = 1 if g2 more recent than g1, 0 otherwise

        # orders_from_tree: order helper dictionary
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
                for (top, topchild, bottoms) in subtrees_snode:
                    self._top_nodes[snode.name].append(top)
                    if topchild:
                        self._top_nodes_with_child[snode.name].append(top)
                    if bottoms:
                        self._bottom_nodes[snode.name].extend(bottoms) 

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
                self._incomparable_nodes.append((gnode.name, descendant.name))
            for descendant in descendants_not_in_species(snode):
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
        # loss variables
        # key = (snode, gnode at top of snode)
        # value = 1 if gnode creates a loss in snode, 0 otherwise
        self._loss_vars = pulp.LpVariable.dicts("loss", list((s,g.name) for (s,g) in loss_keys), 0, 1, pulp.LpBinary) 

        # p_{g,g'}
        # path helper variables
        # key = pair of gene nodes
        # value = 1 if path between gene nodes has at least one duplication, 0 otherwise

        # add path variable for a gnode to itself (used in kappa constraints)
        all_pairs = list(pulp.permutation(all_gnodes, 2))
        for g in all_gnodes:
            all_pairs.append((g,g))
        
        self._path_vars = pulp.LpVariable.dicts("path", list((n1.name,n2.name) for (n1,n2) in all_pairs), 0, 1, pulp.LpBinary) 

        # lambda_sg
        # loss helper variables
        # key = (snode, gnode at top of snode)
        # value = 1 if gnode on same locus as gnode2 mapped to top of snode and gnode2 < gnode, 0 otherwise
        self._lambda_vars = pulp.LpVariable.dicts("lambda", self._loss_vars.keys(), 0, 1, pulp.LpBinary) 
    

        #========================================
        # coalspec variables
    
        coalspec_keys = []
        for snode, tops in self._top_nodes_with_child.iteritems():
            for top in tops:
                coalspec_keys.append((snode, top))

        # c_sg 
        # coalspec variables
        # key = (snode, gene node at top of snode with child)
        # value = 1 if gnode on same locus as another gnode at top of snode with child, 0 otherwise
        self._coalspec_vars = pulp.LpVariable.dicts("coal_spec", list((s,g.name) for (s,g) in coalspec_keys), 0, 1, pulp.LpBinary) 
        
        #========================================
        # coaldup variables

        # omega_{g1,g2}
        # order helper variables used in coaldup constraints
        # key = (gnode1 mapped to bottom of a snode, gnode2 mapped to the bottom of that snode)
        # value = 1 if branch to gnode2 duplicates and branch to gnode1 doesn't and gnode1 is a leaf, 0 otherwise
        
        omega_keys = []
        leaf_species = [snode for snode in self.stree.leaves()]
        for snode in leaf_species:
            for g1 in self._bottom_nodes[snode.name]:
                for g2 in self._gnodes[snode]:
                    if g1 != g2:
                        omega_keys.append((g1, g2))

        self._omega_vars = pulp.LpVariable.dicts("omega", list((n1.name,n2.name) for (n1,n2) in omega_keys), 0, 1, pulp.LpBinary) 


        # kappa_{g1,g2} 
        # coaldup helper variables
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
        # coaldup variables
        # key = gnode
        # value = number of coalescenses due to dup at g
        self._coaldup_vars = pulp.LpVariable.dicts("coal_dup", list(node.name for node in all_gnodes), 0, None, pulp.LpInteger)

   
    #=============================================================================
    # rounding utilities for CPLEX_PY
    def get_binary_dicts(self):
        # every dictionary except self._coal_dup_vars, as it is not made of binary variables
        return [self.dup_vars, self.order_vars, self._loss_vars, self._path_vars, self._lambda_vars, \
            self._coalspec_vars, self._kappa_vars, self._omega_vars]


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
            for i in xrange(1, len(lst)):
                g1 = lst[i]
                j = i-1
                while j >= 0 and lpvars.order_vars[g1.name, lst[j].name].varValue == 1:
                    lst[j+1] = lst[j]
                    j -= 1
                lst[j+1] = g1
            
            # sanity check that all the order variables are satisfied by the order in lst (after the insertion sort)
            for (g1, g2) in list(pulp.combination(lst, 2)):
                assert lpvars.order_vars[g1.name, g2.name].varValue == 1.0, (((g1.name, g2.name),lst), "is not in the correct order \
                    with order value: ", lpvars.order_vars[g1.name, g2.name].varValue)
    
    #========================================
    # put everything together

    labeled_recon = reconlib.LabeledRecon(srecon, lrecon, order)
    #========================================
    # check conversion    
    stree = lpvars.stree.copy()

    coalspec_tuples = [(stree.nodes[snode], gtree.nodes[gnode]) for (snode, gnode), coalspec_var in lpvars._coalspec_vars.iteritems() if coalspec_var.varValue == 1.0]
    coaldup_nodes = [gtree.nodes[g] for g, coaldup_var in lpvars._coaldup_vars.iteritems() if coaldup_var.varValue > 0]
    ncoaldups = int(sum(var.varValue for var in lpvars._coaldup_vars.values()))
    loss_tuples = [(stree.nodes[snode], gtree.nodes[gnode]) for (snode, gnode), loss_var in lpvars._loss_vars.iteritems() if loss_var.varValue == 1.0]
       
    reconlib.init_dup_loss_coal_tree(stree)
    LCT_dups, LCT_losses, LCT_coalspecs, LCT_ncoaldups = get_event_vars_from_LCT(gtree, stree, lpvars, labeled_recon, True)

    if set(dup_nodes) != set(LCT_dups) or set(loss_tuples) != set(LCT_losses) or set(coalspec_tuples) != set(LCT_coalspecs) or ncoaldups != LCT_ncoaldups:
        return None
 

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


def get_event_vars_from_LCT(gtree, stree, lpvars, labeled_recon, conversion_check=False):
    
    # returns event vars from LCT

    dups = []
    losses = [] 
    coalspecs = []
    coaldups = []    

    extra = labeled_recon.get_dict()

    new_srecon = util.mapdict(extra["species_map"], val=lambda snode: stree.nodes[snode.name])
    new_order = util.mapdict(extra["order"], key=lambda snode: stree.nodes[snode.name])
    
    extra = extra.copy()
    
    srecon = new_srecon
    order = new_order

    extra["species_map"] = srecon
    extra["order"] = order

    events = phylo.label_events(gtree, srecon)
    subtrees = reconlib.factor_tree(gtree, stree, srecon, events)

    ncoalspec = 0
    ncoaldup = 0
    
    for snode in stree:
        subtrees_snode = subtrees[snode]
        if len(subtrees_snode) == 0:
            continue

        # count dups
        dup_nodes = reconlib.find_dup_snode(gtree, stree, extra, snode,
                                     subtrees, subtrees_snode, nodefunc=lambda node: node)
        # snode.data["dup"] += ndup_snode
        if len(dup_nodes):
            dups.extend(dup_nodes)

        # count losses
        loss_snodes = reconlib.find_loss_snode(gtree, stree, extra, snode,
                                       subtrees, subtrees_snode, nodefunc=lambda node: node)
        # snode.data["loss"] += nloss_snode
        if len(loss_snodes):
            for loss in loss_snodes:
                losses.append((snode, loss[0]))

        coalspec_lineages = reconlib.find_coal_spec_snode(gtree, stree, extra, snode,
                                       subtrees, subtrees_snode,
                                       nodefunc=lambda node: node,
                                       implied=True)
        
        for lineages in coalspec_lineages:
            assert len(lineages) > 1
            ncoalspec += len(lineages) - 1
            for i in range(1, len(lineages)):
                coalspecs.append((snode, lineages[i].parent))

        coaldup_lineages = reconlib.find_coal_dup_snode(gtree, stree, extra, snode,
                                     subtrees, subtrees_snode,
                                     nodefunc=lambda node: node)        

        for lineages in coaldup_lineages:
            assert len(lineages) > 1
            ncoaldup += len(lineages) - 1
            coaldup = [i for i in lpvars._bottom_nodes[snode.name] if i not in lineages] 
            coaldups.append((snode, lineages))

    # for coal in coaldups:
    #     print "bottom_nodes of " + str(coal[0]) + ": " + str(lpvars._bottom_nodes[coal[0].name])
    #     print "_gnodes of " + str(coal[0]) + ": " + str(lpvars._gnodes[coal[0]])
    #     print "_cnodes of " + str(coal[0]) + ": " + str(lpvars._cnodes[coal[0]])
    #     print type(coal[0])
    


    # print "\ndups: ", dups
    # print "\nlosses: ", losses
    # print "\ncoalspec lineages", coalspecs
    # print "ncoalspecs: ", ncoalspec
    # print "\ncoaldup lineages", coaldups
    # print "ncoaldups: ", ncoaldup

    return dups, losses, coalspecs, ncoaldup