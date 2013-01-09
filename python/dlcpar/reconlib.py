#
# Reconciliation library
# 

# python libraries
import copy
import collections

# rasmus libraries
from rasmus import treelib
from rasmus import util

# compbio libraries
from compbio import phylo

# yjw libraries
from yjw.bio import phyloDLC

# dlcpar libraries
import dlcpar
from dlcpar import common

#=============================================================================
# reconciliation data structures


class LabeledRecon (object):
    """The reconciliation data structure for the DLCpar model

    species_map :   dict with key = gene_tree node, value = species_tree node
    locus_map :     dict with key = gene_tree node, value = locus (int)
    order :         2D dict with key1 = species branch node (snode), key2 = locus,
                      value = ordered list of gene_tree nodes (nodes), where
                              nodes are the internal nodes for branch snode and
                              the parents of these nodes are equal to locus

    The gene_tree should contain all implied speciation nodes.
    """

    def __init__(self, species_map=None, locus_map=None, order=None):
        self.species_map = species_map
        self.locus_map = locus_map
        self.order = order

    def sort_loci(self, gtree):
        """Sorts loci in locus_map.
        
        Loci are sorted using leaf names then by pre-order index.
        """

        # create mapping
        m = {}
        next = 1        
        for node in sorted(gtree.leaves(), key=lambda node: node.name):
            locus = self.locus_map[node]
            if locus not in m:
                m[locus] = next
                next += 1
        for node in gtree.preorder():
            if node.is_leaf():
                continue
            locus = self.locus_map[node]
            if locus not in m:
                m[locus] = next
                next += 1
                
        # remap
        locus_map = {}
        for node, locus in self.locus_map.iteritems():
            locus_map[node] = m[locus]
        order = {}
        for snode, d in self.order_iteritems():
            order[snode] = {}
            for locus, lst in d.iteritems():
                order[snode][m[locus]] = lst
        
        self.locus_map = locus_map
        self.order = order
    

    def get_dict(self):
        return {"species_map": self.species_map,
                "locus_map" : self.locus_map,
                "order" : self.order}


    def write(self, filename, gtree,
              exts={"tree" : ".tree",
                    "recon" : ".recon",
                    "order" : ".order"},
              filenames={}):
        """Write the reconciliation to a file"""
        assert gtree and self.species_map and self.locus_map and (self.order is not None)

        gtree.write(
            filenames.get("tree", filename + exts["tree"]),
            rootData=True)
        
        util.write_delim(
            filenames.get("recon", filename + exts["recon"]),
            [(str(node.name), str(snode.name), self.locus_map[node])
             for node, snode in self.species_map.iteritems()])

        order = {}
        for snode, d in self.order.iteritems():
            for locus, lst in d.iteritems():
                order[snode, locus] = lst
        util.write_dict(
            filenames.get("order", filename + exts["order"]),
            order,
            keyfunc=lambda (snode, locus): "%s\t%s" % (str(snode.name), str(locus)),
            valfunc=lambda lst: ",".join(map(lambda x: str(x.name), lst)))


    def read(self, filename, stree,
             exts={"tree" : ".tree",
                   "recon" : ".recon",
                   "order" : ".order"},
             filenames={}):
        """Read the reconciliation from a file"""

        gtree = treelib.read_tree(
            filenames.get("tree", filename + exts["tree"]))

        self.species_map = {}
        self.locus_map = {}
        for name, sname, locus in util.read_delim(filenames.get("recon", filename + exts["recon"])):
            if name.isdigit(): name = int(name)
            if sname.isdigit(): sname = int(sname)
            assert locus.isdigit()
            locus = int(locus)
            
            node = gtree.nodes[name]
            self.species_map[node] = stree.nodes[sname]
            self.locus_map[node] = locus

        self.order = collections.defaultdict(dict)
        for toks in util.read_delim(filenames.get("order", filename + exts["order"])):
            sname, locus, lst = toks[0], toks[1], toks[2].split(',')
            if sname.isdigit(): sname = int(sname)
            assert locus.isdigit()
            locus = int(locus)
            names = map(lambda x: int(x) if x.isdigit() else x, lst)
            
            snode = stree.nodes[sname]
            nodes = map(lambda x: gtree.nodes[x], names)
            if snode not in self.order:
                self.order[snode] = {}
            self.order[snode][locus] = nodes
        self.order = dict(self.order)

        return gtree, self.get_dict()


#=============================================================================
# convenient input/output


def write_labeled_recon(filename, gtree, extra,
                        exts={"tree" : ".tree",
                              "recon" : ".recon",
                              "order" : ".order"},
                        filenames={}):
    """Writes a labeled reconciliation to files"""
    
    labeled_recon = LabledRecon(extra["species_map"], extra["locus_map"], extra["order"])
    labeled_recon.write(filename, gtree, exts, filenames)


def read_labeled_recon(filename, stree,
                       exts={"tree" : ".tree",
                             "recon" : ".recon",
                             "order" : ".order"},
                       filenames={}):
    """Reads a labeled reconciliation from files"""
    
    labeled_recon = LabeledRecon()
    return labeled_recon.read(filename, stree, exts, filenames)


#=============================================================================
# conversion utilities


def get_tree_depths(tree, distfunc=lambda node: node.dist):
    """Get depth of all nodes (depth = length of path from root to node)"""
    depths = {}
    for node in tree.preorder():
        if not node.parent:
            depths[node] = 0
        else:
            depths[node] = depths[node.parent] + distfunc(node)
    return depths


def recon_to_labeledrecon(coal_tree, recon, stree, gene2species,
                          name_internal="n", locus_mpr=True):
    """Convert from DLCoal to DLCpar reconciliation model

    If locus_mpr is set (default), use MPR from locus_tree to stree.
    """
    
    gene_tree = coal_tree.copy()
    coal_recon = recon.coal_recon
    locus_tree = recon.locus_tree
    if not locus_mpr:
        locus_recon = recon.locus_recon
        daughters = recon.daughters
    else:
        locus_recon = phylo.reconcile(locus_tree, stree, gene2species)
        locus_events = phylo.label_events(locus_tree, locus_recon)
        daughters = filter(lambda node: locus_events[node.parent] == "dup", recon.daughters)

    #========================================
    # find species map

    # find species tree subtree
    substree = treelib.subtree(stree, locus_recon[coal_recon[coal_tree.root]])
    
    # find species map
    species_map = {}
    for node in gene_tree:
        cnode = coal_tree.nodes[node.name]
        lnode = coal_recon[cnode]
        snode = locus_recon[lnode]
        species_map[node] = substree[snode.name]

    # add implied speciation and delay nodes to gene tree
    events = phylo.label_events(gene_tree, species_map)
    added_spec, added_dup, added_delay = add_implied_nodes(gene_tree, substree, species_map, events)

    # rename nodes
    common.rename_nodes(gene_tree, name_internal)

    #========================================
    # helper functions
    
    def walk_up(node):
        if node.name in coal_tree.nodes:
            return coal_tree.nodes[node.name]
        return walk_up(node.parent)
    
    def walk_down(node):
        if node.name in coal_tree.nodes:
            return coal_tree.nodes[node.name]
        assert len(node.children) == 1, (node.name, node.children)
        return walk_down(node.children[0])

    #========================================
    # find locus map
    
    # label loci in locus tree
    loci = {}
    next = 1
    # keep track of duplication ages (measured as dist from leaf since root dist may differ in coal and locus trees)
    locus_times = treelib.get_tree_ages(locus_tree)
    dup_times = {}
    dup_snodes = {}
    for lnode in locus_tree.preorder():
        if not lnode.parent:            # root
            loci[lnode] = next
        elif lnode in daughters:        # duplication
            next += 1
            loci[lnode] = next
            dup_times[next] = locus_times[lnode.parent]
            dup_snodes[next] = locus_recon[lnode.parent]
        else:                           # regular node
            loci[lnode] = loci[lnode.parent]

    # label loci in gene tree
    locus_map = {}
    for node in gene_tree:
        if node.name in coal_tree.nodes:
            # node in coal tree
            cnode = coal_tree.nodes[node.name]
            lnode = coal_recon[cnode]
            locus_map[node] = loci[lnode]
        else:
            # node not in coal tree, so use either parent or child locus
            cnode_up = walk_up(node)
            lnode_up = coal_recon[cnode_up]
            loci_up = loci[lnode_up]

            cnode_down = walk_down(node)
            lnode_down = coal_recon[cnode_down]
            loci_down = loci[lnode_down]
            
            if loci_up == loci_down:
                # parent and child locus match
                locus_map[node] = loci[lnode_up]
            else:
                # determine whether to use parent or child locus
                snode = species_map[node]
                dup_snode = dup_snodes[loci_down]
		if (snode.name == dup_snode.name) or (snode.name in dup_snode.descendant_names()):
                    locus_map[node] = loci_down
                else:
                    locus_map[node] = loci_up
        
    #========================================
    # find order
    
    # find loci that give rise to new loci in each sbranch
    parent_loci = set()
    for node in gene_tree:
        if node.parent:
            locus = locus_map[node]
            plocus = locus_map[node.parent]

            if locus != plocus:
                snode = species_map[node]
                parent_loci.add((snode, plocus))        

    # find order (locus tree and coal tree must use same timescale)
    order = {}
    for node in gene_tree:
        if node.parent:
            snode = species_map[node]
            plocus = locus_map[node.parent]
            
	    if (snode, plocus) in parent_loci:
                order.setdefault(snode, {})
                order[snode].setdefault(plocus, [])
                order[snode][plocus].append(node)

    # find coalescent/duplication times (= negative age) and depths
    coal_times = treelib.get_tree_ages(coal_tree)
    depths = get_tree_depths(gene_tree, distfunc=lambda node: 1)
    def get_time(node):
        if locus_map[node.parent] != locus_map[node]:
            # duplication
            return -dup_times[locus_map[node]], depths[node]
        else:
            # walk up to the nearest node in the coal tree
            # if the node was added (due to spec or dup), it has a single child
            # so it can be placed directly after its parent without affecting the extra lineage count
            if node.name in coal_tree.nodes:
                cnode = coal_tree.nodes[node.name]
            else:
                cnode = walk_up(node)
	    return -coal_times[cnode], depths[node]

    # sort by node times
    # 1) larger age (smaller dist from root) are earlier in sort
    # 2) if equal dist, then smaller depths are earlier in sort
    for snode, d in order.iteritems():
        for locus, lst in d.iteritems():
            lst.sort(key=get_time)

    #========================================
    # put everything together

    return gene_tree, LabeledRecon(species_map, locus_map, order)


def labeledrecon_to_recon(labeled_recon, stree):
    """Convert from DLCpar to DLCoal reconciliation model

    NOTE: This is non-reversible because it produces NON-dated coalescent and locus trees
    """
    
    gene_tree = labeled_recon.gene_tree
    species_map = labeled_recon.species_map
    order = labeled_recon.order

    #========================================
    # coalescent tree
    coal_tree = gene_tree.copy()
    treelib.remove_single_children(coal_tree)

    #========================================
    # locus tree

    # factor gene tree
    events = phylo.label_events(gene_tree, species_map)
    subtrees = factor_tree(gene_tree, stree, species_map, events)

    # build locus tree along each species branch
    locus_tree = treelib.Tree()
    root = treelib.TreeNode(locus_tree.new_name())
    locus_tree.add(root)
    locus_tree.root = root

    coal_recon = {}
    locus_recon = {}
    daughters = []
    for snode in stree.preorder():
        subtrees_snode = subtrees[snode]
        if len(subtrees_snode) == 0:
            continue

##        coal_recon[node] = lnode
##        locus_recon[lnode] = snode

	# duplication
        # speciation

##        for (root, rootchildren, leaves) in subtrees_snode:
##            newnode = treelib.TreeNode(locus_tree.new_name())
    
    locus_events = phylo.label_events(locus_tree, locus_recon)

    #========================================
    # put everything together

    return coal_tree, phyloDLC.Recon(coal_recon, locus_tree, locus_recon, locus_events, daughters)


#=============================================================================
# gene tree factorization functions


def is_full_tree(tree, stree, recon, events):
    """Check that the tree has all implied internal nodes AND no extra nodes"""

    for node in tree:
        if events[node] == "gene":
            continue

        snode = recon[node]
        schildren = snode.children
        nschildren = len(schildren)
        if nschildren != 2 and nschildren != 0:
            raise Exception("Species tree must be binary")
        schildren2 = [recon[child] for child in node.children]

        if events[node] == "spec":
            if len(node.children) == 1:
                # speciation followed by loss
                if schildren[0] != schildren2[0] and \
                   schildren[1] != schildren2[0]:
                    return False
            elif len(node.children) == 2:
                # speciation
                if set(schildren) != set(schildren2):
                    return False
            else:
                raise Exception("Cannot handle non-binary trees")

        elif events[node] == "dup":
            if len(node.children) == 1:
                # extra node
                raise Exception("Tree contains extra node under %s" % node.name)
            elif len(node.children) == 2:
                # duplication
                if not (snode == schildren2[0] == schildren2[1]):
                    return False
            else:
                raise Exception("Cannot handle non-binary trees")

    return True


def add_spec_from_dup_nodes(node, tree, recon, events):
   """Relabel the current speciation node 'node' as a duplication.
      Insert new speciation nodes BELOW gene node 'node'.
      New nodes reconcile to same species node as 'node'.
      Modifies recon and events accordingly.
   """

   assert events[node] == "spec"
   snode = recon[node]
   events[node] = "dup"

   # insert new nodes into tree
   added = []
   for child in list(node.children):
       added.append(phylo.add_spec_node(child, snode, tree, recon, events))
   
   return added


def add_implied_spec_nodes(tree, stree, recon, events):
    """Add speciation nodes to tree that are implied but are not present
    because of gene losses.

    Extends phylo.add_implied_spec_nodes to handle non-MPR.
    Only guaranteed to work for binary trees.
    """

    added_spec = phylo.add_implied_spec_nodes(tree, stree, recon, events)
    
    added_dup = []
    for node in list(tree):
        schildren = [recon[child] for child in node.children]
        if len(schildren) > 1 and len(set(schildren)) == 1 and events[node] != "dup":
            added_dup.extend(add_spec_from_dup_nodes(node, tree, recon, events))

    assert is_full_tree(tree, stree, recon, events)

    return added_spec, added_dup


def add_delay_nodes(node, tree, recon, events):
    """Insert new delay nodes BELOW gene node 'node'.
       New nodes reconcile to same species node as 'node'.
       Modifies recon and events accordingly.
    """
    assert events[node] == "spec"
    snode = recon[node]

    # insert new nodes into tree
    added = []
    for child in list(node.children):
        newnode = phylo.add_spec_node(child, snode, tree, recon, events)
        events[newnode] = "delay"
        added.append(newnode)

    return added


def add_implied_delay_nodes(tree, stree, recon, events):
    """Add nodes to tree after each speciation to allow
       for delay between coalescence and speciation
    """

    added = []
    for node in list(tree):
        if events[node] == "spec":
            added.extend(add_delay_nodes(node, tree, recon, events))
    return added


def add_implied_nodes(tree, stree, recon, events):
    """Wrapper to add implied speciation nodes and delay nodes"""

    added_spec, added_dup = add_implied_spec_nodes(tree, stree, recon, events)
    added_delay = add_implied_delay_nodes(tree, stree, recon, events)
    return added_spec, added_dup, added_delay


def get_subtree(node, snode, recon, events):
    """Returns the leaves of a duplication subtree"""

    leaves = []

    def walk(node):
        if recon[node] != snode:
            return
        if events[node] != "dup":
            leaves.append(node)
        else:
            for child in node.children:
                walk(child)
    walk(node)

    return leaves


def factor_tree(tree, stree, recon, events):
    """Returns subtrees for each species branch

    Output is a dict with key = snode, val = (R(subtree), start, leaves(subtree)),
    where start is the node of the gene tree on which to recur
    (either the root of the subtree or its child).
    """

    # initialize
    subtrees = {}
    for snode in stree:
        subtrees[snode] = []

    # root is a dup
    if events[tree.root] != "spec":
        snode = recon[tree.root]
        subleaves = get_subtree(tree.root, snode, recon, events)
        subtrees[snode].append((tree.root, tree.root, subleaves))

    # find subtrees
    for node in tree:
        if events[node] == "spec":
            snode = recon[node]
            for schild in snode.children:
                nodes2 = [x for x in node.children if recon[x] == schild]
                if len(nodes2) > 0:
                    assert len(nodes2) == 1, (node, nodes2)
                    node2 = nodes2[0]
                    subleaves = get_subtree(node2, schild, recon, events)
                    subtrees[schild].append((node, node2, subleaves))
                else:
                    subtrees[schild].append((node, None, None))

    return subtrees


#=============================================================================    
# event functions


def _subtree_helper_snode(tree, stree, extra, snode, subtrees=None):
    """Returns the subtrees for a species branch

       These subtrees are specified by
         rootchildren : dict of (root, rootchild) pairs
         leaves       : list of leaves
    """

    if subtrees is None:
        subtrees = _subtree_helper(tree, stree, extra)
        subtrees = subtrees[snode]
            
    rootchildren = {}
    all_leaves = []
    for (root, rootchild, leaves) in subtrees:
        rootchildren[root] = rootchild
        if leaves is not None:
            all_leaves.extend(leaves)
    return rootchildren, all_leaves


def _subtree_helper(tree, stree, extra, subtrees=None):
    """Returns a dictionary of subtrees for each species branch"""

    if subtrees is None:
        recon = extra["species_map"]
        events = phylo.label_events(tree, recon)
        subtrees = factor_tree(tree, stree, recon, events)
    return subtrees


def find_dup_snode(tree, stree, extra, snode,
                   subtrees=None, rootchildren=None, leaves=None,
                   nodefunc=lambda node: node):
    """Find duplication nodes (node locus differs from parent node locus)"""

    if (rootchildren is None) and (leaves is None):
        rootchildren, leaves = _subtree_helper_snode(tree, stree, extra, snode, subtrees)
    loci = extra["locus_map"]

    dup_nodes = []
    for root, rootchild in rootchildren.iteritems():
        if not rootchild:
            continue
              
        for node in tree.preorder(rootchild, is_leaf=lambda x: x in leaves):
            if node.parent:
                if loci[nodefunc(node)] != loci[nodefunc(node.parent)]:
                    dup_nodes.append(node)
    return dup_nodes


def count_dup_snode(tree, stree, extra, snode,
                    subtrees=None, rootchildren=None, leaves=None,
                    nodefunc=lambda node: node):
    """Find number of inferred duplications in a species branch"""

    return len(find_dup_snode(tree, stree, extra, snode,
                              subtrees, rootchildren, leaves,
                              nodefunc))


def count_dup(tree, stree, extra,
              subtrees=None,
              nodefunc=lambda node: node):
    """Returns the number of inferred duplications"""

    subtrees = _subtree_helper(tree, stree, extra, subtrees)

    ndup = 0
    for snode in stree:
        ndup += count_dup_snode(tree, stree, extra, snode, subtrees[snode], nodefunc=nodefunc)
    return ndup


def count_loss_snode(tree, stree, extra, snode,
                     subtrees=None, rootchildren=None, leaves=None,
                     nodefunc=lambda node: node):
    """Find number of inferred losses in a species branch"""

    if (rootchildren is None) and (leaves is None):
        rootchildren, leaves = _subtree_helper_snode(tree, stree, extra, snode, subtrees)
    loci = extra["locus_map"]

    all_loci = set()
    leaf_loci = set()
    for root, rootchild in rootchildren.iteritems():
        all_loci.add(loci[nodefunc(root)])
        
        if not rootchild:
            continue

        for node in tree.preorder(rootchild, is_leaf=lambda x: x in leaves):
            all_loci.add(loci[nodefunc(node)])
            if node in leaves:
                leaf_loci.add(loci[nodefunc(node)])
    nloss = len(all_loci.difference(leaf_loci))
    return nloss


def count_loss(tree, stree, extra,
               subtrees=None,
               nodefunc=lambda node: node):
    """Returns the number of inferred losses"""

    subtrees = _subtree_helper(tree, stree, extra, subtrees)

    nloss = 0
    for snode in stree:
        nloss += count_loss_snode(tree, stree, extra, subtrees[snode], nodefunc=nodefunc)
    return nloss


def count_coal_snode_dup(tree, stree, extra, snode,
                         subtrees=None, rootchildren=None, leaves=None,
                         nodefunc=lambda node: node):
    """Returns the number of inferred extra lineages in a species branch
       (at duplication nodes in the locus tree)"""

    if (rootchildren is None) and (leaves is None):
        rootchildren, leaves = _subtree_helper_snode(tree, stree, extra, snode, subtrees)
    loci, order = extra["locus_map"], extra["order"]

    ncoal = 0
    if snode not in order:
        return ncoal
   
    # find "start" branches of each locus
    start = collections.defaultdict(list)    
    # find "parent loci" (loci that give rise to new loci) within the species branch
    dup_nodes = find_dup_snode(tree, stree, extra, snode, subtrees, rootchildren, leaves, nodefunc)
    parent_loci = set()
    for node in dup_nodes:
        if node.parent:
            parent_loci.add(loci[nodefunc(node.parent)])
    # for each locus found, if this locus is a "parent locus",
    # add the children if the dup node is not a leaf
    # (leaves never incur extra lineages in this species branch)
    for node in dup_nodes:
        locus = loci[nodefunc(node)]
        if (node in leaves) or (locus not in parent_loci):
            continue
        for child in node.children:
            start[locus].append(child)
    # for each locus that exists at the top of the species branch,
    # if this locus is a "parent locus",
    # add the child if it exists (assume immediate loss if child does not exist)
    for root, rootchild in rootchildren.iteritems():
        locus = loci[nodefunc(root)]
        if locus not in parent_loci:
            continue
        
        # handle root separately
        if not root.parent:
            for child in root.children:
                start[locus].append(child)
        else:            
            if rootchild:
                start[locus].append(rootchild)
    assert set(start) == set(order[snode]), (dict(start), order[snode])

    # for each locus in the species branch, walk down the subtrees using order
    # to determine the number of extra lineages
    count = {}
    for locus, lst in order[snode].iteritems():
        current = start[locus]
        num_lineages = len(current)
        for next in lst:
            assert num_lineages == len(current), (num_lineages, lst)
            assert next in current, (next, current)
            current.remove(next)
            num_lineages -= 1

            # update lineage count and list of nodes
            next_locus = loci[nodefunc(next)]
            if locus == next_locus:
                # deep coalescence - keep even if next in leaves to allow for delay btwn coalescence and speciation
                for child in next.children:
                    current.append(child)
                    num_lineages += 1
	    else:
                # duplication
		if num_lineages > 1:
                    ncoal += num_lineages - 1

    return ncoal    


def count_coal_snode_spec(tree, stree, extra, snode,
                          subtrees=None, rootchildren=None, leaves=None,
                          nodefunc=lambda node: node,
                          implied=True):
    """Returns the number of inferred extra lineages in a species branch
       (at speciation nodes in the locus tree)"""

    if (rootchildren is None) and (leaves is None):
        rootchildren, leaves = _subtree_helper_snode(tree, stree, extra, snode, subtrees)
    loci, order = extra["locus_map"], extra["order"]

    if not implied:
        raise Exception("not implemented")

    root_loci = {}
    for root, rootchild in rootchildren.iteritems():
        # only count if there is a tree branch in the species branch
        if rootchild:
            locus = loci[nodefunc(root)]
            root_loci.setdefault(locus, 0)
            root_loci[locus] += 1

    ncoal = 0
    for locus, count in root_loci.iteritems():
        ncoal += count-1
    return ncoal


def count_coal_snode(tree, stree, extra, snode,
                     subtrees=None, rootchildren=None, leaves=None,
                     nodefunc=lambda node: node,
                     implied=True):
    
    ncoal_dup = count_coal_snode_dup(tree, stree, extra, snode, subtrees, rootchildren, leaves, nodefunc)
    ncoal_spec = count_coal_snode_spec(tree, stree, extra, snode, subtrees, rootchildren, leaves, nodefunc, implied)
    return ncoal_dup + ncoal_spec
        

def count_coal(tree, stree, extra, subtrees=None, implied=True,
               nodefunc=lambda node: node):
    """Returns the number of inferred extra lineages"""

    subtrees = _subtree_helper(tree, stree, extra, subtrees)

    ncoal = 0
    for snode in stree:
        ncoal += count_coal_snode(tree, stree, extra, snode, subtrees[snode], nodefunc=nodefunc, implied=implied)
    return ncoal



#============================================================================
# duplication loss coal counting
#

init_dup_loss_coal_tree = phyloDLC.init_dup_loss_coal_tree

def count_dup_loss_coal_tree(gene_tree, extra, stree, gene2species,
                             implied=True):
    """count dup loss coal"""

    ndup = 0
    nloss = 0
    ncoal = 0
    nappear = 0

    # use stree to modify internal locus map and order
    new_srecon = util.mapdict(extra["species_map"], val=lambda snode: stree.nodes[snode.name])
    new_order = util.mapdict(extra["order"], key=lambda snode: stree.nodes[snode.name])
    
    srecon = new_srecon
    order = new_order
    extra = extra.copy()
    extra["species_map"] = srecon
    extra["order"] = order
    
    # count appearance
    snode = stree.nodes[srecon[gene_tree.root].name]
    snode.data["appear"] += 1
    nappear += 1

    # factor gene tree
    events = phylo.label_events(gene_tree, srecon)
    subtrees = factor_tree(gene_tree, stree, srecon, events)

    # count events along each species branch
    for snode in stree:
        subtrees_snode = subtrees[snode]
        if len(subtrees_snode) == 0:
            continue

        # count genes
        rootchildren, leaves = _subtree_helper_snode(gene_tree, stree, extra, snode, subtrees_snode)
        if snode.is_leaf():
            snode.data["genes"] += len(leaves)                
        
        # count dups
        ndup_snode = count_dup_snode(gene_tree, stree, extra, snode, subtrees_snode, rootchildren, leaves)
        snode.data["dup"] += ndup_snode
        ndup += ndup_snode

        # count losses
        nloss_snode = count_loss_snode(gene_tree, stree, extra, snode, subtrees_snode, rootchildren, leaves)
        snode.data["loss"] += nloss_snode
        nloss += nloss_snode

        # count deep coalescence (extra lineages)
        ncoal_snode = count_coal_snode(gene_tree, stree, extra, snode, subtrees_snode, rootchildren, leaves, implied=implied)
        snode.data["coal"] += ncoal_snode
        ncoal += ncoal_snode

    return ndup, nloss, ncoal, nappear

count_ancestral_genes = phylo.count_ancestral_genes

def count_dup_loss_coal_trees(gene_trees, extras, stree, gene2species,
                              implied=True):
    """Returns new species tree with dup,loss,coal,appear,genes counts in node's data"""

    stree = stree.copy()
    init_dup_loss_coal_tree(stree)

    for i,gene_tree in enumerate(gene_trees):
        count_dup_loss_coal_tree(gene_tree, extras[i],
	                         stree, gene2species,
                                 implied=implied)
    count_ancestral_genes(stree)
    return stree

