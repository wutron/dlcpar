#!/usr/bin/env python

"""

   Code for the DLC Parsimony Reconciliation
   (duplications, losses, and coalescence)

"""

import random, collections

from rasmus import treelib, util
from compbio import phylo

from yjw import combinatorics

# gtree = G (with implied speciation nodes)
# stree = S
# srecon = M
# lrecon = L


def get_sub_tree(node, snode, recon, events):
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

def factor_tree(gtree, stree, gene2species, recon=None, events=None):
    """Returns subtrees for each species branch

    Output is a dict with key = snode, val = (R(subtree), start, leaves(subtree)),
    where start is the node of the gene tree on which to recur
    (either the root of the subtree or its child).
    """

    if recon is None:
        recon = phylo.reconcile(gtree, stree, gene2species)
    if events is None:
        events = phylo.label_events(gtree, recon)
    
    # initialize
    subtrees = {}
    for snode in stree:
        subtrees[snode] = []

    # root is a dup
    if events[gtree.root] != "spec":
        snode = recon[gtree.root]
        subleaves = get_sub_tree(gtree.root, snode, recon, events)
        subtrees[snode].append((gtree.root, gtree.root, subleaves))

    # find subtrees
    for node in gtree:
        if events[node] == "spec":
            snode = recon[node]
            for schild in snode.children:
                nodes2 = [x for x in node.children if recon[x] == schild]
                if len(nodes2) > 0:
                    assert len(nodes2) == 1, (node, nodes2)
                    node2 = nodes2[0]
                    subleaves = get_sub_tree(node2, schild, recon, events)
                    subtrees[schild].append((node, node2, subleaves))
                else:
                    subtrees[schild].append((node, None, None))

    return subtrees


def find_unique_partitions(leaves, partitions, ids):
    """Returns unique leaf partitions

    Output is a dict with
      key = unique label of leaf nodes,
      val = (labels, mindup)
    where
      labels is the partition label for all nodes in the subtree
      mindup is the min number of dups necessary to attain partition
    
    If multiple partitions have the same mindup, one is chosen randomly.
    """

    # find sort order
    leaf_names = sorted(map(lambda node: node.name, leaves),
                        key=lambda name: ids[name])

    # find unique leaf partitions
    unique = collections.defaultdict(list)
    for p in partitions:
        # unique mapping
        x = map(lambda name: p[name], leaf_names)
        y = []
        m = {}
        next = 0
        for locus in x:
            if locus not in m:
                m[locus] = next
            next += 1
            y.append(m[locus])
        y = tuple(y)

        # internal node mappings
        z = {}
        for (name, locus) in p.iteritems():
            if locus not in m:
                m[locus] = next
            next += 1
            z[name] = m[locus]
        
        # store with dup 
        ndup = len(set(p.values())) - 1
        unique[y].append((z, ndup))

    # out of leaf partitions with mindup,
    # choose a partition for internal nodes randomly
    for y, lst in unique.iteritems():
        mindup = min([x[1] for x in lst])
        lst2 = [x for x in lst if x[1] == mindup]
        unique[y] = random.sample(lst2, 1)[0]
    
##    print leaf_names
##    util.print_dict(unique)
    return unique


def find_locus_partitions_sbranch(sub, ids, is_leaf, K=None, check=True):
    """Finds the valid locus partitions for each subtree in the species branch"""
    
    if K is None:
        K = util.INF
    null_partition = {(): (None, None)}

    # initialize leaf partitions for each subtree
    partitions = {}    
    
    # iterate through subtrees in species branch to find leaf partitions
    # start at rootchild since duplication along (root, rootchild) will be handled when combining subtrees
    for (root, rootchild, leaves) in sub:

        if not rootchild:
            # loss
            print 'T', root, rootchild, leaves
            partitions[root] = null_partition
            continue

        # find maximum number of loci in subtree
        Ksub = len(leaves)
        if K < Ksub:
            Ksub = K

        # initialize partitions storage and rootchild node
        HG = collections.defaultdict(list)
        h = {rootchild.name : 0, "nloci" : 1}
        next = 1
        cur = h[rootchild.name]
        HG[rootchild].append(h)

        print 'T', root, rootchild, leaves, Ksub
        print 'Tnodes', list(gtree.preorder(rootchild, is_leaf=lambda x: x in leaves))

        # recur down subtree to find leaf partitions 
        for node in gtree.preorder(rootchild, is_leaf=lambda x: x in leaves):
            if node in leaves:
                continue

            for h in HG[node]:
                cur = h[node.name]
                children = node.children
                nchildren = len(children)

                if nchildren == 0:
                    raise Exception("invalid number of children: %s" % node.name)
                elif nchildren == 1:
                    child = children[0]
                    h1 = h.copy(); h1[child.name] = cur
                    HG[child].append(h1)
                    if h["nloci"] < Ksub:
                        h2 = h.copy(); h2[child.name] = next; h2["nloci"] += 1
                        next += 1
                        HG[child].append(h2)
                elif nchildren == 2:
                    left, right = children
                    h1 = h.copy(); h1[left.name] = cur;  h1[right.name] = cur
                    HG[left].append(h1)
                    if h["nloci"] < Ksub:
			# daughter lineage matters if subtree has >1 internal node (> 2 leaves)
			# (i.e. if not at rootchild)
			h2 = h.copy(); h2[left.name] = cur;  h2[right.name] = next;   h2["nloci"] += 1
			HG[left].append(h2)
			if node is not rootchild:
                            h3 = h.copy(); h3[left.name] = next; h3[right.name] = cur;    h3["nloci"] += 1
			    HG[left].append(h3)
                        next += 1

                        # duplication along both child lineages incurs a deep coalescence so
			# it is always less parsimonious than dup from parent followed by dup in one child lineage
			# (dup in parent and both children incurs a loss so again is less parsimonious)
			##if h["nloci"] + 1 < Ksub:
                        ##    h4 = h.copy(); h4[left.name] = next; h4[right.name] = next+1; h4["nloci"] += 2
                        ##    next += 1
                        ##    HG[left].append(h4)
                    HG[right] = HG[left]
                else:
                    raise Exception("invalid number of children: %s" % node.name)
            
        # recur up subtree to combine leaf partitions
        IG = collections.defaultdict(list)
        for node in gtree.postorder(rootchild, is_leaf=lambda x: x in leaves):

            # base case
            if node in leaves:
                IG[node] = HG[node]
                continue

            children = node.children
            nchildren = len(children)
        
            if nchildren == 0:
                raise Exception("invalid number of children: %s" % node.name)
            elif nchildren == 1:
                child = children[0]
                IG[node] = IG[child]
            elif nchildren == 2:
                left, right = children
                for hleft in IG[left]:
                    hleftset = set(hleft)
                    for hright in IG[right]:
                        intersect = hleftset.intersection(hright)
                        intersect.remove("nloci")
                        if all([hleft[name] == hright[name] for name in intersect]):
                            h = hleft.copy()
                            h.update(hright)
                            IG[node].append(h)
            else:
                raise Exception("invalid number of children: %s" % node.name)

        print 'IG', rootchild.name, len(IG[rootchild]), IG[rootchild]

        # check size - no easy way to do this with powersets since have to check for dup in both children
        ##if check:
        ##    # Ksub is maximum number of LOCI so there can be [0, Ksub-1] duplications
	##    nleaves = len(leaves)
        ##    nbranches = 2*nleaves - 1
        ##    ##nbranches = nleaves - 1  # also equal to the number of internal nodes
	##    ncombos = combinatorics.num_powerset(nbranches, R=range(Ksub)) - nleaves + 1
	##    assert ncombos == len(IG[rootchild]), (nbranches, Ksub, ncombos, len(IG[rootchild]))
        
        # store using unique leaf labeling
        unique = find_unique_partitions(leaves, IG[rootchild], ids)

        # ensure valid partitions if species branch is leaf branch
        if is_leaf:
            # create leaf pairs
            pairs = []
            for i,leaf1 in enumerate(leaves):
                for leaf2 in leaves[i+1:]:
                    pairs.append((leaf1.name, leaf2.name))

            for p in unique.keys():
               for (leaf1, leaf2) in pairs:
                   if p[leaf1] == p[leaf2]:
                       del unique[p]
        
        # store
        partitions[root] = unique
        
    return partitions


def find_locus_partitions(gtree, stree, gene2species,
                          srecon, sevents=None,
                          check=True):
    """Finds the valid locus partitions"""
    
    # events implied by reconciliation
    if sevents is None:
        sevents = phylo.label_events(gtree, srecon)

    # get node ids
    ids = dict((node.name, gid) for (gid, node) in enumerate(gtree.preorder()))    

    # find maximum number of loci
    recon = phylo.reconcile(gtree, stree, gene2species)
    events = phylo.label_events(gtree, recon)
    K = phylo.count_dup(gtree, events) + 1
    print 'max loci: %d' % K

    # initialize locus partition
    HS = {}

    # find speciation subtrees
    subtrees = factor_tree(gtree, stree, gene2species, srecon, sevents)

    # recur down species tree
    for snode in stree.preorder():
        util.tic(snode)

        # get subtrees in the species branch
        sub = subtrees[snode]

        # get leaf partitions for each subtree in species branch
        partitions = find_locus_partitions_sbranch(sub, ids, snode.is_leaf(), K, check)

        # combine subtrees
        # look at parent partitions (partitions at "top" of sbranch),
        # then use subtrees to determine leaf partitions (partitions at "bottom" of sbranch)
        # HS is a dict with key = snode, val = (partitions, minloci)
        if not snode.parent:
            # root branch
            top = {(0,): None}
        else:
            top = HS[snode.parent]
        
        HS[snode] = {}
        leaf_partitions = [partitions[root].keys() for (root, rootchild, leaves) in sub]
        print 'evolve', top, leaf_partitions
        for partition, cost in top.iteritems():
            # start with partition at top of sbranch
            assert len(partition) == len(leaf_partitions), (len(partition), len(leaf_partitions))
            mult_factor = len(partition)    # dummy to separate partitions in subtrees
            
            for part in combinatorics.combine_lists(*leaf_partitions):
                # find partition at bottom of sbranch
                next = []
                for i, start in enumerate(partition):
                    next.extend([start * mult_factor + leaf for leaf in part[i]])

                # store using unique leaf labeling
		m = {}
		nextid = 0
		y = []
		for locus in next:
		    if locus not in m:
		        m[locus] = nextid
			nextid += 1
		    y.append(m[locus])
		y = tuple(y)
                HS[snode][y] = None
        print 'HS', HS[snode]

        util.toc()           
        
    return HS


if __name__ == "__main__":
##    gtree = treelib.parse_newick("(((a_1,a_2)n3,b2)n2,c2)n1")
##    stree = treelib.parse_newick("((a,b)m2,c)m1")

##    gtree = treelib.parse_newick("((a_1,(b_1,c_1)n3)n2,(b_2,c_2)n4)n1")
    gtree = treelib.parse_newick("(((a_1,(b_1,c_1)n3)n2,(b_2,c_2)n4)n1,a_2)n0")
##    gtree = treelib.parse_newick("((((a_1,(b_1,c_1)n3)n2,(b_2,c_2)n4)n1,a_2)n0,a_3)n10")
    stree = treelib.parse_newick("(a,(b,c)m2)m1")

    gene2species = phylo.make_gene2species([("a_*", "a"), ("b_*", "b"), ("c_*", "c")])

    srecon = phylo.reconcile(gtree, stree, gene2species)
    sevents = phylo.label_events(gtree, srecon)
    phylo.add_implied_spec_nodes(gtree, stree, srecon, sevents)

    treelib.draw_tree_names(gtree, minlen=10, maxlen=10)
    treelib.draw_tree_names(stree, minlen=10, maxlen=10)

##    subtrees = factor_tree(gtree, stree, gene2species, srecon)
##    for snode in stree.preorder():
##        print snode
##        for (root, rootchild, leaves) in subtrees[snode]:
##            print root, rootchild, leaves
##        print
##    print

    util.tic("partitions")
    partitions = find_locus_partitions(gtree, stree, gene2species, srecon)
    util.toc()

