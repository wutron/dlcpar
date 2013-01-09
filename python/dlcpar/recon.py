#!/usr/bin/env python

"""

   Code for the DLC Parsimony Reconciliation
   (duplications, losses, and coalescence)

"""

import sys
import random, collections

from rasmus import treelib, util
from compbio import phylo

from yjw import combinatorics

from dlcpar import common
from dlcpar import reconlib

# gtree = G (with implied speciation nodes)
# stree = S
# srecon = M
# lrecon = L

def dlc_recon(tree, stree, gene2species,
          dupcost=1, losscost=1, coalcost=1, implied=True,
          log=sys.stdout):
    """Perform reconciliation using DLCoal model with parsimony costs"""

    reconer = DLCRecon(tree, stree, gene2species,
                       dupcost=dupcost, losscost=losscost, coalcost=coalcost, implied=implied,
                       log=log)
    return reconer.recon()

    

class DLCRecon(object):

    def __init__(self, gtree, stree, gene2species,
                 dupcost=1, losscost=1, coalcost=1, implied=True,
                 name_internal="n", log=sys.stdout):

        self.gtree = gtree
        common.rename_nodes(self.gtree, name_internal)
        self.stree = stree
        self.gene2species = gene2species
        
        self.dupcost = dupcost
        self.losscost = losscost
        self.coalcost = coalcost
        self.implied = implied

        self.log = log

    #=============================
    # main methods

    def recon(self):
        """Perform reconciliation"""
               
        # find species map
        self._find_species_map()

        # add implied speciation nodes but first start the species tree at the right root
        substree = treelib.subtree(self.stree, self.srecon[self.gtree.root])
        subrecon = util.mapdict(self.srecon, val=lambda snode: substree.nodes[snode.name])
        
        # switch internal storage with subtrees
        self.stree, subtree = substree, self.stree
        self.srecon, subrecon = subrecon, self.srecon

        # add implied speciation nodes
        phylo.add_implied_spec_nodes(self.gtree, self.stree, self.srecon, self.sevents)

        treelib.draw_tree_names(self.gtree, minlen=10, maxlen=10)
        treelib.draw_tree_names(self.stree, minlen=10, maxlen=10)
##        treelib.draw_tree_names(substree, minlen=10, maxlen=10)

        ##    subtrees = self._factor_tree()
        ##    for snode in stree.preorder():
        ##        print snode
        ##        for (root, rootchild, leaves) in subtrees[snode]:
        ##            print root, rootchild, leaves
        ##        print
        ##    print

        self._find_locus_map()

        # revert subtrees
        self.stree, subtree = subtree, self.stree
        self.srecon, subrecon = subrecon, self.srecon

        # convert to LabeledRecon data structure
        labeled_recon = reconlib.LabeledRecon(self.srecon, self.lrecon, self.order)

        return self.gtree, labeled_recon


    def _find_species_map(self):
        self.srecon = phylo.reconcile(self.gtree, self.stree, self.gene2species)
        self.sevents = phylo.label_events(self.gtree, self.srecon)


    def _find_locus_map(self):
        util.tic("find locus map")
        locus_maps = self._enumerate_loci_maps()
        util.toc()

        self._dp_locus_map(locus_maps)

        self.lrecon = dict([(node, 1) for node in self.gtree])
        self.order = dict()
        

    #=============================
    # event/cost methods

    def _count_events(self, loci, snode, subtree):
        rootchildren, leaves = subtree
        extra = {"species_map" : self.srecon, "locus_map" : loci}
        
        ndup = reconlib.count_dup_snode(self.gtree, self.stree, extra, snode, subtree, nodefunc=lambda node: node.name)
        nloss = reconlib.count_loss_snode(self.gtree, self.stree, extra, snode, subtree, use_names=lambda node: node.name)
        ncoal, order = self._count_coal(rootchildren, leaves, loci)
        return ndup, nloss, ncoal, order



    def _count_coal(self, rootchildren, leaves, loci):
        """Find number of inferred extra lineages"""
        
        mincoal = 0
        minorder = None
        
        # try all orderings of internal nodes
        for order in self._enumerate_locus_order(rootchildren, leaves, loci):
            pass
        
        return mincoal, minorder


    def _enumerate_locus_order(self, rootchildren, leaves, loci):
        """Find all internal orderings"""
        
        gtree = self.gtree

        # for each parent locus, find (1) the start and (2) the nodes that have this parent locus
        parent_loci = {}
        for root, rootchild in rootchildren.iteritems():
            # assume immediate loss
            if not rootchild:
                continue

            # handle root
            locus = loci[root.name]
            if locus not in parent_loci:
                parent_loci[locus] = ([], [])
            parent_loci[locus][0].append(root)

            # recur over trees
            for node in gtree.preorder(rootchild, is_leaf=lambda x: x in leaves):
                # loci
                locus = loci[node.name]
                parent = node.parent
                if parent:
                    parent_locus = loci[parent.name]
                    
                # start of parent locus
                if (not parent) or ((parent) and (locus != parent_locus)):
                    if locus not in parent_loci:
                        parent_loci[locus] = ([], [])
                    parent_loci[locus][0].append(node)                        

                # add to parent locus
                if node.parent:
                    parent_loci[parent_locus][1].append(node)

        print 'loci', loci, 'parent_loci', parent_loci

        def walk():
            pass

        for parent_locus, (starts, children) in parent_loci.iteritems():
            pass
           
        return []


    def _compute_cost(self, ndup, nloss, ncoal):
        """Find reconciliation cost"""
        return ndup*self.dupcost + nloss*self.losscost + ncoal*self.coalcost


    #=============================
    # locus methods

    def _find_unique_loci(self, loci, leaves):
        """Returns unique leaf loci"""

        # unique mapping
        x = map(lambda node: loci[node.name], leaves)
        y = []
        m = {}
        next = 1
        for locus in x:
            if locus not in m:
                m[locus] = next
                next += 1
            y.append(m[locus])
        y = tuple(y)

        # internal node mappings
##        z = {}
##        for (name, locus) in loci.iteritems():
##            if locus not in m:
##                m[locus] = next
##                next += 1
##            z[name] = m[locus]
            
        return y, m


    def _evolve_subtree(self, root, leaves, state, start=0, next=1):
        """Given state changes, find the locus at the nodes"""

        gtree = self.gtree
        
        assert next > start, (start, next)
        
        loci = {}
        for node in gtree.preorder(root, is_leaf=lambda x: x in leaves):
            if node is root:
                loci = {node.name : start}
            else:
                if state[node.name]:
                    loci[node.name] = next
                    next += 1
                else:
                    loci[node.name] = loci[node.parent.name]
        return loci


    def _find_locus_states_sbranch(self, sub, ids, is_leaf, K=None, check=True):
        """Finds the valid locus states for each subtree in the species branch"""
        
        gtree = self.gtree
        if K is None:
            K = util.INF

        # storage for states for all subtrees
        # each node of the tree is labeled with T/F corresponding to whether a dup occurred along the branch
        all_states = []

        # iterate through subtrees in species branch to find leaf states
        # start at rootchild since duplication along (root, rootchild) will be handled when combining subtrees
        for (root, rootchild, leaves) in sub:

            if not rootchild:
                # loss
                print 'T', root, rootchild, leaves
                all_states.append([None])
                continue

            # find maximum number of dup in subtree
            # K can be less than Ksub due to speciation nodes 
            Ksub = len(leaves)
            if K < Ksub:
                Ksub = K

            # initialize states storage and rootchild node
            # TODO: more efficient to use matrix to store states
            states_down = collections.defaultdict(list)
            state = {rootchild.name : False}
            states_down[rootchild].append(state)

            print 'T', root, rootchild, leaves, Ksub
            print 'Tnodes', list(gtree.preorder(rootchild, is_leaf=lambda x: x in leaves))

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
                        if ndup < Ksub:
                            s2 = state.copy();  s2[child.name] = True
                            states_down[child].append(s2)
                    elif nchildren == 2:
                        left, right = children
                        s1 = state.copy();      s1[left.name] = False;  s1[right.name] = False
                        states_down[left].append(s1)
                        if ndup < Ksub:
                            # daughter lineage matters
                            s2 = state.copy();  s2[left.name] = False;  s2[right.name] = True
                            states_down[left].append(s2)
                            s3 = state.copy();  s3[left.name] = True;   s3[right.name] = False
                            states_down[left].append(s3)

                            # duplication along both child lineages incurs a deep coalescence so
                            # it is always less parsimonious than dup from parent followed by dup in one child lineage
                            # (dup in parent and both children incurs a loss so again is less parsimonious)
                            ##if ndup + 1 < Ksub:
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

            print 'states_up', rootchild.name, len(states_up[rootchild]), states_up[rootchild]

            # check size - no easy way to do this with powersets since have to check for dup in both children
            ##if check:
            ##    # Ksub is maximum number of LOCI so there can be [0, Ksub-1] duplications
            ##    nleaves = len(leaves)
            ##    nbranches = 2*nleaves - 1
            ##    ##nbranches = nleaves - 1  # also equal to the number of internal nodes
            ##    ncombos = combinatorics.num_powerset(nbranches, R=range(Ksub)) - nleaves + 1
            ##    assert ncombos == len(states_up[rootchild]), (nbranches, Ksub, ncombos, len(states_up[rootchild]))

            # ensure valid states if species branch is leaf branch
            # also, for all species branches, allow duplication along root branch
            states = []       # state for one subtree

            # iterate over states
            for i, state in enumerate(states_up[rootchild][:]):
                valid = True
                
                if is_leaf:
                    # find locus at leaves, then check if valid
                    loci = self._evolve_subtree(rootchild, leaves, state)
                    leaf_loci = [loci[node.name] for node in leaves]
                    valid = len(set(leaf_loci)) == len(leaf_loci)
                    print 'loci', loci, leaf_loci, valid
                        
                if valid:
                    # store original (without duplication along root branch)
                    states.append(state)

                    # allow duplication along root branch
                    if root is not rootchild:
                        ndup = util.counteq(True, state.values())
                        if ndup < Ksub:
                            s1 = state.copy(); s1[rootchild.name] = True
                            states.append(s1)
                        
            # store
            all_states.append(states)
            
        return all_states


    def _enumerate_loci_maps(self, check=True):
        """Finds the valid locus maps"""

        gtree = self.gtree
        stree = self.stree
        gene2species = self.gene2species
        srecon = self.srecon
        sevents = self.sevents
        
        # get node ids
        ids = dict((node.name, gid) for (gid, node) in enumerate(gtree.preorder()))    

        # find maximum number of dup - does not assume that species map is MPR
        recon = phylo.reconcile(gtree, stree, gene2species)
        events = phylo.label_events(gtree, recon)
        K = phylo.count_dup(gtree, events)
        print 'max dup: %d' % K

        # initialize locus partition
        PS = {}
        LS = {}
        if sevents[gtree.root] == "spec":
            sroot = stree.root
        else:
            sroot = None            
        PS[sroot] = {(1,): None}
        LS[sroot] = [gtree.root]

        # find speciation subtrees
        subtrees = reconlib.factor_tree(self.gtree, self.stree,
                                        self.srecon, self.sevents)

        # recur down species tree
        for snode in stree.preorder():
            util.tic(snode)

            # get subtrees in the species branch
            sub = subtrees[snode]

            # nothing happened along the branch if there are no subtrees
            # still have to initialize loci if at root though
            print sub
            if len(sub) == 0:
                util.toc()
                continue

            # find all leaves
            rootchildren = {}
            leaves = []
            for (root, rootchild, myleaves) in sub:
                rootchildren[root] = rootchild
                if myleaves is not None:    # TODO: possible loss!
                    leaves.extend(myleaves)
            leaves.sort(key=lambda node: ids[node.name])
            LS[snode] = leaves

            # sort subtrees by id of root
            sub.sort(key=lambda (root, rootchild, l): ids[root.name])

            # get leaf states for each subtree in species branch
            states = self._find_locus_states_sbranch(sub, ids, snode.is_leaf(), K, check)
            print 'states', states

            # combine subtrees
            # look at parent loci, then use states in subtrees to determine all loci in subtrees
##            if not snode.parent:
##                # root branch
##                top_sbranch = {(1,): None}
##                top_leaves = [gtree.root]
##            else:
            top_sbranch = PS[snode.parent]
            top_leaves = LS[snode.parent]

           # TODO: incorporate K
            PS[snode] = collections.defaultdict(dict)
            print 'evolve', top_sbranch, '-->', states
            for top_loci, _ in top_sbranch.iteritems():
                # start with loci at top of sbranch
                assert len(top_loci) == len(states), (len(top_loci), len(states))

                # initialize next loci with top of sbranch
                init_next = max(top_loci) + 1
                init_loci = {}
                for i, start in enumerate(top_loci):
                    init_loci[top_leaves[i].name] = start
                
                for state in combinatorics.product(*states):
                    next = init_next
                    loci = init_loci.copy()
                    
                    # find next loci using states
    ##                print 'state', state
                    for i, start in enumerate(top_loci):
                        s = state[i]
                        if s is not None:
                            root = top_leaves[i]
                            rootchild = rootchildren[root]

                            # evolve from root to rootchild
    ##                        if root is not rootchild:
                            if s[rootchild.name]:
                                rootchild_loci = next
                                next += 1
                            else:
                                rootchild_loci = start

                            # evolve from rootchild to leaves
    ##                        print 'evolve subtree', root, rootchild, part, rootchild_loci, next
                            l = self._evolve_subtree(rootchild, leaves, s, rootchild_loci, next)
                            assert all([name not in loci for name in l if name != root.name]), (l, loci)
                            loci.update(l)
                            next = max(init_next, max(loci.values())) + 1
    ##                        print 'evolved', start, '-->', l, '-->', l, ':', next
    ##                print 'top_loci', top_loci, 'parts', parts, 'loci', loci

                    # if leaf, see if valid
                    if snode.is_leaf():
                        leaf_loci = [loci[node.name] for node in leaves]
                        if not len(set(leaf_loci)) == len(leaf_loci):   # invalid
                            continue

                    # find cost
                    ndup, nloss, ncoal, order = self._count_events(loci, snode, (rootchildren, leaves))
                    cost = self._compute_cost(ndup, nloss, ncoal)
                    print 'cost', rootchildren, leaves, loci, ndup, nloss, ncoal, cost

                    # (unique) loci at bottom of sbranch
                    bottom_loci, mapping = self._find_unique_loci(loci, leaves)
    ##                print 'bottom_loci', bottom_loci, mapping

                    if top_loci not in PS[snode][bottom_loci]:
                        PS[snode][bottom_loci][top_loci] = []
                    PS[snode][bottom_loci][top_loci].append((loci, mapping, ndup, cost))

            print 'LS', LS[snode.parent] if snode.parent else None, LS[snode]
            print 'PS'
            for bottom_loci, d in PS[snode].iteritems():
                for top_loci, lst in d.iteritems():
                    print bottom_loci, '<--', top_loci, ':', lst

            # TODO: pick optimal cost for each unique top -> bottom
            for bottom_loci, d in PS[snode].iteritems():
                for top_loci, lst in d.iteritems():
                    its, mincost = util.minall(lst, minfunc=lambda (loci, mapping, ndup, cost): cost)

                    # if multiple loci achive the minimum cost,
                    # choose one randomly from the set with min dup (this allows more locis down the tree)
                    if len(its) > 1:
                        its, mindup = util.minall(its, minfunc = lambda (loci, mapping, ndup, cost): ndup)
                    it = random.choice(its)

                    PS[snode][bottom_loci][top_loci] = it

            # unique
            print 'unique'
            for bottom_loci, d in PS[snode].iteritems():
                for top_loci, lst in d.iteritems():
                    print bottom_loci, '<--', top_loci, ':', lst

            util.toc()           
            
        return PS


    def _dp_locus_map(self, locus_map):
        """Find optimal locus map through DP"""
        pass


if __name__ == "__main__":
    gene_tree = treelib.parse_newick("(((a_1,a_2)n3,b_2)n2,c_2)n1")
    species_tree = treelib.parse_newick("((a,b)m2,c)m1")

##    gene_tree = treelib.parse_newick("((a_1,(b_1,c_1)n3)n2,(b_2,c_2)n4)n1")
#####    gene_tree = treelib.parse_newick("(((a_1,(b_1,c_1)n3)n2,(b_2,c_2)n4)n1,a_2)n0")
##    gene_tree = treelib.parse_newick("((((a_1,(b_1,c_1)n3)n2,(b_2,c_2)n4)n1,a_2)n0,a_3)n10")
##    species_tree = treelib.parse_newick("(a,(b,c)m2)m1")
##    species_tree = treelib.parse_newick("((a,(b,c)m2)m1,d)m0")
#####    species_tree = treelib.parse_newick("(a,((b,c)m2,(d,e)m11)m10)m1")

    species_map = phylo.make_gene2species([("a_*", "a"), ("b_*", "b"), ("c_*", "c"), ("d_*", "d")])

    dlc_recon(gene_tree, species_tree, species_map)
    
