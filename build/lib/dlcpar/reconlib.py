#
# Reconciliation library
#

# python libraries
import copy
import collections
import sys

# rasmus libraries
from rasmus import treelib
from rasmus import util

# compbio libraries
from compbio import phylo, phyloDLC

# dlcpar libraries
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
                    note: this excludes leaf nodes of species branch
                          if leaf node has no children or single child
                          and locus of leaf is same as locus of parent
                          (so keep nodes whose branch has a duplication)

    The gene_tree should contain all implied speciation (and delay) nodes.
    """

    def __init__(self, species_map=None, locus_map=None, order=None):
        self.species_map = species_map
        self.locus_map = locus_map
        self.order = order


    def __eq__(self, other):
        """x.__eq__(y) <==> x == y

        NOTE 1: Internal nodes of the coal_tree and stree must be identical!
        NOTE 2: Data are not compared.
        """

        def error(msg):
            print >>sys.stderr, msg
            return False

        # 1) are species maps identical?
        species_map = util.mapdict(self.species_map, key=lambda node: node.name, val=lambda snode: snode.name)
        other_species_map = util.mapdict(other.species_map, key=lambda node: node.name, val=lambda snode: snode.name)
        if species_map != other_species_map:
            return error("species map mismatch")

        # 2) are locus maps identical?
        locus_map = util.mapdict(self.locus_map, key=lambda node: node.name)
        other_locus_map = util.mapdict(other.locus_map, key=lambda node: node.name)

        # 2a) are number of loci identical?
        if len(locus_map) != len(other_locus_map):
            return error("locus map mismatch")

        # 2b) are loci partitions identical?
        # map loci in self to loci in other
        m = {}
        for name, locus in locus_map.iteritems():
            other_locus = other_locus_map[name]
            if locus in m:
                if m[locus] != other_locus:
                    return error("locus map mismatch")
            else:
                m[locus] = other_locus

        # 3) are orders equal?
        order = util.mapdict(self.order, key=lambda snode: snode.name)
        other_order = util.mapdict(other.order, key=lambda snode: snode.name)

        # 3a) are species equal?
        if set(order) != set(other_order):
            return error("order mismatch")

        for sname in order:
            d = util.mapdict(order[sname], key=lambda locus: m[locus])
            other_d = other_order[sname]

            # 3b) are loci equal?
            if set(d) != set(other_d):
                return error("order mismatch")

            for locus in d:
                lst = map(lambda node: node.name, d[locus])
                other_lst = map(lambda node: node.name, other_d[locus])
                if len(lst) != len(other_lst):
                    return error("order mismatch")

                # 3c) are duplications and partitions identical?
                prev = 0
                for i, name in enumerate(lst):
                    pnode = d[locus][i] # cannot get parent from name only
                    if pnode and locus_map[name] != locus_map[pnode.name]:
                        if sorted(lst[prev:i]) != sorted(other_lst[prev:i]):
                            return error("order mismatch")
                        if lst[i] != other_lst[i]:
                            return error("order mismatch")
                        prev = i + 1
                    if sorted(lst[prev:]) != sorted(other_lst[prev:]):
                        return error("order mismatch")

        # everything identical
        return True


    def __ne__(self, other):
        """x.__ne__(y) <==> x != y"""
        return not self.__eq__(other)


    def sort_loci(self, gtree):
        """Sorts loci in locus_map.

        Loci are sorted using leaf names then by pre-order index.
        """

        # create mapping
        m = {}
        next_locus = 1
        for node in sorted(gtree.leaves(), key=lambda node: node.name):
            locus = self.locus_map[node]
            if locus not in m:
                m[locus] = next_locus
                next_locus += 1
        for node in gtree.preorder():
            if node.is_leaf():
                continue
            locus = self.locus_map[node]
            if locus not in m:
                m[locus] = next_locus
                next_locus += 1

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
              filenames={},
              filestreams={}):
        """Write the reconciliation to a file"""
        assert gtree and self.species_map and self.locus_map and (self.order is not None)

        treefile = filestreams.get("tree", filenames.get("tree", filename + exts["tree"]))
        reconfile = filestreams.get("recon", filenames.get("recon", filename + exts["recon"]))
        orderfile = filestreams.get("order", filenames.get("order", filename + exts["order"]))

        gtree.write(treefile,
                    rootData=True)

        util.write_delim(reconfile,
            [(str(node.name), str(snode.name), self.locus_map[node])
             for node, snode in self.species_map.iteritems()])

        order = {}
        for snode, d in self.order.iteritems():
            for locus, lst in d.iteritems():
                order[snode, locus] = lst
        util.write_delim(orderfile,
            [(str(snode.name), str(locus), ",".join(map(lambda x: str(x.name), lst)))
             for (snode, locus), lst in order.iteritems()])


    def read(self, filename, stree,
             exts={"tree" : ".tree",
                   "recon" : ".recon",
                   "order" : ".order"},
             filenames={}, filestreams={}):
        """Read the reconciliation from a file"""

        treefile = filestreams.get("tree", filenames.get("tree", filename + exts["tree"]))
        reconfile = filestreams.get("recon", filenames.get("recon", filename + exts["recon"]))
        orderfile = filestreams.get("order", filenames.get("order", filename + exts["order"]))

        gtree = treelib.read_tree(treefile)

        self.species_map = {}
        self.locus_map = {}
        for name, sname, locus in util.read_delim(reconfile):
            if name.isdigit(): name = int(name)
            if sname.isdigit(): sname = int(sname)
            assert locus.isdigit()
            locus = int(locus)

            node = gtree.nodes[name]
            self.species_map[node] = stree.nodes[sname]
            self.locus_map[node] = locus

        self.order = collections.defaultdict(dict)
        for toks in util.read_delim(orderfile):
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
                        filenames={},
                        filestreams={}):
    """Writes a labeled reconciliation to files"""

    labeled_recon = LabledRecon(extra["species_map"], extra["locus_map"], extra["order"])
    labeled_recon.write(filename, gtree, exts, filenames, filestreams)


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


def recon_to_labeledrecon(coal_tree, recon, stree, gene2species,
                          name_internal="n", locus_mpr=True,
                          delay=True):
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
    added_spec, added_dup, added_delay = add_implied_nodes(gene_tree, substree,
                                                           species_map, events)
    added = added_spec + added_dup + added_delay

    # rename internal nodes
    common.rename_nodes(gene_tree, name_internal)

    # revert to use input species tree
    species_map = util.mapdict(species_map, val=lambda snode: stree.nodes[snode.name])

    #========================================
    # find locus map

    # label loci in locus tree
    loci = {}
    locus = 1
    dup_snodes = {}     # key = daughter locus, val = snode with duplication
    daughter_locus = {} # key = lnode mapped to duplication, val = daughter locus
    for lnode in locus_tree.preorder():
        pnode = lnode.parent

        if not pnode:               # root
            loci[lnode] = locus
        elif lnode in daughters:    # duplication
            locus += 1
            loci[lnode] = locus

            snode = locus_recon[pnode]
            dup_snodes[locus] = snode
            daughter_locus[pnode] = locus
        else:                       # regular node
            loci[lnode] = loci[pnode]

    # label loci in gene tree
    locus_map = {}
    for cnode in coal_tree:
        node = gene_tree.nodes[cnode.name]
        lnode = coal_recon[cnode]
        locus_map[node] = loci[lnode]
    for node in added:
        # node not in coal tree (because added implied node)
        # use either parent or child locus
        assert node not in locus_map, node.name

        ptr = node
        while ptr not in locus_map:
            ptr = ptr.parent
        locus_up = locus_map[ptr]

        ptr = node
        while ptr not in locus_map:
            assert len(ptr.children) == 1, (ptr.name, ptr.children)
            ptr = ptr.children[0]
        locus_down = locus_map[ptr]

        if locus_up == locus_down:
            # parent and child locus match
            locus_map[node] = locus_up
        else:
            # determine whether to use parent or child locus
            snode = species_map[node]
            dup_snode = dup_snodes[locus_down]
            if (snode.name == dup_snode.name) or (snode.name in dup_snode.descendant_names()):
                locus_map[node] = locus_down
            else:
                locus_map[node] = locus_up

    #========================================
    # find order

    # find order of daughter loci for each locus
    dup_order = {}      # key1 = snode, key2 = plocus, val = lst of daughter loci
    for pnode in locus_tree.preorder():
        if pnode in daughter_locus:
            snode = locus_recon[pnode]
            plocus = loci[pnode]
            locus = daughter_locus[pnode]

            dup_order.setdefault(snode, {})
            dup_order[snode].setdefault(plocus, [])
            dup_order[snode][plocus].append(locus)

    # find loci that give rise to new loci in each sbranch
    parent_loci = set()
    locus_to_node = {}  # key = locus, val = duplicated node with locus
    for node in gene_tree:
        pnode = node.parent
        if pnode:
            locus = locus_map[node]
            plocus = locus_map[pnode]

            if locus != plocus:
                snode = species_map[node]
                parent_loci.add((snode, plocus))
                locus_to_node[locus] = node

    # find nodes for each species and locus
    order = {}
    times = {} # for sorting witin partitions
    for i, node in enumerate(gene_tree.preorder()):
        pnode = node.parent
        times[node] = i

        if pnode:
            snode = species_map[node]
            locus = locus_map[node]
            plocus = locus_map[pnode]

            if (snode, plocus) in parent_loci:
                # skip if same locus as parent and leaf node
                if locus == plocus and (node.is_leaf() or \
                    (len(node.children) == 1 and all([snode != species_map[child] for child in node.children]))):
                    continue

                order.setdefault(snode, {})
                order[snode].setdefault(plocus, [])
                order[snode][plocus].append(node)

    # find duplication order and partitions
    # partition0   dup0   partition1   dup1   partition2 ...
    for snode, d in order.iteritems():
        for plocus, lst in d.iteritems():
            duplicated_loci = dup_order[snode][plocus]
            duporder = map(lambda locus: locus_to_node[locus], duplicated_loci)

            npartitions = len(duporder) + 1
            partitions = [set() for _ in xrange(npartitions)]

            for node in lst:
                if locus_map[node] != plocus:
                    # duplication
                    assert node in duporder, (snode.name, plocus, node.name, duporder)
                    continue

                # walk to nearest locus tree node
                ptr = node
                while ptr.name not in coal_tree.nodes:
                    ptr = ptr.parent
                cnode_up = coal_tree.nodes[ptr.name]
                lnode_up = coal_recon[cnode_up]

                ptr = node
                while ptr.name not in coal_tree.nodes:
                    assert len(ptr.children) == 1, (ptr.name, ptr.children)
                    ptr = ptr.children[0]
                cnode_down = coal_tree.nodes[ptr.name]
                lnode_down = coal_recon[cnode_down]

                assert lnode_up == lnode_down, (node.name, lnode_up, lnode_down)
                lnode = lnode_up

                # find partition
                if lnode in daughter_locus:
                    # lnode mapped to duplication, put in partition before duplication
                    locus = daughter_locus[lnode]
                    ndx = duplicated_loci.index(locus)
                else:
                    # walk up locus tree to duplication
                    while lnode is not None:
                        if locus_recon[lnode] != snode:
                            # lnode in different species, put in first partition
                            ndx = 0
                            break
                        if lnode in daughter_locus:
                            # lnode mapped below duplication, put in partition after duplication
                            locus = daughter_locus[lnode]
                            ndx = duplicated_loci.index(locus) + 1
                            break
                        lnode = lnode.parent

                partitions[ndx].add(node)

            # put duplication order and partitions together
            # also order partitions using temporal constraints
            lst2 = []
            for i in xrange(npartitions-1):
                lst2.extend(sorted(partitions[i], key=lambda node: times[node]))
                lst2.append(duporder[i])
            lst2.extend(sorted(partitions[i+1], key=lambda node: times[node]))

            assert set(lst) == set(lst2), (snode.name, plocus, duporder, lst, lst2)
            d[plocus] = lst2

    #========================================
    # try to remove implied delay nodes

    if not delay:
        for node in added_delay:
            pnode = node.parent
            if locus_map[node] == locus_map[pnode]:
                # delay nodes are leaves of sbranch and so are never in order
                # if the node has the same locus as its parent
                #snode = species_map[node]
                #plocus = locus_map[pnode]
                #if snode in order and plocus in order[snode]:
                #    order[snode][plocus].remove(node)
                del locus_map[node]
                del species_map[node]
                phylo.remove_spec_node(node, gene_tree)
            else:
                raise Exception("Cannot remove implied delay node")

    #========================================
    # put everything together

    return gene_tree, LabeledRecon(species_map, locus_map, order)


def labeledrecon_to_recon(gene_tree, labeled_recon, stree,
                          name_internal="n"):
    """Convert from DLCpar to DLCoal reconciliation model

    NOTE: This is non-reversible because it produces NON-dated coalescent and locus trees
    """

    locus_map = labeled_recon.locus_map
    species_map = labeled_recon.species_map
    order = labeled_recon.order

    # utility function to find longest common substring in a list of strings
    # used to find locus name from a list of coalescent (gene) names
    def get_locus_name(genenames):
        # if only one gene, use gene name as locus name
        if len(genenames) == 1:
            return genenames[0]

        shortest = min(genenames, key=len)

        # find longest common prefix
        prefix = ''
        for i, current_char in enumerate(shortest):
            if any(genename[i] != current_char for genename in genenames):
                break
            prefix = prefix + current_char

        # find longest common suffix
        suffix = ''
        n = len(shortest)-1
        for i, current_char in enumerate(reversed(shortest)):
            if any(genename[n-i] != current_char for genename in genenames):
                break
            suffix = current_char + suffix

        # locus name is two strings combined, with "__" replaced by "_"
        if prefix == '':
            result = suffix
        elif suffix == '':
            result = prefix
        elif prefix[-1] == '_' and suffix[0] == '_':
            result = prefix + suffix[1:]
        else:
            result = prefix + suffix

        # strip leading and trailing "_"
        assert len(result) != 0
        result = result.strip("_")

        return result

    # coalescent tree equals gene tree
    coal_tree = gene_tree.copy()

    # factor gene tree
    events = phylo.label_events(gene_tree, species_map)
    subtrees = factor_tree(gene_tree, stree, species_map, events)

    # 2D dict to keep track of gene names
    # genenames[snode][locus] = list of genes in species and locus
    genenames = {}
    for snode in stree:
        genenames[snode] = collections.defaultdict(list)
    for leaf in gene_tree.leaves():
        genenames[species_map[leaf]][locus_map[leaf]].append(leaf.name)

    # 2D dict to keep track of locus tree nodes by hashing by speciation node and locus
    # locus_tree_map[snode][locus] = list of nodes (sorted from oldest to most recent)
    locus_tree_map = {}
    for snode in stree:
        locus_tree_map[snode] = {}

    # initialize locus tree, coal/locus recon, and daughters
    locus_tree = treelib.Tree()
    coal_recon = {}
    locus_recon = {}
    locus_events = {}
    daughters = []

    # initialize root of locus tree
    root = treelib.TreeNode(locus_tree.new_name())
    locus_tree.add(root)
    locus_tree.root = root
    sroot = species_map[gene_tree.root]
    locus = locus_map[gene_tree.root]
    coal_recon[coal_tree.root] = root
    locus_recon[root] = sroot
    locus_tree_map[sroot][locus] = [root]

    # build locus tree along each species branch
    for snode in stree.preorder(sroot):
        subtrees_snode = subtrees[snode]

        # skip if no branches in this species branch
        if len(subtrees_snode) == 0:
            continue

        # build locus tree
        # 1) speciation
        if snode.parent:
            for (root, rootchild, leaves) in subtrees_snode:
                if rootchild:
                    locus = locus_map[root]     # use root locus!

                    # create new locus tree node in this species branch
                    if locus not in locus_tree_map[snode]:
                        old_node = locus_tree_map[snode.parent][locus][-1]

                        new_node = treelib.TreeNode(locus_tree.new_name())
                        locus_tree.add_child(old_node, new_node)
                        locus_recon[new_node] = snode

                        # update event
                        locus_events[old_node] = "spec"

                        # update nodes in locus tree
                        locus_tree_map[snode][locus] = [new_node]

                    # update coal_recon
                    cnode = coal_tree.nodes[rootchild.name]
                    lnode = locus_tree_map[snode][locus][-1]
                    coal_recon[cnode] = lnode

        # 2) duplication
        if snode in order:
            # may have to reorder loci (in case of multiple duplications)
            queue = collections.deque(order[snode].keys())
            while len(queue) > 0:
                plocus = queue.popleft()
                if plocus not in locus_tree_map[snode]:
                    # parent locus not yet created, punt to future
                    queue.append(plocus)
                    continue

                # handle this ordered list
                lst =  order[snode][plocus]
                for gnode in lst:
                    locus = locus_map[gnode]
                    cnode = coal_tree.nodes[gnode.name]

                    if locus != plocus:     # duplication
                        old_node = locus_tree_map[snode][plocus][-1]

                        # child of duplication in mother locus
                        new_node1 = treelib.TreeNode(locus_tree.new_name())
                        locus_tree.add_child(old_node, new_node1)
                        locus_recon[new_node1] = snode

                        # child of duplication in daughter locus
                        new_node2 = treelib.TreeNode(locus_tree.new_name())
                        locus_tree.add_child(old_node, new_node2)
                        locus_recon[new_node2] = snode
                        daughters.append(new_node2)
                        coal_recon[cnode] = new_node2

                        # update event
                        locus_events[old_node] = "dup"

                        # update nodes in locus tree
                        locus_tree_map[snode][plocus].append(new_node1)
                        locus_tree_map[snode][locus] = [new_node2]

                    else:                   # deep coalescence
                        # update coal_recon
                        lnode = locus_tree_map[snode][locus][-1]
                        coal_recon[cnode] = lnode

        # reconcile remaining coal tree nodes to locus tree
        for (root, rootchild, leaves) in subtrees_snode:
            if rootchild:
                for gnode in gene_tree.preorder(rootchild, is_leaf=lambda x: x in leaves):
                    cnode = coal_tree.nodes[gnode.name]
                    if cnode not in coal_recon:
                        locus = locus_map[gnode]
                        lnode = locus_tree_map[snode][locus][-1]
                        coal_recon[cnode] = lnode

        # tidy up if at an extant species
        if snode.is_leaf():
            for locus, nodes in locus_tree_map[snode].iteritems():
                names = genenames[snode][locus] # list of gene names
                lnode = nodes[-1]

                # find locusname as common substring of list of gene names
                # then relabel genes in locus tree and update event
                locus_name = get_locus_name(names)
                locus_tree.rename(lnode.name, locus_name)
                locus_events[lnode] = "gene"

                # reconcile genes (genes in coal tree reconcile to genes in locus tree)
                # possible mismatch due to genes having an internal ordering even though all exist to present time
                # [could also do a new round of "speciation" at bottom of extant species branches,
                # but this introduces single children nodes that would just be removed anyway]
                cnodes = [coal_tree.nodes[name] for name in names]
                for cnode in cnodes:
                    coal_recon[cnode] = lnode

    # rename internal nodes
    common.rename_nodes(locus_tree, name_internal)

    # simplify coal_tree (and reconciliations)
    removed = treelib.remove_single_children(coal_tree)
    for cnode in removed:
        del coal_recon[cnode]

    # simplify locus_tree (and reconciliations + daughters)
    removed = treelib.remove_single_children(locus_tree)
    for cnode, lnode in coal_recon.items():
        if lnode in removed:
            # reconciliation updates to first child that is not removed
            new_lnode = lnode
            while new_lnode in removed:
                new_lnode = new_lnode.children[0]
            coal_recon[cnode] = new_lnode
    for lnode in removed:
        del locus_recon[lnode]
        del locus_events[lnode]
    for ndx, lnode in enumerate(daughters):
        if lnode in removed:
            # daughter updates to first child that is not removed
            new_lnode = lnode
            while new_lnode in removed:
                new_lnode = new_lnode.children[0]
            daughters[ndx] = new_lnode

##    locus_events = phylo.label_events(locus_tree, locus_recon)
    assert all([lnode in locus_events for lnode in locus_tree])

    #========================================
    # put everything together

    return coal_tree, phyloDLC.Recon(coal_recon, locus_tree, locus_recon, locus_events, daughters)


#=============================================================================
# gene tree factorization functions


def is_full_tree(tree, stree, recon, events):
    """Checks that the tree has all implied internal nodes AND no extra nodes

    Does NOT handle delay nodes
    """

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
                    print >>sys.stderr, "Reconciliation mismatch under speciation-loss node %s" % node.name
                    return False
            elif len(node.children) == 2:
                # speciation
                if set(schildren) != set(schildren2):
                    print >>sys.stderr, "Reconciliation mismatch under speciation node %s" % node.name
                    return False
            else:
                raise Exception("Cannot handle non-binary trees")

        elif events[node] == "dup":
            if len(node.children) == 1:
                # extra node
                print >>sys.stderr, "Single child under duplication node %s" % node.name
                return False
            elif len(node.children) == 2:
                # duplication
                if not (snode == schildren2[0] == schildren2[1]):
                    print >>sys.stderr, "Reconciliation mismatch under duplication node %s" % node.name
                    return False
            else:
                raise Exception("Cannot handle non-binary trees")

        else:
            raise Exception("Invalid event %s at node %s" % (events[node], node.name))

    return True


def add_spec_from_dup_nodes(node, tree, recon, events):
   """
   Relabel the current speciation node 'node' as a duplication.
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
    """
    Add speciation nodes to tree that are implied but are not present because of gene losses.

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
    """
    Insert new delay nodes BELOW gene node 'node'.
    New nodes reconcile to same species node as 'node'.
    Modifies recon and events accordingly.

    TODO: same as add_spec_from_dup_nodes
    """

    assert events[node] == "spec"
    snode = recon[node]
    events[node] = "dup"    # move current node to be internal to species branch

    # insert new nodes into tree
    added_delay = []
    for child in list(node.children):
        added_delay.append(phylo.add_spec_node(child, snode, tree, recon, events))

    return added_delay


def add_implied_delay_nodes(tree, stree, recon, events):
    """
    Add nodes to tree after each speciation to allow
    for delay between coalescence and speciation
    """

    added = []
    for node in list(tree):
        if events[node] == "spec" and len(node.children) > 1:
            added.extend(add_delay_nodes(node, tree, recon, events))
    return added


def add_implied_nodes(tree, stree, recon, events, delay=True):
    """Wrapper to add implied speciation nodes and delay nodes

    If 'delay' is set, then add delay nodes as well.
    """

    added_spec, added_dup = add_implied_spec_nodes(tree, stree, recon, events)
    if delay:
        added_delay = add_implied_delay_nodes(tree, stree, recon, events)
    else:
        added_delay = []
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

    Output is a dict with key = snode, val = (root(subtree), start, leaves(subtree)),
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
    """Returns the subtrees for a species branch"""

    if subtrees is None:
        subtrees = _subtree_helper(tree, stree, extra)
    return subtrees[snode]


def _subtree_helper(tree, stree, extra,
                    subtrees=None):
    """Returns a dictionary of subtrees for each species branch"""

    if subtrees is None:
        recon = extra["species_map"]
        events = phylo.label_events(tree, recon)
        subtrees = factor_tree(tree, stree, recon, events)
    return subtrees


def find_dup_snode(tree, stree, extra, snode,
                   subtrees=None, subtrees_snode=None,
                   nodefunc=lambda node: node):
    """Find duplication nodes (node locus differs from parent node locus)"""

    if subtrees_snode is None:
        subtrees_snode = _subtree_helper_snode(tree, stree, extra, snode, subtrees)
    lrecon = extra["locus_map"]

    dup_nodes = []
    for (root, rootchild, leaves) in subtrees_snode:
        if not rootchild:
            continue

        for node in tree.preorder(rootchild, is_leaf=lambda x: x in leaves):
            pnode = node.parent
            if pnode and lrecon[nodefunc(node)] != lrecon[nodefunc(pnode)]:
                dup_nodes.append(node)
    return dup_nodes


def count_dup_snode(tree, stree, extra, snode,
                    subtrees=None, subtrees_snode=None,
                    nodefunc=lambda node: node):
    """Find number of inferred duplications in a species branch"""

    return len(find_dup_snode(tree, stree, extra, snode,
                              subtrees, subtrees_snode,
                              nodefunc))


def count_dup(tree, stree, extra,
              subtrees=None,
              nodefunc=lambda node: node):
    """Returns the number of inferred duplications"""

    subtrees = _subtree_helper(tree, stree, extra, subtrees)

    ndup = 0
    for snode in stree:
        ndup += count_dup_snode(tree, stree, extra, snode,
                                subtrees, subtrees[snode],
                                nodefunc=nodefunc)
    return ndup


def count_loss_snode(tree, stree, extra, snode,
                     subtrees=None, subtrees_snode=None,
                     nodefunc=lambda node: node):
    """Find number of inferred losses in a species branch"""

    if subtrees_snode is None:
        subtrees_snode = _subtree_helper_snode(tree, stree, extra, snode, subtrees)
    lrecon = extra["locus_map"]

    all_loci = set()
    leaf_loci = set()
    for (root, rootchild, leaves) in subtrees_snode:
        all_loci.add(lrecon[nodefunc(root)])

        if not rootchild:
            continue

        for node in tree.preorder(rootchild, is_leaf=lambda x: x in leaves):
            all_loci.add(lrecon[nodefunc(node)])
            if node in leaves:
                leaf_loci.add(lrecon[nodefunc(node)])
    nloss = len(all_loci.difference(leaf_loci))
    return nloss


def count_loss(tree, stree, extra,
               subtrees=None,
               nodefunc=lambda node: node):
    """Returns the number of inferred losses"""

    subtrees = _subtree_helper(tree, stree, extra, subtrees)

    nloss = 0
    for snode in stree:
        nloss += count_loss_snode(tree, stree, extra, snode,
                                  subtrees, subtrees[snode],
                                  nodefunc=nodefunc)
    return nloss


def count_coal_snode_dup(tree, stree, extra, snode,
                         subtrees=None, subtrees_snode=None,
                         nodefunc=lambda node: node):
    """Returns the number of inferred extra lineages in a species branch
       (at duplication nodes in the locus tree)"""

    if subtrees_snode is None:
        subtrees_snode = _subtree_helper_snode(tree, stree, extra, snode, subtrees)
    srecon, lrecon, order = extra["species_map"], extra["locus_map"], extra["order"]

    ncoal = 0
    if snode not in order:
        return ncoal
    order = order[snode]

    # find all leaves
    all_leaves = []
    for (root, rootchild, leaves) in subtrees_snode:
        if leaves is not None:
            all_leaves.extend(leaves)

    # find "start" branches of each locus
    start = collections.defaultdict(list)
    # find "parent loci" (loci that give rise to new loci) within the species branch
    dup_nodes = find_dup_snode(tree, stree, extra, snode,
                               subtrees, subtrees_snode,
                               nodefunc)
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
        if (node in all_leaves) or (locus not in parent_loci):
            continue
        for child in node.children:
            start[locus].append(child)
    # for each locus that exists at the top of the species branch,
    # if this locus is a "parent locus",
    # add the child if it exists (assume immediate loss if child does not exist)
    for (root, rootchild, leaves) in subtrees_snode:
        locus = lrecon[nodefunc(root)]
        if locus not in parent_loci:
            continue

        # handle root separately
        if not root.parent:
            for child in root.children:       # add all children of root
                if srecon[child] is snode:    # ensure node in sbranch
                    start[locus].append(child)
        else:
            if rootchild:
                start[locus].append(rootchild)
    assert set(start) == set(order), (dict(start), order)

    # for each locus in the species branch, walk down the subtrees using order
    # to determine the number of extra lineages
    for plocus, nodes in order.iteritems():
        current = start[plocus]
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
                # note: keep even if next_node in leaves to allow for delay btwn coalescence and speciation
                #       this special case may not be necessary since leaf nodes no longer in order
                for child in next_node.children:
                    current.append(child)
                    num_lineages += 1
            else:
                # duplication
                if num_lineages > 1:
                    ncoal += num_lineages - 1

    return ncoal


def count_coal_snode_spec(tree, stree, extra, snode,
                          subtrees=None, subtrees_snode=None,
                          nodefunc=lambda node: node,
                          implied=True):
    """Returns the number of inferred extra lineages in a species branch
       (at speciation nodes in the locus tree)"""

    if subtrees_snode is None:
        subtrees_snode = _subtree_helper_snode(tree, stree, extra, snode, subtrees)
    lrecon = extra["locus_map"]

    if not implied:
        raise Exception("not implemented")

    root_loci = {}
    for (root, rootchild, leaves) in subtrees_snode:
        # only count if there is a tree branch in the species branch
        # TODO: if root == rootchild (which only occurs if at root of gene tree
        #       and this root is a duplication), current implementation
        #       inaccurately counts branch(root) as a root lineage,
        #       but this makes no difference since there is a single subtree
        #       so we would have num_lineages = 1 and ncoal = 0
        if rootchild:
            locus = lrecon[nodefunc(root)]
            root_loci.setdefault(locus, 0)
            root_loci[locus] += 1

    ncoal = 0
    for locus, num_lineages in root_loci.iteritems():
        ncoal += num_lineages - 1
    return ncoal


def count_coal_snode(tree, stree, extra, snode,
                     subtrees=None, subtrees_snode=None,
                     nodefunc=lambda node: node,
                     implied=True):
    """Returns the number of inferred extra lineages in a species branch"""

    ncoal_dup = count_coal_snode_dup(tree, stree, extra, snode,
                                     subtrees, subtrees_snode,
                                     nodefunc=nodefunc)
    ncoal_spec = count_coal_snode_spec(tree, stree, extra, snode,
                                       subtrees, subtrees_snode,
                                       nodefunc=nodefunc,
                                       implied=implied)
    return ncoal_dup + ncoal_spec


def count_coal(tree, stree, extra, subtrees=None, implied=True,
               nodefunc=lambda node: node):
    """Returns the number of inferred extra lineages"""

    subtrees = _subtree_helper(tree, stree, extra, subtrees)

    ncoal = 0
    for snode in stree:
        ncoal += count_coal_snode(tree, stree, extra, snode,
                                  subtrees, subtrees[snode],
                                  nodefunc=nodefunc,
                                  implied=implied)
    return ncoal



#============================================================================
# duplication loss coal counting
#

init_dup_loss_coal_tree = phyloDLC.init_dup_loss_coal_tree

def count_dup_loss_coal_tree(gene_tree, extra, stree, gene2species,
                             gene2locus=None, implied=True):
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
        if snode.is_leaf():
            all_leaves = set()
            for (root, rootchild, leaves) in subtrees_snode:
                if leaves is not None:
                    if gene2locus:    # if multiple samples, map genes to locus
                        leaves = set([gene2locus(leaf.name) for leaf in leaves])
                    all_leaves.update(leaves)
            snode.data["genes"] += len(all_leaves)

        # count dups
        ndup_snode = count_dup_snode(gene_tree, stree, extra, snode,
                                     subtrees, subtrees_snode)
        snode.data["dup"] += ndup_snode
        ndup += ndup_snode

        # count losses
        nloss_snode = count_loss_snode(gene_tree, stree, extra, snode,
                                       subtrees, subtrees_snode)
        snode.data["loss"] += nloss_snode
        nloss += nloss_snode

        # count deep coalescence (extra lineages)
        ncoal_snode = count_coal_snode(gene_tree, stree, extra, snode,
                                       subtrees, subtrees_snode,
                                       implied=implied)
        snode.data["coal"] += ncoal_snode
        ncoal += ncoal_snode

    return ndup, nloss, ncoal, nappear

count_ancestral_genes = phylo.count_ancestral_genes

def count_dup_loss_coal_trees(gene_trees, extras, stree, gene2species,
                              gene2locus=None, implied=True):
    """Returns new species tree with dup,loss,coal,appear,genes counts in node's data"""

    stree = stree.copy()
    init_dup_loss_coal_tree(stree)

    for i,gene_tree in enumerate(gene_trees):
        count_dup_loss_coal_tree(gene_tree, extras[i],
                                 stree, gene2species,
                                 gene2locus, implied=implied)
    count_ancestral_genes(stree)
    return stree

