#
# Reconciliation library for dup/loss/deep coalescence
# Based on three-tree model of DLCoal
#

# python libraries
import sys
import copy

# rasmus libraries
from rasmus import util
from rasmus import treelib

# compbio libraries
from compbio import phylo, coal

#=============================================================================
# Reconciliation Data Structures

class Recon (object):
    """The reconciliation data structure for the DLCoal model"""

    def __init__(self, coal_recon=None,
                 locus_tree=None, locus_recon=None, locus_events=None,
                 daughters=None, data=None):
        if coal_recon:
            self.coal_recon = coal_recon
        else:
            self.coal_recon = {}

        if locus_tree:
            self.locus_tree = locus_tree
        else:
            self.locus_tree = {}

        if locus_recon:
            self.locus_recon = locus_recon
        else:
            self.locus_recon = {}

        if locus_events:
            self.locus_events = locus_events
        else:
            self.locus_events = {}

        if daughters:
            self.daughters = daughters
        else:
            self.daughters = set()

        if data:
            self.data = data
        else:
            self.data = {}


    def copy(self):
        return Recon(self.coal_recon,
                     self.locus_tree, self.locus_recon, self.locus_events,
                     self.daughters, data=copy.deepcopy(self.data))


    def get_dict(self):
        return {"coal_recon": self.coal_recon,
                "locus_tree": self.locus_tree,
                "locus_recon": self.locus_recon,
                "locus_events": self.locus_events,
                "daughters": self.daughters}


    def __repr__(self):
        return repr({"coal_recon": [(x.name, y.name) for x,y in
                                    self.coal_recon.iteritems()],
                     "locus_tree": self.locus_tree.get_one_line_newick(
                         root_data=True),
                     "locus_top": phylo.hash_tree(self.locus_tree),
                     "locus_recon": [(x.name, y.name) for x,y in
                                     self.locus_recon.iteritems()],
                     "locus_events": [(x.name, y) for x,y in
                                      self.locus_events.iteritems()],
                     "daughters": [x.name for x in self.daughters],
                     "data": self.data})


    def __eq__(self, other):
        """x.__eq__(y) <==> x==y

        NOTE 1: Only the locus tree node names are allowed to change.
                (Internal nodes of coal_tree and stree must be identical!)
        NOTE 2: Data are not compared.
        """

        def error(msg):
            print >>sys.stderr, msg
            return False

        # are locus_trees identical?
        if phylo.hash_tree(self.locus_tree) != phylo.hash_tree(other.locus_tree):
            return error("locus_tree mismatch")

        # map nodes using leaf names -- TODO : more efficient to use hash_tree_names to map?
        def get_leaf_dct(tree):
            """return dict with key=leaves, val=node"""
            leaves = {}
            for node in tree.postorder():
                if node.is_leaf():
                    leaves[node] = [node.name]
                else:
                    leaves[node] = []
                    for child in node.children:
                        leaves[node].extend(leaves[child])

            dct = {}
            for node in tree:
                dct[tuple(sorted(leaves[node]))] = node
            return dct
        def get_map(tree1, tree2):
            """remap tree1 nodes into tree2 nodes"""
            m = {}
            tree1_dict = get_leaf_dct(tree1)
            tree2_dict = get_leaf_dct(tree2)
            for leaves, node in tree1_dict.iteritems():
                m[node] = tree2_dict[leaves]
            return m
        locus_map = get_map(self.locus_tree, other.locus_tree)

        # are locus_recon identical?
        locus_recon = util.mapdict(self.locus_recon, key=lambda lnode: locus_map[lnode].name, val=lambda snode: snode.name)
        other_locus_recon = util.mapdict(other.locus_recon, key=lambda lnode: lnode.name, val=lambda snode: snode.name)
        if locus_recon != other_locus_recon:
            return error("locus_recon mismatch")

        # are locus_events identical?
        locus_events = util.mapdict(self.locus_events, key=lambda lnode: locus_map[lnode].name)
        other_locus_events = util.mapdict(other.locus_events, key=lambda lnode: lnode.name)
        if locus_events != other_locus_events:
            return error("locus_events mismatch")

        # are daughters identical?
        daughters = set([locus_map[lnode].name for lnode in self.daughters])
        other_daughters = set([lnode.name for lnode in other.daughters])
        if daughters != other_daughters:
            return error("daughters mismatch")

        # are coal_recon identical?
        coal_recon = util.mapdict(self.coal_recon, key=lambda node: node.name, val=lambda lnode: locus_map[lnode].name)
        other_coal_recon = util.mapdict(other.coal_recon, key=lambda node: node.name, val=lambda lnode: lnode.name)
        if coal_recon != other_coal_recon:
            return error("coal_recon mismatch")

        # everything identical
        return True


    def __ne__(self, other):
        """x.__ne__(y) <==> x != y"""
        return not self.__eq__(other)


    def write(self, filename, coal_tree,
              exts={"coal_tree": ".coal.tree",
                    "coal_recon": ".coal.recon",
                    "locus_tree": ".locus.tree",
                    "locus_recon": ".locus.recon",
                    "daughters": ".daughters"
                    },
              filenames={},
              filestreams={}):
        """Writes a reconciled gene tree to files"""

        assert coal_tree and self.coal_recon and \
               self.locus_tree and self.locus_recon and self.locus_events and \
               (self.daughters is not None)

        # coal tree and recon
        coal_tree.write(
            filestreams.get("coal_tree", filenames.get("coal_tree", filename + exts["coal_tree"])),
            rootData=True)
        phylo.write_recon_events(
            filestreams.get("coal_recon", filenames.get("coal_recon", filename + exts["coal_recon"])),
            self.coal_recon, noevent="none")

        # locus tree and recon
        self.locus_tree.write(
            filestreams.get("locus_tree", filenames.get("locus_tree", filename + exts["locus_tree"])),
            rootData=True)
        phylo.write_recon_events(
            filestreams.get("locus_recon", filenames.get("locus_recon", filename + exts["locus_recon"])),
            self.locus_recon, self.locus_events)

        # daughters
        util.write_list(
            filestreams.get("daughters", filenames.get("daughters", filename + exts["daughters"])),
            [x.name for x in self.daughters])


    def read(self, filename, stree,
             exts={"coal_tree": ".coal.tree",
                   "coal_recon": ".coal.recon",
                   "locus_tree": ".locus.tree",
                   "locus_recon": ".locus.recon",
                   "daughters": ".daughters"
                   },
             filenames={},
             check=True):
        """Reads a reconciled gene tree from files"""

        # trees
        coal_tree = treelib.read_tree(
            filenames.get("coal_tree", filename + exts["coal_tree"]))
        self.locus_tree = treelib.read_tree(
            filenames.get("locus_tree", filename + exts["locus_tree"]))

        # recons
        self.coal_recon, junk = phylo.read_recon_events(
            filenames.get("coal_recon", filename + exts["coal_recon"]),
            coal_tree, self.locus_tree)
        self.locus_recon, self.locus_events = phylo.read_recon_events(
            filenames.get("locus_recon", filename + exts["locus_recon"]),
            self.locus_tree, stree)

        self.daughters = set(
            self.locus_tree.nodes[x] for x in util.read_strings(
            filenames.get("daughters", filename + exts["daughters"])))

        assert (not check) or (check and self.is_valid(coal_tree))

        return coal_tree, self.get_dict()


    def is_valid(self, coal_tree):
        if not assert_daughters(self.locus_events, self.daughters):
            print >>sys.stderr, "Locus events and daughters do not match"
            return False
        if not assert_bounded_coal(coal_tree, self.coal_recon, self.locus_tree, self.daughters):
            print >>sys.stderr, "Bounded coalescent violated"
            return False
        return True


#=============================================================================
# Reconciliation Functions

find_loss = phylo.find_loss
count_dup = phylo.count_dup
count_loss = phylo.count_loss

def assert_daughters(locus_events, daughters):
    """Ensure duplications and daughters match"""

    dup_nodes = []
    for node, event in locus_events.iteritems():
        if event == "dup":
            dup_nodes.append(node)

    if len(dup_nodes) != len(daughters):
        return False

    for node in dup_nodes:
        assert len(node.children) == 2
        if (node.children[0] not in daughters) and \
           (node.children[1] not in daughters):
            return False

    return True


def assert_bounded_coal_lca(coal_tree, coal_recon, locus_tree, daughters):
    """Ensure bounded coalescent for daughter lineages (using lca mapping).
       Consider the lca (in the coal tree) of the subtrees (of the daughter lineage)
       If this lca reconciles above the daughter branch, a violation occurred"""

    for dnode in daughters:
        cnodes = map(lambda name: coal_tree.nodes[name], dnode.leaf_names())
        lca = treelib.lca(cnodes)
        lnode = coal_recon[lca]

        def walk(node):
            if node is dnode:
                return False
            else:
                if node.parent:
                    return walk(node.parent)
            return True
        violated = walk(lnode)

        if violated:
            return False
    return True


def assert_bounded_coal_lineages(coal_tree, coal_recon, locus_tree, daughters):
    """Ensure bounded coalescent for daughter lineages (using lineage counts)."""

    lineages = coal.count_lineages_per_branch(coal_tree, coal_recon, locus_tree)

    for dnode in daughters:
        if lineages[dnode][1] != 1:
            return False
    return True


assert_bounded_coal = assert_bounded_coal_lineages


def count_coal(tree, recon, stree):
    """Count the number of (deep) coalescences in a gene tree"""
    counts = coal.count_lineages_per_branch(tree, recon, stree)
    ncoal = sum(max(x[1]-1, 0) for x in counts.itervalues())
    return ncoal

def count_dup_loss_coal(coal_tree, extra, stree, implied=True):
    """Returns the number of duplications + transfers + losses in a gene tree"""

    locus_tree = extra["locus_tree"]
    locus_recon = extra["locus_recon"]
    locus_events = extra["locus_events"]
    coal_recon = extra["coal_recon"]

    ndup = count_dup(locus_tree, locus_events)
    nloss = count_loss(locus_tree, stree, locus_recon)

    if implied:
        # add implied speciation nodes if desired
        # this must be added AFTER counting dups and losses since it affects loss inference
        added = phylo.add_implied_spec_nodes(locus_tree, stree, locus_recon, locus_events)
    ncoal = count_coal(coal_tree, coal_recon, locus_tree)
    if implied:
        phylo.remove_implied_spec_nodes(locus_tree, added, locus_recon, locus_events)

    return ndup + nloss + ncoal


#=============================================================================
# Convenient Reconciliation Input/Output

def write_dlcoal_recon(filename, coal_tree, extra,
                       exts={"coal_tree": ".coal.tree",
                             "coal_recon": ".coal.recon",
                             "locus_tree": ".locus.tree",
                             "locus_recon": ".locus.recon",
                             "daughters": ".daughters"
                            },
                       filenames={},
                       filestreams={}):
    """Writes a reconciled gene tree to files"""

    recon = Recon(extra["coal_recon"], extra["locus_tree"], extra["locus_recon"], extra["locus_events"],
                  extra["daughters"])
    recon.write(filename, coal_tree, exts, filenames, filestreams)

def read_dlcoal_recon(filename, stree,
                      exts={"coal_tree": ".coal.tree",
                            "coal_recon": ".coal.recon",
                            "locus_tree": ".locus.tree",
                            "locus_recon": ".locus.recon",
                            "daughters": ".daughters"
                           },
                      filenames={},
                      check=True):
    """Reads a reconciled gene tree from files"""

    recon = Recon()
    return recon.read(filename, stree,
                      exts, filenames,
                      check=check)

#============================================================================
# duplication loss coal counting
#

def init_dup_loss_coal_tree(stree):
    """initalize counts to zero"""

    def walk(node):
        node.data['dup'] = 0
        node.data['loss'] = 0
        node.data['coal'] = 0
        node.data['appear'] = 0
        node.data['genes'] = 0
        node.recurse(walk)
    walk(stree.root)

def count_dup_loss_coal_tree(coal_tree, extra, stree, gene2species,
                             implied=True, locus_mpr=True):
    """count dup loss coal"""

    if not locus_mpr:
        raise Exception("not implemented")

    # TODO: use locus_recon and locus_events rather than MPR
    #       (currently, phylo.py reconciliation functions fail for non-MPR)
    locus_tree = extra["locus_tree"]
    locus_recon = phylo.reconcile(locus_tree, stree, gene2species)
    locus_events = phylo.label_events(locus_tree, locus_recon)
    coal_recon = extra["coal_recon"]

    ndup, nloss, nappear = phylo.count_dup_loss_tree(locus_tree, stree, gene2species,
                                                     locus_recon, locus_events)

    # add implied speciation nodes if desired
    # this must be added AFTER counting dups and losses since it affects loss inference
    if implied:
        added = phylo.add_implied_spec_nodes(locus_tree, stree, locus_recon, locus_events)

    # count coals
    ncoal = 0
    counts = coal.count_lineages_per_branch(coal_tree, coal_recon, locus_tree)
    for lnode, (count_bot, count_top) in counts.iteritems():
        n = max(count_top-1, 0)
        locus_recon[lnode].data['coal'] += n
        ncoal += n

    if implied:
        phylo.remove_implied_spec_nodes(locus_tree, added, locus_recon, locus_events)

    return ndup, nloss, ncoal, nappear

count_ancestral_genes = phylo.count_ancestral_genes

def count_dup_loss_coal_trees(coal_trees, extras, stree, gene2species,
                              implied=True, locus_mpr=True):
    """Returns new species tree with dup,loss,coal,appear,genes counts in node's data"""

    if not locus_mpr:
        raise Exception("not implemented")

    stree = stree.copy()
    init_dup_loss_coal_tree(stree)

    for i,coal_tree in enumerate(coal_trees):
        # copy locus_recon - must do since stree has been copied
        # can skip since we assume MPR rather than using locus_recon
##        locus_recon = extras[i]["locus_recon"]
##        locus_recon_copy = {}
##        for lnode, snode in locus_recon.iteritems():
##            locus_recon_copy[lnode] = stree.nodes[snode.name]
##        extras[i]["locus_recon"] = locus_recon_copy

        count_dup_loss_coal_tree(coal_tree, extras[i],
                                 stree, gene2species,
                                 implied=implied)
    count_ancestral_genes(stree)
    return stree

