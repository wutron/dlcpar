"""
phylo.py
Phylogeny classes and functions
"""

# python libraries
import sys

# rasmus libraries
from rasmus import stats
from rasmus import treelib
from rasmus import util


#=============================================================================
# counting functions

def num_rooted_trees(nleaves):
    """Return number of rooted trees

    (2n-3)!! = (2n-3)!/[2^(n-2)*(n-2)!] for n >= 2
    """
    assert nleaves >= 1
    if nleaves < 2:
        return 1
    return stats.factorial(2*nleaves-3)/(2**(nleaves-2) * stats.factorial(nleaves-2))

def num_unrooted_trees(nleaves):
    """Return number of unrooted trees

    (2n-5)!! = (2n-5)!/[2^(n-3)*(n-3)!] for n >= 3
    """
    assert nleaves >= 1
    if nleaves < 3:
        return 1
    return stats.factorial(2*nleaves-5)/(2**(nleaves-3) * stats.factorial(nleaves-3))


#=============================================================================
# gene to species mapping functions


def gene2species(genename):
    """Default gene2species mapping"""
    return genename


def make_gene2species(maps):
    """Returns a function that maps gene names to species names

    maps -- a list of tuples [(gene_pattern, species_name), ... ]
    """

    # find exact matches and expressions
    exacts = {}
    exps = []
    for mapping in maps:
        if "*" not in mapping[0]:
            exacts[mapping[0]] = mapping[1]
        else:
            exps.append(mapping)

    # create mapping function
    def gene2species_helper(gene):
        """helper function"""
        # eval expressions first in order of appearance
        for exp, species in exps:
            if exp[-1] == "*":
                if gene.startswith(exp[:-1]):
                    return species
            elif exp[0] == "*":
                if gene.endswith(exp[1:]):
                    return species

        if gene in exacts:
            return exacts[gene]

        raise Exception("Cannot map gene '%s' to any species" % gene)
    return gene2species_helper


def read_gene2species(*filenames):
    """
    Reads a gene2species file

    Returns a function that will map gene names to species names.
    """

    for filename in filenames:
        maps = []
        maps.extend(util.read_delim(util.skip_comments(
            util.open_stream(filename))))
    return make_gene2species(maps)


#=============================================================================
# reconciliation functions

def reconcile(gtree, stree, gene2species=gene2species):
    """
    Returns a reconciliation dict for a gene tree 'gtree' and species tree 'stree'
    """

    recon = {}

    # determine the preorder traversal of the stree
    order = {}
    def walk(node):
        """helper function"""
        order[node] = len(order)
        node.recurse(walk)
    walk(stree.root)

    # label gene leaves with their species
    for node in gtree.leaves():
        recon[node] = stree.nodes[gene2species(node.name)]

    # recur through gene tree
    def walk2(node):
        """helper function"""
        node.recurse(walk2)

        if not node.is_leaf():
            # this node's species is lca of children species
            child_snodes = [recon[child] for child in node.children]
            recon[node] = _reconcile_lca(order, child_snodes)
    walk2(gtree.root)

    return recon


def _reconcile_lca(order, nodes):
    """Helper function for reconcile"""

    # handle simple and complex cases
    if len(nodes) == 1:
        return nodes[0]
    if len(nodes) > 2:
        return treelib.lca(nodes)

    # 2 node case
    node1, node2 = nodes
    index1 = order[node1]
    index2 = order[node2]

    while index1 != index2:
        if index1 > index2:
            node1 = node1.parent
            index1 = order[node1]
        else:
            node2 = node2.parent
            index2 = order[node2]
    return node1


def reconcile_node(node, stree, recon):
    """Reconcile a single gene node to a species node"""
    snodes = [recon[x] for x in node.children]
    return treelib.lca(snodes)


def assert_recon(tree, stree, recon):
    """Assert that a reconciliation is valid"""

    def below(node1, node2):
        """Return True if node1 is below node2"""
        while node1:
            if node1 == node2:
                return True
            node1 = node1.parent
        return False

    for node in tree:
        # Every node in gene tree should be in reconciliation
        assert node in recon

        # Every node should map to a species node equal to or
        # below their parent's mapping
        if node.parent:
            snode = recon[node]
            parent_snode = recon[node.parent]
            assert below(snode, parent_snode)


def label_events(gtree, recon):
    """Returns a dict with gene node keys and values indicating
    'gene', 'spec', or 'dup'"""
    events = {}

    def walk(node):
        """helper function"""
        events[node] = label_events_node(node, recon)
        node.recurse(walk)
    walk(gtree.root)

    return events


def label_events_node(node, recon):
    """Returns event as 'gene', 'spec', or 'dup'"""
    if not node.is_leaf():
        child_snodes = [recon[child] for child in node.children]
        if recon[node] in child_snodes:
            return "dup"
        return "spec"
    return "gene"


def find_loss_node(node, recon):
    """Finds the loss events for a branch in a reconciled gene tree"""
    loss = []

    # if not parent, then no losses
    if not node.parent:
        return loss

    # determine starting and ending species
    sstart = recon[node]
    send = recon[node.parent]

    # determine species path of this gene branch (node, node.parent)
    ptr = sstart
    spath = []
    while ptr != send:
        spath.append(ptr)
        ptr = ptr.parent

    # determine whether node.parent is a dup
    # if so, send (species end) is part of species path
    if label_events_node(node.parent, recon) == "dup":
        spath.append(send)

    # go up species path (skip starting species)
    # every node on the list is at least one loss
    for i, snode in enumerate(spath[1:]):
        for schild in snode.children:
            if schild != spath[i]:
                loss.append([node, schild])

    return loss


def find_loss_under_node(node, recon):
    """Finds the loss events for a branch in a reconciled gene tree"""
    loss = []
    snodes = {}
    internal = {}
    species1 = recon[node]

    # walk from child species to parent species
    for child in node.children:
        ptr = recon[child]
        snodes[ptr] = 1
        while ptr != species1:
            ptr = ptr.parent
            snodes[ptr] = 1
            internal[ptr] = 1

    # foreach internal node in partial speciation tree, all children
    # not in speciation are loss events
    for i in internal:
        for child in i.children:
            if child not in snodes:
                loss.append([node, child])
    return loss


def find_loss(gtree, recon, node=None):
    """Returns a list of gene losses in a gene tree

    # TODO: generalize to non-MPR recon
            (in particular, to handle duplication followed immediately by loss)
    """
    loss = []

    def walk(node):
        """helper function"""
        loss.extend(find_loss_node(node, recon))

        # add losses (for non-MPR)
        #snode = recon[node]
        #child_snodes = set()
        #for child in node.children:
        #    child_snodes.update([recon[child]] + recon[child].descendants())
        #
        #for child_snode in snode.children:
        #    if child_snode not in child_snodes:
        #        loss.append([node, child_snode])

        node.recurse(walk)
    if node:
        walk(node)
    else:
        walk(gtree.root)

    return loss


def count_dup(gtree, events, node=None):
    """Returns the number of duplications in a gene tree

    TODO: generalize to non-MPR events
          (in particular, to handle duplication followed immediately by loss)
    """
    var = {"dups": 0}

    def walk(node):
        """helper function"""
        if events[node] == "dup":
            var["dups"] += len(node.children) - 1
        node.recurse(walk)
    if node:
        walk(node)
    else:
        walk(gtree.root)

    return var["dups"]


def count_loss(gtree, recon, node=None):
    """Returns the number of losses in a gene tree"""
    return len(find_loss(gtree, recon, node))


def count_dup_loss(gtree, stree, recon, events=None):
    """Returns the number of duplications + losses in a gene tree"""
    if events is None:
        events = label_events(gtree, recon)

    nloss = count_loss(gtree, recon)
    ndups = count_dup(gtree, events)
    return nloss + ndups


def find_species_roots(tree, stree, recon):
    """Find speciation nodes in the gene tree that reconcile to the
       species tree root"""

    roots = []
    def walk(node):
        """helper function"""
        found = False
        for child in node.children:
            found = walk(child) or found
        if not found and recon[node] == stree.root:
            roots.append(node)
            found = True
        return found
    walk(tree.root)
    return roots


def find_orthologs(gtree, stree, recon, events=None,
                   counts=True, species_branch=False):
    """Find all ortholog pairs within a gene tree"""

    if events is None:
        events = label_events(gtree, recon)
    orths = []

    for node, event in events.items():
        if event == "spec":
            leavesmat = [x.leaves() for x in node.children]
            sp_counts = [util.hist_dict([recon[i] for i in row])
                         for row in leavesmat]

            for i in xrange(len(leavesmat)):
                for j in xrange(i+1, len(leavesmat)):
                    for gene1 in leavesmat[i]:
                        for gene2 in leavesmat[j]:
                            if gene1.name > gene2.name:
                                g1, g2 = gene2, gene1
                                a, b = j, i
                            else:
                                g1, g2 = gene1, gene2
                                a, b = i, j

                            orth = [g1.name, g2.name]
                            if counts:
                                orth.extend([sp_counts[a][recon[g1]],
                                             sp_counts[b][recon[g2]]])
                            if species_branch:
                                orth.append(recon[node])
                            orths.append(tuple(orth))

    return orths


def find_paralogs(gtree, stree, recon, events=None,
                  counts=True, species_branch=False,
                  split=True):
    """Find all paralog pairs within a gene tree

    Same as find_orthologs but looks for event == "dup".
    """

    if events is None:
        events = label_events(gtree, recon)
    paralogs = []

    for node, event in events.items():
        if event == "dup":
            # ignore pre-root duplications
            if split and recon[node] == stree.root:
                continue

            leavesmat = [x.leaves() for x in node.children]
            sp_counts = [util.hist_dict([recon[i] for i in row])
                         for row in leavesmat]

            for i in xrange(len(leavesmat)):
                for j in xrange(i+1, len(leavesmat)):
                    for gene1 in leavesmat[i]:
                        for gene2 in leavesmat[j]:
                            if gene1.name > gene2.name:
                                g1, g2 = gene2, gene1
                                a, b = j, i
                            else:
                                g1, g2 = gene1, gene2
                                a, b = i, j

                            paralog = [g1.name, g2.name]
                            if counts:
                                paralog.extend([sp_counts[a][recon[g1]],
                                                sp_counts[b][recon[g2]]])
                            if species_branch:
                                paralog.append(recon[node])
                            paralogs.append(tuple(paralog))

    return paralogs


def subset_recon(tree, recon, events=None):
    """Ensure the reconciliation only refers to nodes in tree"""

    # get all nodes that are walkable
    nodes = set(tree.postorder())
    for node in list(recon):
        if node not in nodes:
            del recon[node]
    if events:
        for node in list(events):
            if node not in nodes:
                del events[node]


#=============================================================================
# reconciliation input/output

def write_recon(filename, recon):
    """Write a reconciliation to a file"""
    util.write_delim(filename, [(str(a.name), str(b.name))
                                for a, b in recon.items()])


def read_recon(filename, tree1, tree2):
    """Read a reconciliation from a file"""
    recon = {}
    for a, b in util.read_delim(filename):
        if a.isdigit(): a = int(a)
        if b.isdigit(): b = int(b)
        recon[tree1.nodes[a]] = tree2.nodes[b]
    return recon


def write_events(filename, events):
    """Write events data structure to file"""
    util.write_delim(filename, [(str(a.name), b)
                                for a, b in events.items()])


def read_events(filename, tree):
    """Read events data structure from file"""
    events = {}
    for name, event in util.read_delim(filename):
        if name.isdigit(): name = int(name)
        events[tree.nodes[name]] = event
    return events


def write_recon_events(filename, recon, events=None, noevent=""):
    """Write a reconciliation and events to a file"""

    if events is None:
        events = dict.fromkeys(recon.keys(), noevent)

    util.write_delim(filename, [(str(a.name), str(b.name), events[a])
                                for a, b in recon.items()])


def read_recon_events(filename, tree1, tree2):
    """Read a reconciliation and events data structure from file"""

    recon = {}
    events = {}
    for a, b, event in util.read_delim(filename):
        if a.isdigit(): a = int(a)
        if b.isdigit(): b = int(b)
        node1 = tree1.nodes[a]
        recon[node1] = tree2.nodes[b]
        events[node1] = event
    return recon, events


#============================================================================
# duplication loss counting

def init_dup_loss_tree(stree):
    """initialize counts to zero"""

    for node in stree:
        node.data['dup'] = 0
        node.data['loss'] = 0
        node.data['appear'] = 0
        node.data['genes'] = 0


def count_dup_loss_tree(tree, stree, gene2species, recon=None, events=None):
    """count dup loss

    TODO: generalize to non-MPR recon/events
          (in particular, to handle duplication followed immediately by loss)
    """

    if recon is None:
        recon = reconcile(tree, stree, gene2species)
    if events is None:
        events = label_events(tree, recon)
    losses = find_loss(tree, recon)

    dup = 0
    loss = 0
    appear = 0

    # count appearance
    recon[tree.root].data["appear"] += 1
    appear += 1

    # count dups
    for node, event in events.iteritems():
        if event == "dup":
            recon[node].data['dup'] += 1
            dup += 1
        elif event == "gene":
            recon[node].data['genes'] += 1

    # count losses
    for _, snode in losses: # gnode, snode
        snode.data['loss'] += 1
        loss += 1

    return dup, loss, appear


def count_ancestral_genes(stree):
    """count ancestral genes"""
    def walk(node):
        """helper function"""
        if not node.is_leaf():
            counts = []
            for child in node.children:
                walk(child)
                counts.append(child.data['genes']
                              - child.data['appear']
                              - child.data['dup']
                              + child.data['loss'])
            assert util.equal(* counts), (node.name, str(counts))
            node.data['genes'] = counts[0]
    walk(stree.root)


def count_dup_loss_trees(trees, stree, gene2species):
    """
    Returns new species tree with dup,loss,appear,genes counts in node's data
    """

    stree = stree.copy()
    init_dup_loss_tree(stree)

    for tree in trees:
        count_dup_loss_tree(tree, stree, gene2species)
    count_ancestral_genes(stree)

    return stree


def write_event_tree(stree, out=sys.stdout):
    """
    Draw tree with event counts
    """
    labels = {}
    for name, node in stree.nodes.iteritems():
        labels[name] = "[%s]\nD=%d,L=%d;\nG=%d;" % \
                       (str(name),
                        node.data['dup'], node.data['loss'],
                        node.data['genes'])

    treelib.draw_tree(stree, labels=labels, minlen=15, spacing=4,
                      label_offset=-3,
                      out=out)


#=============================================================================
# duplication consistency

def dup_consistency(tree, recon, events):
    """
    Calculate duplication consistency scores for a reconcilied tree

    See Vilella2009 (Ensembl)
    """

    if len(tree.leaves()) == 1:
        return {}

    spset = {}
    def walk(node):
        """helper function"""
        for child in node.children:
            walk(child)
        if node.is_leaf():
            spset[node] = set([recon[node]])
        elif len(node.children) == 1:
            pass
        elif len(node.children) == 2:
            spset[node] = (spset[node.children[0]] |
                           spset[node.children[1]])
        else:
            raise Exception("too many children (%d)" % len(node.children))
    walk(tree.root)

    conf = {}
    for node in tree:
        if events[node] == "dup":
            conf[node] = (len(spset[node.children[0]] &
                              spset[node.children[1]]) /
                          float(len(spset[node])))

    return conf


def fix_ils_errors(events, dupcons, new_copy=True):
    """
    Relabels "dup" nodes as "spec" if the duplication consistency
    score for the node is 0, i.e. looks like ILS.
    """
    if new_copy:
        events = events.copy()

    for node, score in dupcons.iteritems():
        if score == 0:
            events[node] = "spec"
    return events


#=============================================================================
# tree rooting

def recon_root(gtree, stree, gene2species=gene2species,
               rootby="duploss", new_copy=True,
               keep_name=False, return_cost=False,
               dupcost=1, losscost=1):
    """
    Reroot a tree by minimizing the number of duplications/losses/both

    Note that rootby trumps dupcost/losscost.

    gtree -- gene tree to reroot
    stree -- species tree
    gene2species -- mapping from gene names to species names
    rootby -- method to root by ("dup", "loss", "duploss")
    returnCost -- if True, also return event cost of rerooted gene tree
    dupcost -- cost of gene duplication
    losscost -- cost of gene loss
    keepName -- if True, reuse existing root name for new root node
    """
    # assert valid inputs
    assert rootby in ["dup", "loss", "duploss"], "unknown rootby value '%s'" % rootby
    assert dupcost >= 0 and losscost >= 0

    # make a consistent unrooted copy of gene tree
    if new_copy:
        gtree = gtree.copy()

    if len(gtree.leaves()) == 2:
        if return_cost:
            recon = reconcile(gtree, stree, gene2species)
            events = label_events(gtree, recon)
            cost = 0
            if rootby in ["dup", "duploss"] and dupcost != 0:
                cost += count_dup(gtree, events) * dupcost
            if rootby in ["loss", "duploss"] and losscost != 0:
                cost += count_loss(gtree, recon) * losscost
            return gtree, cost
        else:
            return gtree

    if keep_name:
        oldroot = gtree.root.name

    treelib.unroot(gtree, new_copy=False)
    treelib.reroot(gtree,
                   gtree.nodes[sorted(gtree.leaf_names())[0]].parent.name,
                   on_branch=False, new_copy=False)

    # make recon root consistent for rerooting tree of the same names
    # TODO: there is the possibility of ties, they are currently broken
    # arbitrarily.  In order to make comparison of reconRooted trees with
    # same gene names accurate, hashOrdering must be done, for now.
    hash_order_tree(gtree, gene2species)

    # get list of edges to root on
    edges = []
    def walk(node):
        """helper function"""
        edges.append((node, node.parent))
        if not node.is_leaf():
            node.recurse(walk)
            edges.append((node, node.parent))
    for child in gtree.root.children:
        walk(child)

    # try initial root and recon
    treelib.reroot(gtree, edges[0][0].name, new_copy=False)
    if keep_name:
        gtree.rename(gtree.root.name, oldroot)
    recon = reconcile(gtree, stree, gene2species)
    events = label_events(gtree, recon)

    # find reconciliation that minimizes dup/loss
    minroot = edges[0]
    rootedge = sorted(edges[0])
    cost = 0

    if rootby in ["dup", "duploss"] and dupcost != 0:
        cost += count_dup(gtree, events) * dupcost
    if rootby in ["loss", "duploss"] and losscost != 0:
        cost += count_loss(gtree, recon) * losscost
    mincost = cost

    # try rooting on everything
    edge = None
    for edge in edges[1:]:
        if sorted(edge) == rootedge:
            continue
        rootedge = sorted(edge)

        node1, node2 = edge
        if node1.parent != node2:
            node1, node2 = node2, node1
        assert node1.parent == node2, "%s %s" % (node1.name, node2.name)

        # uncount cost
        if rootby in ["dup", "duploss"] and dupcost != 0:
            if events[gtree.root] == "dup":
                cost -= dupcost
            if events[node2] == "dup":
                cost -= dupcost
        if rootby in ["loss", "duploss"] and losscost != 0:
            cost -= len(find_loss_under_node(gtree.root, recon)) * losscost
            cost -= len(find_loss_under_node(node2, recon)) * losscost

        # new root and recon
        treelib.reroot(gtree, node1.name, new_copy=False, keep_name=keep_name)

        recon[node2] = reconcile_node(node2, stree, recon)
        recon[gtree.root] = reconcile_node(gtree.root, stree, recon)
        events[node2] = label_events_node(node2, recon)
        events[gtree.root] = label_events_node(gtree.root, recon)

        if rootby in ["dup", "duploss"] and dupcost != 0:
            if events[gtree.root] == "dup":
                cost += dupcost
            if events[node2] == "dup":
                cost += dupcost
        if rootby in ["loss", "duploss"] and losscost != 0:
            cost += len(find_loss_under_node(gtree.root, recon)) * losscost
            cost += len(find_loss_under_node(node2, recon)) * losscost

        # keep track of min cost
        if cost < mincost:
            mincost = cost
            minroot = edge

    # root tree by minroot
    if edge != minroot:
        node1, node2 = minroot
        if node1.parent != node2:
            node1, node2 = node2, node1
        assert node1.parent == node2

        treelib.reroot(gtree, node1.name, new_copy=False, keep_name=keep_name)

    if return_cost:
        return gtree, mincost
    return gtree


def midroot_recon(tree, stree, recon, events, params, generate):
    """Reroot a tree at midpoint"""

    node1, node2 = tree.root.children

    specs1 = []
    specs2 = []

    # find nearest specs/genes
    def walk(node, specs):
        """helper function"""
        if events[node] == "dup":
            for child in node.children:
                walk(child, specs)
        else:
            specs.append(node)
    #walk(node1, specs1)
    #walk(node2, specs2)
    specs1 = node1.leaves()
    specs2 = node2.leaves()

    def get_dists(start, end):
        """helper function"""
        exp_dist = 0
        obs_dist = 0

        sstart = recon[start]
        send = recon[end]
        while sstart != send:
            exp_dist += params[sstart.name][0]
            sstart = sstart.parent

        while start != end:
            obs_dist += start.dist
            start = start.parent

        return exp_dist, obs_dist / generate

    diffs1 = []
    for spec in specs1:
        if events[tree.root] == "spec":
            exp_dist1, obs_dist1 = get_dists(spec, tree.root)
        else:
            exp_dist1, obs_dist1 = get_dists(spec, node1)
        diffs1.append(obs_dist1 - exp_dist1)

    diffs2 = []
    for spec in specs2:
        if events[tree.root] == "spec":
            exp_dist2, obs_dist2 = get_dists(spec, tree.root)
        else:
            exp_dist2, obs_dist2 = get_dists(spec, node2)
        diffs2.append(obs_dist2 - exp_dist2)

    totdist = (node1.dist + node2.dist) / generate

    left = node1.dist - stats.mean(diffs1)
    right = totdist - node2.dist + stats.mean(diffs2)

    #print diffs1, diffs2
    #print stats.mean(diffs1), stats.mean(diffs2)

    mid = util.clamp((left + right) / 2.0, 0, totdist)

    node1.dist = mid * generate
    node2.dist = (totdist - mid) * generate


def stree2gtree(stree, genes, gene2species):
    """Create a gene tree with the same topology as the species tree"""
    tree = stree.copy()

    for gene in genes:
        tree.rename(gene2species(gene), gene)
    return tree


#=============================================================================
# relationships
# encoded using gene name tuples


def get_gene_dups(tree, events):
    """Returns duplications as gene name tuples"""
    return set(tuple(sorted([tuple(sorted(child.leaf_names()))
                             for child in node.children]))
               for node, kind in events.iteritems()
               if kind == "dup")

def get_speciations(tree, events):
    """Returns speciations as gene name tuples"""
    return set(tuple(sorted([tuple(sorted(child.leaf_names()))
                             for child in node.children]))
               for node, kind in events.iteritems()
               if kind == "spec")


def get_gene_losses(tree, recon):
    """Returns losses as gene name, species name tuples"""
    return set((loss[0].name, loss[1].name)
               for loss in find_loss(tree, recon))


def get_orthologs(tree, events):
    """Returns orthologs as gene name pairs"""

    specs = [sorted([sorted(child.leaf_names())
                     for child in node.children])
             for node in events
             if events[node] == "spec"]

    return set(tuple(sorted((a, b)))
               for x in specs
               for a in x[0]
               for b in x[1])


#=============================================================================
# Tree hashing

def hash_tree_compose(child_hashes, node=None):
    """Compose hash using leaf names"""
    return "(%s)" % ",".join(child_hashes)


def hash_tree_compose_names(child_hashes, node=None):
    """Compose hash using leaf and internal names"""
    return "(%s)%s" % (",".join(child_hashes), node.name)


def hash_tree(tree, smap=lambda x: x, compose=hash_tree_compose):
    """Hash tree using using leaf names"""
    def walk(node):
        """helper function"""
        if node.is_leaf():
            return smap(node.name)

        child_hashes = map(walk, node.children)
        child_hashes.sort()
        return compose(child_hashes, node)

    if isinstance(tree, treelib.Tree) or hasattr(tree, "root"):
        return walk(tree.root)
    elif isinstance(tree, treelib.TreeNode):
        return walk(tree)
    else:
        raise Exception("Expected Tree object")


def hash_tree_names(tree, smap=lambda x: x, compose=hash_tree_compose_names):
    """Hash tree using leaf and internal names"""
    return hash_tree(tree, smap, compose)


def hash_order_tree(tree, smap=lambda x: x, compose=hash_tree_compose):
    """Reorder tree using leaf names"""

    def walk(node):
        """helper function"""
        if node.is_leaf():
            return smap(node.name)

        child_hashes = map(walk, node.children)
        ind = util.sort_index(child_hashes)
        child_hashes = [child_hashes[i] for i in ind]
        node.children = [node.children[i] for i in ind]
        return compose(child_hashes)

    walk(tree.root)


#=============================================================================
# add implied speciation nodes to a gene tree

def add_spec_node(node, snode, tree, recon, events):
    """
    insert new speciation node above gene node 'node' from gene tree 'tree'

    new node reconciles to species node 'snode'.  Modifies recon and events
    accordingly
    """

    newnode = treelib.TreeNode(tree.new_name())
    parent = node.parent

    # find index of node in parent's children
    nodei = parent.children.index(node)

    # insert new node into tree
    tree.add_child(parent, newnode)
    parent.children[nodei] = newnode
    parent.children.pop()
    tree.add_child(newnode, node)

    # add recon and events info
    recon[newnode] = snode
    events[newnode] = "spec"

    return newnode


def remove_spec_node(node, tree, recon=None, events=None):
    """
    removes speciation node 'node' from gene tree 'tree'

    Modifies recon and events accordingly
    """
    assert len(node.children) == 1
    parent = node.parent
    child = node.children[0]

    # remove node from tree - handle root node specially
    if parent is None:
        tree.root = child
        tree.root.parent = None
        node.children = []
    else:
        nodei = parent.children.index(node)
        parent.children[nodei] = child
        child.parent = parent
        node.parent = None
        node.children = []
    tree.remove(node)

    # remove recon and events info
    if recon:
        del recon[node]
    if events:
        del events[node]


def add_implied_spec_nodes(tree, stree, recon, events):
    """
    adds speciation nodes to tree that are implied but are not present
    because of gene losses
    """

    added_nodes = []

    for node in list(tree):
        # process this node and the branch above it

        # handle root node specially
        if node.parent is None:
            # ensure root of gene tree properly reconciles to
            # root of species tree
            if recon[node] == stree.root:
                continue
            tree.root = treelib.TreeNode(tree.new_name())
            tree.add_child(tree.root, node)
            recon[tree.root] = stree.root
            events[tree.root] = "spec"
            added_nodes.append(tree.root)

        # determine starting and ending species
        sstart = recon[node]
        send = recon[node.parent]

        # the species path is too short to have implied speciations
        if sstart == send:
            continue

        parent = node.parent

        # determine species path of this gene branch (node, node->parent)
        snode = sstart.parent

        while snode != send:
            added_nodes.append(add_spec_node(node, snode, tree, recon, events))
            node = node.parent
            snode = snode.parent


        # determine whether node.parent is a dup
        # if so, send (a.k.a. species end) is part of species path
        if events[parent] == "dup":
            added_nodes.append(add_spec_node(node, send, tree, recon, events))

    return added_nodes


def remove_implied_spec_nodes(tree, added_nodes, recon=None, events=None):
    """
    removes speciation nodes from tree
    """
    for node in added_nodes:
        remove_spec_node(node, tree, recon, events)


#=============================================================================
# reconciliation rearrangements

def change_recon_up(recon, node, events=None):
    """
    Move the mapping of a node up one branch
    """

    if events is not None and events[node] == "spec":
        # promote speciation to duplication
        # R'(v) = e(R(u))
        events[node] = "dup"
    else:
        # R'(v) = p(R(u))
        recon[node] = recon[node].parent


def change_recon_down(recon, node, schild, events=None):
    """
    Move the mapping of a node down one branch
    """

    if events is not None and recon[node] == schild:
        events[node] = "spec"
    else:
        recon[node] = schild


def can_change_recon_up(recon, node, events=None):
    """Returns True is recon can remap node one 'step' up"""

    if events is not None and events[node] == "spec" and not node.is_leaf():
        # promote speciation to duplication
        return True

    # move duplication up one branch
    rnode = recon[node]
    prnode = rnode.parent

    # rearrangement is valid if
    return (not node.is_leaf() and
            prnode is not None and #  1. there is parent sp. branch
            (node.parent is None or # 2. no parent to restrict move
             rnode != recon[node.parent] # 3. not already matching parent
            )
           )


def enum_recon(tree, stree, depth=None,
               step=0, preorder=None,
               recon=None, events=None,
               gene2species=None):
    """
    Enumerate reconciliations between a gene tree and a species tree
    """

    if recon is None:
        recon = reconcile(tree, stree, gene2species)
        events = label_events(tree, recon)

    if preorder is None:
        preorder = list(tree.preorder())

    # yield current recon
    yield recon, events

    if depth is None or depth > 0:
        for i in xrange(step, len(preorder)):
            node = preorder[i]
            if can_change_recon_up(recon, node, events):
                schild = recon[node]
                change_recon_up(recon, node, events)

                # recurse
                depth2 = depth - 1 if depth is not None else None
                for r, e in enum_recon(tree, stree, depth2,
                                       i, preorder,
                                       recon, events):
                    yield r, e

                change_recon_down(recon, node, schild, events)


#=============================================================================
# Phylogenetic reconstruction: Neighbor-Joining

def neighborjoin(distmat, genes, usertree=None):
    """Neighbor joining algorithm"""

    tree = treelib.Tree()
    leaves = {}    # key = (gene1, gene2)
    dists = util.Dict(dim=2)
    restdists = {}

    # initialize distances
    for i in xrange(len(genes)):
        r = 0
        for j in xrange(len(genes)):
            dists[genes[i]][genes[j]] = distmat[i][j]
            r += distmat[i][j]
        restdists[genes[i]] = r / (len(genes) - 2)

    # initialize leaves
    for gene in genes:
        tree.add(treelib.TreeNode(gene))
        leaves[gene] = 1

    # if usertree is given, determine merging order
    merges = []
    newnames = {}
    if usertree != None:
        def walk(node):
            if not node.is_leaf():
                assert len(node.children) == 2, \
                    Exception("usertree is not binary")

                for child in node:
                    walk(child)
                merges.append(node)
                newnames[node] = len(merges)
            else:
                newnames[node] = node.name
        walk(usertree.root)
        merges.reverse()

    # join loop
    while len(leaves) > 2:
        # search for closest genes
        if not usertree:
            low = util.INF
            lowpair = (None, None)
            leaveslst = leaves.keys()

            for i in xrange(len(leaves)):
                for j in xrange(i+1, len(leaves)):
                    gene1, gene2 = leaveslst[i], leaveslst[j]
                    dist = dists[gene1][gene2] - restdists[gene1] \
                                               - restdists[gene2]
                    if dist < low:
                        low = dist
                        lowpair = (gene1, gene2)
        else:
            node = merges.pop()
            lowpair = (newnames[node.children[0]],
                       newnames[node.children[1]])

        # join gene1 and gene2
        gene1, gene2 = lowpair
        parent = treelib.TreeNode(tree.new_name())
        tree.add_child(parent, tree.nodes[gene1])
        tree.add_child(parent, tree.nodes[gene2])

        # set distances
        tree.nodes[gene1].dist = (dists[gene1][gene2] + restdists[gene1] -
                                  restdists[gene2]) / 2.0
        tree.nodes[gene2].dist = dists[gene1][gene2] - tree.nodes[gene1].dist

        # gene1 and gene2 are no longer leaves
        del leaves[gene1]
        del leaves[gene2]

        gene3 = parent.name
        r = 0
        for gene in leaves:
            dists[gene3][gene] = (dists[gene1][gene] + dists[gene2][gene] -
                                  dists[gene1][gene2]) / 2.0
            dists[gene][gene3] = dists[gene3][gene]
            r += dists[gene3][gene]
        leaves[gene3] = 1

        if len(leaves) > 2:
            restdists[gene3] = r / (len(leaves) - 2)

    # join the last two genes into a tribranch
    gene1, gene2 = leaves.keys() # pylint: disable=unbalanced-tuple-unpacking
    if not isinstance(gene1, int):
        gene1, gene2 = gene2, gene1
    tree.add_child(tree.nodes[gene1], tree.nodes[gene2])
    tree.nodes[gene2].dist = dists[gene1][gene2]
    tree.root = tree.nodes[gene1]

    # root tree according to usertree
    if usertree != None and treelib.is_rooted(usertree):
        roots = set([newnames[usertree.root.children[0]],
                     newnames[usertree.root.children[1]]])
        newroot = None
        for child in tree.root.children:
            if child.name in roots:
                newroot = child

        assert newroot != None

        treelib.reroot(tree, newroot.name, new_copy=False)

    return tree


#=============================================================================

def tree2distmat(tree, leaves):
    """Returns pair-wise distances between leaves of a tree"""

    # TODO: not implemented efficiently
    mat = []
    for i in xrange(len(leaves)):
        mat.append([])
        for j in xrange(len(leaves)):
            mat[-1].append(treelib.find_dist(tree, leaves[i], leaves[j]))

    return mat


#============================================================================
# branch splits

def find_splits(tree, rooted=False):
    """Find branch splits for a tree

    If 'rooted' is True, then orient splits based on rooting
    """

    all_leaves = set(tree.leaf_names())
    nall_leaves = len(all_leaves)

    # find descendants
    descendants = {}
    def walk(node):
        """helper function"""
        if node.is_leaf():
            s = descendants[node] = set([node.name])
        else:
            s = set()
            for child in node.children:
                s.update(walk(child))
            descendants[node] = s
        return s
    for child in tree.root.children:
        walk(child)

    # left child's descendants immediately defines
    # right child's descendants (by complement)
    if len(tree.root.children) == 2:
        # in order to work with rooted, be consistent about which descendents
        # to keep
        a, b = tree.root.children
        if len(descendants[a]) < len(descendants[b]):
            del descendants[a]
        elif len(descendants[b]) < len(descendants[a]):
            del descendants[b]
        elif min(descendants[a]) < min(descendants[b]):
            del descendants[a]
        else:
            del descendants[b]

    # build splits list
    splits = []
    for leaves in descendants.itervalues():
        if len(leaves) > 1 and \
           (rooted or len(leaves) < nall_leaves - 1):
            set1 = tuple(sorted(leaves))
            set2 = tuple(sorted(all_leaves - leaves))
            if not rooted:
                if len(set1) > len(set2) or \
                   (len(set1) == len(set2) and min(set1) > min(set2)):
                    set1, set2 = set2, set1
            splits.append((set1, set2))

    return splits


def robinson_foulds_error(tree1, tree2, rooted=False):
    """Returns RF error

    This definition of RF error is the fraction of branches in larger
    tree that are not present in the smaller tree.

    Of course, trees can be the same size as well.
    """
    splits1 = find_splits(tree1, rooted=rooted)
    splits2 = find_splits(tree2, rooted=rooted)

    overlap = set(splits1) & set(splits2)

    #assert len(splits1) == len(splits2)

    denom = float(max(len(splits1), len(splits2)))

    if denom == 0.0:
        return 0.0
    return 1 - (len(overlap) / denom)


#=============================================================================
# bootstrap

def add_bootstraps(tree, trees, rooted=False):
    """Add bootstrap support to tree"""

    # get bootstrap counts
    ntrees = 0
    split_counts = util.Dict(dim=1, default=0)
    for gtree in trees:
        ntrees += 1
        for split in find_splits(gtree, rooted=rooted):
            split_counts[split] += 1

    counts = {}
    for split, count in split_counts.iteritems():
        counts[split[0]] = count
        counts[split[1]] = count

    # add bootstrap support to tree
    def walk(node):
        """helper function"""
        if node.is_leaf():
            s = set([node.name])
        else:
            s = set()
            for child in node.children:
                s.update(walk(child))
            node.data["boot"] = counts.get(tuple(sorted(s)), 0)/float(ntrees)
        return s
    for child in tree.root.children:
        walk(child)

    if rooted:
        if tree.root.children[0].is_leaf() or \
           tree.root.children[1].is_leaf():
            tree.root.children[0].data["boot"] = 0
            tree.root.children[1].data["boot"] = 0
        else:
            b = (tree.root.children[0].data["boot"] + \
                 tree.root.children[1].data["boot"]) / 2.0
            tree.root.children[0].data["boot"] = b
            tree.root.children[1].data["boot"] = b

    return tree


#=============================================================================
# consensus methods

def _consensus_special(trees, rooted=False):
    """Handle special cases for consensus tree"""

    nleaves = len(trees[0].leaves())
    tree = treelib.Tree()

    # handle special cases
    if not rooted and nleaves == 3:
        leaves = trees[0].leaf_names()
        root = tree.make_root()
        n = tree.add_child(root, treelib.TreeNode(tree.new_name()))
        tree.add_child(n, treelib.TreeNode(leaves[0]))
        tree.add_child(n, treelib.TreeNode(leaves[1]))
        tree.add_child(root, treelib.TreeNode(leaves[2]))
        return tree

    elif nleaves == 2:
        leaves = trees[0].leaf_names()
        root = tree.make_root()
        tree.add_child(root, treelib.TreeNode(leaves[0]))
        tree.add_child(root, treelib.TreeNode(leaves[1]))
        return tree

    return None


def consensus_majority_rule(trees, extended=True, rooted=False):
    """Performs majority rule on a set of trees

    extended -- if True, performs the extended majority rule
    rooted   -- if True, assumes trees are rooted
    """

    nleaves = len(trees[0].leaves())
    ntrees = len(trees)

    # handle special cases
    contree = _consensus_special(trees, rooted)
    if contree is not None:
        return contree

    # consensus tree
    contree = treelib.Tree()

    # count all splits
    split_counts = util.Dict(dim=1, default=0)
    for tree in trees:
        for split in find_splits(tree, rooted):
            split_counts[split] += 1
    contree.nextname = max(tree.nextname for tree in trees)

    #util.print_dict(split_counts)

    # choose splits
    pick_splits = 0
    rank_splits = split_counts.items()
    rank_splits.sort(key=lambda x: x[1], reverse=True)

    # add splits to the contree in increasing frequency
    for split, count in rank_splits:
        if not extended and count <= ntrees / 2.0:
            continue

        # choose split if it is compatiable
        if _add_split_to_tree(contree, split, count / float(ntrees), rooted):
            pick_splits += 1

        # stop if enough splits are choosen
        if ((rooted and pick_splits >= nleaves - 2) \
            or (not rooted and pick_splits >= nleaves - 3)):
            break

    # add remaining leaves and remove clade data
    _post_process_split_tree(contree)

    return contree


def splits2tree(splits, rooted=False):
    """Builds a tree from a set of splits

    Silently reject splits that are in conflict.  Process splits in order.

    splits -- iterable of splits
    rooted -- if True treat splits as rooted/polarized
    """

    tree = treelib.Tree()
    for split in splits:
        _add_split_to_tree(tree, split, 1.0, rooted)
    _post_process_split_tree(tree)
    return tree


def _add_split_to_tree(tree, split, count, rooted=False):
    """Add split to tree"""

    split = (set(split[0]), set(split[1]))

    # init first split
    if len(tree) == 0:
        root = tree.make_root()
        root.data["leaves"] = split[0] | split[1]

        if len(split[0]) == 1:
            node = tree.add_child(root, treelib.TreeNode(list(split[0])[0]))
            node.data["leaves"] = split[0]
            node.data["boot"] = count
        else:
            node = tree.add_child(root, treelib.TreeNode(tree.new_name()))
            node.data["leaves"] = split[0]
            node.data["boot"] = count

        if len(split[1]) == 1:
            node = tree.add_child(root, treelib.TreeNode(list(split[1])[0]))
            node.data["leaves"] = split[1]
            node.data["boot"] = count

        return True

    def walk(node, clade):
        """helper function"""
        if node.is_leaf():
            # make new child
            if len(clade) == 1:
                name = list(clade)[0]
            else:
                name = tree.new_name()
            child = tree.add_child(node, treelib.TreeNode(name))
            child.data["leaves"] = clade
            child.data["boot"] = count
            return True

        # which children intersect this clade?
        intersects = []
        for child in node:
            leaves = child.data["leaves"]
            intersect = clade & leaves

            if len(clade) == len(intersect):
                if len(intersect) < len(leaves):
                    # subset, recurse
                    return walk(child, clade)
                # len(intersect) == len(leaves), split is already present
                return True

            elif not intersect: # len(intersect) == 0
                continue

            elif len(intersect) == len(leaves):
                # len(clade) > len(leaves)
                # superset
                intersects.append(child)
            else:
                # conflict
                return False

        # insert new node
        if len(clade) == 1:
            name = list(clade)[0]
        else:
            name = tree.new_name()
        new_node = tree.add_child(node, treelib.TreeNode(name))
        new_node.data["leaves"] = clade
        new_node.data["boot"] = count
        for child in intersects:
            tree.remove(child)
            tree.add_child(new_node, child)

        return True

    # try to place split into tree
    if rooted:
        walk(tree.root, split[0])
    else:
        if walk(tree.root, split[0]):
            return True
        return walk(tree.root, split[1])

    # split is in conflict
    return False


def _post_process_split_tree(tree):
    """Post-process a tree built from splits"""

    for node in list(tree):
        if len(node.data["leaves"]) > 1:
            for leaf_name in node.data["leaves"]:
                for child in node:
                    if leaf_name in child.data.get("leaves", ()):
                        break
                else:
                    child = tree.add_child(node, treelib.TreeNode(leaf_name))
        else:
            assert node.name == list(node.data["leaves"])[0], node.name

    # remove leaf data and set root
    for node in tree:
        if "leaves" in node.data:
            del node.data["leaves"]


#=============================================================================
# multifucating trees

def ensure_binary_tree(tree):
    """Arbitrarily expand multifurcating nodes"""

    # first tree just rerooting root branch
    if len(tree.root.children) > 2:
        treelib.reroot(tree, tree.root.children[0].name, new_copy=False)

    multibranches = [node for node in tree
                     if len(node.children) > 2]

    for node in multibranches:
        children = list(node.children)

        # remove children
        for child in children:
            tree.remove(child)

        # add back in binary
        while len(children) > 2:
            left = children.pop()
            right = children.pop()
            newnode = treelib.TreeNode(tree.new_name())
            newnode.data['boot'] = 0
            tree.add_child(newnode, left)
            tree.add_child(newnode, right)
            children.append(newnode)

        # add last two to original node
        tree.add_child(node, children.pop())
        tree.add_child(node, children.pop())
