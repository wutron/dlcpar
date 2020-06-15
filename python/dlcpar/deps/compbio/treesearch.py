"""
treesearch.py
Propose trees based on rearrangement operations
"""

# python libraries
import math
import random
from abc import ABCMeta, abstractmethod

# rasmus libraries
from rasmus import stats
from rasmus import treelib
from rasmus import util

# compbio libraries
from compbio import phylo


#=============================================================================
# local rearrangements


def perform_nni(tree, node1, node2, change=0, rooted=True):
    r"""Propose new tree using Nearest Neighbor Interchange (NNI)

    Branch for NNI is specified by giving its two incident nodes (node1 and
    node2).  Change specifies which  subtree of node1 will be swapped with
    the uncle.  See figure below.

         node2
        /     \
      uncle    node1
               /  \
         child[0]  child[1]

    special case with rooted branch and rooted=False:

              node2
             /     \
        node2'      node1
       /     \     /     \
      uncle   * child[0] child[1]
    """

    if node1.parent != node2:
        node1, node2 = node2, node1

    # try to see if edge is one branch (not root edge)
    if not rooted and treelib.is_rooted(tree) and \
       node2 == tree.root:
        # special case of specifying root edge
        if node2.children[0] == node1:
            node2 = node2.children[1]
        else:
            node2 = node2.children[0]

        # edge is not an internal edge, give up
        if len(node2.children) < 2:
            return

    if node1.parent == node2.parent == tree.root:
        uncle = 0

        if len(node2.children[0].children) < 2 and \
           len(node2.children[1].children) < 2:
            # can't do NNI on this branch
            return
    else:
        assert node1.parent == node2

        # find uncle
        uncle = 0
        if node2.children[uncle] == node1:
            uncle = 1

    # swap parent pointers
    node1.children[change].parent = node2
    node2.children[uncle].parent = node1

    # swap child pointers
    node2.children[uncle], node1.children[change] = \
        node1.children[change], node2.children[uncle]


def propose_random_nni(tree):
    """Propose random NNI rearrangement"""

    nodes = tree.nodes.values()

    # find edges for NNI
    while True:
        node1 = random.sample(nodes, 1)[0]
        if not node1.is_leaf() and node1.parent is not None:
            break

    node2 = node1.parent
    #a = node1.children[random.randint(0, 1)]
    #b = node2.children[1] if node2.children[0] == node1 else node2.children[0]
    #assert a.parent.parent == b.parent

    return node1, node2, random.randint(0, 1)


def perform_spr(tree, subtree, newpos):
    r"""Propose new tree using Subtree Pruning and Regrafting (SPR)

    a = subtree
    e = newpos

    BEFORE
             ...
         f         d
        /           \
       c             e
      / \           ...
     a   b
    ... ...

    AFTER
            ...
        f         d
       /           \
      b             c
     ...           / \
                  a   e
                 ... ...

    Requirements:
    1. a (subtree) is not root or children of root
    2. e (newpos) is not root, a, descendant of a, c (parent of a), or
       b (sibling of a)
    3. tree is binary
    """

    # TODO: check requirements

    a = subtree
    e = newpos

    c = a.parent
    f = c.parent
    bi = 1 if c.children[0] == a else 0
    b = c.children[bi]
    ci = 0  if f.children[0] == c else 1
    d = e.parent
    ei = 0 if d.children[0] == e else 1

    d.children[ei] = c
    c.children[bi] = e
    f.children[ci] = b
    b.parent = f
    c.parent = d
    e.parent = c


def propose_random_spr(tree):
    r"""Propose random SPR rearrangement

    What if e == f  (also equivalent to NNI) this is OK

    BEFORE
           d
          / \
         e  ...
        / \
       c  ...
      / \
     a   b
    ... ...

    AFTER
           d
          / \
         c
        / \
       a   e
      ... / \
         b  ...
        ...

    What if d == f  (also equivalent to NNI) this is OK

    BEFORE
         f
        / \
       c   e
      / \  ...
     a   b
    ... ...

    AFTER
         f
        / \
       b   c
      ... / \
         a   e
        ... ...

    Requirements:
    1. a (subtree) is not root or children of root
    2. e (newpos) is not root, a, descendant of a, c (parent of a), or
       b (sibling of a)
    3. tree is binary
    """

    # TODO: check requirements

    assert len(tree.nodes) >= 5, "Tree is too small"

    # find subtree (a) to cut off (any node that is not root or child of root)
    nodes = tree.nodes.values()
    while True:
        a = random.sample(nodes, 1)[0]
        if (a.parent is not None and a.parent.parent is not None):
            break
    subtree = a

    # find sibling (b) of a
    c = a.parent
    bi = 1 if c.children[0] == a else 0
    b = c.children[bi]

    # choose newpos (e)
    e = None
    while True:
        e = random.sample(nodes, 1)[0]

        # test if e is a valid choice
        if e.parent is None or e == a or e == c or e == b:
            continue

        # also test if e is a descendent of a
        under_a = False
        ptr = e.parent
        while ptr is not None:
            if ptr == a:
                under_a = True
                break
            ptr = ptr.parent
        if under_a:
            continue

        break
    newpos = e

    return subtree, newpos


#=============================================================================
# tree seach

class TreeSearch(object):
    """Abstract class for proposing trees"""

    __metaclass__ = ABCMeta


    def __init__(self, tree):
        self.tree = tree


    def __iter__(self):
        return self


    def set_tree(self, tree):
        """Set current tree"""
        self.tree = tree


    def get_tree(self):
        """Get current tree"""
        return self.tree


    @abstractmethod
    def propose(self):
        """Propose tree"""
        raise NotImplementedError("subclasses must override propose()")


    @abstractmethod
    def revert(self):
        """Revert search"""
        raise NotImplementedError("subclasses must override revert()")


    def reset(self):
        """Reset search"""
        pass


    def next(self):
        """Return next tree"""
        return self.propose()


class TreeSearchNni(TreeSearch):
    """Propose trees using NNI"""

    def __init__(self, tree):
        TreeSearch.__init__(self, tree)
        self.set_tree(tree)
        self.node1 = None
        self.node2 = None
        self.child = None


    def set_tree(self, tree):
        """Set tree"""
        self.tree = tree
        self.reset()


    def propose(self):
        """Propose tree"""
        self.node1, self.node2, self.child = propose_random_nni(self.tree)
        perform_nni(self.tree, self.node1, self.node2, self.child)
        return self.tree


    def revert(self):
        """Revert search"""
        if self.node1 is not None:
            perform_nni(self.tree, self.node1, self.node2, self.child)
        return self.tree


    def reset(self):
        """Reset search"""
        self.node1 = None
        self.node2 = None
        self.child = None


class TreeSearchSpr(TreeSearch):
    """Propose trees using SPR"""

    def __init__(self, tree):
        TreeSearch.__init__(self, tree)
        self.set_tree(tree)
        self.node1 = None
        self.node2 = None


    def set_tree(self, tree):
        """Set tree"""
        self.tree = tree
        self.reset()


    def propose(self):
        """Propose tree"""

        # choose SPR move
        self.node1, node3 = propose_random_spr(self.tree)

        # remember sibling of node1
        p = self.node1.parent
        self.node2 = (p.children[1] if p.children[0] == self.node1
                      else p.children[0])

        # perform SPR move
        perform_spr(self.tree, self.node1, node3)
        return self.tree


    def revert(self):
        """Revert search"""
        if self.node1 is not None:
            perform_spr(self.tree, self.node1, self.node2)
        return self.tree


    def reset(self):
        """Reset search"""
        self.node1 = None
        self.node2 = None


class TreeSearchMix(TreeSearch):
    """Propose trees using mix of operations"""

    def __init__(self, tree):
        TreeSearch.__init__(self, tree)
        self.total_weight = 0.0
        self.last_propose = 0
        self.methods = []
        self.set_tree(tree)


    def set_tree(self, tree):
        """Set tree"""
        self.tree = tree
        for method in self.methods:
            method[0].set_tree(tree)


    def add_proposer(self, proposer, weight):
        """Add proposer with weight"""
        self.total_weight += weight
        self.methods.append((proposer, weight))


    def propose(self):
        """Propose tree"""
        # randomly choose method
        choice = random.random() * self.total_weight
        s = self.methods[0][1]
        i = 0
        while i < len(self.methods)-1 and s < choice:
            i += 1
            s += self.methods[i][1]

        # make proposal
        self.last_propose = i
        self.tree = self.methods[i][0].propose()
        return self.tree


    def revert(self):
        """Revert search"""
        self.tree = self.methods[self.last_propose][0].revert()
        return self.tree


    def reset(self):
        """Reset search"""
        for method in self.methods:
            method[0].reset()


class TreeSearchUnique(TreeSearch):
    """Propose unique tree topologies"""

    def __init__(self, tree, search, tree_hash=None, maxtries=5,
                 auto_add=True):
        TreeSearch.__init__(self, tree)
        self.search = search
        self.seen = set()
        self._tree_hash = tree_hash if tree_hash else phylo.hash_tree
        self.maxtries = maxtries
        self.auto_add = auto_add
        self.set_tree(tree)


    def set_tree(self, tree):
        """Set tree"""
        self.tree = tree
        self.search.set_tree(tree)


    def propose(self):
        """Propose tree"""
        for i in xrange(self.maxtries):
            if i > 0:
                self.search.revert()
            tree = self.search.propose()
            top = self._tree_hash(tree)
            if top not in self.seen:
                #util.logger("tried", i, len(self.seen))
                break
        else:
            #util.logger("maxtries", len(self.seen))
            pass

        if self.auto_add:
            self.seen.add(top)
        self.tree = tree
        return self.tree


    def revert(self):
        """Revert search"""
        self.tree = self.search.revert()
        return self.tree


    def reset(self):
        """Reset search"""
        self.seen.clear()
        self.search.reset()


    def add_seen(self, tree):
        """Add tree to set of trees proposed"""
        top = self._tree_hash(tree)
        self.seen.add(top)


class TreeSearchPrescreen(TreeSearch):
    """Propose trees with prescreening"""

    def __init__(self, tree, search, prescreen, poolsize):
        TreeSearch.__init__(self, tree)
        self.search = TreeSearchUnique(tree, search, auto_add=False)
        self.prescreen = prescreen
        self.poolsize = poolsize
        self.oldtree = None
        self.set_tree(tree)


    def set_tree(self, tree):
        """Set tree"""
        self.tree = tree
        self.search.set_tree(tree)


    def propose(self):
        """Propose tree"""

        # save old topology
        self.oldtree = self.tree.copy()

        pool = []
        best_score = self.prescreen(self.tree)
        total = -util.INF

        # TODO: add unique filter

        # make many subproposals
        self.search.reset()
        for _ in xrange(self.poolsize):
            self.search.propose()
            score = self.prescreen(self.tree)
            tree = self.tree.copy()

            # save tree and logl
            pool.append((tree, score))
            total = stats.logadd(total, score)

            if score > best_score:
                # make more proposals off this one
                best_score = score
            else:
                self.search.revert()

        # propose one of the subproposals
        choice = random.random()
        partsum = -util.INF

        for tree, score in pool:
            partsum = stats.logadd(partsum, score)
            if choice < math.exp(partsum - total):
                # propose tree i
                treelib.set_tree_topology(self.tree, tree)
                break

        self.search.add_seen(self.tree)


    def revert(self):
        """Revert search"""
        if self.oldtree:
            treelib.set_tree_topology(self.tree, self.oldtree)


    def reset(self):
        """Reset search"""
        self.oldtree = None
        self.search.reset()
