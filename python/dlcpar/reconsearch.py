"""
reconsearch.py
Library to solve DLC MPR Problem using heuristic search
"""

# based on dlcoal recon.py

# python libraries
import random
import sys

# dlcpar libraries
from dlcpar import common
from dlcpar import constants

# rasmus libraries
from rasmus import stats
from rasmus import timer
from rasmus import treelib
from rasmus import util

# compbio libraries
from compbio import coal
from compbio import phylo
from compbio import phyloDLC
from compbio import treesearch

#=============================================================================
# reconciliation

def dlc_recon(tree, stree, gene2species,
              dupcost=constants.DEFAULT_DUP_COST,
              losscost=constants.DEFAULT_LOSS_COST,
              coalcost=constants.DEFAULT_COAL_COST,
              implied=True,
              nsearch=1000, nprescreen=20, nconverge=None,
              search=None,
              init_locus_tree=None,
              log=sys.stdout):
    """Perform reconciliation using heuristic search"""
    if search is None:
        search = lambda tree: LocusTreeSearch(tree, stree, gene2species,
                                              dupcost, losscost,
                                              nprescreen=nprescreen)

    reconer = DLCRecon(tree, stree, gene2species,
                       dupcost=dupcost, losscost=losscost, coalcost=coalcost,
                       implied=implied,
                       init_locus_tree=init_locus_tree,
                       log=log)
    reconer.set_proposer(DLCReconProposer(
        tree, stree, gene2species, search=search))
    return_val, mincost, runtime = reconer.recon(nsearch, nconverge)
    return return_val.get_dict(), mincost, runtime


class DLCRecon(object):
    """Reconciliation class using heuristic search"""

    def __init__(self, tree, stree, gene2species,
                 dupcost=constants.DEFAULT_DUP_COST,
                 losscost=constants.DEFAULT_LOSS_COST,
                 coalcost=constants.DEFAULT_COAL_COST,
                 implied=True,
                 init_locus_tree=None,
                 name_internal="n", log=sys.stdout):

        # rename input tree nodes
        common.rename_nodes(tree, name_internal)

        self.coal_tree = tree
        self.stree = stree
        self.gene2species = gene2species

        self.dupcost = dupcost
        self.losscost = losscost
        self.coalcost = coalcost
        self.implied = implied

        self.name_internal = name_internal
        self.log_stream = log
        self.init_locus_tree = init_locus_tree \
                               if init_locus_tree else tree.copy()

        self.proposer = DLCReconProposer(tree, stree, gene2species)
        self.log = timer.Timer(log)

        # these attributes are assigned when performing reconciliation using self.recon()
        self.maxrecon = None
        self.mincost = util.INF


    def set_proposer(self, proposer):
        """Set proposal algorithm"""
        self.proposer = proposer


    def recon(self, nsearch=1000, nconverge=None):
        """Perform reconciliation"""

        self.log.start("Reconciling")

        # initialize
        self.init_search()
        proposal = self.proposer.init_proposal()
        self.maxrecon = proposal.copy()

        # keep track of convergence
        mincost = util.INF
        if nconverge:
            iconverge = 0

        # search
        for _ in xrange(nsearch):
            # evaluate cost of proposal
            cost = self.eval_proposal(proposal)
##            util.print_dict(proposal.data)
##            print '\t'.join(map(lambda key: str(proposal.data[key]),
##                                ("cost", "ndup", "nloss", "ncoal")))

            # update maxrecon based on accepting / rejecting proposal
            self.eval_search(cost, proposal)

            # stop if converged
            # why not check accept? because can toggle between
            # multiple optimal solutions with the same cost
            if cost < mincost:
                mincost = cost
                if nconverge:
                    iconverge = 0
            else:
                if nconverge:
                    iconverge += 1
                    if iconverge == nconverge:
                        break

            # make new proposal
            proposal = self.proposer.next_proposal()

        # rename locus tree nodes
        common.rename_nodes(self.maxrecon.locus_tree, self.name_internal)

        # calculate runtime
        runtime = self.log.stop()

        return self.maxrecon, mincost, runtime


    def init_search(self):
        """Initialize new search"""

        # init locus tree as congruent to coal tree
        # equivalent to assuming no ILS
        #self.proposer.set_locus_tree(self.coal_tree.copy())
        self.proposer.set_locus_tree(self.init_locus_tree.copy())

        self.mincost = util.INF
        self.maxrecon = None


    def eval_proposal(self, proposal):
        """Compute cost of proposal"""

        if not phyloDLC.assert_daughters(proposal.locus_events, proposal.daughters):
            # ensure locus events (duplications) and daughters match
            ndup, nloss, ncoal = None, None, None
            dupcost, losscost, coalcost = util.INF, util.INF, util.INF
        else:
            # find dup cost
            if self.dupcost == 0:
                ndup = None
                dupcost = 0
            else:
                ndup = phylo.count_dup(proposal.locus_tree, proposal.locus_events)
                dupcost = ndup * self.dupcost

            # find loss cost
            if self.losscost == 0:
                nloss = None
                losscost = 0
            else:
                nloss = phylo.count_loss(proposal.locus_tree, proposal.locus_recon)
                losscost = nloss * self.losscost

            # find coal cost (first ensure bounded coalescent is satisfied
            # - should always be true based on how daughters are proposed)
            phyloDLC.assert_bounded_coal(self.coal_tree, proposal.coal_recon,
                                         proposal.locus_tree, proposal.daughters)
            if self.coalcost == 0:
                ncoal = None
                coalcost = 0
            else:
                # add implied speciation nodes if desired
                # this must be added AFTER counting dups and losses since it affects loss inference
                if self.implied:
                    added = phylo.add_implied_spec_nodes(proposal.locus_tree, self.stree,
                                                         proposal.locus_recon,
                                                         proposal.locus_events)

                ncoal = phyloDLC.count_coal(self.coal_tree,
                                            proposal.coal_recon,
                                            proposal.locus_tree)
                coalcost = ncoal * self.coalcost

                if self.implied:
                    phylo.remove_implied_spec_nodes(proposal.locus_tree,
                                                    added,
                                                    proposal.locus_recon,
                                                    proposal.locus_events)

        # total cost
        cost = dupcost + losscost + coalcost

        # logging info
        info = {}
        info["ndup"] = ndup
        info["nloss"] = nloss
        info["ncoal"] = ncoal
        info["cost"] = cost
        proposal.data = info

        return cost


    def eval_search(self, cost, proposal):
        """Evaluate proposal based on cost"""

        self.log_proposal(proposal)

        if cost < self.mincost:
            accept = True
        elif cost == self.mincost and random.random() < 0.5:
            accept = True
        else:
            accept = False

        if accept:
            self.mincost = cost
            self.maxrecon = proposal.copy()
            # search with a new copy
            self.proposer.accept()
        else:
            self.proposer.reject()


    def log_proposal(self, proposal):
        """Log proposal"""
        self.log_stream.write(repr(proposal) + "\n")
        self.log_stream.flush()



class DLCReconProposer(object):
    """Propose three-tree reconciliation"""

    def __init__(self, coal_tree, stree, gene2species,
                 search=treesearch.TreeSearchNni):
        self._coal_tree = coal_tree
        self._stree = stree
        self._gene2species = gene2species
        self._locus_search = search(None)

        self._accept_locus = False

        self._recon = None


    def set_locus_tree(self, locus_tree):
        """Set locus tree"""
        self._locus_search.set_tree(locus_tree)


    def init_proposal(self):
        """Get first proposal"""

        if self._locus_search.get_tree() is None:
            self._locus_search.set_tree(self._coal_tree.copy())
        self._recon = self._recon_lca(self._locus_search.get_tree().copy())

        return self._recon


    def next_proposal(self):
        """Get next proposal"""

        if len(self._locus_search.get_tree().leaves()) <= 2:
            return self._recon

        # if locus_tree has not yet been accepted, then revert it
        if not self._accept_locus:
            self._locus_search.revert()

        # propose new locus_tree
        self._locus_search.propose()
        self._accept_locus = False
        locus_tree = self._locus_search.get_tree().copy()

        # TODO: make recon root optional
        phylo.recon_root(locus_tree, self._stree,
                         self._gene2species,
                         new_copy=False)
        common.rename_nodes(locus_tree)

        # propose remaining parts of dlcoal recon
        self._recon = self._recon_lca(locus_tree)

        return self._recon


    def _recon_lca(self, locus_tree):
        """Reconcile using LCA"""
        # get locus tree, and LCA (MPR) locus_recon
        locus_recon = phylo.reconcile(locus_tree, self._stree,
                                      self._gene2species)
        locus_events = phylo.label_events(locus_tree, locus_recon)

        # propose LCA (MPR) coal_recon
        coal_recon = phylo.reconcile(self._coal_tree,
                                     locus_tree, lambda x: x)

        # propose daughters
        daughters = self._propose_daughters(
            self._coal_tree, coal_recon,
            locus_tree, locus_recon, locus_events)

        return phyloDLC.Recon(coal_recon, locus_tree, locus_recon, locus_events,
                              daughters)


    def _propose_daughters(self, coal_tree, coal_recon,
                           locus_tree, locus_recon, locus_events):
        """Propose daughters"""
        return _propose_daughters(coal_tree, coal_recon,
                                  locus_tree, locus_events)


    def accept(self):
        """Accept proposal"""
        self._accept_locus = True


    def reject(self):
        """Reject proposal"""
        pass


def _propose_daughters(coal_tree, coal_recon, locus_tree, locus_events):
    """Propose daughters"""

    lineages = coal.count_lineages_per_branch(coal_tree, coal_recon, locus_tree)
    daughters = set()

    for node, event in locus_events.iteritems():
        if event == "dup":
            # choose one of the children of node to be a daughter
            children = [child for child in node.children
                        if lineages[child][1] == 1]
            if children: # len(children) > 0
                daughters.add(children[stats.sample([1] * len(children))])

    return daughters


#=============================================================================
# tree search

class LocusTreeSearchPrescreen(treesearch.TreeSearchPrescreen):
    """Propose trees with prescreening"""

    def propose(self):
        """Propose tree"""

        # save old topology
        self.oldtree = self.tree.copy()

        pool = []
        best_score = self.prescreen(self.tree)

        # TODO: add unique filter

        # make many subproposals
        self.search.reset()
        for _ in xrange(self.poolsize):
            self.search.propose()
            score = self.prescreen(self.tree)
            tree = self.tree.copy()

            # save tree and score
            pool.append((tree, score))

            if score < best_score:
                # make more proposals off this one
                best_score = score
            elif score == best_score:
                # randomly decide whether to revert
                if random.random() < 0.5:
                    self.search.revert()
            else:
                self.search.revert()

        # propose one of the subproposals
        trees, _ = util.minall(pool, keyfunc=lambda it: it[0], minfunc=lambda it: it[1])
        tree = random.choice(trees)
        treelib.set_tree_topology(self.tree, tree)

        self.search.add_seen(self.tree)


class LocusTreeSearch(treesearch.TreeSearch):
    """Propose locus trees"""

    def __init__(self, tree, stree, gene2species,
                 dupcost, losscost,
                 tree_hash=None, nprescreen=20, weight=.2):
        treesearch.TreeSearch.__init__(self, tree)

        self.stree = stree
        self.gene2species = gene2species

        self.dupcost = dupcost
        self.losscost = losscost

        #self.search = UniqueTreeSearch(tree, treesearch.TreeSearchNni(tree),
        #                               tree_hash)
        #self.search = UniqueTreeSearch(tree, treesearch.TreeSearchSpr(tree),
        #                               tree_hash)

        mix = treesearch.TreeSearchMix(tree)
        mix.add_proposer(treesearch.TreeSearchNni(tree), .4)
        mix.add_proposer(treesearch.TreeSearchSpr(tree), .6)
        #self.search = treesearch.TreeSearchUnique(tree, mix, tree_hash)

        if nprescreen == 1:
            self.search = mix
        else:
            prescreen = LocusTreeSearchPrescreen(tree, mix,
                                                 self.prescreen,
                                                 poolsize=nprescreen)

            mix2 = treesearch.TreeSearchMix(tree)
            mix2.add_proposer(prescreen, 1.0-weight)
            mix2.add_proposer(mix, weight)

            self.search = mix2

    def set_tree(self, tree):
        """Set current tree"""
        self.tree = tree
        self.search.set_tree(tree)

    def reset(self):
        """Reset search"""
        self.search.reset()

    def propose(self):
        """Propose tree"""
        self.search.propose()
        return self.tree

    def revert(self):
        """Revert search"""
        self.search.revert()
        return self.tree

    def prescreen(self, tree):
        """Prescreen search based on dup-loss cost"""
        recon = phylo.reconcile(tree, self.stree, self.gene2species)
        events = phylo.label_events(tree, recon)

        if self.dupcost == 0:
            dupcost = 0
        else:
            ndup = phylo.count_dup(tree, events)
            dupcost = ndup * self.dupcost

        if self.losscost == 0:
            losscost = 0
        else:
            nloss = phylo.count_loss(tree, recon)
            losscost = nloss * self.losscost

        return dupcost + losscost
