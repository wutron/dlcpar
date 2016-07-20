# based on dlcoal recon.py
# relies on heuristic search

# python libraries
import sys, copy
import random

# dlcpar libraries
from dlcpar import common

# rasmus libraries
from rasmus import util, stats, treelib

# compbio libraries
from compbio import phylo
from compbio import coal
from compbio import phyloDLC

#=============================================================================
# reconciliation

def dlc_recon(tree, stree, gene2species,
              dupcost=1, losscost=1, coalcost=1, implied=True,
              nsearch=1000, nprescreen=20,
              search=None,
              init_locus_tree=None,
              log=sys.stdout):
    """Perform reconciliation using DLCoal model with parsimony costs"""
    if search is None:
        search = lambda tree: LocusTreeSearch(tree, stree, gene2species,
                                              dupcost, losscost,
                                              nprescreen=nprescreen)

    reconer = DLCRecon(tree, stree, gene2species,
                       dupcost=dupcost, losscost=losscost, coalcost=coalcost, implied=implied,
                       init_locus_tree=init_locus_tree,
                       log=log)
    reconer.set_proposer(DLCReconProposer(
        tree, stree, gene2species, search=search))
    return reconer.recon(nsearch).get_dict()


class DLCRecon (object):

    def __init__(self, tree, stree, gene2species,
                 dupcost=1, losscost=1, coalcost=1, implied=True,
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


    def set_proposer(self, proposer):
        """Set the proposal algorithm"""
        self.proposer = proposer


    def set_log(self, log):
        self.log_stream = log
        

    def recon(self, nsearch=1000):
        """Perform reconciliation"""
        
        self.init_search()
        proposal = self.proposer.init_proposal()
        self.maxrecon = proposal.copy()
        for i in xrange(nsearch):
##            if i % 10 == 0:
##                print "search", i

##            util.tic("eval")
            cost = self.eval_proposal(proposal)
##            util.print_dict(proposal.data)
##            print '\t'.join(map(lambda key: str(proposal.data[key]),
##                                ("cost", "ndup", "nloss", "ncoal")))
##            util.toc()

##            util.tic("prop")
            self.eval_search(cost, proposal)
            proposal = self.proposer.next_proposal()
##            util.toc()
        
        # rename locus tree nodes
        common.rename_nodes(self.maxrecon.locus_tree, self.name_internal)
        
        return self.maxrecon


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
                nloss = phylo.count_loss(proposal.locus_tree, self.stree, proposal.locus_recon)
                losscost = nloss * self.losscost

            # find coal cost (first ensure bounded coalescent is satisfied - should always be true based on how daughters are proposed)
            phyloDLC.assert_bounded_coal(self.coal_tree, proposal.coal_recon, proposal.locus_tree, proposal.daughters)
            if self.coalcost == 0:
                ncoal = None
                coalcost = 0
            else:
                # add implied speciation nodes if desired
                # this must be added AFTER counting dups and losses since it affects loss inference
                if self.implied:
                    added = phylo.add_implied_spec_nodes(proposal.locus_tree, self.stree,
                                                         proposal.locus_recon, proposal.locus_events)

                ncoal = phyloDLC.count_coal(self.coal_tree, proposal.coal_recon, proposal.locus_tree)
                coalcost = ncoal * self.coalcost

                if self.implied:
                    phylo.remove_implied_spec_nodes(proposal.locus_tree, added,
                                                    proposal.locus_recon, proposal.locus_events)

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
        """Evaluate a proposal for search"""
        
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
        self.log_stream.write(repr(proposal) + "\n")
        self.log_stream.flush()



class DLCReconProposer (object):

    def __init__(self, coal_tree, stree, gene2species,
                 search=phylo.TreeSearchNni):
        self._coal_tree = coal_tree
        self._stree = stree
        self._gene2species = gene2species
        self._locus_search = search(None)

        self._accept_locus = False

        self._recon = None
        

    def set_locus_tree(self, locus_tree):
        self._locus_search.set_tree(locus_tree)
        

    def init_proposal(self):
        """Get first proposal"""

        if self._locus_search.get_tree() is None:
            self._locus_search.set_tree(self._coal_tree.copy())
        self._recon = self._recon_lca(self._locus_search.get_tree().copy())
        
        return self._recon


    def next_proposal(self):

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
                         newCopy=False)
        common.rename_nodes(locus_tree)

        # propose remaining parts of dlcoal recon
        self._recon = self._recon_lca(locus_tree)

        return self._recon


    def _recon_lca(self, locus_tree):
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
        return propose_daughters(coal_tree, coal_recon,
                                 locus_tree, locus_events)


    def accept(self):
        self._accept_locus = True


    def reject(self):
        pass


def propose_daughters(coal_tree, coal_recon, locus_tree, locus_events):

    lineages = coal.count_lineages_per_branch(coal_tree, coal_recon, locus_tree)
    daughters = set()

    for node, event in locus_events.iteritems():
        if event == "dup":
            # choose one of the children of node to be a daughter
            children = [child for child in node.children
                        if lineages[child][1] == 1]
            if len(children) > 0:
                daughters.add(children[stats.sample([1] * len(children))])

    return daughters


#=============================================================================
# tree search

class LocusTreeSearchPrescreen (phylo.TreeSearchPrescreen):

    def propose(self):

        # save old topology
        self.oldtree = self.tree.copy()

        pool = []
        best_score = self.prescreen(self.tree)
        total = -util.INF

        # TODO: add unique filter

        # make many subproposals
        self.search.reset()
        for i in xrange(self.poolsize):
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
        trees, minscore = util.minall(pool, keyfunc=lambda it: it[0], minfunc=lambda it: it[1])
        tree = random.choice(trees)
        treelib.set_tree_topology(self.tree, tree)

        self.search.add_seen(self.tree)


class LocusTreeSearch (phylo.TreeSearch):

    def __init__(self, tree, stree, gene2species,
                 dupcost, losscost,
                 tree_hash=None, nprescreen=20, weight=.2):
        phylo.TreeSearch.__init__(self, tree)

        self.stree = stree
        self.gene2species = gene2species

        self.dupcost = dupcost
        self.losscost = losscost
        
        #self.search = UniqueTreeSearch(tree, phylo.TreeSearchNni(tree),
        #                               tree_hash)
        #self.search = UniqueTreeSearch(tree, phylo.TreeSearchSpr(tree),
        #                               tree_hash)

        mix = phylo.TreeSearchMix(tree)
        mix.add_proposer(phylo.TreeSearchNni(tree), .4)
        mix.add_proposer(phylo.TreeSearchSpr(tree), .6)        
        #self.search = phylo.TreeSearchUnique(tree, mix, tree_hash)

        if nprescreen == 1:
            self.search = mix
        else:
            prescreen = LocusTreeSearchPrescreen(tree, mix,
                                                 self.prescreen,
                                                 poolsize=nprescreen)
            
            mix2 = phylo.TreeSearchMix(tree)
            mix2.add_proposer(prescreen, 1.0-weight)
            mix2.add_proposer(mix, weight)

            self.search = mix2
        
    def set_tree(self, tree):
        self.tree = tree
        self.search.set_tree(tree)

    def reset(self):
        self.search.reset()

    def propose(self):
        self.search.propose()
        return self.tree
        
    def revert(self):
        self.search.revert()
        return self.tree

    def prescreen(self, tree):
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
            nloss = phylo.count_loss(tree, self.stree, recon)
            losscost = nloss * self.losscost

        return dupcost + losscost
    

