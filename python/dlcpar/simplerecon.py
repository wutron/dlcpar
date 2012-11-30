# based on dlcoal recon.py

import sys, copy
import random

import dlcpar
from dlcpar import common

from rasmus import util, stats, treelib
from compbio import phylo, coal
from yjw.bio import phyloDLC

#=============================================================================
# reconciliation

def dlc_recon(tree, stree, gene2species,
              dupcost=1, losscost=1, coalcost=1,
              nsearch=1000,
              search=None,
              init_locus_tree=None,
              log=sys.stdout):
    """Perform reconciliation using DLCoal model with parsimony costs"""
    if search is None:
        search = lambda tree: LocusTreeSearch(tree, stree, gene2species)

    reconer = DLCRecon(tree, stree, gene2species,
                       dupcost=dupcost, losscost=losscost, coalcost=coalcost,
                       init_locus_tree=init_locus_tree,
                       log=log)
    reconer.set_proposer(DLCReconProposer(
        tree, stree, gene2species, search=search))
    return reconer.recon(nsearch).get_dict()


class DLCRecon (object):

    def __init__(self, tree, stree, gene2species,
                 dupcost=1, losscost=1, coalcost=1,
                 init_locus_tree=None,
                 name_internal="n", log=sys.stdout):

        self.coal_tree = tree
        self.stree = stree
        self.gene2species = gene2species

        self.dupcost = dupcost
        self.losscost = losscost
        self.coalcost = coalcost
        
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
            if i % 10 == 0:
                print "search", i

            util.tic("eval")
            cost = self.eval_proposal(proposal)
            util.print_dict(proposal.data)
            util.toc()

            util.tic("prop")
            self.eval_search(cost, proposal)
            proposal = self.proposer.next_proposal()
            util.toc()
        
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


    def next_proposal(self):
        """Returns next proposal"""
        self.proposal.next_proposal()


    def eval_proposal(self, proposal):
        """Compute cost of proposal"""

        # find dup/loss cost
        ndup = phylo.count_dup(proposal.locus_tree, proposal.locus_events)
        nloss = phylo.count_loss(proposal.locus_tree, self.stree, proposal.locus_recon)

        # find coal cost
	if phyloDLC.assert_bounded_coal(self.coal_tree, proposal.coal_recon, proposal.daughters):
	    ncoal = phyloDLC.count_coal(self.coal_tree, proposal.coal_recon, proposal.locus_tree)
	else:
	    ncoal = util.INF

        # total cost
        cost = ndup*self.dupcost + nloss*self.losscost + ncoal*self.coalcost

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
            self.proposer.accept()
        else:
            self.proposer.reject()


    def log_proposal(self, proposal):
        self.log_stream.write(repr(proposal) + "\n")
        self.log_stream.flush()



class DLCReconProposer (object):

    def __init__(self, coal_tree, stree, gene2species,
                 search=phylo.TreeSearchNni,
                 num_coal_recons=1): # DEBUG
        self._coal_tree = coal_tree
        self._stree = stree
        self._gene2species = gene2species
        self._locus_search = search(None)

        # coal recon search
        self._num_coal_recons = num_coal_recons
        self._i_coal_recons = 0
        self._coal_recon_enum = None
        self._coal_recon_depth = 2
        self._accept_locus = False

        self._recon = None
        

    def set_locus_tree(self, locus_tree):
        self._locus_search.set_tree(locus_tree)

    def init_proposal(self):
        """Get first proposal"""

        if self._locus_search.get_tree() is None:
            self._locus_search.set_tree(self._coal_tree.copy())
        self._i_coal_recons = 0
        self._recon = self._recon_lca(self._locus_search.get_tree().copy())
        
        return self._recon


    def next_proposal(self):

        if len(self._locus_search.get_tree().leaves()) <= 2:
            return self._recon
        
        if self._i_coal_recons >= self._num_coal_recons:
            # propose new locus_tree
            
            # if locus_tree has not yet been accepted, then revert it
            if not self._accept_locus:
                self._locus_search.revert()
                
            self._locus_search.propose()
            self._accept_locus = False
            self._i_coal_recons = 0
            locus_tree = self._locus_search.get_tree().copy()
            
            # TODO: make recon root optional
            phylo.recon_root(locus_tree, self._stree,
                             self._gene2species,
                             newCopy=False)
            common.rename_nodes(locus_tree)

            # propose remaining parts of dlcoal recon
            self._recon = self._recon_lca(locus_tree)
        else:
            # modify coal_recon
            
            try:
                self._i_coal_recons += 1
                self._coal_recon_enum.next()
            except StopIteration:
                self._i_coal_recon = self._num_coal_recons
                return self.next_proposal()

        return self._recon


    def _recon_lca(self, locus_tree):
        # get locus tree, and LCA locus_recon
        locus_recon = phylo.reconcile(locus_tree, self._stree,
                                      self._gene2species)
        locus_events = phylo.label_events(locus_tree, locus_recon)

        # propose LCA coal_recon
        coal_recon = phylo.reconcile(self._coal_tree,
                                     locus_tree, lambda x: x)

        # propose daughters
        daughters = self._propose_daughters(
            self._coal_tree, coal_recon,
            locus_tree, locus_recon, locus_events)


        self._coal_recon_enum = phylo.enum_recon(
            self._coal_tree, locus_tree,
            recon=coal_recon,
            depth=self._coal_recon_depth)


        return Recon(coal_recon, locus_tree, locus_recon, locus_events,
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


class Recon (object):
    """
    The reconciliation datastructure for the DLCoal model
    """
    
    def __init__(self, coal_recon, locus_tree, locus_recon, locus_events,
                 daughters, data=None):
        self.coal_recon = coal_recon
        self.locus_tree = locus_tree
        self.locus_recon = locus_recon
        self.locus_events = locus_events
        self.daughters = daughters

        if data is None:
            self.data = {}
        else:
            self.data = data

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

#=============================================================================
# tree search

class LocusTreeSearch (phylo.TreeSearch):

    def __init__(self, tree, stree, gene2species, 
                 tree_hash=None):
        phylo.TreeSearch.__init__(self, tree)

        self.stree = stree
        self.gene2species = gene2species
        
        #self.search = UniqueTreeSearch(tree, phylo.TreeSearchNni(tree),
        #                               tree_hash)
        #self.search = UniqueTreeSearch(tree, phylo.TreeSearchSpr(tree),
        #                               tree_hash)

        mix = phylo.TreeSearchMix(tree)
        mix.add_proposer(phylo.TreeSearchNni(tree), .4)
        mix.add_proposer(phylo.TreeSearchSpr(tree), .6)        
        #self.search = phylo.TreeSearchUnique(tree, mix, tree_hash)

        self.search = mix
        
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

