"""
   Code for the DLC ILP Reconciliation
   (duplications, losses, and coalescence)
"""

from compbio import phylo
from dlcpar import common, reconlib
from dlcpar.recon import *
from dlcpar import ilprecon_solver
from rasmus import treelib, util
from recon import _random_choice
# integer linear programming
from pulp import *


# The following attributes in DLCLPRecon correspond to variables described DLCLP paper
# gtree = T_g (but with implied speciation nodes)
# stree = T_s
# srecon = M
# lrecon = L


def ilp_recon(tree, stree, gene2species,
              dupcost=1, losscost=1, coalcost=1, coaldupcost=None,
              delay=True, log=sys.stdout):
    """Perform reconciliation using DLCoal model with parsimony costs"""

    reconer = DLCLPRecon(tree, stree, gene2species,
                         dupcost=dupcost, losscost=losscost, coalcost=coalcost, coaldupcost=coaldupcost,
                         delay=delay, log=log)
    return reconer.recon()


class DLCLPRecon(DLCRecon):

    def __init__(self, gtree, stree, gene2species,
                 dupcost=1, losscost=1, coalcost=1, coaldupcost=None,
                 delay=True, name_internal="n", log=sys.stdout):

        # rename gene tree nodes
        common.rename_nodes(gtree, name_internal)

        self.gtree = gtree
        self.stree = stree
        self.gene2species = gene2species

        assert (dupcost >= 0) and (losscost >= 0) and (coalcost >= 0) and (coaldupcost >= 0), (dupcost, losscost, coalcost, coaldupcost)
        self.dupcost = dupcost
        self.losscost = losscost
        self.coalcost = coalcost  # actually coalspeccost, using coalcost for backwards compatibility
        self.coaldupcost = coaldupcost if coaldupcost is not None else coalcost

        self.delay = delay


        self.name_internal = name_internal
        self.log = util.Timer(log)


    def recon(self):
        """Perform reconciliation

        assigns the following attributes:
            self.ilp
            self.srecon
            self.lrecon
            self.order
            self.nsoln
            self.cost

        Returns:
            a tuple with the gene tree, a LabeledRecon, the runtime, and the cost
        """

        self.log.start("Reconciling")

        # log input gene and species trees
        # log_tree(self.gtree, self.log, func=treelib.draw_tree_names)
        # log_tree(self.stree, self.log, func=treelib.draw_tree_names)

        self._infer_species_map()

        # add implied speciation nodes but first start the species tree at the right root
        substree = treelib.subtree(self.stree, self.srecon[self.gtree.root])
        subsrecon = util.mapdict(self.srecon, val=lambda snode: substree.nodes[snode.name])

        # switch internal storage with subtrees
        self.stree, substree = substree, self.stree
        self.srecon, subsrecon = subsrecon, self.srecon

        # add implied nodes (standard speciation, speciation from duplication, delay nodes)
        # then relabel events (so that factor_tree works)
        reconlib.add_implied_nodes(self.gtree, self.stree, self.srecon, self.sevents, delay=self.delay)

        self.sevents = phylo.label_events(self.gtree, self.srecon)
        common.rename_nodes(self.gtree, self.name_internal)

        # log gene tree (with species map)
        # log_tree(self.gtree, self.log, func=draw_tree_srecon, srecon=self.srecon)

        # find subtrees
        self.subtrees = reconlib.factor_tree(self.gtree, self.stree, self.srecon, self.sevents)

        ilp, lpvars = ilprecon_solver.solve(self.gtree, self.srecon, self.subtrees,
                                            self.dupcost, self.losscost, self.coalcost, self.coaldupcost)

        dup_nodes = [g for g, dup_var in lpvars.dup_vars.iteritems() if dup_var.varValue == 1.0]

        self.cost = value(ilp.objective)

        self.lrecon = self._infer_locus_map(dup_nodes)

        # optimum_order = self._infer_opt_order(dup_nodes)

        # log gene tree (with species map and locus map)
        # log_tree(self.gtree, self.log, func=draw_tree_recon, srecon=self.srecon, lrecon=self.lrecon)

        # revert to use input species tree
        self.stree = substree
        self.srecon = util.mapdict(self.srecon, val=lambda snode: self.stree.nodes[snode.name])
        # TODO make order
        self.order = util.mapdict({}, key=lambda snode: self.stree.nodes[snode.name])

        labeled_recon = reconlib.LabeledRecon(self.srecon, self.lrecon, self.order)

        draw_tree_recon(self.gtree, self.srecon, self.lrecon, minlen= 15, maxlen = 15)

        runtime = self.log.stop()

        return self.gtree, labeled_recon, runtime, self.cost

    def _infer_locus_map(self, dup_nodes):
        """
        Infer (and assign) locus map

        :param list dup_nodes: gene tree nodes whose parents are mapped to different loci
        """
        locus = 1
        lrecon = {}

        for gnode in self.gtree.preorder():
            # if it's the root, it has the first locus
            if gnode == self.gtree.root:
                lrecon[gnode] = locus
            # if it has a dup, it has a different locus
            elif gnode in dup_nodes:
                locus += 1
                lrecon[gnode] = locus
            # if there's no dup, it's the same as the parent
            else:
                lrecon[gnode] = lrecon[gnode.parent]

        return lrecon
