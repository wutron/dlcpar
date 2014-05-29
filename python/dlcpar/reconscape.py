"""

   Code for the DLC Parsimony Reconciliation Landscape

"""

# python libraries
import sys
import collections
import re

# geometry libraries
from shapely import geometry
from shapely.wkt import dumps, loads

# rasmus libraries
from rasmus import treelib, util
from compbio import phylo

# yjw libraries
from yjw import combinatorics

# dlcpar libraries
from dlcpar import common
from dlcpar import reconlib

from dlcpar.recon import *
from dlcpar.countvector import *

#==========================================================
# globals

DEFAULT_RANGE = (0.2, 10)
DEFAULT_RANGE_STR = '-'.join(map(str, DEFAULT_RANGE))

#==========================================================
# reconciliation

def dlc_recon(tree, stree, gene2species, gene2locus=None,
              duprange=DEFAULT_RANGE, lossrange=DEFAULT_RANGE,
              implied=True, delay=True,
              max_loci=INF, max_dups=INF, max_losses=INF,
              log=sys.stdout):
    """Perform reconciliation using DLCoal model with parsimony costs"""

    reconer = DLCRecon(tree, stree, gene2species, gene2locus,
                       duprange=duprange, lossrange=lossrange,
                       max_loci=max_loci, max_dups=max_dups, max_losses=max_losses,
                       log=log)
    return reconer.recon()
    

class DLCRecon(DLCRecon):

    def __init__(self, gtree, stree, gene2species, gene2locus=None,
                 duprange=DEFAULT_RANGE, lossrange=DEFAULT_RANGE,
                 max_loci=INF, max_dups=INF, max_losses=INF,
                 name_internal="n", log=sys.stdout):

        # rename gene tree nodes
        common.rename_nodes(gtree, name_internal)
        
        self.gtree = gtree
        self.stree = stree
        self.gene2species = gene2species
        self.gene2locus = gene2locus

        dup_min, dup_max = duprange
        loss_min, loss_max = lossrange
        assert (dup_min > 0) and (dup_max > 0) and (dup_min < dup_max) and \
               (loss_min > 0) and (loss_max > 0) and (loss_min < loss_max)
        self.duprange = duprange
        self.lossrange = lossrange
        
        self.implied = True
        self.delay = False
        self.prescreen = False
       
        assert (max_loci > 0) and (max_dups > 0) and (max_losses > 0)
        self.max_loci = max_loci
        self.max_dups = max_dups
        self.max_losses = max_losses

        self.name_internal = name_internal
        self.log = util.Timer(log)

    #=============================
    # main methods

    def recon(self):
        """Perform reconciliation"""

        self.log.start("Reconciling")

        # log input gene and species trees
        self.log.log("gene tree\n")
        log_tree(self.gtree, self.log, func=treelib.draw_tree_names)
        self.log.log("species tree\n")
        log_tree(self.stree, self.log, func=treelib.draw_tree_names)
               
        # infer species map
        self._infer_species_map()
        self.log.log("\n\n")

        # add implied speciation nodes but first start the species tree at the right root
        substree = treelib.subtree(self.stree, self.srecon[self.gtree.root])
        subrecon = util.mapdict(self.srecon, val=lambda snode: substree.nodes[snode.name])
        
        # switch internal storage with subtrees
        self.stree, subtree = substree, self.stree
        self.srecon, subrecon = subrecon, self.srecon

        # add implied nodes (standard speciation, speciation from duplication, delay nodes)
        # then relabel events (so that factor_tree works)
        reconlib.add_implied_nodes(self.gtree, self.stree, self.srecon, self.sevents, delay=self.delay)
        self.sevents = phylo.label_events(self.gtree, self.srecon)
        common.rename_nodes(self.gtree, self.name_internal)

        # log gene tree (with species map)
        self.log.log("gene tree (with species map)\n")
        log_tree(self.gtree, self.log, func=draw_tree_srecon, srecon=self.srecon)

        # infer locus map
        self._infer_locus_map()
        self.log.log("\n\n")        

        self.log.stop()
        
        return self.count_vectors


    #=============================
    # event/cost methods -- these operate at the species branch level

    def _count_events(self, lrecon, subtrees, nodefunc=lambda node: node.name,
                      all_leaves=None,
                      max_dups=INF, max_losses=INF,
                      min_cvs=None):
        """
        Count number of dup, loss, coal events

        Returns
        - number of duplications
        - number of losses
        - minimum number of extra lineages over all internal node orderings
        - the optimal internal node ordering
        """
        
        extra = {"species_map" : self.srecon, "locus_map" : lrecon}

        # defaults
        ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln = INF, INF, INF, INF, {}, INF

        # duplications
##        ndup = reconlib.count_dup_snode(self.gtree, self.stree, extra, snode=None,
##                                        subtrees_snode=subtrees,
##                                        nodefunc=nodefunc)
        dup_nodes = reconlib.find_dup_snode(self.gtree, self.stree, extra, snode=None,
                                            subtrees_snode=subtrees,
                                            nodefunc=nodefunc)
        ndup = len(dup_nodes)
        if ndup > max_dups:     # skip rest if exceed max_dups
            return ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln

        # losses
        nloss = reconlib.count_loss_snode(self.gtree, self.stree, extra, snode=None,
                                          subtrees_snode=subtrees,
                                          nodefunc=nodefunc)
        if nloss > max_losses:  # skip rest if exceed max_losses
            return ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln
        
        # extra lineages at speciations
        ncoal_spec = reconlib.count_coal_snode_spec(self.gtree, self.stree, extra, snode=None,
                                                    subtrees_snode=subtrees,
                                                    nodefunc=nodefunc,
                                                    implied=self.implied)
        if (min_cvs is not None) and is_maximal(CountVector(ndup, nloss, ncoal_spec), min_cvs):  # skip rest if already not Pareto-optimal
            return ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln
        
        # extra lineages at duplications
        ncoal_dup, order, nsoln = self._count_min_coal_dup(lrecon, subtrees, nodefunc=nodefunc,
                                                           dup_nodes=dup_nodes, all_leaves=all_leaves)

        return ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln



    #=============================
    # locus partition methods (used by _enumerate_locus_maps)
    # partition structure: key1 = snode, key2 = bottom_loci, key3 = top_loci, value = CountVectorSet
    
    def _initialize_partitions_sroot(self, partitions, bottom_loci, top_loci):
        cvs = CountVectorSet([ZERO_COUNT_VECTOR])
        partitions[bottom_loci][top_loci] = cvs


    def _find_optimal_cost(self, partitions, bottom_loci, top_loci,
                           lrecon, subtrees, leaves,
                           max_dups, max_losses):
        mincvs = None
        if (bottom_loci in partitions) and (top_loci in partitions[bottom_loci]):
            mincvs = partitions[bottom_loci][top_loci]

        ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln = \
              self._count_events(lrecon, subtrees, all_leaves=leaves,
                                 max_dups=max_dups, max_losses=max_losses,
                                 min_cvs=mincvs)
        ncoal = ncoal_spec + ncoal_dup

        return ndup, nloss, ncoal_spec, ncoal_dup, ncoal, order, nsoln


    def _update_partitions(self, partitions, bottom_loci, top_loci,
                           lrecon, order,
                           ndup, nloss, ncoal, nsoln):
        if top_loci not in partitions[bottom_loci]:
            partitions[bottom_loci][top_loci] = CountVectorSet()
        cv = CountVector(ndup, nloss, ncoal, nsoln)
        partitions[bottom_loci][top_loci].add(cv)
        

    def _filter_partitions(self, partitions):
        self.log.log("optimal count vectors")
        for bottom_loci, d in partitions.iteritems():
            for top_loci, cvs in d.iteritems():
                # filter down to set of Pareto-optimal count vectors
                new_cvs = cvs.pareto_filter(self.duprange, self.lossrange)
                partitions[bottom_loci][top_loci] = new_cvs

                # log
                self.log.log("\t%s -> %s" % (top_loci, bottom_loci))
                for cv in new_cvs:
                    self.log.log("\t\tvector: %s" % cv)
                self.log.log()


    #=============================
    # DP table methods (used by _infer_opt_locus_map)

    def _dp_table(self, locus_maps, subtrees):
        # locus_maps is a multi-dimensional dict with the structure
        # key1 = snode, key2 = bottom_loci, key3 = top_loci, value = CountVectorSet

        stree = self.stree

        # dynamic programming storage
        F = {}      # key1 = snode, key2 = top_loci, val = cost-to-go

        # recur up species tree to determine (optimal) cost-to-go
        for snode in stree.postorder():
            self.log.start("Working on snode %s" % snode.name)
            F[snode] = {}

            if snode not in locus_maps:
                # nothing in this sbranch
                F[snode][()] = CountVectorSet([ZERO_COUNT_VECTOR])
                self.log.log("Empty sbranch")
                self.log.stop()
                continue

            # get stored values for the sbranch
            locus_maps_snode = locus_maps[snode]
            subtrees_snode = subtrees[snode]

            if snode.is_leaf():
                # leaf base case
                for bottom_loci, d in locus_maps_snode.iteritems():
                    for top_loci, cvs in d.iteritems():
                        assert top_loci not in F[snode]
                        F[snode][top_loci] = cvs
            else:
                if len(snode.children) != 2:
                    raise Exception("non-binary species tree")

                # find cost-to-go of assigning top_loci to top of sbranch
                # = cost of top_loci to bottom_loci along sbranch
                #   + cost of bottom_loci at top of left child
                #   + cost of bottom_loci at top of right child
                sleft, sright = snode.children
                costs = collections.defaultdict(CountVectorSet) # separate CostVectorSet for each assignment of top_loci for this sbranch
                for bottom_loci, d in locus_maps_snode.iteritems():
                    # find cost-to-go in children
                    # locus assignment may have been removed due to search heuristics
                    if bottom_loci in F[sleft]:
                        cvs_left = F[sleft][bottom_loci]
                    else:
                        cvs_left = CountVectorSet([MAX_COUNT_VECTOR])
                    if bottom_loci in F[sright]:
                        cvs_right = F[sright][bottom_loci]
                    else:
                        cvs_right = CountVectorSet([MAX_COUNT_VECTOR])
                    children_cvs = cvs_left * cvs_right

                    # add cost in this sbranch
                    for top_loci, cvs in d.iteritems():
                        cvs_to_go = cvs * children_cvs
                        costs[top_loci].update(cvs_to_go)
                        
                # for each assignment of top_loci to top of sbranch,
                # filter down to set of Pareto-optimal count vectors
                for top_loci, cvs in costs.iteritems():
                    F[snode][top_loci] = cvs.pareto_filter(self.duprange, self.lossrange)

            self.log.log("DP table")
            for top_loci, cvs in F[snode].iteritems():
                self.log.log("\t%s :" % str(top_loci))
                for cv in cvs:
                    self.log.log("\t\tvector: %s" % cv)
            self.log.stop()

        return F

            
    def _dp_terminate(self, F):
        stree = self.stree
        
        assert len(F[stree.root]) == 1, F[stree.root]
        cvs = F[stree.root].values()[0]
        self.log.log("")
        self.log.log("Optimal count vectors:")
        for cv in cvs:
            self.log.log("\tvector: %s" % cv)
        self.log.log("")

        self.count_vectors = cvs


    def _dp_traceback(self, locus_maps, subtrees, F):
        pass


#==========================================================
# regions

def get_regions(cvs, duprange, lossrange, restrict=True):
    """
    Return dictionary mapping cost vectors to regions.

    cvs       -- Pareto-optimal CostVectorSet
    duprange  -- (low, high) relative to unit cost of (deep) coalescence
    lossrange -- (low, high) relative to unit cost of (deep) coalescence
    restrict  -- Keep only regions in the specified range
    """

    regions = {}

    dup_min, dup_max = duprange
    loss_min, loss_max = lossrange

    # bounding box
    boundingbox = [(dup_min, loss_min), (dup_min, loss_max),
                   (dup_max, loss_max), (dup_max, loss_min)]
    bb = geometry.Polygon(boundingbox)

    # find region for each event count vector
    # Let x and y denote the dup and loss costs (relative to deep coalescence cost).
    # A cost vector cv has cost = cv.d * x + cv.l * y + cv.c.
    # For cv1 to be optimal, then for each cv2 (s.t. cv2 != cv1), it must be that
    #     cv1.d * x + cv1.l * y + cv1.c <= cv2.d * x + cv2.l * y + cv2.c
    # that is
    #     y <= x * (cv2.d - cv1.d)/(cv1.l - cv2.l) + (cv2.c - cv1.c)/(cv1.l - cv2.l).
    # For the special case in which cv1.l == cv2.l,
    #     x <= (cv2.c - cv1.c) / (cv1.d - cv2.d)

    for cv1 in cvs:
        region = bb
        for cv2 in cvs:
            if cv2 == cv1:
                continue    # skip comparison

            if cv1.l == cv2.l:  # vertical line defines the half-space
                xint = 1.0*(cv2.c - cv1.c)/(cv1.d - cv2.d)  # x-intercept
                
                # the variable "below" is True iff half-space is to left of line
                below = cv1.d - cv2.d > 0
                
                # test if half-space is to left or right of the line
                if below:      # below
                    lowestx = min(xint, dup_min)
                    poly = [(xint, loss_min), (xint, loss_max),
                            (lowestx, loss_max), (lowestx, loss_min)]
                else:          # above
                    highestx = max(xint, dup_max)
                    poly = [(xint, loss_min), (xint, loss_max),
                            (highestx, loss_max), (highestx, loss_min)]    
            else:
                m = 1.0*(cv2.d - cv1.d)/(cv1.l - cv2.l)  # slope
                b = 1.0*(cv2.c - cv1.c)/(cv1.l - cv2.l)  # y-intercept
                
                # the variable "below" is True iff half-space is below line
                below = cv1.l - cv2.l > 0
                
                if m == 0:      # horizontal line defines the half-space
                    # test if half-space is below or above the line
                    if below:  # below
                        lowesty = min(b, loss_min)
                        poly = [(dup_min, b), (dup_max, b),
                                (dup_max, lowesty), (dup_min, lowesty)]
                    else:      # above
                        highesty = max(b, loss_max)
                        poly = [(dup_min, b), (dup_max, b),
                                (dup_max, highesty), (dup_min, highesty)]
                else:           # "slanted" line defines the half-space
                    # y-coord of intersection with left/right edge of boundingbox 
                    lefty = m * dup_min + b
                    righty = m * dup_max + b
                    
                    # test if half-space is below or above the line
                    if below:  # below
                        lowesty = min(loss_min, lefty, righty)
                        poly = [(dup_min, lefty), (dup_max, righty),
                                (dup_max, lowesty), (dup_min, lowesty)]
                    else:      # above
                        highesty = max(loss_max, lefty, righty)
                        poly = [(dup_min, lefty), (dup_max, righty),
                                (dup_max, highesty), (dup_min, highesty)]
                
            # update region
            P = geometry.Polygon(poly)
            region = region.intersection(P)

        # keep region?
        if restrict:
            if (region.is_empty) or (not region.intersects(bb)):
                continue           
        regions[cv1] = region

    # sort
    regions = collections.OrderedDict(sorted(regions.items(),
                                             key=lambda it: it[0],
                                             cmp=CountVector.lex))

    return regions


def read_regions(filename):
    regions = collections.OrderedDict()
    
    for i, toks in enumerate(util.DelimReader(filename)):
        if i == 0:
            duprange, lossrange = map(float, toks[:2]), map(float, toks[2:])
            continue

        cv_str, coords_str, area_str = toks
        cv = parse_count_vector(cv_str)
        region = loads(coords_str)
        area = float(area_str)
        regions[cv] = region
                    
    return regions, duprange, lossrange


def write_regions(filename, regions, duprange, lossrange):
    out = util.open_stream(filename, 'w')
    print >>out, '\t'.join(map(str, duprange + lossrange))
    for cv, region in regions.iteritems():
        coords = None; area = None
        if isinstance(region, geometry.Polygon):                                              # non-degenerate
            coords = list(region.exterior.coords)
            area = region.area
        elif isinstance(region, geometry.LineString) or isinstance(region, geometry.Point):   # degenerate
            coords = list(region.coords)
            area = region.area
        else:                                                                                 # non-degenerate (collection)
            coords = []; area = 0
            for r in region:
                if isinstance(r, geometry.Polygon):                                           # non-degenerate
                    c = list(r.exterior.coords)
                elif isinstance(r, geometry.LineString) or isinstance(r, geometry.Point):     # degenerate
                    c = list(r.coords)
                else:
                    raise Exception("cost vector (%s) has invalid subregion (%s)" % (cv, type(r)))
                coords.append(c)
                area += r.area
                
        coords = dumps(region)
        toks = (cv, coords, area)
        print >>out, '\t'.join(map(str, toks))
    out.close()

#==========================================================
# visualization

def draw_landscape(regions, duprange, lossrange,
                   filename=sys.stdout,
                   log=False):
    """
    Plot the landscape.
    
    The x-axis represents dup cost and the y-axis represents the loss cost
    (both relative to unit cost of deep coalescence).
    """
    
    from rasmus import plotting
    import matplotlib.pyplot as plt

    dup_min, dup_max = duprange
    loss_min, loss_max = lossrange

    # axes
    if log:
        plt.xscale('log')
        plt.yscale('log')
    plt.axis('scaled')
    plt.axis([dup_min, dup_max, loss_min, loss_max])
    plt.xlabel("relative duplication cost")
    plt.ylabel("relative loss cost")
    
    # color map
    n = len(regions)
    colormap = plotting.rainbow_color_map(low=0, high=n-1)

    # plot regions
    for i, (cv, region) in enumerate(regions.iteritems()):
       
        # output
        label = str(cv)
        color = colormap.get(i)

        if isinstance(region, geometry.Polygon):        # non-degenerate
            coords = list(region.exterior.coords)
            plt.gca().add_patch(plt.Polygon(coords,
                                            label=label, color=color))
        elif isinstance(region, geometry.LineString):   # degenerate
            coords = list(region.coords)
            plt.plot([coords[0][0], coords[1][0]],
                     [coords[0][1], coords[1][1]],
                     linewidth=4,
                     label=label, color=color)
        elif isinstance(region, geometry.Point):        # degenerate
            coords = list(region.coords)
            plt.plot(coords[0][0], coords[0][1],
                     'o', markersize=4,
                     label=label, color=color)
        else:                                           # non-degenerate (collection)
            for r in region:
                if isinstance(r, geometry.Polygon):         # non-degenerate
                    coords = list(r.exterior.coords)
                    plt.gca().add_patch(plt.Polygon(coords,
                                                    label=label, color=color))
                elif isinstance(r, geometry.LineString):    # degenerate
                    coords = list(r.coords)
                    plt.plot([coords[0][0], coords[1][0]],
                             [coords[0][1], coords[1][1]],
                             linewidth=4,
                             label=label, color=color)
                elif isinstance(r, geometry.Point):         # degenerate
                    coords = list(r.coords)
                    plt.plot(coords[0][0], coords[0][1],
                             'o', markersize=4,
                             label=label, color=color)
                else:
                     raise Exception("cost vector (%s) has invalid subregion (%s)" % (cv, type(r)))
    
    # legend
    leg = plt.legend()
    for i in range(len(leg.legendHandles)):  # adjust legend marker thickness
        leg.legendHandles[i].set_linewidth(2.0)
    
    if filename:
        plt.savefig(filename, format="pdf")
    else:
        plt.show()
