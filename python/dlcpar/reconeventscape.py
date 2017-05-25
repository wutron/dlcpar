"""

   Code for the DLC Parsimony Reconciliation Landscape

"""

# python libraries
import sys
import collections
import itertools
from collections import Counter

# geometry libraries
from shapely import geometry
from shapely.ops import cascaded_union
from shapely.wkt import dumps, loads
from fractions import Fraction

# rasmus libraries
from rasmus import treelib, util
from compbio import phylo

# dlcpar libraries
from dlcpar import common
from dlcpar import reconlib

# output
import csv

from dlcpar.recon import *
from dlcpar.countvector import *

#==========================================================
# globals

DEFAULT_RANGE = (0.2, 5)
DEFAULT_RANGE_STR = '-'.join(map(str, DEFAULT_RANGE))

#==========================================================
# reconciliation

def dlceventscape_recon(tree, stree, gene2species, gene2locus=None,
                   duprange=DEFAULT_RANGE, lossrange=DEFAULT_RANGE,
                   implied=True, delay=True,
                   max_loci=INF, max_dups=INF, max_losses=INF,
                   log=sys.stdout, allow_both=False,
                   intersect=False, outfile="data.csv"):
    """Perform reconciliation using DLCoal model with parsimony costs"""

    reconer = DLCScapeRecon(tree, stree, gene2species, gene2locus,
                       duprange=duprange, lossrange=lossrange,
                       max_loci=max_loci, max_dups=max_dups, max_losses=max_losses,
                       log=log, allow_both=allow_both, intersect=intersect, outfile=outfile)
    return reconer.recon()


class DLCScapeRecon(DLCRecon):

    def __init__(self, gtree, stree, gene2species, gene2locus=None,
                 duprange=DEFAULT_RANGE, lossrange=DEFAULT_RANGE,
                 max_loci=INF, max_dups=INF, max_losses=INF,
                 name_internal="n", log=sys.stdout, allow_both=False,
                 intersect=False, outfile="data.csv"):

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
        self.allow_both = allow_both

        self.intersect = intersect
        self.outfile = outfile

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

        self.log.stop()

        self.output()

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
        events = Counter()
        # TODO:can we take in the events for ndup, nloss.

        # duplications
##        ndup = reconlib.count_dup_snode(self.gtree, self.stree, extra, snode=None,
##                                        subtrees_snode=subtrees,
##                                        nodefunc=nodefunc)
        dup_nodes = reconlib.find_dup_snode(self.gtree, self.stree, extra, snode=None,
                                            subtrees_snode=subtrees,
                                            nodefunc=nodefunc)
        # make the dup events
        for node in dup_nodes:
            events[("D", node)] = 1

        ndup = len(dup_nodes)
        if ndup > max_dups:     # skip rest if exceed max_dups
            return ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events

        # losses
        nloss, losses = reconlib.count_loss_snode(self.gtree, self.stree, extra, snode=None,
                                          subtrees_snode=subtrees,
                                          nodefunc=nodefunc)
        # make the loss events
        for loss in losses:
            event = ["L"]
            event.extend(loss)
            events[tuple(event)] = 1

        if nloss > max_losses:  # skip rest if exceed max_losses
            return ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events

        # extra lineages at speciations
        ncoal_spec, coal_lineages = reconlib.count_coal_snode_spec(self.gtree, self.stree, extra, snode=None,
                                                    subtrees_snode=subtrees,
                                                    nodefunc=nodefunc,
                                                    implied=self.implied)
        # make the coalescence events
        for lineage in coal_lineages:
            event = ["C"]
            event.extend(lineage)
            events[tuple(event)] = 1

        if (min_cvs is not None) and is_maximal(CountVector(ndup, nloss, ncoal_spec), min_cvs):  # skip rest if already not Pareto-optimal
            return ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events

        # extra lineages at duplications
        ncoal_dup, order, nsoln = self._count_min_coal_dup(lrecon, subtrees, nodefunc=nodefunc,
                                                           dup_nodes=dup_nodes, all_leaves=all_leaves)
        # make the speciation events
        nspec, speciation = reconlib.count_spec_snode(self.gtree,self.stree, extra, snode=None,
                                                    subtrees_snode=subtrees,
                                                    nodefunc=nodefunc)
        for locus, lineages in speciation.iteritems():
            event = ["S"]
            event.extend(lineages)
            events[tuple(event)] = 1

        return ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events



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

        ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events = \
              self._count_events(lrecon, subtrees, all_leaves=leaves,
                                 max_dups=max_dups, max_losses=max_losses,
                                 min_cvs=mincvs)
        # count coalescences due to duplication
        ncoal = ncoal_spec + ncoal_dup
        #ncoal = ncoal_spec 
        # don't count multiple orderings as different solutions
        #nsoln = 1


        return ndup, nloss, ncoal_spec, ncoal_dup, ncoal, order, nsoln, events


    def _update_partitions(self, partitions, bottom_loci, top_loci,
                           lrecon, order,
                           ndup, nloss, ncoal_spec, ncoal_dup, nsoln, events):
        if top_loci not in partitions[bottom_loci]:
            partitions[bottom_loci][top_loci] = CountVectorSet()
        #rmawhorter: here if you want, you can disregard the cost of coalescence due to duplication
        ncoal = ncoal_spec + ncoal_dup
        cv = CountVector(ndup, nloss, ncoal, nsoln, events)
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
                costs = collections.defaultdict(CountVectorSet) # separate CountVectorSet for each assignment of top_loci for this sbranch
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
        
        #TODO: F still has some pareto-optimal vectors that don't have a region of optimality
        # after this returns the dp table. Either costvector._filter isn't working properly
        # or something else strange is going on
        return F


    def _dp_terminate(self, F):
        stree = self.stree

        assert len(F[stree.root]) == 1, F[stree.root]
        cvs = F[stree.root].values()[0]
        for cv in cvs:
            print cv

        self.log.log("")
        self.log.log("Optimal count vectors:")
        for cv in cvs:
            self.log.log("\tvector: %s" % cv)
        self.log.log("")

        self.count_vectors = cvs


    def _dp_traceback(self, locus_maps, subtrees, F):
        pass


    def _union(self, cvs):  
        out = defaultdict(list)
        for cv in cvs:
            for event, count in cv.events.iteritems():
                #use every event
                out[cv].append(event)
        return out

    def _intersect(self, cvs):
        out = defaultdict(list)
        for cv in cvs:
            for event, count in cv.events.iteritems():
                #only accept events that occur in every MPR
                if count == cv.count:
                    out[cv].append(event)
        return out

    def output(self, regions=None):
        #TODO: might also want to compute the percentage of area of the reconscape
        # an event is present in, by calculating the area of the regions it is present in.
        # a better metric might be to weight each area's contribution by the frequency of the
        # event in that area.
        event_dict = {}
        if self.intersect:
            event_dict = self._intersect(self.count_vectors)
        else:
            event_dict = self._union(self.count_vectors)
        event_counts = Counter()

        # determine the number of regions that each event was in
        for cv, event_list in event_dict.iteritems():
                count_events = Counter(event_list)
                event_counts += count_events
        #self.outfile
        # open the output file
        ofile = open(self.outfile, "wb")
        writer = csv.writer(ofile, delimiter = ",")
        # write each vector with its associated events (union or intersection)
        for cv in self.count_vectors:
            l = [cv.d, cv.l, cv.c, cv.count]
            l.extend(event_dict[cv])
            writer.writerow(l)
        # write the events, in order of how many regions they appear in
        nregions = 0
        line = []
        for eventcount in event_counts.most_common():
            if eventcount[1] != nregions:
                writer.writerow(line)
                line = []
                nregions = eventcount[1]
                line.append(nregions)
            line.append(eventcount[0])
        # compute the fraction of area of the reconscape that an event is present in
        # weighted by the frequency of that event in each region
        if regions:
            # area of the reconscape
            #maybe it's better to just add the areas??
            recon_area = (self.duprange(1) - self.duprange(0)) * (self.lossrange[1] - self.lossrange[0])
            recon_shape_area = 0
            for cv, poly in regions.iteritems():
                recon_shape_area += poly.area
            #the shape area should agree with the actual area
            assert abs(recon_area - recon_shape_area) < 0.01, (recon_area, recon_shape_area)
            #TODO: finish this up


#==========================================================
# regions

def get_regions(cvs, duprange, lossrange, restrict=True,
                log=sys.stdout):
    """
    Return dictionary mapping count vectors to regions.

    cvs       -- Pareto-optimal CountVectorSet
    duprange  -- (low, high) relative to unit cost of (deep) coalescence
    lossrange -- (low, high) relative to unit cost of (deep) coalescence
    restrict  -- Keep only regions in the specified range
    """

    log = util.Timer(log)
    log.start("Finding regions")

    def Frac(numerator=0, denominator=None):
        return Fraction(numerator, denominator).limit_denominator()

    dup_min, dup_max = map(Frac, duprange)
    loss_min, loss_max = map(Frac, lossrange)

    # bounding box
    bb = geometry.box(dup_min, loss_min, dup_max, loss_max)

    # empty
    EMPTY = geometry.GeometryCollection()

    # find region for each event count vector
    # Let x and y denote the dup and loss costs (relative to deep coalescence cost).
    # A cost vector cv has cost = cv.d * x + cv.l * y + cv.c.
    # For cv1 to be optimal, then for each cv2 (s.t. cv2 != cv1), it must be that
    #     cv1.d * x + cv1.l * y + cv1.c <= cv2.d * x + cv2.l * y + cv2.c
    # that is
    #     y <= x * (cv2.d - cv1.d)/(cv1.l - cv2.l) + (cv2.c - cv1.c)/(cv1.l - cv2.l).
    # For the special case in which cv1.l == cv2.l,
    #     x <= (cv2.c - cv1.c) / (cv1.d - cv2.d)

    # find all lines
    lines_diag, lines_horiz, lines_vert = set(), set([loss_min, loss_max]), set([dup_min, dup_max])
    for cv1 in cvs:
        for cv2 in cvs:
            if cv2 == cv1:
                continue    # skip comparison

            if cv1.l == cv2.l:  # vertical line defines the half-space
                xint = Frac(cv2.c - cv1.c, cv1.d - cv2.d)  # x-intercept
                lines_vert.add(xint)
            else:
                m = Frac(cv2.d - cv1.d, cv1.l - cv2.l)  # slope
                b = Frac(cv2.c - cv1.c, cv1.l - cv2.l)  # y-intercept
                if m == 0:      # horizontal line defines the half-space
                    lines_horiz.add(b)
                else:           # "slanted" line defines the half-space
                    lines_diag.add((m, b))
    lines_diag, lines_horiz, lines_vert = map(list, (lines_diag, lines_horiz, lines_vert))

    # find all intersection points
    points = set()
    for (m1, b1), (m2, b2) in itertools.combinations(lines_diag, r=2):
        if m1 == m2:
            continue    # parallel lines
        # y1 = m1 * x + b1 = m2 * x + b2 = y2, that is, x = (b2 - b1) / (m1 - m2)
        x = Frac(b2 - b1, m1 - m2)
        y1 = m1 * x + b1
        y2 = m2 * x + b2
        assert y1 == y2, (y1, y2)
        points.add((x, y1))
    for (m, b), x in itertools.product(lines_diag, lines_horiz):
        y = m * x + b
        points.add((x, y))
    for (m, b), y in itertools.product(lines_diag, lines_vert):
        x = Frac(y - b, m)
        points.add((x, y))
    for x, y in itertools.product(lines_vert, lines_horiz):
        points.add((x, y))
    points = filter(lambda pt: bb.intersects(geometry.Point(pt)), points)

    def find_closest_point(point):
        min_points, min_dist = util.minall(points, minfunc=lambda pt: geometry.Point(point).distance(geometry.Point(pt)))
        if min_dist > 1e-10:
            return None     # no point found (probably a point collinear with neighbors)
        assert len(min_points) == 1, (point, min_points, min_dist)
        return min_points[0]

    # find regions
    log.start("Mapping polygons")
    regions = {}
    for cv1 in cvs:
        region = bb
        for cv2 in cvs:
            if cv2 == cv1:
                continue    # skip comparison

            if cv1.l == cv2.l:  # vertical line defines the half-space
                xint = Frac(cv2.c - cv1.c, cv1.d - cv2.d)  # x-intercept

                # the variable "below" is True iff half-space is to left of line
                below = cv1.d - cv2.d > 0

                # test if half-space is to left or right of the line
                if below:      # below
                    lowestx = min(xint, dup_min)
                    poly = geometry.box(lowestx, loss_min, xint, loss_max)
                else:          # above
                    highestx = max(xint, dup_max)
                    poly = geometry.box(xint, loss_min, highestx, loss_max)
            else:
                m = Frac(cv2.d - cv1.d, cv1.l - cv2.l)  # slope
                b = Frac(cv2.c - cv1.c, cv1.l - cv2.l)  # y-intercept

                # the variable "below" is True iff half-space is below line
                below = cv1.l - cv2.l > 0

                if m == 0:      # horizontal line defines the half-space
                    # test if half-space is below or above the line
                    if below:  # below
                        lowesty = min(b, loss_min)
                        poly = geometry.box(dup_min, lowesty, dup_max, b)
                    else:      # above
                        highesty = max(b, loss_max)
                        poly = geometry.box(dup_min, b, dup_max, highesty)
                else:           # "slanted" line defines the half-space
                    # y-coord of intersection with left/right edge of boundingbox
                    lefty = m * dup_min + b
                    righty = m * dup_max + b

                    # test if half-space is below or above the line
                    if below:  # below
                        lowesty = min(loss_min, lefty, righty)
                        coords = [(dup_min, lowesty), (dup_max, lowesty),
                                  (dup_max, righty), (dup_min, lefty)]
                    else:      # above
                        highesty = max(loss_max, lefty, righty)
                        coords = [(dup_min, lefty), (dup_max, righty),
                                  (dup_max, highesty), (dup_min, highesty)]
                    poly = geometry.Polygon(coords)

            # update region
            keep = False
            if region.intersects(poly):
                region = region.intersection(poly)

                # valid?
                if not region.is_empty:
                    # correction
                    try:
                        filtered_region = filter(lambda r: isinstance(r, geometry.Polygon), region)
                        assert len(filtered_region) == 1
                        region = filtered_region[0]
                    except:
                        pass

                    # finesse coordinates (due to floating-point approximations)
                    if isinstance(region, geometry.Polygon):
                        coords = list(region.exterior.coords)[:-1]

                        # find closest coordinates
                        coords = map(find_closest_point, coords)
                        coords = filter(lambda pt: pt is not None, coords)

                        # collapse coordinates
                        new_coords = [coords[0]]
                        for i in range(1, len(coords)):
                            p1, p2 = coords[i-1], coords[i]
                            if not geometry.Point(p1).almost_equals(geometry.Point(p2)):
                                new_coords.append(p2)
                        coords = new_coords

                        # make new region
                        if len(coords) > 2:
                            keep = True
                            region = geometry.Polygon(coords)
            if not keep:
                region = EMPTY
                break

        # keep only polygons
        if (not region.is_empty) and isinstance(region, geometry.Polygon):
            regions[cv1] = region
    log.stop()

    # go back and get lines and points
    log.start("Mapping lines and points")
    new_regions = collections.defaultdict(set)
    for cv, region in regions.iteritems():
        new_regions[cv].add(region)     # polygons
        coords = map(find_closest_point, region.exterior.coords)
        coords = filter(lambda pt: pt is not None, coords)
        for i in range(1, len(coords)): # lines
            x1, y1 = coords[i-1]
            min_cvs1, min_cost1 = util.minall(cvs, minfunc=lambda cv: cv.d * x1 + cv.l * y1 + cv.c)
            x2, y2 = coords[i]
            min_cvs2, min_cost2 = util.minall(cvs, minfunc=lambda cv: cv.d * x2 + cv.l * y2 + cv.c)
            min_cvs = set(min_cvs1) & set(min_cvs2)
            poly = geometry.LineString(coords[i-1:i+1])
            for min_cv in min_cvs:
                new_regions[min_cv].add(poly)
        for x, y in coords:             # points
            min_cvs, min_cost = util.minall(cvs, minfunc=lambda cv: cv.d * x + cv.l * y + cv.c)
            poly = geometry.Point(x, y)
            for min_cv in min_cvs:
                new_regions[min_cv].add(poly)
    for cv, geoms in new_regions.iteritems():
        region = cascaded_union(geoms)
        assert not region.is_empty, cv
        assert isinstance(region, geometry.Polygon) or \
               isinstance(region, geometry.LineString) or \
               isinstance(region, geometry.Point), \
               (cv, dumps(region))
        new_regions[cv] = region
    log.stop()

    # empty regions
    if not restrict:
        for cv in cvs:
            if cv not in regions:
                new_regions[cv] = EMPTY

    # sort
    regions = collections.OrderedDict(sorted(new_regions.items(),
                                             key=lambda it: it[0],
                                             cmp=CountVector.lex))

    log.stop()

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
        else:
            raise Exception("count vector (%s) has invalid region (%s)" % (cv, dumps(region)))

        coords = dumps(region)
        toks = (cv, coords, area)
        print >>out, '\t'.join(map(str, toks))
    out.close()

#==========================================================
# visualization

def draw_landscape(regions, duprange, lossrange,
                   filename=None,
                   log=False):
    """
    Plot the landscape.

    The x-axis represents dup cost and the y-axis represents the loss cost
    (both relative to unit cost of deep coalescence).
    """

    from rasmus import plotting
    import matplotlib.pyplot as plt

    # axes
    dup_min, dup_max = duprange
    loss_min, loss_max = lossrange
    if log:
        plt.xscale('log')
        plt.yscale('log')
    plt.axis('scaled')
    plt.axis([dup_min, dup_max, loss_min, loss_max])
    plt.xlabel("relative duplication cost")
    plt.ylabel("relative loss cost")

    # legend
    legend_handles = []
    legend_labels = []

    # color map
    n = len(regions)
    colormap = plotting.rainbow_color_map(low=0, high=n-1)

    # plot regions
    for i, (cv, region) in enumerate(regions.iteritems()):
        # output
        label = str((cv.d, cv.l, cv.c)) + ":" + str(cv.count)
        color = colormap.get(i)

        if isinstance(region, geometry.Polygon):        # non-degenerate
            coords = list(region.exterior.coords)
            h = plt.gca().add_patch(plt.Polygon(coords, color=color))
        elif isinstance(region, geometry.LineString):   # degenerate
            coords = list(region.coords)
            h, = plt.plot([coords[0][0], coords[1][0]],
                          [coords[0][1], coords[1][1]],
                          linewidth=4, color=color)
        elif isinstance(region, geometry.Point):        # degenerate
            coords = list(region.coords)
            h, = plt.plot(coords[0][0], coords[0][1],
                          'o', markersize=4, color=color)
        else:                                           # non-degenerate (collection)
            raise Exception("count vector (%s) has invalid region (%s)" % (cv, dumps(region)))

        legend_handles.append(h)
        legend_labels.append(label)

    # legend
    leg = plt.legend(legend_handles, legend_labels, numpoints=1,
                     fontsize=10, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    # save
    if filename:
        plt.savefig(filename, format="pdf",
                    bbox_extra_artists=(leg,), bbox_inches='tight')
    else:
        plt.show()
    plt.close()
