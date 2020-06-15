"""
reconscape.py
Library to solve for DLC MPR Landscapes using DP
"""

# python libraries
import sys
import collections
import itertools

# geometry libraries
from fractions import Fraction
from shapely import geometry
from shapely.ops import cascaded_union
from shapely.wkt import dumps, loads

# rasmus libraries
from rasmus import timer
from rasmus import treelib
from rasmus import util
from compbio import phylo

# dlcpar libraries
from dlcpar import common
from dlcpar import constants
from dlcpar import reconlib
from dlcpar import countvector as cvlib
from dlcpar import recondp

#==========================================================
# globals

INF = recondp.INF
INIT_LOCUS = recondp.INIT_LOCUS

#==========================================================
# reconciliation

def dlcscape_recon(tree, stree, gene2species, gene2locus=None,
                   duprange=constants.DEFAULT_DUP_RANGE,
                   lossrange=constants.DEFAULT_LOSS_RANGE,
                   max_loci=INF, max_dups=INF, max_losses=INF,
                   allow_both=False,
                   compute_events=False,
                   log=sys.stdout):
    """Find reconciliation landscape with parsimony costs through Pareto optimality"""

    reconer = DLCScapeRecon(tree, stree, gene2species, gene2locus,
                            duprange=duprange, lossrange=lossrange,
                            max_loci=max_loci, max_dups=max_dups, max_losses=max_losses,
                            allow_both=allow_both,
                            compute_events=compute_events,
                            log=log)
    return reconer.recon()


class DLCScapeRecon(recondp.DLCRecon):
    """Reconciliation class using Pareto optimality"""

    def __init__(self, gtree, stree, gene2species, gene2locus=None,
                 duprange=constants.DEFAULT_DUP_RANGE,
                 lossrange=constants.DEFAULT_LOSS_RANGE,
                 max_loci=INF, max_dups=INF, max_losses=INF,
                 allow_both=False,
                 compute_events=False,
                 name_internal="n", log=sys.stdout):
        recondp.DLCRecon.__init__(self, gtree, stree, gene2species, gene2locus,
                                  max_loci=max_loci, max_dups=max_dups, max_losses=max_losses,
                                  allow_both=allow_both,
                                  name_internal=name_internal, log=log)

        dup_min, dup_max = duprange
        loss_min, loss_max = lossrange
        assert (dup_min > 0) and (dup_max > 0) and (dup_min < dup_max)
        assert (loss_min > 0) and (loss_max > 0) and (loss_min < loss_max)
        self.duprange = duprange
        self.lossrange = lossrange

        # compute_events is True to return events
        # False returns a CountVector where events is an empty Counter
        self.compute_events = compute_events

        # these attributes are assigned when performing reconciliation using self.recon()
        # all locus maps tiles - used for reconciliation graph
        self.locus_maps = None
        self.count_vectors = None


    #=============================
    # main methods

    def recon(self):
        """Perform reconciliation"""

        # check feasibility
        self.log.start("Checking feasibility")
        feasible = self._are_constraints_consistent()
        self.log.stop()
        if not feasible:
            self.log.log("tree not feasible")
            return self.gtree, None

        self.log.start("Reconciling")

        # log input gene and species trees
        self.log.log("gene tree\n")
        reconlib.log_tree(self.gtree, self.log, func=treelib.draw_tree_names)
        self.log.log("species tree\n")
        reconlib.log_tree(self.stree, self.log, func=treelib.draw_tree_names)

        # infer species map
        self._infer_species_map()
        self.log.log("\n\n")

        # add implied speciation nodes but first start the species tree at the right root
        substree = treelib.subtree(self.stree, self.srecon[self.gtree.root])
        subrecon = util.mapdict(self.srecon, val=lambda snode: substree.nodes[snode.name])

        # switch internal storage with subtrees
        self.stree = substree
        self.srecon = subrecon

        # add implied nodes (standard speciation, speciation from duplication, delay nodes)
        # then relabel events (so that factor_tree works)
        sevents = phylo.label_events(self.gtree, self.srecon)
        reconlib.add_implied_nodes(self.gtree, self.stree, self.srecon, sevents, delay=self.delay)
        common.rename_nodes(self.gtree, self.name_internal)

        # log gene tree (with species map)
        self.log.log("gene tree (with species map)\n")
        reconlib.log_tree(self.gtree, self.log, func=reconlib.draw_tree_recon,
                          srecon=self.srecon)

        # infer locus map
        # sets self.count_vectors and self.locus_maps
        self.locus_maps = self._infer_locus_map()

        # calculate runtime
        runtime = self.log.stop()

        return self.count_vectors, runtime


    #=============================
    # event/cost methods -- these operate at the species branch level

    def _count_events(self, lrecon, subtrees, nodefunc=lambda node: node.name,
                      all_bottoms=None,
                      max_dups=INF, max_losses=INF,
                      snode=None):
        """
        Count number of dup, loss, coal events

        Returns list of tuples, where each element of list corresponds to one reconciliation.
        See also DLCRecon._count_events.
        """

        extra = {"species_map" : self.srecon, "locus_map" : lrecon}

        # defaults
        ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln = INF, INF, INF, INF, {}, INF
        events = collections.Counter()

        # duplications
        dup_nodes = reconlib.find_dup_snode(self.gtree, self.stree, extra, snode=None,
                                            subtrees_snode=subtrees,
                                            nodefunc=nodefunc)
        ndup = len(dup_nodes)
        if ndup > max_dups:     # skip rest if exceed max_dups
            return [(ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events)]

        # losses
        losses = reconlib.find_loss_snode(self.gtree, self.stree, extra, snode=None,
                                          subtrees_snode=subtrees,
                                          nodefunc=nodefunc)
        nloss = len(losses)
        if nloss > max_losses:  # skip rest if exceed max_losses
            return [(ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events)]

        # extra lineages at speciations
        coal_specs = reconlib.find_coal_spec_snode(self.gtree, self.stree, extra, snode=None,
                                                   subtrees_snode=subtrees,
                                                   nodefunc=nodefunc,
                                                   implied=self.implied)
        ncoal_spec = 0
        for lineages in coal_specs:
            ncoal_spec += len(lineages) - 1

        # extra lineages at duplications
        start = self._find_locus_orders_start(lrecon, subtrees, nodefunc=nodefunc,
                                              dup_nodes=dup_nodes, all_bottoms=all_bottoms)
        orders, nsoln = self._find_locus_orders(lrecon, subtrees, start, nodefunc=nodefunc,
                                                dup_nodes=dup_nodes, all_bottoms=all_bottoms)

        # calculate ncoal_dup using a random order (all of these orders are optimal)
        order = {}
        for locus, orderings in orders.iteritems():
            order[locus] = common.random_choice(orderings)
        ncoal_dup = self._count_coal_dup(lrecon, order, start, nodefunc=nodefunc)

        #====================
        # make events

        if not self.compute_events:
            return [(ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events)]

        # speciations
        specs = reconlib.find_spec_snode(self.gtree, self.stree, extra, snode=snode,
                                         subtrees_snode=subtrees,
                                         nodefunc=nodefunc)
        spec_events = self._make_spec_events(specs, snode)

        # losses
        loss_events = self._make_loss_events(losses, snode)

        # extra lineages at speciations
        coal_spec_events = self._make_coal_spec_events(coal_specs, snode)

        # merge events so far
        events = spec_events + loss_events + coal_spec_events

        # duplications and coalescences due to duplication
        # do last because depends on partial order and we want to reuse other events
        solns = []
        order_events = self._make_order_events(lrecon, orders, start, snode,
                                               nodefunc=nodefunc)
        for order, events_for_order in order_events:
            # merge dup and coal_dup events with rest of events
            # each partial order contributes one solution
            merged_events = events.copy() + events_for_order.copy()
            soln = (ndup, nloss, ncoal_spec, ncoal_dup, order, 1, merged_events)
            solns.append(soln)
        return solns


    def _make_spec_events(self, specs, snode):
        """Make spec events"""

        srecon = self.srecon

        events = collections.Counter()
        if snode.is_leaf():
            assert len(specs) == 0, "speciations for leaf species not allowed"
            return events

        assert len(snode.children) == 2, "species tree not binary"
        left_snode, right_snode = snode.children

        for nodes in specs:
            lefts = set()
            rights = set()
            for node in nodes:
                for child in node.children:
                    if srecon[child] == left_snode:
                        assert child not in lefts
                        lefts.add(child)
                    elif srecon[child] == right_snode:
                        assert child not in rights
                        rights.add(child)
                    else:
                        raise Exception("invalid child species")
            spec = ["S", tuple(sorted(nodes)), tuple(sorted(lefts)), tuple(sorted(rights)), snode]
            events[tuple(spec)] = 1
        return events


    def _make_loss_events(self, losses, snode):
        """Make loss events"""

        srecon = self.srecon

        events = collections.Counter()
        for nodes in losses:
            survived = set()
            for node in nodes:
                for child in node.children:
                    if srecon[child] != snode:
                        assert child not in survived
                        survived.add(child)
            loss = ["L", tuple(sorted(nodes)), tuple(sorted(survived)), snode]
            events[tuple(loss)] = 1
        return events


    def _make_coal_spec_events(self, coal_specs, snode):
        """Make coal_spec events"""

        events = collections.Counter()
        for lineages in coal_specs:
            assert len(lineages) == len(set(lineages))
            coal_spec = ["C", tuple(sorted(lineages)), snode]
            events[tuple(coal_spec)] = 1
        return events


    def _make_order_events(self, lrecon, orders, start, snode,
                           nodefunc=lambda node: node.name):
        """Make dup events and coal_dup events based on set of dup nodes and all orders

        Return list of (dup_order, events), where events is a Counter for the duplication order.
        """

        # skip if no order - order is empty, and events counter is empty
        if not orders: # len(orders) == 0
            return [({}, collections.Counter())]

        # generate list of dictionaries for all possible combinations of optimal locus orders
        all_opt_orders = []
        ranges = [range(len(orders[locus])) for locus in orders] # length of node list for each locus
        for choice in itertools.product(*ranges): # Cartesian product
            order = {}
            for index, locus in enumerate(orders):
                order[locus] = orders[locus][choice[index]] # new dict with choices made for ordering
            all_opt_orders.append(order)

        # make events
        events = []    # list of (dup_order, events_for_this_order)
        for order in all_opt_orders:
            events_for_order = collections.Counter()

            # get contemporary lineages at duplications
            coals = self._find_contemporary_lineages(lrecon, order, start, nodefunc=nodefunc)

            # dup events
            for dup_node, lineages in coals.iteritems():
                assert len(lineages) == len(set(lineages))

                # left: lineage in new locus
                lefts = (dup_node,)

                # right: contemporary lineages in parent locus
                rights = set(lineages)

                # put left and right together
                dup = ["D", dup_node, tuple(sorted(lefts)), tuple(sorted(rights)), snode]
                events_for_order[tuple(dup)] = 1

            # coalescence due to duplication events
            for lineages in coals.itervalues():
                if len(lineages) > 1:
                    coal_dup = ["K", tuple(sorted(lineages)), snode]
                    events_for_order[tuple(coal_dup)] = 1

            # put together
            events.append((order, events_for_order))

        # return list of possible event sets
        return events


    #=============================
    # locus tiles methods (used by _enumerate_locus_maps)
    # tiles structure:
    #     key1 = snode, key2 = (bottom_loci, top_loci)
    #     value = CountVectorSet

    def _initialize_tiles_sroot(self, tiles):
        """Initialize tiles for species root"""
        top_loci = bottom_loci = (INIT_LOCUS,)
        cvs = cvlib.CountVectorSet([cvlib.ZERO_COUNT_VECTOR])
        tiles[bottom_loci, top_loci] = cvs


    def _find_optimal_solns(self, tiles, bottom_loci, top_loci,
                            lrecon, subtrees, bottoms=None,
                            max_dups=INF, max_losses=INF,
                            snode=None):
        """Find solutions for tiles"""
        # each soln has format (ndup, nloss, ncoal_spec, ncoal_dup, order, nsoln, events)
        solns = self._count_events(lrecon, subtrees, all_bottoms=bottoms,
                                   max_dups=max_dups, max_losses=max_losses,
                                   snode=snode)
        return solns


    def _update_tiles(self, tiles, bottom_loci, top_loci,
                      lrecon, order,
                      ndup, nloss, ncoalspec, ncoaldup, nsoln, events):
        """Update tiles to keep those with Pareto-optimal count vectors"""

        key = (bottom_loci, top_loci)
        if key not in tiles:
            tiles[key] = cvlib.CountVectorSet()

        # create new CountVector
        ncoal = ncoalspec + ncoaldup    # can ignore cost of coal due to dup if desired
        cv = cvlib.CountVector(ndup, nloss, ncoal, nsoln, events)
        tiles[key].add(cv)

        # filter cvs to keep it pareto-optimal: updates min_cvs so that
        # count_events uses best cvs possible when deciding when to skip work
        tiles[key] = tiles[key].pareto_filter(self.duprange, self.lossrange)


    def _choose_tile(self, tiles):
        """Choose single tile with Pareto-optimal count vector"""

        self.log.log("optimal count vectors")

        for (bottom_loci, top_loci), cvs in tiles.iteritems():
            new_cvs = cvs.pareto_filter(self.duprange, self.lossrange)
            tiles[bottom_loci, top_loci] = new_cvs

            # log
            self.log.log("\t%s -> %s" % (top_loci, bottom_loci))
            for cv in new_cvs:
                self.log.log("\t\tvector: %s" % cv)
            self.log.log()


    #=============================
    # DP table methods (used by _infer_opt_locus_map)

    def _dp_compute(self, locus_maps):
        """Compute DP table"""
        # locus_maps is a multi-dimensional dict with the structure
        # key1 = snode, key2 = (bottom_loci, top_loci), value = CountVectorSet

        stree = self.stree

        # dynamic programming storage
        dp_table = {}    # key1 = snode, key2 = top_loci, val = cvs-to-go

        # recur up species tree to determine (optimal) cvs-to-go
        for snode in stree.postorder():
            self.log.start("Working on snode %s" % snode.name)
            dp_table[snode] = {}

            if snode not in locus_maps:
                # nothing in this sbranch
                dp_table[snode][()] = cvlib.CountVectorSet([cvlib.ZERO_COUNT_VECTOR])
                self.log.log("Empty sbranch")
                self.log.stop()
                continue

            # get stored values for the sbranch
            locus_maps_snode = locus_maps[snode]

            if snode.is_leaf():
                # leaf base case
                for (bottom_loci, top_loci), cvs in locus_maps_snode.iteritems():
                    assert top_loci not in dp_table[snode]
                    dp_table[snode][top_loci] = cvs
            else:
                if len(snode.children) != 2:
                    raise Exception("non-binary species tree")

                # find cost-to-go of assigning top_loci to top of sbranch
                # = cost of top_loci to bottom_loci along sbranch
                #   + cost of bottom_loci at top of left child
                #   + cost of bottom_loci at top of right child
                sleft, sright = snode.children
                costs = collections.defaultdict(cvlib.CountVectorSet) # key = top_loci, val = cvs

                for (bottom_loci, top_loci), cvs in locus_maps_snode.iteritems():
                    # find cost-to-go in children
                    # locus assignment may have been removed due to search heuristics
                    cvs_left = dp_table[sleft].get(bottom_loci,
                                                   cvlib.CountVectorSet([cvlib.MAX_COUNT_VECTOR]))
                    cvs_right = dp_table[sright].get(bottom_loci,
                                                     cvlib.CountVectorSet([cvlib.MAX_COUNT_VECTOR]))
                    children_cvs = cvs_left * cvs_right

                    # add cost in this sbranch
                    cvs_to_go = cvs * children_cvs
                    costs[top_loci] = costs[top_loci].merge(cvs_to_go)

                # for each assignment of top_loci to top of sbranch,
                # filter down to set of Pareto-optimal count vectors
                for top_loci, cvs in costs.iteritems():
                    dp_table[snode][top_loci] = cvs.pareto_filter(self.duprange, self.lossrange)

            self.log.log("DP table")
            for top_loci, cvs in dp_table[snode].iteritems():
                self.log.log("\t%s :" % str(top_loci))
                for cv in cvs:
                    self.log.log("\t\tvector: %s" % cv)
            self.log.stop()

        # TODO: dp_table has some Pareto-optimal vectors
        # such that for all possible dup/coal cost and loss/coal cost,
        # the vector never has minimum cost
        return dp_table


    def _dp_terminate(self, dp_table):
        """Terminate DP table to find optimal cost vector"""
        stree = self.stree

        # not necessary since cost along root sbranch already determined,
        # and by design, dp_table[sroot] is always assigned locus = INIT_LOCUS
        assert len(dp_table[stree.root]) == 1, dp_table[stree.root]
        cvs = dp_table[stree.root].values()[0]

        self.log.log("")
        self.log.log("Optimal count vectors:")
        for cv in cvs:
            self.log.log("\tvector: %s" % cv)
        self.log.log("")

        self.count_vectors = cvs


    def _dp_traceback(self, locus_maps, subtrees, dp_table):
        """Traceback through DP table (null)"""
        pass


#==========================================================
# regions

def _get_frac(numerator=0, denominator=None):
    """Fraction representation"""
    return Fraction(numerator, denominator).limit_denominator()


def get_regions(cvs, duprange, lossrange, restrict=True,
                log=sys.stdout):
    """
    Return dictionary mapping count vectors to regions.

    cvs       -- Pareto-optimal CountVectorSet
    duprange  -- (low, high) relative to unit cost of (deep) coalescence
    lossrange -- (low, high) relative to unit cost of (deep) coalescence
    restrict  -- Keep only regions in the specified range
    """

    log = timer.Timer(log)
    log.start("Finding regions")

    dup_min, dup_max = map(_get_frac, duprange)
    loss_min, loss_max = map(_get_frac, lossrange)

    # bounding box
    bb = geometry.box(dup_min, loss_min, dup_max, loss_max)

    # empty
    empty = geometry.GeometryCollection()

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
                xint = _get_frac(cv2.c - cv1.c, cv1.d - cv2.d)  # x-intercept
                lines_vert.add(xint)
            else:
                m = _get_frac(cv2.d - cv1.d, cv1.l - cv2.l)  # slope
                b = _get_frac(cv2.c - cv1.c, cv1.l - cv2.l)  # y-intercept
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
        x = _get_frac(b2 - b1, m1 - m2)
        y1 = m1 * x + b1
        y2 = m2 * x + b2
        assert y1 == y2, (y1, y2)
        points.add((x, y1))
    for (m, b), x in itertools.product(lines_diag, lines_horiz):
        y = m * x + b
        points.add((x, y))
    for (m, b), y in itertools.product(lines_diag, lines_vert):
        x = _get_frac(y - b, m)
        points.add((x, y))
    for x, y in itertools.product(lines_vert, lines_horiz):
        points.add((x, y))
    points = [pt for pt in points if bb.intersects(geometry.Point(pt))]


    def find_closest_point(point):
        """find closest point"""
        geom_point = geometry.Point(point)
        min_points, min_dist = \
            util.minall(points, minfunc=lambda pt: geom_point.distance(geometry.Point(pt)))
        if min_dist > 1e-10:
            return None     # no point found (probably a point co-linear with neighbors)
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
                xint = _get_frac(cv2.c - cv1.c, cv1.d - cv2.d)  # x-intercept

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
                m = _get_frac(cv2.d - cv1.d, cv1.l - cv2.l)  # slope
                b = _get_frac(cv2.c - cv1.c, cv1.l - cv2.l)  # y-intercept

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
                        filtered_regions = [r for r in region if isinstance(r, geometry.Polygon)]
                        assert len(filtered_regions) == 1
                        region = filtered_regions[0]
                    except:
                        pass

                    # finesse coordinates (due to floating-point approximations)
                    if isinstance(region, geometry.Polygon):
                        coords = list(region.exterior.coords)[:-1]

                        # find closest coordinates
                        coords = map(find_closest_point, coords)
                        coords = [pt for pt in coords if pt is not None]

                        # collapse coordinates
                        new_coords = [coords[0]]
                        for i in xrange(1, len(coords)):
                            p1, p2 = coords[i-1], coords[i]
                            if not geometry.Point(p1).almost_equals(geometry.Point(p2)):
                                new_coords.append(p2)
                        coords = new_coords

                        # make new region
                        if len(coords) > 2:
                            keep = True
                            region = geometry.Polygon(coords)
            if not keep:
                region = empty
                break

        # keep only polygons
        if (not region.is_empty) and isinstance(region, geometry.Polygon):
            regions[cv1] = region
    log.stop()

    # go back and get lines and points
    # Note that new_regions could have value = set, but shapely no longer supports Polygon hashing
    log.start("Mapping lines and points")
    new_regions = collections.defaultdict(list)
    for cv, region in regions.iteritems():
        if region not in new_regions[cv]:
            new_regions[cv].append(region)
        coords = map(find_closest_point, region.exterior.coords)
        coords = [pt for pt in coords if pt is not None]
        for i in range(1, len(coords)): # lines
            x1, y1 = coords[i-1]
            min_cvs1, _ = util.minall(cvs, minfunc=lambda cv: cv.d * x1 + cv.l * y1 + cv.c)
            x2, y2 = coords[i]
            min_cvs2, _ = util.minall(cvs, minfunc=lambda cv: cv.d * x2 + cv.l * y2 + cv.c)
            min_cvs = set(min_cvs1) & set(min_cvs2)
            poly = geometry.LineString(coords[i-1:i+1])
            for min_cv in min_cvs:
                if poly not in new_regions[min_cv]:
                    new_regions[min_cv].append(poly)
        for x, y in coords:             # points
            min_cvs, _ = util.minall(cvs, minfunc=lambda cv: cv.d * x + cv.l * y + cv.c)
            poly = geometry.Point(x, y)
            for min_cv in min_cvs:
                if poly not in new_regions[min_cv]:
                    new_regions[min_cv].append(poly)
    for cv, geoms in new_regions.iteritems():
        region = cascaded_union(geoms)
        assert not region.is_empty, cv
        assert isinstance(region, (geometry.Polygon, geometry.LineString, geometry.Point)), \
               (cv, dumps(region))
        new_regions[cv] = region
    log.stop()

    # empty regions
    if not restrict:
        for cv in cvs:
            if cv not in regions:
                new_regions[cv] = empty

    # sort
    regions = collections.OrderedDict(sorted(new_regions.items(),
                                             key=lambda it: it[0],
                                             cmp=cvlib.CountVector.lex))

    log.stop()

    return regions


def read_regions(filename):
    """Read regions from file"""
    regions = collections.OrderedDict()

    for i, toks in enumerate(util.DelimReader(filename)):
        if i == 0:
            duprange, lossrange = map(float, toks[:2]), map(float, toks[2:])
            continue

        cv_str, coords_str, _ = toks # (cv_str, coords_str, area_str)
        cv = cvlib.parse_count_vector(cv_str)
        region = loads(coords_str)
        #area = float(area_str)
        regions[cv] = region

    return regions, duprange, lossrange


def write_regions(filename, regions, duprange, lossrange, close=False):
    """Write regions to file in tab-delimited format

    duprange_low    duprange_high    lossrange_low    lossrange_high
    count_vector_string    polygon_region    area
    ...
    """
    out = util.open_stream(filename, 'w')
    print >>out, '\t'.join(map(str, duprange + lossrange))

    for cv, region in regions.iteritems():
        coords = None
        area = None
        if isinstance(region, geometry.Polygon):                          # non-degenerate
            coords = list(region.exterior.coords)
            area = region.area
        elif isinstance(region, (geometry.LineString, geometry.Point)):   # degenerate
            coords = list(region.coords)
            area = region.area
        else:
            raise Exception("count vector (%s) has invalid region (%s)" % (cv, dumps(region)))

        coords = dumps(region)
        toks = (cv.to_string(), coords, area) # do not write out events
        print >>out, '\t'.join(map(str, toks))

    if close:
        out.close()


#==========================================================
# events

def write_events(filename, regions, intersect, close=False):
    """Write events to file in tab-delimited format

    (header row) Count Vector    Events
    (data row)   count_vector_string    event string (one per column)
    """

    # create CountVectorSet from regions
    cvs = cvlib.CountVectorSet()
    for cv in regions.iterkeys():
        cvs.add(cv)

    # event: key1 = (d,l,c)
    #        key2 = (unformatted) event
    #        val = count of event across all MPRs with that cost vector.
    if intersect:
        events = cvs.intersect_events()
    else:
        events = cvs.union_events()

    # format events, and calculate number of regions each event appears in
    event_dict = {}
    event_region_counts = collections.Counter() # key = formatted event (without counts)
                                                # value = number of regions that event appears in
    for cv, d in events.iteritems():
        formatted_events = [format_event(event, d[event]) for event in d]
        event_dict[cv] = formatted_events

        formatted_events = [format_event(event) for event in d] # ignores counts
        cv_formatted_events = collections.Counter(formatted_events)
        event_region_counts.update(cv_formatted_events)

    # open output file
    out = util.open_stream(filename, 'w')
    print >>out, '\t'.join(["Count Vector", "Events"])

    # write each vector with its associated events (union or intersection)
    for cv in cvs:
        cv_key = cv.to_tuple(count=False)
        print >>out, '\t'.join([cv.to_string()] + sorted(event_dict[cv_key]))
    print >>out, ''

    # write number of regions with events that occur in that number of regions
    print >>out, '\t'.join(["# Regions", "Events"])
    num_regions = sorted(set(event_region_counts.values()), reverse=True)
    for nr in num_regions:
        # get events
        events = []
        for (event, count) in event_region_counts.iteritems():
            if count == nr:
                events.append(event)
        assert events # len(events) > 0

        # write events
        print >>out, '\t'.join([str(nr)] + sorted(events))

    if close:
        out.close()


def format_event(event, count=None, sep=' '):
    """
    Return string representation of event using only leaves of gene tree

    Converts from internal event representation, which uses internal gene tree nodes
    (including implied speciation nodes).  Outputs events in similar format as tree-relations.

    S,D,L uses subtrees of the locus tree. C,K uses subtrees of the gene tree.
    All leaves are sorted within their own list.
    speciation:  S    leaves on one subtree
                      leaves on other subtree
                      species above the speciation
    duplication: D    leaves on subtree with daughter locus
                      leaves on subtree with parent locus
                      species where the duplication occurred
    loss:        L    leaves on subtree that is not lost
                      species where the locus was lost
    coal_spec:   C    tuples of leaves for each subtree which did not coalesce
                      species in which lineages could have coalesced
    coal_dup:    K    tuples of leaves for each subtree which did not coalesce at a duplication
                      species with the duplication
    """

    event_type = event[0]
    snode = event[-1]
    info = event[1:-1]

    # new event is same kind of event (S, D, L, C, K) as old one
    new_event = [event_type]

    # speciation
    if event_type == 'S':
        _, lefts, rights = info # (nodes, lefts, rights)

        # get leaves
        left_leaves = []
        right_leaves = []
        for left in lefts:
            left_leaves.extend(left.leaf_names())
        for right in rights:
            right_leaves.extend(right.leaf_names())

        # sort leaves, then flip left and right if needed
        left_leaves.sort()
        right_leaves.sort()
        if right_leaves < left_leaves:
            left_leaves, right_leaves = right_leaves, left_leaves

        # create string
        new_event.append('(' + sep.join(left_leaves) + ')' \
                         + sep \
                         + '(' + sep.join(right_leaves) + ')')

    # duplication
    elif event_type == 'D':
        _, lefts, rights = info # (dup_node, lefts, rights)

        # get leaves
        assert len(lefts) == 1, left # only one (starting lineage) after dup node
        left_leaves = lefts[0].leaf_names()
        right_leaves = reduce(lambda x, y: x+y, [node.leaf_names() for node in rights])

        # sort leaves
        left_leaves.sort()
        right_leaves.sort()

        # create string
        new_event.append('(' + sep.join(left_leaves) + ')' \
                         + sep \
                         + '(' + sep.join(right_leaves) + ')')

    # loss
    elif event_type == 'L':
        _, survived = info # (nodes, survived)

        # get leaves
        leaves = []
        for other in survived:
            leaves.extend(other.leaf_names())

        # sort leaves
        leaves.sort()

        # create string
        new_event.append('(' + sep.join(leaves) + ')')

    # coalescence - the gene leaves which are children of lineages on the same locus
    elif event_type == 'C' or event_type == 'K':
        nodes = info[0]

        # get leaves for each subtree
        leaf_strs = []
        for node in nodes:
            leaf_strs.append('(' + sep.join(sorted(node.leaf_names())) + ')')

        # create string
        new_event.append('(' + sep.join(sorted(leaf_strs)) + ')')

    else:
        raise Exception("invalid event type: %s" % str(event))

    # species node of new event is same as old event
    new_event.append(str(snode.name))
    event_str = sep.join(new_event)
    if count:
        event_str += ":" + str(count)

    return event_str


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

    from rasmus import colors
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
    colormap = colors.rainbow_color_map(low=0, high=n-1)

    # plot regions
    for i, (cv, region) in enumerate(regions.iteritems()):
        # output
        label = cv.to_string()
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
        plt.savefig(filename,
                    bbox_extra_artists=(leg,), bbox_inches='tight')
    else:
        plt.show()
    plt.close()
