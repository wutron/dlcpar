"""
countvector.py
Event count vectors
"""

# python libraries
import re
from collections import defaultdict
from collections import Counter

# rasmus libraries
from rasmus import util

#==========================================================
# event count vectors

class CountVector(object):
    """
    A count vector is in the following format:
    d(uplication), l(oss), c(oalescence)

    Additionally, this class keeps track of the count,
    (number of reconciliations with the count vector).
    """

    def __init__(self, d, l, c, count=1, events=None):
        self.d = d
        self.l = l
        self.c = c
        self.count = count

        # self.events is a list of Counters with key = event, value = count
        # where each Counter counts the events for one or more partial order(s).
        if events is None:
            events = Counter()
        self.events = [events]

    def __add__(self, other):
        """Add two CountVectors.

        Each CountVector tracks events for a sub-MPR for part of the species tree.
        (1) The # of solutions multiplies because you may use any MPR for either subtree
            to create a valid MPR for both parts of the species tree.
        (2) The frequency of each event increases by the # of solutions that include that event,
            which is the # of possible sub-MPRs for the other part of the species tree.
        """
        # flatten the event lists to guarantee attribute events is a single-element list
        s = self.flatten()
        o = other.flatten()

        # create the events
        d = Counter()
        for sevent, scount in s.events[0].iteritems():
            d[sevent] = scount * other.count
        for oevent, ocount in o.events[0].iteritems():
            # oevent cannot already be in d if we add CountVectors for different species
            assert oevent not in d
            d[oevent] = ocount * self.count

        return CountVector(self.d + other.d,
                           self.l + other.l,
                           self.c + other.c,
                           self.count * other.count,
                           d)

    def __repr__(self):
        """String representation with events"""
        return "<%s,%s,%s>:%s:%s" % (self.d, self.l, self.c, self.count, self.events)


    def to_string(self):
        """String representation without events"""
        return "<%s,%s,%s>:%s" % (self.d, self.l, self.c, self.count)


    def __eq__(self, other):
        return (self.d == other.d) and (self.l == other.l) and (self.c == other.c)


    def __lt__(self, other):
        return ((self.d <= other.d)
                and (self.l <= other.l)
                and (self.c <= other.c)
                and ((self.d < other.d)
                     or (self.l < other.l)
                     or (self.c < other.c)
                    )
               )


    def __lte__(self, other):
        return self.__lt__(other) or self.__eq__(other)


    def lex(self, other):
        """lexicographic sorting"""
        if self.__eq__(other):
            return 0

        if self.d < other.d:
            return -1
        if self.d == other.d and self.l < other.l:
            return -1
        if self.d == other.d and self.l == other.l and self.c < other.c:
            return -1

        return 1


    def to_tuple(self, count=False):
        """tuple representation"""
        if count:
            return (self.d, self.l, self.c, self.count)
        return (self.d, self.l, self.c)


    def flatten(self):
        """Return a new CountVector in which events have been added
        over all partial orders to yield a frequency.
        """
        # unions the keys, and sums the counts for matching keys
        sum_events = reduce(lambda x, y: x + y, self.events)
        return CountVector(self.d, self.l, self.c, self.count, sum_events)


def parse_count_vector(string):
    """Parse count vector string"""
    # TODO: currently ignores events
    pattern = r"^<(\d+),(\d+),(\d+)>:(\d+)$"
    m = re.match(pattern, string)
    if not m:
        raise Exception("invalid CountVector string: %s" % string)
    vals = map(int, m.groups())
    result = CountVector(*vals)
    return result


class CountVectorSet(object):
    """
    Dictionary of CountVectors, where key = (d,l,c) and value = CountVector.
    """

    def __init__(self, iterable=None):
        self.dict = dict()
        if iterable is not None:
            for v in iterable:
                k = v.to_tuple()
                if k in self.dict:
                    raise Exception("duplicate tuple '%s' in CountVectorSet" % str(k))
                self.dict[k] = v


    def __iter__(self):
        return self.dict.itervalues()


    def add(self, v):
        """Add CV v to this CVS."""
        k = v.to_tuple()
        if k not in self.dict:
            self.dict[k] = v
        else:
            self.dict[k].count += v.count
            # If the two vectors match, they are both a set of solutions for the same tile,
            # so merge their events into a single set of solutions.
            # There could be duplicate events in events after this call.
            self.dict[k].events.extend(v.events)


    def merge(self, other):
        """Return a new CountVectorSet that is the merge of this CVS
        with CVS other."""
        # flatten each vector in self first
        fself = self.flatten_set()

        # now add each flattened element of other to self
        for v in other:
            fself.add(v.flatten())
        return fself


    def __mul__(self, other):
        """Return a new CountVectorSet computed from
        (1) taking Cartesian product of self and other, then
        (2) converting each result into single cost vector
            by adding two cost vectors.
        """

        # create cvs with flattened events first
        fself = self.flatten_set()
        fother = other.flatten_set()

        # now make the cartesian product
        result = CountVectorSet()
        for v in fself:
            for w in fother:
                result.add(v + w)
        return result


    def _filter(self, duprange, lossrange):
        """Return a new CountVectorSet in which cost vectors
        that cannot be optimal in given cost range are removed."""

        dup_min, dup_max = duprange
        loss_min, loss_max = lossrange

        # compute lowest upper bound
        lub = min([v.d * dup_max + v.l * loss_max + v.c for v in self])

        # filter list
        result = CountVectorSet()
        for v in self:
            mincost = v.d * dup_min + v.l * loss_min + v.c
            if mincost <= lub:
                result.add(v)

        return result


    def pareto_filter(self, duprange, lossrange):
        """Return a new CountVectorSet consisting only of
        Pareto-optimal cost vectors."""

        lst = self._filter(duprange, lossrange).dict.values()
        #lst = self.values()
        lst.sort(cmp=CountVector.lex)

        lex_lst = [lst[0]]
        for i in xrange(1, len(lst)):
            predecessor, current = lst[i-1], lst[i]
            if predecessor.d < current.d or predecessor.l < current.l:
                lex_lst.append(current)

        result = CountVectorSet()
        for v in lex_lst:
            if is_minimal(v, lst):
                result.add(v)

        return result


    def flatten_set(self):
        """Return a new CountVectorSet where the event lists
        have been flattened to contain only one counter.
        """
        result = CountVectorSet()
        for v in self:
            result.add(v.flatten())
        return result


    def __combine_events(self, intersect=True):
        """Return a dictionary of events, where key1 = (d,l,c), key2 = event,
        val = count of event across all MPRs with that cost vector.

        intersect indicates event appears in all MPRs, and
        union indicates events appears in at least one MPR.
        """
        result = defaultdict(dict)
        for v in self:
            fv = v.flatten() # flatten count vector to combine event counts
            key = v.to_tuple()
            for event_dict in fv.events:
                for event, count in event_dict.iteritems():
                    if (not intersect) or (intersect and count == v.count):
                        result[key][event] = count
        return result


    def union_events(self):
        """See __combine_events."""
        return self.__combine_events(intersect=False)


    def intersect_events(self):
        """See __combine_events."""
        return self.__combine_events(intersect=True)


def is_minimal(v, cvs):
    """Return True if CountVector v is smaller than all (non-equal) cost vectors
    in CountVectorSet cvs."""
    for w in cvs:
        if w < v:
            return False
    return True


#==========================================================
# globals

INF = util.INF

ZERO_COUNT_VECTOR = CountVector(0, 0, 0, 1)
MAX_COUNT_VECTOR = CountVector(INF, INF, INF, INF)
