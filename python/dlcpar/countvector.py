"""

   Event Count Vectors

"""

# python libraries
import re

# rasmus libraries
from rasmus import util
from collections import defaultdict
from collections import Counter

#==========================================================
# event count vectors

class CountVector(object):
    """
    A count vector is in the following format:
    d(uplication), l(oss), c(oalescence)

    Additionally, this class keeps track of the count,
    (number of reconciliations with the count vector).
    """

    def __init__(self, d, l, c, count=1, events=[Counter()]):
        self.d = d
        self.l = l
        self.c = c
        self.count = count
        # events is a list of possible sets of events for a given tile
        self.events = [events]
        

    def __add__(self, other):
        # flatten the event lists
        s = self.flatten()
        o = other.flatten()
        d = Counter()
        for sevent, scount in s.events[0].iteritems():
            d[sevent] = scount * other.count
        for oevent, ocount in o.events[0].iteritems():
            d[oevent] = ocount * self.count

        return CountVector(self.d + other.d,
                           self.l + other.l,
                           self.c + other.c,
                           self.count * other.count,
                           d
                           )

    def __repr__(self):
        return "<%s,%s,%s>:%s:%s" % (self.d, self.l, self.c, self.count, self.events)

    def __eq__(self, other):
        return (self.d == other.d) and (self.l == other.l) and (self.c == other.c)

    def __lt__(self, other):
        return (self.d <= other.d) and \
               (self.l <= other.l) and \
               (self.c <= other.c) and \
               ((self.d < other.d) or (self.l < other.l) or (self.c < other.c))

    def __lte__(self, other):
        return self.__lt__(other) or self.__eq__(other)

    def __mul__(self, cvs):
        result = CountVectorSet()
        for v in cvs:
            result.add(self + v)
        return result

    def lex(self, other):
        # lexicographic sorting
        if self.__eq__(other): return 0
        if self.d < other.d: return -1
        elif self.d == other.d and self.l < other.l: return -1
        elif self.d == other.d and self.l == other.l and self.c < other.c: return -1
        else: return 1

    def to_tuple(self, count=False):
        if count:
            return (self.d, self.l, self.c, self.count)
        else:
            return (self.d, self.l, self.c)

    def flatten(self):
        sevents = reduce(lambda x,y: x + y, self.events)
        return CountVector(self.d, self.l, self.c, self.count, sevents)

def parse_count_vector(string):
    pattern = "^<(\d+),(\d+),(\d+)>:(\d+)$"
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
        k = v.to_tuple()
        if k not in self.dict:
            self.dict[k] = v
        else:
            self.dict[k].count += v.count
            # add means there's multiple possible events
            self.dict[k].events.extend(v.events)

    def update(self, other):
        # flatten each vector in self first
        fself = self.flatten_set()
        # now add each flattened element of other to self
        for v in other:
            fself.add(v.flatten())
        return fself

    def __mul__(self, other):
        """
        Returns new CountVectorSet computed from (1) taking Cartesian product of self and other, then
        (2) converting each result into single cost vector by adding two cost vectors.
        """
        # create cvs with flattened events first
        fself = self.flatten_set()
        fother = other.flatten_set()
        # now make the cartesian product
        result = CountVectorSet()
        for v in fself:
            for w in fother:
                result.add(v + w)
        return result#.flatten_set()

    def _filter(self, duprange, lossrange):
        """Returns new CountVectorSet in which cost vectors that cannot be optimal in given cost range are removed"""

        dup_min, dup_max = duprange
        loss_min, loss_max = lossrange

        # compute lowest upper bound
        LUB = min(map(lambda v: v.d * dup_max + v.l * loss_max + v.c, self))

        # filter list
        result = CountVectorSet()
        for v in self:
            mincost = v.d * dup_min + v.l * loss_min + v.c
            if mincost <= LUB:
                result.add(v)

        return result

    def pareto_filter(self, duprange, lossrange):
        """Returns new CountVectorSet consisting only of Pareto optimal cost vectors"""

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
        """Returns a new CountVectorSet where the event lists have been flattened to contain only one counter"""
        result = CountVectorSet()
        for v in self:
            result.add(v.flatten())
        return result

    def union_events(self):
        """Returns a dictionary of events, where the key is a CV and the value is a list of all events
        appearing in any MPR with that cost"""
        result = defaultdict(list)
        # flatten self for no duplicates
        for v in self:
            fv = v.flatten()
            for event_dict in fv.events:
                for event, count in event_dict.iteritems():
                    result[v].append(event)
        return result


    def intersect_events(self):
        """Returns a dictionary of events, where the key is a CV and the value is a list of all events
        appearing in ALL MPRs with that cost"""
        result = defaultdict(list)
        # flatten self to know the true event counts
        for v in self:
            fv = v.flatten()
            for event_dict in fv.events:
                for event, count in event_dict.iteritems():
                    if count == v.count:
                        result[v].append(event)
        return result

def is_minimal(v, cvs):
    """Returns True if CountVector v is smaller than all (non-equal) cost vectors in CountVectorSet cvs"""
    for w in cvs:
        if w < v:
            return False
    return True

def is_maximal(v, cvs):
    """Returns True if CountVector v is at least as large as all (non-equal) cost vectors in CountVectorSet cvs"""
    for w in cvs:
        if w > v:
            return False
    return True

def is_maximal_lte(v, cvs):
    """Returns True if CountVector v is larger than all (non-equal) cost vectors in CountVectorSet cvs"""
    for w in cvs:
        # not using v <= w because it's broken somehow
        if v < w or v == w:
            return False
    return True

#==========================================================
# globals

INF = util.INF

ZERO_COUNT_VECTOR = CountVector(0, 0, 0, 1)
MAX_COUNT_VECTOR = CountVector(INF, INF, INF, INF)
