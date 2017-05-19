"""

   Event Count Vectors

"""

# python libraries
import re

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

    def __init__(self, d, l, c, count=1):
        self.d = d
        self.l = l
        self.c = c
        self.count = count

    def __add__(self, other):
        return CountVector(self.d + other.d,
                           self.l + other.l,
                           self.c + other.c,
                           self.count * other.count)

    def __repr__(self):
        return "<%s,%s,%s>:%s" % (self.d, self.l, self.c, self.count)

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

    def update(self, other):
        for v in other:
            self.add(v)

    def __mul__(self, other):
        """
        Returns new CountVectorSet computed from (1) taking Cartesian product of self and other, then
        (2) converting each result into single cost vector by adding two cost vectors.
        """
        result = CountVectorSet()
        for v in self:
            for w in other:
                result.add(v + w)
        return result

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

def is_minimal(v, cvs):
    """Returns True if CountVector v is smaller than all (non-equal) cost vectors in CountVectorSet cvs"""
    for w in cvs:
        if w < v:
            return False
    return True

def is_maximal(v, cvs):
    """Returns True is CountVector v is larger than all (non-equal) cost vectors in CountVectorSet cvs"""
    for w in cvs:
        if w > v:
            return False
    return True

#==========================================================
# globals

INF = util.INF

ZERO_COUNT_VECTOR = CountVector(0, 0, 0, 1)
MAX_COUNT_VECTOR = CountVector(INF, INF, INF, INF)
