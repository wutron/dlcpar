from itertools import *
import math,random

"""
See also the Combinatorial Functions available in SAGE (software)
http://www.sagemath.org/doc/reference/sage/combinat/combinat.html
"""


def factorial(n):
    """Iterative factorial function"""
    f = 1
    while n >= 2:
        f,n = f*n, n-1
    return f


"""
The next four function are in itertools starting from Python 2.6.
"""
           
# Actually a class method of chain
def from_iterable(iterables):
    """
    Generate elements from the first iterable until it is exhausted,
    then proceeds to the next iterable, until all of the iterables are exhausted.
    Used for treating consecutive sequences as a single sequence.
    
    chain.from_iterable(['ABC', 'DEF']) --> A B C D E F
    """
    for it in iterables:
        for element in it:
            yield element


def product_old(*args, **kwds):
    """
    Generate Cartesian product of input iterables.

    product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111

    Note: This function is deprecated.  Use product instead.
    """
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)


def product(*args, **kwds):
    """
    Generate Cartesian product of input iterables.

    product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111

    Alternative fast implementation of product for python < 2.6
    http://stackoverflow.com/questions/1681269/how-code-a-function-similar-to-itertools-product-in-python-2-5
    """
    def cycle(sequence, uplevel):
        while True:
            vals = uplevel.next()  # advance upper level, raises if done
            it = iter(sequence)    # (re-)start iteration of current level
            try:
                while True: yield vals + (it.next(),)
            except StopIteration:
                pass

    step = iter(((),))
    pools = map(tuple, args) * kwds.get('repeat', 1)
    for pool in pools:
        step = cycle(pool, step)   # build stack of iterators
    return step


def permutations(iterable, r=None):
    """
    Generate r-length tuples, all possible orderings, no repeated elements.

    permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
    permutations(range(3)) --> 012 021 102 120 201 210
    """
    pool = tuple(iterable)
    n = len(pool)
    r = n if r is None else r
    if r > n:
        return
    indices = range(n)
    cycles = range(n, n-r, -1)
    yield tuple(pool[i] for i in indices[:r])
    while n:
        for i in reversed(range(r)):
            cycles[i] -= 1
            if cycles[i] == 0:
                indices[i:] = indices[i+1:] + indices[i:i+1]
                cycles[i] = n - i
            else:
                j = cycles[i]
                indices[i], indices[-j] = indices[-j], indices[i]
                yield tuple(pool[i] for i in indices[:r])
                break
        else:
            return


def combinations(iterable, r=None):
    """
    Generate r-length tuples, in sorted order, no repeated elements.
    
    combinations('ABCD', 2) --> AB AC AD BC BD CD
    combinations(range(4), 3) --> 012 013 023 123
    """
    pool = tuple(iterable)
    n = len(pool)
    r = n if r is None else r
    if r > n:
        return
    indices = range(r)
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)


"""
Recipes from http://docs.python.org/library/itertools.html
"""
def powerset(iterable, R=None):
    """
    For the set S = iterable, generate the set of all subsets of S, including the empty set and S itself.
    The number of subsets should be in R, where R = None means generate all powersets.

    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    powerset([1,2,3],[1,2]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3)
    """
    s = list(iterable)
    R = range(len(s)+1) if R is None else R
    return from_iterable(combinations(s, r) for r in R)


def combinations_with_replacement(iterable, r):
    """
    Generate r-length tuples, in sorted order, with repeated elements sampled from iterable.
    
    combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC"
    number items returned: (n+r-1)! / r! / (n-1)!
    """
    pool = tuple(iterable)
    n = len(pool)
    if not n and r:
        return
    indices = [0] * r
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != n - 1:
                break
        else:
            return
        indices[i:] = [indices[i] + 1] * (r - i)
        yield tuple(pool[i] for i in indices)


def roundrobin(*iterables):
    """
    Iterate through the iterables, selecting one element from each until exhausted.
    
    roundrobin('ABC', 'D', 'EF') --> A D E B F C
    Recipe credited to George Sakkis
    """
    pending = len(iterables)
    nexts = cycle(iter(it).next for it in iterables)
    while pending:
        try:
            for next in nexts:
                yield next()
        except StopIteration:
            pending -= 1
            nexts = cycle(islice(nexts, pending))


def random_product(*args, **kwds):
    "Random selection from itertools.product(*args, **kwds)"
    pools = map(tuple, args) * kwds.get('repeat', 1)
    return tuple(random.choice(pool) for pool in pools)


def random_permutation(iterable, r=None):
    "Random selection from itertools.permutations(iterable, r)"
    pool = tuple(iterable)
    r = len(pool) if r is None else r
    return tuple(random.sample(pool, r))


def random_combination(iterable, r):
    "Random selection from itertools.combinations(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(xrange(n), r))
    return tuple(pool[i] for i in indices)


def random_combination_with_replacement(iterable, r):
    "Random selection from itertools.combinations_with_replacement(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.randrange(n) for i in xrange(r))
    return tuple(pool[i] for i in indices)



"""
New functions
"""
def num_combinations(n,r):
    """Returns the number of combinations"""
    # n!/(r!*(n-r)!) if 0 <= r <= n, or 0 if r > n
    if r>=0 and r<=n:
        return factorial(n)/(factorial(r)*factorial(n-r))
    if r>n:
        return 0
    return


def num_permutations(n,r=None):
    """Returns the number of permutations"""
    # n!/(n-r)! if 0 <= r <= n, or 0 if r > n
    r = n if r is None else r
    if r>=0 and r<=n:
        return factorial(n)/factorial(n-r)
    if r>n:
        return 0
    return


def num_powerset(n,R=None):
    """Returns the number of powersets"""
    # sum_{r in R}{num_combinations(n,r)}
    R = range(n+1) if R is None else R
    return sum([num_combinations(n,r) for r in R]) 


def num_powerset_orderedsubsets(n,R=None):
    """Returns the number of powersets with ordered subsets"""
    # sum_{r in R}{num_permutations(iterable,r)}
    R = range(n+1) if R is None else R
    return sum([num_permutations(n,r) for r in R])


def powerset_orderedsubsets(iterable,R=None):
    """
    For set S = iterable, generate the set of all ordered subsets of S, including the empty set and S itself.
    The number of subsets should be in R, where R = None means generate all powersets with ordered subsets.
    
    powerset_orderedsubsets([1,2,3]) --> () (1,) (2,) (3,)
                                         (1,2) (1,3) (2,1) (2,3) (3,1) (3,2)
                                         (1,2,3) (1,3,2) (2,1,3) (2,3,1) (3,1,2) (3,2,1)
    powerset_orderedsubsets([1,2,3],[1,2]) --> (1,) (2,) (3,)
                                               (1,2) (1,3) (2,1) (2,1) (3,1) (3,2)
    """
    s = list(iterable)
    R = range(len(s)+1) if R is None else R
    return from_iterable(permutations(s, r) for r in R)



def make_list_with_counts(iterable,counts):
    """
    Generate a list using the elements from iterable and the specified counts.

    make_list_with_counts('ABC',[0,1,2]) --> BCC
    """
    # repeat element in iterable for the specified count
    for i in xrange(len(iterable)):
        for r in repeat(iterable[i],counts[i]):
            yield r


combine_lists = product


"""
Partition Functions
"""
def partition(iterable, chain=chain, map=map):
    """
    Generate distinct partitions of a sequence.
    Recipe from Python Codebook (http://code.activestate.com/recipes/576795)
    """
##    s = iterable if hasattr(iterable, '__getslice__') else tuple(iterable)
##    n = len(s)
##    first, middle, last = [0], range(1, n), [n]
##    getslice = s.__getslice__
##    return [map(getslice, chain(first, div), chain(div, last))
##            for i in range(n) for div in combinations(middle, i)]

    # edited to be generator and return tuples
    s = iterable if hasattr(iterable, '__getslice__') else tuple(iterable)
    n = len(s)
    first, middle, last = [0], range(1, n), [n]
    getslice = s.__getslice__
    for i in range(n):
        for div in combinations(middle, i):
            yield map(getslice, chain(first, div), chain(div, last))



"""Old Code"""
##def filter_partitions(partitions):
##    """
##    Filter partitions so that elements are treated as unique based on their value.
##    Elements are treated as unique based on their value,
##    NOT based on their location (unlike other combinatoric functions).
##    """
##    for partition in partitions[:]:
##        keep = True
##        for s in partition:
##            if len(set(s)) != len(s):
##                keep = False
##
##        if keep:
##            yield(partition)


# import combinatorics from other files, 
# so that only this file needs to be included
try:
    from yjw.integerpartitions import *
except ImportError:
    try:
        from integerpartitions import *
    except ImportError:
        pass

try:
    from yjw.setpartitions import *
except ImportError:
    try:
        from setpartitions import *
    except ImportError:
        pass
