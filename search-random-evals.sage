import random
from tools.exhaust import RegeneratingRSSchemeFinder


def random_combination(iterable, r):
    """
    Random selection from itertools.combinations(iterable, r).
    <https://docs.python.org/2/library/itertools.html#recipes>
    """
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(xrange(n), r))
    return tuple(pool[i] for i in indices)


def main():
    F.<g> = GF(256)
    count = 0
    best_count = 0
    best_bw = Infinity
    while True:
        finder = RegeneratingRSSchemeFinder(
            F=F,
            t=2,
            n=5,
            k=3,
            evals=random_combination(F, 5),
        )
        bw = finder.exhaust(
            good_enough=16,
            search=RegeneratingRSSchemeFinder.Search.LINES,
            outf='output/output-%d.tex' % count,
        )
        if bw < best_bw:
            best_bw = bw
            best_count = count
        count += 1
        print 'Best evals: index %d (total bw %d)' % (best_count, best_bw)
        print '=' * 40
        print

if __name__ == '__main__':
    main()
