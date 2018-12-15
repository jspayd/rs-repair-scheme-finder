from tools.exhaust import RegeneratingRSSchemeFinder

def main()
    F.<g> = GF(256^2)
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
            special=True,
            outf='output/output-%d.tex' % count,
        )
        if bw < best_bw:
            best_bw = bw
            best_count = count
        count += 1
        print 'Best evals: index %d (total bw %d)' % (best_count, best_bw)
        print '=' * 40
        print
