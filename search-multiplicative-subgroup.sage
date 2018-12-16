from tools.exhaust import RegeneratingRSSchemeFinder

def main():
    F.<g> = GF(256)
    finder = RegeneratingRSSchemeFinder(
        F=F,
        t=2,
        n=5,
        k=3,
        evals=[x for x in F if x^5 == 1],
    )
    finder.exhaust(
        good_enough=16,
        search=RegeneratingRSSchemeFinder.Search.LINES,
        outf='output/multiplicative-subgroup.tex',
    )

if __name__ == '__main__':
    main()
