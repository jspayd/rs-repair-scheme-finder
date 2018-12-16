from tools.exhaust import RegeneratingRSSchemeFinder

def main():
    F.<g> = GF(256)
    finder = RegeneratingRSSchemeFinder(
        F=F,
        t=2,
        n=5,
        k=3,
    )
    finder.exhaust(
        good_enough=16,
        search=RegeneratingRSSchemeFinder.Search.LINES,
        outf='output/default-points.tex',
    )

if __name__ == '__main__':
    main()
