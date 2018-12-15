from tools.exhaust import RegeneratingRSSchemeFinder

def main():
    F.<g> = GF(256^2)
    finder = RegeneratingRSSchemeFinder(
        F=F,
        t=2,
        n=5,
        k=3,
    )
    finder.exhaust(
        good_enough=32,
        special=True,
        outf='output/large-field.tex',
    )

if __name__ == '__main__':
    main()
