"""
This file finds a linear repair scheme for the (10,14)-GRS code implemented at
https://github.com/facebookarchive/hadoop-20/tree/master/src/contrib/raid/src/java/org/apache/hadoop/raid

This code is a GRS code with evaluation points A = {1, a, a^2, ..., a^13},
where a is a primitive element of GF(2^8).

By our work, it suffices to find, for each alpha^* = a^{i*} in A, a set of
two cubic polynomials over GF(2^4), so that p_1(alpha), p_2(alpha) are
linearly dependent over GF(2^4) for many alpha \neq alpha^*, but
p_1(alpha*), p_2(alpha*) are linearly independent.

This file exhaustively searches over cubic polynomials over GF(2^8) which
have three roots in A to find such a scheme.

   -- Mary Wootters, September 2015
"""

from itertools import combinations, product
import sys

PADDING = 10

# Constants related to particular repair scheme
n = 5
k = 3
t = 2
subfield = 16
bits = log(subfield, 2)


# Set up our field and polynomial ring:
# a is a primitive element for F.
F.<a> = GF(subfield^t)
R.<X> = PolynomialRing(F)


def lin_ind_over_B(x, y):
    """
    Returns true if x and y are linearly independent over B. This is a pretty
    easy test :)
    """
    return x^(subfield-1) != y^(subfield-1)


def rank_over_B(vals):
    """
    Given a set of values in F, returns the rank of that set over B.
    """
    eq_classes = []
    for val in vals:
        if val == 0:
            continue
        placed = False
        for eq_class in eq_classes:
            if not lin_ind_over_B(val, eq_class[0]):
                eq_class.append(val)
                placed = True
        if not placed:
            eq_classes.append([val])
    return len(eq_classes)


def check_scheme_over_B(evals, P, istar):
    """
    Returns None if the scheme corresponding to P will fail to recover
    f(alpha^*) otherwise returns the number of bits needed to recover
    f(alpha^*) with this scheme.

    P should consist of t polynomials in R; istar is an index so that alpha^* =
    evals[istar].
    """
    if len(P) != t:
        print 'P should consist of only', t, 'polynomials!'
        return None

    s = [p(evals[istar]) for p in P]
    if rank_over_B(s) != t:
        return None
    ret = 0
    for i in range(n):
        if i != istar:
            t_ = [p(evals[i]) for p in P]
            ret += rank_over_B(t_)
    ret = ret * bits
    return ret


def combinations_count(p, r):
    result = 1
    for x in range(p - r + 1, p + 1):
        result *= x
    result /= factorial(r)
    return result


def product_count(p, r):
    return p^r


def exhaust_over_B(evals=[a^i for i in range(n)],
                   good_enough=64,
                   istar=0,
                   factored=False):
    """
    For a fixed alpha^* = evals[istar], exhaust over all linear repair schemes
    over B. If factored is True, this is limited to the repair schemes where
    the polynomials each have n-k-1 roots in the evaluation set. This is pretty
    restrictive, but already it's good enough to get nontrivial results.

    goodEnough is a threshold: stop searching if you find a scheme with that
    many bits or fewer.
    """
    print 'Searching for polynomials for recovering evaluation at %s.' % (
        evals[istar],
    )
    best_bw = Infinity
    best_polys = []
    count = 0
    polys = None
    poly_count = None
    if factored:
        polys = combinations(evals, n-k-1)
        poly_count = combinations_count(len(evals), n-k-1)
    else:
        polys = product(F, repeat=n-k)
        poly_count = product_count(len(F), n-k)
    schemes = combinations(polys, t)
    num_schemes = combinations_count(poly_count, t)
    for T in schemes:
        count += 1
        print ('\rChecking scheme %d/%d (%.1f%%). '
               'Best scheme (bw %s): %s.%s'
              ) % (
            count,
            num_schemes,
            count * 100.0 / num_schemes,
            best_bw,
            best_polys,
            ' ' * PADDING,
        ),
        sys.stdout.flush()
        P = None
        if factored:
            P = [prod([(X + x) for x in T_cur]) for T_cur in T]
        else:
            P = [sum([(X^j * x) for x, j in zip(T[j], range(n-k))])
                 for j in range(len(T))
                ]
        bw = check_scheme_over_B(evals, P, istar)
        if bw is not None:
            if bw < best_bw:
                best_bw = bw
                best_polys = P

        if best_bw <= good_enough:
            print
            return best_bw, best_polys
    print
    return best_bw, best_polys


def exhaust_more(evals=[a^i for i in range(n)],
                 good_enough=64,
                 outf='output.tex',
                 factored=False):
    """
    Run exhaust_over_B for each alpha^*, and write the results to outf (as
    LaTeX code).
    """
    print ('Exhausting over polynomials with coefficients in F_%d for t=%d, '
           'n=%d, k=%d, and evaluation points %s.\n'
          ) % (subfield^t, t, n, k, evals)
    if factored:
        print 'Only checking polynomials with %d roots in %s.' % (n-k-1, evals)

    F = open(outf, 'w')
    F.write('\\documentclass{article}\n\n')
    F.write('\\begin{document}')
    F.write('\\begin{tabular}{|c|cc|c|}\n')
    F.write('\\hline')
    F.write(('$ \\alpha^*$   & Polynomials & &Bandwidth (in bits) for $ '
             '\\alpha^*$   \\\\ \n'
            ))
    F.write('\\hline')
    for i in range(n):
        bw, polys = exhaust_over_B(evals, good_enough, i, factored)
        F.write('$ \\zeta^{' + evals[i]._log_repr() + '}$ & ')
        print i, bw
        for p in polys:
            poly_str = None
            if factored:
                poly_str = pretty_print_factored(p)
            else:
                poly_str = pretty_print(p)
            print poly_str
            F.write('$  ' + poly_str + '$ & ')
        print '===================='
        F.write(str(bw) + '\\\\ \\hline \n')
    F.write('\\end{tabular}\n')
    F.write('\\end{document}')
    F.close()


def pretty_print(p):
    """
    Print out a polynomial as nice LaTeX code.
    """
    c = p.coefficients(sparse=False)
    ret = ''
    for i in range(len(c)):
        if c[i] == 0:
            continue
        elif c[i] == 1:
            pass
        else:
            ret += '\\zeta^{' + c[i]._log_repr() + '}'
        ret += '\\,\\mathbf{X}^{' + str(i) + '}'
        if i < len(c) - 1:
            ret += ' + '
    return ret


def pretty_print_factored(p):
    """
    Print out a polynomial with roots in r as nice LaTeX code.
    """
    rts = p.roots()
    ret = ''
    for rt, multiplicity in rts:
        if multiplicity != 1:
            print 'Now that\'s strange'
            return None
        ret += '(\\mathbf{X} + \\zeta^{' + rt._log_repr() + '})'
    return ret


def main():
    exhaust_more(good_enough=16, factored=True)


if __name__ == '__main__':
    main()
