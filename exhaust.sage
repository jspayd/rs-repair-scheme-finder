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


def combinations_count(p, r):
    """
    Returns the number of combinations of length r drawn from p items.
    """
    result = 1
    for x in range(p - r + 1, p + 1):
        result *= x
    result /= factorial(r)
    return result


def product_count(p, r):
    """
    Returns the number of elements in the r-fold Cartesian product of a set of
    cardinality p with itself.
    """
    return p^r


def pretty_power(p):
    """
    Given some p in F, return a nice LaTeX string representation of that
    number.
    """
    if p == 0:
        return '0'
    elif p == 1:
        return '1'
    exp = int(p._log_repr())
    if exp == 1:
        return '\\zeta'
    else:
        return '\\zeta^{%d}' % exp


def pretty_print(p):
    """
    Print out a polynomial as nice LaTeX code.
    """
    c = p.coefficients(sparse=False)
    ret = ''
    for i in range(len(c)):
        printed = False
        if c[i] == 0:
            continue
        elif c[i] != 1:
            ret += pretty_power(c[i])
            printed = True
        if i == 0 and not printed:
            ret += '\\,1'
        elif i == 1:
            ret += '\\,\\mathbf{X}'
        elif i != 0:
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
        zeta = None
        exp = int(rt._log_repr())
        if p == 1:
            zeta = '1'
        elif exp == 1:
            zeta = '\\zeta'
        else:
            zeta = '\\zeta^{' + rt._log_repr() + '}'
        ret += '(\\mathbf{X} + ' + zeta + ')'
    return ret


class RegeneratingRSSchemeFinder:
    """
    A class for finding regenerating code repair schemes for Reed-Solomon
    codes.
    """

    def __init__(self, F, t, n, k, evals=None):
        """
        Creates a RegeneratingRSSchemeFinder for a code over field F with
        symbol size t, message length k, and codeword length n.
        """
        self.F = F
        self.a = F.gen()
        R.<X> = PolynomialRing(F)
        self.R = R
        self.X = X
        self.t = t
        self.n = n
        self.k = k
        self.B_size = len(F) ^ (1 / t)
        self.B_bits = log(self.B_size, 2)
        if evals is None:
            self.evals=[self.a^i for i in range(n)]
        else:
            self.evals = evals


    def __lin_ind_over_B(self, x, y):
        """
        Returns true if x and y are linearly independent over B. This is a
        pretty easy test :)
        """
        return x^(self.B_size-1) != y^(self.B_size-1)


    def __rank_over_B(self, vals):
        """
        Given a set of values in F, returns the rank of that set over B.
        """
        eq_classes = []
        for val in vals:
            if val == 0:
                continue
            placed = False
            for eq_class in eq_classes:
                if not self.__lin_ind_over_B(val, eq_class[0]):
                    eq_class.append(val)
                    placed = True
            if not placed:
                eq_classes.append([val])
        return len(eq_classes)


    def __check_scheme_over_B(self, P, istar):
        """
        Returns None if the scheme corresponding to P will fail to recover
        f(alpha^*) otherwise returns the number of bits needed to recover
        f(alpha^*) with this scheme.

        P should consist of t polynomials in R; istar is an index so that
        alpha^* = evals[istar].
        """
        if len(P) != self.t:
            print 'P should consist of only', t, 'polynomials!'
            return None

        s = [p(self.evals[istar]) for p in P]
        if self.__rank_over_B(s) != self.t:
            return None
        ret = 0
        for i in range(self.n):
            if i != istar:
                vals = [p(self.evals[i]) for p in P]
                ret += self.__rank_over_B(vals)
        ret = ret * self.B_bits
        return ret


    def __exhaust_over_B(self,
                         good_enough=64,
                         istar=0,
                         factored=False):
        """
        For a fixed alpha^* = evals[istar], exhaust over all linear repair
        schemes over B. If factored is True, this is limited to the repair
        schemes where the polynomials each have n-k-1 roots in the evaluation
        set. This is pretty restrictive, but already it's good enough to get
        nontrivial results.

        goodEnough is a threshold: stop searching if you find a scheme with
        that many bits or fewer.
        """
        print 'Searching for polynomials for recovering evaluation at %s.' % (
            self.evals[istar],
        )
        best_bw = Infinity
        best_polys = []
        count = 0
        polys = None
        poly_count = None
        if factored:
            polys = combinations(self.evals, self.n-self.k-1)
            poly_count = combinations_count(len(self.evals), self.n-self.k-1)
        else:
            polys = product(self.F, repeat=self.n-self.k)
            poly_count = product_count(len(self.F), self.n-self.k)
        schemes = combinations(polys, self.t)
        num_schemes = combinations_count(poly_count, self.t)
        for T in schemes:
            count += 1
            P = None
            if factored:
                P = [prod([(self.X + x) for x in T_cur]) for T_cur in T]
            else:
                P = [sum([(self.X^j * x) for x, j
                     in zip(T[j], range(self.n-self.k))])
                     for j in range(len(T))
                    ]
            bw = self.__check_scheme_over_B(P, istar)
            if bw is not None:
                if bw < best_bw:
                    best_bw = bw
                    best_polys = P

            if count % 1000 == 0:
                print ('Checked scheme %d/%d (%.1f%%). '
                       'Best scheme (bw %s): %s.'
                      ) % (
                    count,
                    num_schemes,
                    count * 100.0 / num_schemes,
                    best_bw,
                    best_polys,
                )

            if best_bw <= good_enough:
                print
                return best_bw, best_polys
        print
        return best_bw, best_polys


    def exhaust(self,
                good_enough=64,
                outf='output.tex',
                factored=False):
        """
        Run __exhaust_over_B for each alpha^*, and write the results to outf
        (as LaTeX code).
        """
        print ('Exhausting over polynomials with coefficients in F_%d for '
               't=%d, n=%d, k=%d, and evaluation points %s.'
              ) % (len(self.F), self.t, self.n, self.k, self.evals)
        if factored:
            print 'Only checking polynomials with %d roots in %s.' % (
                self.n-self.k-1,
                self.evals,
            )
        print

        F = open(outf, 'w')
        F.write('\\documentclass{article}\n\n')
        F.write('\\begin{document}')
        F.write('\\begin{tabular}{|c|cc|c|}\n')
        F.write('\\hline')
        F.write(('$ \\alpha^*$   & Polynomials & &Bandwidth (in bits) for $ '
                 '\\alpha^*$   \\\\ \n'
                ))
        F.write('\\hline')
        for i in range(self.n):
            bw, polys = self.__exhaust_over_B(good_enough, i, factored)
            F.write('$ ' + pretty_power(self.evals[i]) + '$ & ')
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


def main():
    finder = RegeneratingRSSchemeFinder(
        F=GF(256),
        t=2,
        n=5,
        k=3,
    )
    finder.exhaust(
        good_enough=20,
        factored=False,
    )


if __name__ == '__main__':
    main()
