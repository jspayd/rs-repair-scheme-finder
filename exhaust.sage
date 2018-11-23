###############################################################################
##
## This file finds a linear repair scheme for the (10,14)-GRS code implemented
## at https://github.com/facebookarchive/hadoop-20/tree/master/src/contrib/raid/src/java/org/apache/hadoop/raid
##
## This code is a GRS code with evaluation points A = {1, a, a^2, ..., a^13},
## where a is a primitive element of GF(2^8).
##
## By our work, it suffices to find, for each alpha^* = a^{i*} in A, a set of
## two cubic polynomials over GF(2^4), so that p_1(alpha), p_2(alpha) are
## linearly dependent over GF(2^4) for many alpha \neq alpha^*, but
## p_1(alpha*), p_2(alpha*) are linearly independent.
##
## This file exhaustively searches over cubic polynomials over GF(2^8) which
## have three roots in A to find such a scheme.
##
##    -- Mary Wootters, September 2015
##
###############################################################################

n = 5
k = 3
t = 2
subfield = 16

# Set up our field and polynomial ring:
# a is a primitive element for F.
F.<a> = GF(subfield^t)
R.<X> = PolynomialRing(F)

#
# returns true if x and y are linearly independent over B.  This is a pretty
# easy test :)
#
def linIndOverB(x, y):
    return x^(subfield-1) == y^(subfield-1)

#
# given a set of values in F, returns the rank of that set over B.
#
def rankOverB(vals):
    eq_classes = []
    for val in vals:
        if val == 0: continue
        placed = False
        for eq_class in eq_classes:
            if not linIndOverB(val, eq_class[0]):
                eq_class.append(val)
                placed = True
        if not placed:
            eq_classes.append([val])
    return len(eq_classes)

#
# returns None if the scheme corresponding to P will fail to recover f(alpha^*)
# otherwise returns the number of bits needed to recover f(alpha^*) with this
# scheme
#
# P should consist of t polynomials in R; istar is an index so that alpha^* =
# evals[istar].
#
def checkSchemeOverB(evals, P, istar):
    if len(P) != t:
        print 'P should consist of only', t, 'polynomials!'
        return None

    s = [ p(evals[istar]) for p in P ]
    if rankOverB(s) != t:
        return None
    ret = 0
    for i in range(n):
        if i != istar:
            t_ = [ p(evals[i]) for p in P ]
            ret += rankOverB(t_)
    ret = ret * log(subfield, 2)
    return ret

#
# For a fixed alpha^* = evals[istar], exhaust over all linear repair schemes
# over B. If factored is True, this is limited to the repair schemes where the
# polynomials each have n-k-1 roots in the evaluation set. This is pretty
# restrictive, but already it's good enough to get nontrivial results.
#
# goodEnough is a threshold: stop searching if you find a scheme with that many
# bits or fewer.
#
def exhaust_over_B(evals=[ a^i for i in range(n) ], goodEnough=64, istar=0,
                     factored=False):
    bestBW = Infinity
    bestPolys = []
    count = 0
    drawFrom = None
    indices_list = None
    subsets = None
    if factored:
        subsets = Subsets(evals, n-k-1)
    else:
        subsets = Subsets(F.list() * (n-k), n-k, submultiset=True)
    # Combinations of Subsets convert Subsets to lists, which is very slow, so
    # we do combinations of indices.
    indices_list = Combinations(range(len(subsets)), t)
    for indices in indices_list:
        count += 1
        T = [ subsets[index] for index in indices ]
        P = None
        if factored:
            P = [ prod([ (X + x) for x in T[j] ]) for j in range(len(T)) ]
        else:
            P = [ sum([ (X^j * x) for x, j in zip(T[j], range(n-k)) ])
                    for j in range(len(T)) ]
        bw = checkSchemeOverB(evals, P, istar)
        if bw != None:
            if bw < bestBW:
                bestBW = bw
                bestPolys = P
            # if bw <= goodEnough:
            #    print 'SCHEME WITH BW=', bw, ':', p0, p1

        if count % 1000 == 0:
            print 'Check', count, 'best is', bestBW, 'with', bestPolys
        if bestBW <= goodEnough:
            return bestBW, bestPolys
    return bestBW, bestPolys

#
# Run exhaust_over_B for each alpha^*, and write the results to outf (as
# LaTeX code).
#
def exhaustMore(evals=[ a^i for i in range(n) ], goodEnough=64,
                outf='output.tex', factored=False):
    print ('Exhausting over polynomials with coefficients in F_%d for t=%d, '
           'n=%d, k=%d, and evaluation points %s.'
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
        bw, polys = exhaust_over_B(evals, goodEnough, i, factored)
        F.write('$ \\zeta^{' + evals[i]._log_repr() +  '}$ & ')
        print i, bw
        for p in polys:
            polyStr = None
            if factored:
                polyStr = prettyPrintFactored(p)
            else:
                polyStr = prettyPrint(p)
            print polyStr
            F.write('$  ' + polyStr + '$ & ')
        print '===================='
        F.write(str(bw) + '\\\\ \hline \n')
    F.write('\\end{tabular}\n')
    F.write('\\end{document}')
    F.close()

#
# print out a polynomial as nice latex code.
#
def prettyPrint(p):
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
        if i < len(c) -1:
            ret +=  ' + '
    return ret

#
# print out a polynomial with roots in r as nice latex code
#
def prettyPrintFactored(p):
    rts = p.roots()
    ret = ''
    for rt, multiplicity in rts:
        if multiplicity != 1:
            print 'Now that\'s strange'
            return None
        ret += '(\\mathbf{X} + \\zeta^{' + rt._log_repr() + '})'
    return ret

def main():
    exhaustMore(goodEnough=16, factored=True)

if __name__ == '__main__':
    main()
