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

# Set up our field and polynomial ring:
# a is a primitive element for F; in this version of sage, it is a root of 
# x^8 + x^4 + x^3 + x^2 + 1
F.<a> = GF(2^8)
R.<X> = PolynomialRing(F)

#
# returns true if x and y are linearly independent over B = GF(2^4).  This is a
# pretty easy test :)
#
def linIndOverB(x,y):
    return x^(15) == y^(15)

#
# returns None if the scheme corresponding to P will fail to recover f(alpha^*)
# otherwise returns the number of bits needed to recover f(alpha^*) with this
# scheme
#
# P should consist of two polynomials in R; istar is an index so that alpha^* =
# a^(istar).
#
# Currently only works for B = F_16 and F = F_256
#
def checkSchemeOverB(evals, P, istar):
    if len(P) != t:
        print 'P should consist of only', t, 'polynomials!'
        return None
    s0 = P[0](evals[istar])
    s1 = P[1](evals[istar])
    if (not linIndOverB(s0, s1)) or s0 == 0 or s1 == 0:
        return None
    ret = 0
    for i in range(n):
        if i != istar:
            bw = 0
            t0 = P[0](evals[i])
            t1 = P[1](evals[i])
            if t0 == 0 and t1 != 0:
                bw += 1
            elif t0 != 0 and t1 == 0:
                bw += 1
            elif t0 == 0 and t1 == 0:
                bw += 0
            elif not linIndOverB(t0, t1):
                bw += 1
            else:
                bw += 2
            ret += bw
    ret = ret*4
    return ret

#
# For a fixed alpha^* = a^(istar), exhaust over all linear repair schemes over
# B = GF(16) so that the two polynomials p0 and p1 each have three roots in the
# evaluation set. This is pretty restrictive, but already it's good enough to
# get nontrivial results.
#
# goodEnough is a threshold: stop searching if you find a scheme with that many
# bits or fewer.
#
def exhaust_over_F16(evals=[ a^i for i in range(n) ], goodEnough=64, istar=0,
                     factored=False):
    bestBW = Infinity
    bestPolys = []
    count = 0
    drawFrom = None
    if factored:
        drawFrom = evals
    else:
        drawFrom = F
    subsets = Subsets(drawFrom, n-k)
    for y in range(len(subsets)):
        for z in range(y+1, len(subsets)):
            T = [ subsets[y], subsets[z] ]
            count += 1
            if factored:
                p0 = prod([ (X + x) for x in T[0] ])
                p1 = prod([ (X + x) for x in T[1] ])
            else:
                p0 = sum([ (X^j * x) for x, j in zip(T[0], range(n-k)) ])
                p1 = sum([ (X^j * x) for x, j in zip(T[1], range(n-k)) ])
            P = [ p0, p1 ]
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
# Run exhaust_over_F16 for each alpha^*, and write the results to outf (as
# LaTeX code).
#
def exhaustMore(evals=[ a^i for i in range(n) ], goodEnough=64,
                outf='output.txt', factored=False):
    F = open(outf, 'w')
    F.write('\\begin{tabular}{|c|cc|c|}\n') 
    F.write(('$ \\alpha^*$   & Polynomials & &Bandwidth (in bits) for $ '
                 '\\alpha^*$   \\\\ \n'
                ))
    F.write('\\hline')
    for i in range(n):
        bw, polys = exhaust_over_F16(evals, goodEnough, i, factored)
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
    F.write('\end{tabular}\n')
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
        ret += '\\,\\b X^{' + str(i) + '}'
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
        ret += '(\\b X + \\zeta^{' + rt._log_repr() + '})'
    return ret

def main():
    exhaustMore(goodEnough=23, factored=False)

if __name__ == '__main__':
    main()
