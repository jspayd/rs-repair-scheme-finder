###############################################################################
##
## This file finds a linear repair scheme for the (10,14)-GRS code implemented
## at https://github.com/facebookarchive/hadoop-20/tree/master/src/contrib/raid/src/java/org/apache/hadoop/raid
##
## This code is a GRS code with evaluation points A = {a, a^2, ..., a^13},
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
    if x^(15) == y^(15):
        return False
    return True

#
# returns None if the scheme corresponding to P will fail to recover f(alpha^*)
# otherwise returns the number of bits needed to recover f(alpha^*) with this
# scheme
#
# P should consist of two polynomials in R; istar is an index so that alpha^* =
# a^(istar).
#
def checkSchemeOverB(P, istar):
    if len(P) != 2:
        print('P should consist of only two polynomials!')
        return None
    s0 = P[0](a^(istar))
    s1 = P[1](a^(istar))
    if (not linIndOverB(s0, s1)) or s0 == 0 or s1 == 0:
        return None
    ret = 0
    for i in range(n):
        if i != istar:
            bw = 0
            t0 = P[0](a^i)
            t1 = P[1](a^i)
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
def exhaust_over_F16(goodEnough=64, istar=0):
    bestBW = Infinity
    bestPolys = []
    count = 0
    S = [ a^i for i in range(1,n) ]
    for T0 in Subsets(F, n-k):
        for T1 in Subsets(F, n-k):
            count += 1
            p0 = sum([ (X^j * x) for x, j in zip(T0, range(n-k)) ])
            p1 = sum([ (X^j * x) for x, j in zip(T1, range(n-k)) ])
            P = [p0, p1]
            bw = checkSchemeOverB(P, istar)
            if bw != None:
                if bw < bestBW:
                    bestBW = bw
                    bestPolys = P
                if bw <= goodEnough:
                    print 'SCHEME WITH BW=',bw,':',p0,p1

            # if count%1000 == 0:
            #     print "Check",count,"best is",bestBW,"with",bestPolys
            if bestBW <= goodEnough:
                return bestBW, bestPolys
    return bestBW, bestPolys
#
# Run exhaust_over_F16 for each alpha^*, and write the results to outf (as
# LaTeX code).
#
def exhaustMore(goodEnough=64, outf='output.txt'):
    F = open(outf, 'w')
    F.write('\\begin{tabular}{|c|cc|c|}\n') 
    F.write(('$ \\alpha^*$   & Polynomials & &Bandwidth (in bits) for $ '
                 '\\alpha^*$   \\\\ \n'
                ))
    F.write('\\hline')
    for i in range(n):
        bw, polys = exhaust_over_F16( goodEnough, i )
        F.write('$ \\zeta^{' + str(i) +  '}$ & ')
        print i, bw
        for p in polys:
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
    exhaustMore(goodEnough=23)

if __name__ == '__main__':
    main()
