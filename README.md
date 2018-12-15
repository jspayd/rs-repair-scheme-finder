This program finds linear repair schemes for Reed-Solomon codes. It is based on
code by Venkatesan Guruswami and Mary Wootters for finding a scheme for the
(14,10)-GRS code used in [Facebook's Hadoop analytics
cluster](https://github.com/facebookarchive/hadoop-20/tree/master/src/contrib/raid/src/java/org/apache/hadoop/raid).
For a code over finite field *F* with message length *k*, codeword length *n*,
and evaluation points *A*, by the [work of Guruswami and
Wootters](https://arxiv.org/abs/1509.04764v2), it suffices to find, for each
*a\** in A, a set *P* of *t* (*n*-*k*-1)-degree polynomials over a subfield *B*
of *F* where *t* is the degree of *F* over *B*, so that the set {*p*(*a*) : *p*
in *P*} is low rank over *B* for all *a* in *A* \ {*a\**}, and full rank for
*a* = *a\**. This file has three kinds of searches:
  1) exhausting over all (*n*-*k*-1)-degree polynomials over *F*.
  2) exhausting over all polynomials with (*n*-*k*-1) roots in *F*.
  3) a special search designed specifically for the case that (*n*-*k*) = 2 and
     *t* = 2 that searches over sets of polynomials guaranteed to have linearly
     dependent evaluations for at least two evaluation points.
