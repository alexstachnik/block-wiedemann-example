import Poly
import Mat
import pdb

#Convert a polynomial with matrix coefficients to a matrix with
#polynomial entries
def polyMatToMatPoly(Ms,field):
    ring=Poly.Poly(field)
    m=Ms[0].m
    n=Ms[0].n
    A=Mat.Mat(ring,m,n)
    for k in xrange(len(Ms)):
        for i in xrange(m):
            for j in xrange(n):
                d=A.getElt(i,j)
                d=ring.timesXPower(d,1)
                e=ring.fromScalar(Ms[k].getElt(i,j))
                d=ring.add(d,e)
                A.setElt(i,j,d)
    return A
