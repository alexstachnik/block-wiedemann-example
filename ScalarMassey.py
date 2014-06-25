
import Poly

# Implementation of the Berlekamp-Massey algorithm for finding the
# min-poly of a linearly generated sequence of scalars values

def disc(poly,seq,n,polyRing,field):
    total=field.zero()
    for i in xrange(polyRing.degree(poly)+1):
        d=field.mul(polyRing.getElt(poly,i),seq[n-i])
        total=field.add(total,d)
    return total

def scalarMassey(seq,field):
    polyRing=Poly.Poly(field)
    C=[field.one()]
    B=[field.one()]
    x=1
    L=0
    b=field.one()
    N=0
    for N in xrange(len(seq)):
        d=disc(C,seq,N,polyRing,field)
        if (field.isZero(d)):
            x=x+1
        else:
            if 2*L > N:
                mB=polyRing.scalarMul(B,field.div(d,b))
                C=polyRing.sub(C,polyRing.timesXPower(mB,x))
                x=x+1
            else:
                mB=polyRing.scalarMul(B,field.div(d,b))
                B=C
                C=polyRing.sub(C,polyRing.timesXPower(mB,x))
                L=N+1-L
                b=d
                x=1
    C.reverse()
    return C


