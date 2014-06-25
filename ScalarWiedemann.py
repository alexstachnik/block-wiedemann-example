
import Mat
import ScalarMassey
import Poly

# Wiedemann's original paper describes two algorithms, named here
# wiedemann1 and wiedemann2. I think only wiedemann2 actually works

# Ref: Solving sparse linear equations over finite fields,
# D. H. Wiedemann, IEEE Trans. Inform. Theory (1986)

def genSeq(L,u,A,b):
    retVal=[b]
    L=max(L,0)
    for i in xrange(L):
        retVal.append(A.apply(retVal[-1]))
    for i in xrange(L+1):
        retVal[i]=u.apply(retVal[i]).getElt(0,0)
    return retVal

def wiedemann2Wrapper(A,field):
    n=A.n
    while True:
        v=Mat.Mat(field,n,1)
        v.randomize()
        fAv=wiedemann2(A,v,(2*n)+1,field)
        if isMinPoly(fAv,A):
            return fAv

def wiedemann2(A,v,d,field):
    polyRing=Poly.Poly(field)
    n=A.n
    ident=Mat.ident(field,n)
    u=Mat.Mat(field,1,n)
    u.randomize()
    matSeq=[v]
    for i in xrange(2*A.n-1):
        matSeq.append(A.apply(matSeq[-1]))
    for i in xrange(2*A.n):
        matSeq[i]=u.apply(matSeq[i]).getElt(0,0)
    S=ScalarMassey.scalarMassey(matSeq,field)
    if d==polyRing.degree(S):
        return S
    else:
        fA=polyRing.apply(S,
                          lambda a,b:a.add(b),
                          lambda a,b:a.apply(b),
                          lambda a,b:b.scalarMul(a),
                          ident,
                          A)
        vPrime=fA.apply(v)
        if vPrime.isZero():
            return S
        else:
            fPrime=wiedemann2(A,vPrime,d-polyRing.degree(S),field)
            return polyRing.mul(fPrime,S)

def isMinPoly(p,A):
    ring=A.ring
    polyRing=Poly.Poly(ring)
    ident=Mat.ident(ring,A.n)
    pA=polyRing.apply(p,
                      lambda a,b:a.add(b),
                      lambda a,b:a.apply(b),
                      lambda a,b:b.scalarMul(a),
                      ident,
                      A)
    return pA

def wiedemann1(A,b,field):
    b0=b
    k=0
    n=A.n
    y=Mat.Mat(field,n,1)
    d=0
    polyRing=Poly.Poly(field)
    ident=Mat.ident(field,n)
    while not field.isZero(b.nnz()):
        u=Mat.Mat(field,1,n)
        u.randomize()
        S=genSeq(2*(50+n-d),u,A,b)
        f=ScalarMassey.scalarMassey(S,field)
        fA=polyRing.apply(polyRing.polyMinus(f),
                          lambda a,b:a.add(b),
                          lambda a,b:a.apply(b),
                          lambda a,b:b.scalarMul(a),
                          ident,
                          A)
        fAPrime=polyRing.apply(f,
                          lambda a,b:a.add(b),
                          lambda a,b:a.apply(b),
                          lambda a,b:b.scalarMul(a),
                          ident,
                          A)
        y=y.add(fA.apply(b))
        y=y.scalarMul(field.sub(field.zero(),field.inv(f[0])))
        b=b0.sub(A.apply(y))
        d=d+polyRing.degree(f)
    return y
