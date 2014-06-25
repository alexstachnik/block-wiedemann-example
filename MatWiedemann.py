
import Mat
import Poly
from PolyAlg import *
import pdb

#Block Wiedemann implementation. Computes the matrix minimal
#polynomial from M given random U and V

#Ref: Solving Homogeneous Linear Equations Over GF(2) Via Block
#Wiedemann Algorithm, Don Coppersmith, Mathematics of Computations
#(1994)

def fixedMatWiedemannWrapper(M,b,U,V):
    field=M.ring
    n=M.n
    delta=4*n
    matSeq=[V]
    for i in xrange(delta):
        matSeq.append(M.apply(matSeq[-1]))
    for i in xrange(delta):
        matSeq[i]=U.apply(matSeq[i])
    return matWiedemann(matSeq,delta,b,field)

def matWiedemannWrapper(M,b):
    field=M.ring
    n=M.n
    delta=2*n+1
    U=Mat.Mat(field,b,n)
    U.randomize()
    V=Mat.Mat(field,n,b)
    V.randomize()
    matSeq=[V]
    for i in xrange(delta):
        matSeq.append(M.apply(matSeq[-1]))
    for i in xrange(delta):
        matSeq[i]=U.apply(matSeq[i])
    return matWiedemann(matSeq,delta,b,field)

def coeff(j,f,a,l,mu,n,field):
    tot=field.zero()
    for nu in xrange(n):
        for k in xrange(j+1):
            fElt=f[k].getElt(l,nu)
            aElt=a[j-k].getElt(nu,mu)
            tot=field.add(tot,field.mul(fElt,aElt))
    return tot

def coeffMat(j,f,a,n,field):
    result=Mat.Mat(field,n,2*n)
    for l in xrange(n):
        for mu in xrange(2*n):
            result.setElt(l,mu,coeff(j,a,f,l,mu,n,field))
    return result

def matWiedemann(Ms,delta,n,field):
    ring=Poly.Poly(field)
    f=[]
    for i in xrange(n):
        fi=Mat.Mat(field,n,2*n)
        f.append(fi)
    for i in xrange(n):
        f[0].setElt(i,i,field.one())
    d=([0]*n)+([1]*n)
    t=-1
    beta=1
    sigma=0
    mu=0
    while beta<delta-(sigma+mu+1):
        t=t+1
        CapDelta=coeffMat(t,f,Ms,n,field)
        (tau,d)=auxGauss(CapDelta,d)
        for i in xrange(n):
            d[n+i]=d[n+i]+1
        beta=d[n]
        for i in xrange(n,2*n):
            if d[i]<beta:
                beta=d[i]
        mu=d[0]
        for i in xrange(n):
            if d[i]>mu:
                mu=d[i]
        sigma=0
        for i in xrange(n):
            sigma=sigma+d[i]
        if sigma >= delta+1:
            print "Insufficient bound"
            return 0
        #f_i=f_i*tau
        #mulin last n columns of f by z
        #  -> shift last n columns up one mat
        for i in xrange(len(f)):
            f[i]=f[i].apply(tau)
        f.append(Mat.Mat(field,n,2*n))
        for k in xrange(len(f)-1,0,-1):
            for i in xrange(n):
                for j in xrange(n,2*n):
                    f[k].setElt(i,j,f[k-1].getElt(i,j))
        for i in xrange(n):
            for j in xrange(n,2*n):
                f[0].setElt(i,j,field.zero())
    result=[]
    for k in xrange(mu+1):
        result.append(Mat.Mat(field,n,n))
    for j in xrange(n):
        for k in xrange(mu+1):
            for i in xrange(n):
                result[d[j]-k].setElt(i,j,f[k].getElt(i,j))
    return result


def auxGauss(Delta,degs):
    ring=Delta.ring
    #Delta is n*2n
    n=Delta.m
    tau=Mat.ident(ring,2*n)
    gamma=set(xrange(n))
    for i in xrange(n):
        pi=set([j for j in gamma if not ring.isZero(Delta.getElt(i,j))])
        pi.add(n+i)
        l=n+i
        for j in pi:
            if degs[j]<=degs[l]:
                if degs[j] == degs[l]:
                    if l<n+i and j<l:
                        l=j
                else:
                    l=j
        pi.remove(l)
        pivel=Delta.getElt(i,l)
        for j in pi:
            if l == n+i:
                d=ring.neg(Delta.getElt(i,j))
                d=ring.div(d,pivel)
                Mat.colAdd(Delta,j,n+i,d)
                Mat.colAdd(tau,j,n+i,d)
            if l < n+i:
                if j < n+i:
                    d=ring.neg(Delta.getElt(i,j))
                    d=ring.div(d,pivel)
                    Mat.colAdd(Delta,j,l,d)
                    Mat.colAdd(tau,j,l,d)
                if j == n+i:
                    auxel=Delta.getElt(i,n+i)
                    if not ring.isZero(auxel):
                        d=ring.neg(auxel)
                        d=ring.div(d,pivel)
                        Mat.colAdd(Delta,n+i,l,d)
                        Mat.colAdd(tau,n+i,l,d)
                        Mat.colSwap(Delta,l,n+i)
                        Mat.colSwap(tau,l,n+i)
                        temp=degs[l]
                        degs[l]=degs[n+i]
                        degs[n+i]=temp
                    if ring.isZero(auxel):
                        Mat.colAdd(tau,n+i,l,ring.one())
                        Mat.colAdd(Delta,n+i,l,ring.one())
                        temp=degs[l]
                        degs[l]=degs[n+i]
                        degs[n+i]=temp
                        gamma.remove(l)
    return (tau,degs)
