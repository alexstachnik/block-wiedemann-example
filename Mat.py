import itertools
import Algorithms
import Modular

# Matrix over a ring
class Mat(object):
    def __init__(self,ring,m,n):
        self.ring=ring
        self.m=m
        self.n=n
        self.elts = [ring.zero() for i in xrange(m*n)]

    def __str__(self):
        retVal=[]
        for i in xrange(self.m):
            retVal.append('\t'.join([self.ring.toStr(self.getElt(i,j)) for j in xrange(self.n)]))
        return '\n'.join(retVal)

    def sage(self):
        retVal='['
        rows=[]
        for i in xrange(self.m):
            rows.append('['+','.join([self.ring.toStr(self.getElt(i,j)) for j in xrange(self.n)])+']')
        retVal=retVal+','.join(rows)+']'
        return retVal

    def matrixMarket(self):
        retVal=[]
        for i in xrange(self.m):
            for j in xrange(self.n):
                retVal.append(self.ring.toStr(self.getElt(i,j)))
        return '\n'.join(retVal)

    def nnz(self):
        count = 0
        for i in xrange(self.m):
            for j in xrange(self.n):
                if not self.ring.isZero(self.getElt(i,j)):
                    count = count + 1
        return count

    def write(self,file):
        file.write(str(self.m)+' '+str(self.n)+' '+str(self.nnz())+'\n')
        for i in xrange(self.m):
            for j in xrange(self.n):
                if not self.ring.isZero(self.getElt(i,j)):
                    file.write(str(i)+' '+str(j)+' '+
                               self.ring.toStr(self.getElt(i,j))+'\n')

    def getElt(self,i,j):
        return self.elts[i*self.n+j]

    def setElt(self,i,j,d):
        self.elts[i*self.n+j]=d

    def scalarMul(self,d):
        result=Mat(self.ring,self.m,self.n)
        for i in xrange(self.m):
            for j in xrange(self.n):
                result.setElt(i,j,self.ring.mul(d,self.getElt(i,j)))
        return result

    def det(self):
        (L,U,P)=LUP(self)
        if permMatSign(P)>0:
            total=self.ring.one()
        else:
            total=self.ring.sub(self.ring.zero(),self.ring.one())
        for i in xrange(self.n):
            total=self.ring.mul(U.getElt(i,i),total)
        return total

    def crappyDet(self):
        sum=self.ring.zero()
        for perm in itertools.permutations(range(self.n)):
            product=self.ring.one()
            i=0
            for j in perm:
                product=self.ring.mul(product,self.getElt(i,j))
                i=i+1
            if Algorithms.permSign(perm)<0:
                sum=self.ring.sub(sum,product)
            else:
                sum=self.ring.add(sum,product)
        return sum

    def apply(self,rhs):
        result=Mat(self.ring,self.m,rhs.n)
        for i in xrange(self.m):
            for j in xrange(rhs.n):
                tot=self.ring.zero()
                for k in xrange(self.n):
                    d=self.ring.mul(self.getElt(i,k),rhs.getElt(k,j))
                    tot=self.ring.add(tot,d)
                result.setElt(i,j,tot)
        return result

    def add(self,rhs):
        result = Mat(self.ring,self.m,self.n)
        for i in xrange(self.m):
            for j in xrange(self.n):
                result.setElt(i,j,self.ring.add(self.getElt(i,j),rhs.getElt(i,j)))
        return result

    def sub(self,rhs):
        result = Mat(self.ring,self.m,self.n)
        for i in xrange(self.m):
            for j in xrange(self.n):
                result.setElt(i,j,self.ring.sub(self.getElt(i,j),rhs.getElt(i,j)))
        return result

    def randomize(self):
        for i in xrange(self.m):
            for j in xrange(self.n):
                self.setElt(i,j,self.ring.random())

    def copy(self):
        newMat = Mat(self.ring,self.m,self.n)
        for i in xrange(self.m):
            for j in xrange(self.n):
                newMat.setElt(i,j,self.getElt(i,j))
        return newMat

    def transpose(self):
        for i in xrange(1,self.m):
            for j in xrange(i):
                d=self.getElt(j,i)
                self.setElt(j,i,self.getElt(i,j))
                self.setElt(i,j,d)

    def equal(self,rhs):
        if self.m != rhs.m:
            return False
        if self.n != rhs.n:
            return False
        for i in xrange(self.m):
            for j in xrange(self.n):
                if not self.ring.areEqual(self.getElt(i,j),rhs.getElt(i,j)):
                    return False
        return True

    def invert(self):
        n=self.n
        (L,U,P)=LUP(self)
        for i in xrange(n):
            if self.ring.isZero(U.getElt(i,i)):
                raise Exception("singular")
        I=ident(self.ring,n)
        LInv=Mat(self.ring,n,n)
        UInv=Mat(self.ring,n,n)
        for j in xrange(n):
            x=ViewMat(LInv,0,j,n,1)
            b=ViewMat(I,0,j,n,1)
            solveBotTriangular(L,b,x)
        for j in xrange(n):
            x=ViewMat(UInv,0,j,n,1)
            b=ViewMat(I,0,j,n,1)
            solveTopTriangular(U,b,x)
        return UInv.apply(LInv.apply(P))

    def isZero(self):
        for i in xrange(self.m):
            for j in xrange(self.n):
                if not self.ring.isZero(self.getElt(i,j)):
                    return False
        return True

#Elementary row/column operations:

# A_a <- A_a + dA_b
def colAdd(A,a,b,d):
    ring = A.ring
    m=A.m
    for i in xrange(m):
        e=ring.mul(d,A.getElt(i,b))
        e=ring.add(A.getElt(i,a),e)
        A.setElt(i,a,e)

# A_a <- dA_a
def colMul(A,a,d):
    ring = A.ring
    m=A.m
    for i in xrange(m):
        A.setElt(i,a,ring.mul(d,A.getElt(i,a)))

def colSwap(A,a,b):
    m=A.m
    for i in xrange(m):
        d=A.getElt(i,a)
        A.setElt(i,a,A.getElt(i,b))
        A.setElt(i,b,d)

#Compute the sign of a permutation matrix
def permMatSign(P):
    IntRing=Modular.Integer()
    v=Mat(IntRing,P.n,1)
    for i in xrange(P.n):
        v.setElt(i,0,IntRing.make(i))
    permVec=P.apply(v)
    perm=[]
    for i in xrange(permVec.m):
        perm.append(permVec.getElt(i,0))
    return Algorithms.permSign(perm)

def solveBotTriangular(L,b,x):
    n=L.n
    ring=L.ring
    for i in xrange(n):
        x.setElt(i,0,ring.zero())
    for i in xrange(n):
        rowSum=b.getElt(i,0)
        for j in xrange(i):
            d=ring.mul(L.getElt(i,j),x.getElt(j,0))
            rowSum=ring.sub(rowSum,d)
        d=L.getElt(i,i)
        if not ring.isZero(d):
            x.setElt(i,0,ring.div(rowSum,d))

def solveTopTriangular(U,b,x):
    n=U.n
    ring=U.ring
    for i in xrange(n):
        x.setElt(i,0,ring.zero())
    for i in xrange(n):
        rowSum=b.getElt(n-1-i,0)
        for j in xrange(i):
            d=ring.mul(U.getElt(n-1-i,n-1-j),x.getElt(n-1-j,0))
            rowSum=ring.sub(rowSum,d)
        d=U.getElt(n-1-i,n-1-i)
        if not ring.isZero(d):
            x.setElt(n-1-i,0,ring.div(rowSum,d))

#Access the elements of a submatrix, acts like a new (m-i)X(n-j)
#matrix from the perspective of indexing but all setElt operations
#effect the underlying matrix
class ViewMat(Mat):
    def __init__(self,mat,i,j,m,n):
        self.ring=mat.ring
        self.baseM=mat.m
        self.baseN=mat.n
        self.m=m
        self.n=n
        self.rowOffset=i
        self.colOffset=j
        self.elts=mat.elts
    def getElt(self,i,j):
        return self.elts[(i+self.rowOffset)*self.baseN+j+self.colOffset]
    def setElt(self,i,j,d):
        self.elts[(i+self.rowOffset)*self.baseN+j+self.colOffset]=d

#Return the nXn identity matrix
def ident(field,n):
    A = Mat(field,n,n)
    for i in xrange(n):
        A.setElt(i,i,field.one())
    return A

def readMat(file,ring):
    [m,n,l]=file.readline().split()
    m=int(m)
    n=int(n)
    l=int(l)
    newMat = Mat(ring,m,n)
    for k in xrange(l):
        [i,j,x]=file.readline().split()
        i=int(i)
        j=int(j)
        d=ring.make(x)
        newMat.setElt(i,j,d)
    return newMat

def swap(d,i,j):
    x=d[i]
    d[i]=d[j]
    d[j]=x

#LUP Decomposition of a matrix
def LUP(mat):
    n=mat.n
    ring=mat.ring
    L=ident(ring,n)
    U=Mat(ring,n,n)
    P=Mat(ring,n,n)
    A=mat.copy()
    perm={}
    for i in xrange(n):
        perm[i]=i
    for p in xrange(n-1):
        for j in xrange(p+1,n):
            if (ring.isZero(A.getElt(perm[p],p))):
                swap(perm,j,p)
        if ring.isZero(A.getElt(perm[p],p)):
            continue
        for k in xrange(p+1,n):
            d=ring.div(A.getElt(perm[k],p),A.getElt(perm[p],p))
            A.setElt(perm[k],p,d)
            for i in xrange(p+1,n):
                d=A.getElt(perm[k],i)
                d=ring.sub(d,ring.mul(A.getElt(perm[k],p),A.getElt(perm[p],i)))
                A.setElt(perm[k],i,d)
    for i in xrange(n):
        P.setElt(i,perm[i],ring.one())
    for i in xrange(n):
        for j in xrange(i,n):
            d=A.getElt(perm[i],j)
            U.setElt(i,j,d)
    for i in xrange(1,n):
        for j in xrange(i):
            d=A.getElt(perm[i],j)
            L.setElt(i,j,d)
    return (L,U,P)
