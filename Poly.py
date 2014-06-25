#Generates a polynomial ring over an existing ring
class Poly:
    def __init__(self,ring):
        self.ring=ring

    def toStr(self,x):
        strElts=[]
        if len(x)==0 or (len(x)==1 and self.ring.isZero(x[0])):
            return "0"
        for i in xrange(len(x)-1,1,-1):
            if not self.ring.isZero(x[i]):
                if self.ring.isOne(x[i]):
                    strElts.append('x^'+str(i))
                else:
                    strElts.append(str(x[i])+'*x^'+str(i))
        if not (len(x)<2 or self.ring.isZero(x[1])):
            if not self.ring.isOne(x[1]):
                strElts.append(str(x[1])+'*x')
            else:
                strElts.append('x')
        if not (len(x)<1 or self.ring.isZero(x[0])):
            strElts.append(str(x[0]))
        return " + ".join(strElts)

    def degree(self,x):
        return len(x)-1

    def apply(self,poly,resultAdd,indetMul,coeffMul,one,x):
        lastX=x
        total=coeffMul(self.getElt(poly,0),one)
        for i in xrange(1,len(poly)):
            total=resultAdd(total,coeffMul(self.getElt(poly,i),lastX))
            lastX=indetMul(x,lastX)
        return total

    def polyMinus(self,x):
        return x[1:]

    def add(self,lhs,rhs):
        d=max(len(lhs),len(rhs))
        retVal=[self.ring.zero() for i in xrange(d)]
        for i in xrange(d):
            retVal[i]=self.ring.add(self.getElt(lhs,i),self.getElt(rhs,i))
        self.trim(retVal)
        return retVal

    def sub(self,lhs,rhs):
        d=max(len(lhs),len(rhs))
        retVal=[self.ring.zero() for i in xrange(d)]
        for i in xrange(d):
            retVal[i]=self.ring.sub(self.getElt(lhs,i),self.getElt(rhs,i))
        self.trim(retVal)
        return retVal

    def scalarMul(self,x,a):
        return [self.ring.mul(a,c) for c in x]

    def timesXPower(self,x,n):
        zeros=[self.ring.zero() for i in xrange(n)]
        zeros.extend(x)
        return zeros

    def fromScalar(self,n):
        return [n]

    def zero(self):
        return []

    def one(self):
        return [1]

    def isZero(self,x):
        return x==0

    def mul(self,lhs,rhs):
        result = []
        for i in xrange(len(lhs)):
            x=self.scalarMul(rhs,lhs[i])
            result=self.add(result,self.timesXPower(x,i))
        return result

    def getElt(self,x,i):
        if i >= len(x):
            return self.ring.zero()
        else:
            return x[i]

    def trim(self,x):
        for i in xrange(len(x)-1,-1,-1):
            if self.ring.isZero(x[i]):
                x.pop()
            else:
                break

    def setElt(self,i,x):
        if len(x) <= i:
            x.extend([self.ring.zero() for z in xrange(1+i-len(x))])
        self.elts[i]=x
        self.trim()

    def normalize(self,x):
        d0 = x[0]
        for i in xrange(len(x)):
            x[i]=self.ring.div(x[i],d0)

