import random
import Algorithms

#Represents GF(p) for prime p
class Modular:
    def __init__(self,q):
        self.q=q

    def add(self,lhs,rhs):
        return (lhs+rhs)%self.q

    def neg(self,x):
        return self.q-x

    def sub(self,lhs,rhs):
        return (self.q+lhs-rhs)%self.q

    def mul(self,lhs,rhs):
        return (lhs*rhs)%self.q

    def make(self,x):
        return int(self.q+x)%self.q

    def inv(self,x):
        if (x==0):
            raise Exception("div by zero")
        (gcd,a,b)=Algorithms.euclid(x,self.q)
        return self.make(a)

    def div(self,lhs,rhs):
        return self.mul(lhs,self.inv(rhs))

    def zero(self):
        return 0

    def isZero(self,x):
        return x == 0

    def one(self):
        return 1

    def isOne(self,x):
        return x == 1

    def random(self):
        return random.randrange(0,self.q)

    def toStr(self,x):
        return str(x)

    def areEqual(self,x,y):
        return x==y

class Integer:
    def add(self,lhs,rhs):
        return lhs+rhs

    def neg(self,x):
        return -x

    def sub(self,lhs,rhs):
        return lhs-rhs

    def mul(self,lhs,rhs):
        return lhs*rhs

    def make(self,x):
        return int(x)

    def inv(self,x):
        raise Exception("Not a divison ring")

    def div(self,lhs,rhs):
        raise Exception("Not a divison ring")

    def zero(self):
        return 0

    def isZero(self,x):
        return x == 0

    def one(self):
        return 1

    def isOne(self,x):
        return x == 1

    def random(self):
        return random.randrange(0,9999999)

    def toStr(self,x):
        return str(x)

    def areEqual(self,x,y):
        return x==y
