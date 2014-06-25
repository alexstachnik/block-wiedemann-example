
import sys

from Algorithms import *
from Mat import *
from Modular import *
from Poly import *
from ScalarMassey import *
from MatWiedemann import *
from ScalarWiedemann import *
from PolyAlg import *

#Generates a random non-singular matrix and finds its min-poly by
#block-wiedemann and by scalar-wiedemann

field = Modular(101)
#field=Modular(536870909)
ring = Poly.Poly(field)

n=20
A=Mat.Mat(field,n,n)
A.randomize()
while A.det() == 0:
    A.randomize()
print A.det()
b=5
U=Mat.Mat(field,b,n)
V=Mat.Mat(field,n,b)
U.randomize()
V.randomize()

print A.sage()
print "-"
f=fixedMatWiedemannWrapper(A,b,U,V)
f.reverse()

matPoly=polyMatToMatPoly(f,field)
det=matPoly.crappyDet()
det.reverse()
ring.normalize(det)
print det

print ""
otherPoly=wiedemann2Wrapper(A,field)
otherPoly.reverse()
print otherPoly
