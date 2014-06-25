
#Misc algorithms

def euclid(a,b):
    rOld = a
    rNew = b

    sOld = 1
    sNew = 0

    tOld = 0
    tNew = 1

    while (rNew != 0):
        q=rOld/rNew
        (rNew,rOld)=(rOld-q*rNew,rNew)
        (sNew,sOld)=(sOld-q*sNew,sNew)
        (tNew,tOld)=(tOld-q*tNew,tNew)

    return (rOld,sOld,tOld)

def cycles(x):
    visited={}
    cycles=[]
    for i in x:
        cur=i
        cycle=[]
        while not cur in visited:
            visited[cur]=True
            cycle.append(cur)
            cur=x[cur]
        if cycle:
            cycles.append(cycle)
    return cycles

def permSign(x):
    cycleList=cycles(x)
    sign=1
    for cycle in cycleList:
        if (len(cycle)%2)==0:
            sign=-sign
    return sign
