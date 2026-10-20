import json, itertools, os
from multiprocessing import Pool
import scanMB
jobs=[(5.,1.52,sense,s,M,b) for sense in ("retro25","pro25") for M in (0.1,0.3,0.5,0.7,0.95) for b in (2.5,3.0) for s in (0,1,2) if b>=1.52+0.95*(M/0.95)**0.27]
def f(j):
    P,a,sense,s,M,b=j
    return scanMB.run(P,a,sense,s,M,b)
if __name__=="__main__":
    with Pool(2) as p: res=p.map(f,jobs,chunksize=2)
    json.dump(res,open("grid2.json","w")); print("ok")
