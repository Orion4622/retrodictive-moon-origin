import json, itertools, sys, os
from multiprocessing import Pool
import scanMB
RB=lambda M: 0.95*(M/0.95)**0.27
cases=[(5.,1.52),(6.,1.70)]
Ms=[0.1,0.2,0.3,0.5,0.7,0.95]
bs=[1.8,2.0,2.2,2.5,2.75,3.0,3.5]
jobs=[]
for (P,a),M,b,s in itertools.product(cases,Ms,bs,(0,1,2)):
    if b < (a+RB(M)):   # contact avec la proto-Terre : exclu
        continue
    jobs.append((P,a,"retro25",s,M,b))
def f(j):
    P,a,sense,s,M,b=j
    try: return scanMB.run(P,a,sense,s,M,b)
    except Exception as e: return dict(P=P,a=a,sense=sense,seed=s,muB=M,b=b,err=str(e))
if __name__=="__main__":
    print(len(jobs),"jobs",flush=True)
    with Pool(int(os.environ.get("NP","4"))) as p:
        res=p.map(f,jobs,chunksize=2)
    json.dump(res,open("grid_MB.json","w"))
    print("done",flush=True)
