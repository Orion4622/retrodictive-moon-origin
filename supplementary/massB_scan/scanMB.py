import sys, os, json, time, numpy as np
os.environ["RMAX"]="5.0"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phaseB_band as pb
from phaseB_band import Config, integrate, classify, GM, RE, AROCHE
def make_disc(cfg,n,beta,seed,sig,rout):
    rng=np.random.default_rng(seed); sb=np.sin(np.radians(beta))
    lat=np.arcsin(rng.uniform(-sb,sb,n)); lon=rng.uniform(0,2*np.pi,n)
    p=sig+2.0
    u=rng.uniform(0,1,n); a0,a1=AROCHE**p,rout**p
    r=(a0+u*(a1-a0))**(1/p)*RE
    pos=np.stack([r*np.cos(lat)*np.cos(lon),r*np.cos(lat)*np.sin(lon),r*np.sin(lat)],1)
    vphi=np.sqrt(GM/r**3)*r*np.cos(lat)
    vel=np.stack([-vphi*np.sin(lon),vphi*np.cos(lon),np.zeros(n)],1)
    return pos,vel
def run(P,a,sense,seed,muB,b,sig=1,rout=5.0,N=700,vinf=2.0e3):
    Om=2*np.pi/(P*3600); wt=Om/np.sqrt(GM/(a*RE)**3); pb.J2=0.048*(wt/0.97)**2
    kw=dict(iB_deg=25.,retrograde=True) if sense=="retro25" else dict(iB_deg=25.,retrograde=False)
    cfg=Config(wtilde=wt,a=a,muB=muB,b=b,vinf=vinf,**kw)
    pos,vel=make_disc(cfg,N,26.,seed,sig,rout)
    out=integrate(cfg,pos,vel,dt=25.); c=classify(cfg,out); cat=c["cat"]; lig=cat=="liee"
    eps=c["eps"]; j=c["j"]; aa=-GM/(2*eps[lig]); e=np.sqrt(np.clip(1+2*eps[lig]*j[lig]**2/GM**2,0,None)); peri=aa*(1-e)/RE
    beyond=lig.copy(); beyond[lig]=peri>=AROCHE
    return dict(P=P,a=a,sense=sense,seed=seed,muB=muB,b=b,abs=float(np.mean(cat=="retombee_vol")),lig=float(lig.mean()),B=float(np.mean(cat=="emportee_B")),esc=float(np.mean(cat=="echappee")),
                beyond=float(beyond.mean()),tilt=float(np.median(c["tilt"][lig])) if lig.any() else None)
if __name__=="__main__":
    t=time.time(); r=run(5.,1.52,"retro25",0,0.95,3.0); print(r, f"{time.time()-t:.1f}s")
