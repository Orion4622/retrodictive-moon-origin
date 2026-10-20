"""
Phase B avec reservoir en BANDE de basses latitudes (+-beta de part et d'autre de l'equateur), Terre OBLATE
(equipotentielle de Roche, Rpol/Req ~ 0.6 ; import os
J2 = 0.048
RMAX = float(os.environ.get("RMAX", str(AROCHE)))   # bord externe du reservoir (R_terre) ; defaut = aRoche dans le champ), sans ligne equatoriale.
Reutilise Config / B_state / tidal_field / accel_B / surface_radius de stripping.py (lecture seule).
Masses egales ; lat ~ cos(lat) (couche spherique, uniforme en sin lat) ; r uniforme de rs(lat) a aRoche (Sigma ~ 1/r).
vmode="hydro"  : vitesse azimutale = Omega_K(r_sph) * R_cyl  (rotation d'ensemble, soutien vertical par la pression)
vmode="kepler" : vitesse azimutale = sqrt(GM/r_sph)  (orbites circulaires inclinees de lat)
"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stripping import Config, B_state, tidal_field, accel_B, surface_radius, G, ME, RE, GM, AROCHE

import os
J2 = 0.048
RMAX = float(os.environ.get("RMAX", str(AROCHE)))   # bord externe du reservoir (R_terre) ; defaut = aRoche
LAT = np.radians(np.linspace(-90, 90, 361))

def surf_table(cfg): return surface_radius(cfg, LAT)

def make_band(cfg, n=1000, beta_deg=26., seed=0, vmode="hydro"):
    rng = np.random.default_rng(seed)
    sb = np.sin(np.radians(beta_deg))
    lat = np.arcsin(rng.uniform(-sb, sb, n)); lon = rng.uniform(0, 2*np.pi, n)
    rs = np.interp(lat, LAT, surf_table(cfg))
    r = rs + (RMAX*RE - rs) * rng.uniform(0, 1, n)
    pos = np.stack([r*np.cos(lat)*np.cos(lon), r*np.cos(lat)*np.sin(lon), r*np.sin(lat)], 1)
    Rc = r*np.cos(lat)
    vphi = np.sqrt(GM/r**3)*Rc if vmode == "hydro" else np.sqrt(GM/r)
    vel = np.stack([-vphi*np.sin(lon), vphi*np.cos(lon), np.zeros(n)], 1)
    return pos, vel, lat

def a_J2(cfg, r):
    rn = np.linalg.norm(r, axis=1, keepdims=True); z2 = (r[:, 2:3]/rn)**2
    f = -1.5*J2*GM*cfg.Req**2/rn**5
    a = f*r*(1 - 5*z2); a[:, 2:3] = (f*r[:, 2:3]*(3 - 5*z2))
    return a

def integrate(cfg, pos, vel, dt=20., r_stop=8., useJ2=True):
    tab = surf_table(cfg)
    rB, vB = B_state(cfg, r_stop)
    sB = float(np.sign(np.dot(np.cross(rB, vB), [0, 0, 1])))
    r, v = pos.copy(), vel.copy(); alive = np.ones(len(r), bool); fell = np.zeros(len(r), bool)
    def acc(rr, rBc):
        a = -GM*rr/np.linalg.norm(rr, axis=1, keepdims=True)**3 + tidal_field(cfg, rr, rBc)
        return a + a_J2(cfg, rr) if useJ2 else a
    past = False; n = 0
    while True:
        a = alive; rr, vv = r[a], v[a]
        k1r, k1v = vv, acc(rr, rB)
        rB2 = rB + .5*dt*vB
        k2r, k2v = vv+.5*dt*k1v, acc(rr+.5*dt*k1r, rB2)
        k3r, k3v = vv+.5*dt*k2v, acc(rr+.5*dt*k2r, rB2)
        vB2 = vB + .5*dt*accel_B(cfg, rB); rB3 = rB + dt*vB2
        k4r, k4v = vv+dt*k3v, acc(rr+dt*k3r, rB3)
        r[a] = rr + dt/6*(k1r+2*k2r+2*k3r+k4r); v[a] = vv + dt/6*(k1v+2*k2v+2*k3v+k4v)
        k1r_, k1v_ = vB, accel_B(cfg, rB); k2r_, k2v_ = vB+.5*dt*k1v_, accel_B(cfg, rB+.5*dt*k1r_)
        k3r_, k3v_ = vB+.5*dt*k2v_, accel_B(cfg, rB+.5*dt*k2r_); k4r_, k4v_ = vB+dt*k3v_, accel_B(cfg, rB+dt*k3r_)
        rB = rB + dt/6*(k1r_+2*k2r_+2*k3r_+k4r_); vB = vB + dt/6*(k1v_+2*k2v_+2*k3v_+k4v_)
        n += 1
        idx = np.where(a)[0]; rn = np.linalg.norm(r[idx], axis=1)
        lat = np.arcsin(np.clip(r[idx, 2]/rn, -1, 1))
        hit = rn < np.interp(lat, LAT, tab)          # surface OBLATE = magma ocean : absorbe
        alive[idx[hit]] = False; fell[idx[hit]] = True
        if not past and np.dot(rB, vB) > 0: past = True
        if past and np.linalg.norm(rB) > r_stop*cfg.qB: break
        if n > 400000: break
    return dict(r=r, v=v, alive=alive, fell=fell, rB=rB, vB=vB, sB=sB, nstep=n)

def classify(cfg, out):
    r, v = out["r"], out["v"]; rn = np.linalg.norm(r, axis=1)
    eps = 0.5*np.sum(v**2, 1) - GM/rn
    jv = np.cross(r, v); j = np.linalg.norm(jv, axis=1)
    dr, dv = r-out["rB"], v-out["vB"]
    epsB = 0.5*np.sum(dv**2, 1) - G*cfg.MB/np.linalg.norm(dr, axis=1)
    cat = np.empty(len(r), dtype=object)
    bound = (eps < 0) & out["alive"]
    cat[:] = "echappee"
    cat[bound] = "liee"
    cat[(epsB < 0) & (eps >= 0) & out["alive"]] = "emportee_B"
    cat[out["fell"]] = "retombee_vol"
    rcz = jv[:, 2]**2/GM/RE
    tilt = np.degrees(np.arccos(np.clip(jv[:, 2]/np.maximum(j, 1e-30), -1, 1)))
    return dict(cat=cat, eps=eps, j=j, jz=jv[:, 2], rcz=rcz, tilt=tilt, rc=j**2/GM/RE)

if __name__ == "__main__":
    import itertools
    for case, kw in [("retro25", dict(iB_deg=25., retrograde=True)), ("retro0", dict(iB_deg=0., retrograde=True)),
                     ("pro25", dict(iB_deg=25.))]:
        cfg = Config(**kw)
        for beta, vm in itertools.product((5., 15., 26.), ("hydro", "kepler")):
            pos, vel, lat = make_band(cfg, 1000, beta, 0, vm)
            out = integrate(cfg, pos, vel); c = classify(cfg, out)
            N = len(pos); cat = c["cat"]
            lig = cat == "liee"
            win = lig & (c["rcz"] >= AROCHE)
            print(f"{case:8s} beta={beta:4.0f} {vm:6s} sB={out['sB']:+.0f} retombee_vol={np.mean(cat=='retombee_vol'):.3f} "
                  f"liee={lig.mean():.3f} (jz>0: {np.mean(lig&(c['jz']>0)):.3f}) emportee_B={np.mean(cat=='emportee_B'):.3f} "
                  f"echappee={np.mean(cat=='echappee'):.3f} fenetre_jz={win.mean():.3f} "
                  f"tilt_med_liees={np.median(c['tilt'][lig]) if lig.any() else float('nan'):.1f}", flush=True)
