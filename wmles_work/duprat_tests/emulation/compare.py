"""Compare kernel outputs against the CPU reference."""
import sys
import numpy as np


def load(p, n=3000):
    rows = [l.split() for l in open(p)]
    comp = np.array([[float(x) for x in r] for r in rows[:2*n]])
    upd = np.array([[float(x) for x in r] for r in rows[2*n:]])
    return comp, upd


ref_c, ref_u = load(sys.argv[1])
inp = np.loadtxt('inputs.txt', skiprows=1)
for p in sys.argv[2:]:
    c, u = load(p)
    tau_r = np.linalg.norm(ref_c[:, :3], axis=1)
    scale = np.maximum(tau_r, 1e-30)
    dtau = np.linalg.norm(c[:, :3] - ref_c[:, :3], axis=1) / scale
    dal = np.abs(c[:, 3] - ref_c[:, 3]) / np.maximum(np.abs(ref_c[:, 3]), 1e-30)
    du = np.abs(u - ref_u) / np.maximum(np.abs(ref_u), 1e-30)
    # reversed nodes: stress along the wall-parallel velocity
    vel = inp[:, 0:3]; nrm = inp[:, 3:6]
    vpar = vel - np.sum(vel*nrm, axis=1)[:, None]*nrm
    rev = np.sum(ref_c[:3000, :3]*vpar, axis=1) > 0
    rev_k = np.sum(c[:3000, :3]*vpar, axis=1) > 0
    print(f'{p}: max rel |dtau| {dtau.max():.2e} (median {np.median(dtau):.1e}), '
          f'max rel dalpha {dal.max():.2e}, max rel d(dpds) {du.max():.2e}, '
          f'reversed nodes ref {rev.sum()} kernel {rev_k.sum()} '
          f'mismatch {(rev != rev_k).sum()}, nan {np.isnan(c).sum() + np.isnan(u).sum()}')
