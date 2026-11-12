"""Cross-check the profiles of ph_stats.py against pySEMTools probes."""
import sys
import numpy as np
from mpi4py import MPI
from pysemtools.datatypes.msh import Mesh
from pysemtools.datatypes.field import FieldRegistry
from pysemtools.io.ppymech.neksuite import pynekread
from pysemtools.interpolation.probes import Probes

fname, csv, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
comm = MPI.COMM_WORLD
msh = Mesh(comm, create_connectivity=False)
fld = FieldRegistry(comm)
pynekread(fname, comm, data_dtype=np.double, msh=msh, fld=fld)
prof = np.loadtxt(csv, delimiter=',', skiprows=1)
# The probes need 3D elements: extrude the 2D elements into one layer of
# thickness 1 with lz = lx GLL planes, the fields constant in z.
lx = msh.x.shape[-1]
zg = 0.5*(np.polynomial.legendre.Legendre.basis(lx - 1).deriv().roots())
zg = np.concatenate([[-0.5], np.sort(zg.real), [0.5]]) + 0.5


def ext(a):
    return np.repeat(a, lx, axis=1)


z3 = np.broadcast_to(zg[None, :, None, None], ext(msh.x).shape).copy()
msh3 = Mesh(comm, x=ext(msh.x), y=ext(msh.y), z=z3, create_connectivity=False)
pts = np.column_stack([prof[:, 0], prof[:, 1], np.full(len(prof), 0.5)])
pr = Probes(comm, probes=pts, msh=msh3, write_coords=False,
            output_fname=outdir + '/interpolated_fields.csv')
reg = {k: ext(v) for k, v in fld.registry.items()}
u, v = reg['u'], reg['v']
upup = reg['s0'] - u**2
upvp = reg['s3'] - u*v
pr.interpolate_from_field_list(0.0, [u, upup, upvp], comm, write_data=False)
res = pr.interpolated_fields[:, 1:]
ok = ~np.isnan(prof[:, 2])
for k, (name, col) in enumerate([('u', 2), ("u'u'", 4), ("u'v'", 7)]):
    d = np.abs(res[ok, k] - prof[ok, col])
    j = np.argmax(d)
    print(f"{name:5s} median {np.median(d):.1e}, worst at x={prof[ok][j, 0]:.2f} y={prof[ok][j, 1]:.4f} (hill {prof[ok][prof[ok][:, 0] == prof[ok][j, 0]][0, 1]:.4f}); 99th pct {np.percentile(d, 99):.1e}")
    print(f"{name:5s} max |pySEMTools - ph_stats| = {d.max():.2e} "
          f"(max |value| {np.abs(prof[ok, col]).max():.3e}), "
          f"points compared {ok.sum()} of {len(ok)}")
