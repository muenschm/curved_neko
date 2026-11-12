"""Cross-check read_fld and ph_stats against pySEMTools on one stats file."""
import sys
import numpy as np
from mpi4py import MPI
from pysemtools.datatypes.msh import Mesh
from pysemtools.datatypes.field import FieldRegistry
from pysemtools.io.ppymech.neksuite import pynekread

sys.path.insert(0, sys.argv[2])
from wall_dpds_check import read_fld   # noqa: E402

fname = sys.argv[1]
comm = MPI.COMM_WORLD
msh = Mesh(comm, create_connectivity=False)
fld = FieldRegistry(comm)
pynekread(fname, comm, data_dtype=np.double, msh=msh, fld=fld)
print('pySEMTools keys:', list(fld.registry.keys()), 'x shape', msh.x.shape)
r = read_fld(fname)
# pySEMTools arrays are (nel, lz, ly, lx); ours (nel, ly, lx)
def sq(a):
    return a.reshape(a.shape[0], a.shape[-2], a.shape[-1])
pairs = [('x', sq(msh.x), r['X'][:, 0]), ('y', sq(msh.y), r['X'][:, 1]),
         ('u', sq(fld.registry['u']), r['U'][:, 0]),
         ('v', sq(fld.registry['v']), r['U'][:, 1]),
         ('p', sq(fld.registry['p']), r['P']),
         ('t (= pp)', sq(fld.registry['t']), r['T'])]
for k in range(r['S'].shape[1]):
    pairs.append((f's{k} ', sq(fld.registry[f's{k}']), r['S'][:, k]))
for name, a, b in pairs:
    print(f'{name:9s} max |pySEMTools - read_fld| = {np.max(np.abs(a - b)):.2e}'
          f'  (max |value| {np.max(np.abs(a)):.3e})')
