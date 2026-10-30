# Tests of the Duprat kernels

Checks for the branch `feature/duprat-local-dpds` of muenschm/neko: the
reversed wall shear stress, sign(tau_w) in Eq. (5) of Duprat et al. (2011),
and the optional `stokes_limit`.

- `duprat_reference.py`: SciPy reference values for
  `tests/unit/duprat/test_duprat.pf` (adaptive quadrature of U*, brentq for
  u_tau, branch selection as in the kernels). Run with `python3 -I`.
- `branch_scan.py`: tau along the flow as a function of the sampled velocity
  for hill-like parameters. It shows that the stress passes zero continuously
  at the separation velocity u_sep, and that each branch has a single root.
- `emulation/`: the GPU kernels on the host, compared with the CPU kernel on
  3000 random wall nodes (FPG, APG, reversed, no flow direction).
  - `gen_inputs.py` writes `inputs.txt`.
  - `drv_cpu.f90`: CPU kernel, linked against libneko (dp and sp builds).
  - `emu_cuda.cpp`: CUDA and HIP kernel headers compiled with g++
    (`-DKERNEL='"<path>/duprat_kernel.h"'`), double and float.
  - `emu_metal.cpp` with `fakeinc/metal_stdlib`: the Metal kernels, float.
  - `ocl_host.c`: the OpenCL kernels on PoCL; `real` replaced by
    `double`/`float` with sed, as `cltostring.sh` does.
  - `compare.py REF OUT...`: largest relative differences and the number of
    nodes on the reversed branch.

Results on 30 Oct 2026 (relative to the CPU kernel of the same precision):

| backend | dp tau | sp tau | reversed nodes (of 3000) |
|---|---|---|---|
| CUDA (emulated) | 3.1e-8 | 7.7e-5 | 547, identical |
| HIP (emulated) | 3.1e-8 | 7.7e-5 | 547, identical |
| Metal (emulated, float) | - | 7.7e-5 | 547, identical |
| OpenCL (PoCL) | 3.1e-8 | 2.0e-4 | 547, identical |

The dp differences are at the solver tolerance (1e-8). No NaN.
