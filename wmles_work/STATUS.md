# Status of the work

## Upstream pull requests

- Reichardt: [ExtremeFLOW/neko#2884](https://github.com/ExtremeFLOW/neko/pull/2884),
  branch `feature/reichardt-wall-model` (head 7e9dc34).
- Duprat: [ExtremeFLOW/neko#2888](https://github.com/ExtremeFLOW/neko/pull/2888),
  branch `feature/duprat-wall-model` (head 479331c).
- Whichever is merged second needs a merge of `develop`, and
  `WALLM_KNOWN_TYPES` then has 7 entries.

## In progress: Duprat local-mode improvements

Branch `feature/duprat-local-dpds` in muenschm/neko, based on
`feature/duprat-wall-model`. It becomes a follow-up PR once #2888 is merged.

Commits:
- 50b9589: (A) the pressure gradient is averaged across element boundaries
  (gather-scatter and `coef%mult`) before sampling; (B) the sampled gradient
  vector is filtered in time, then projected on the current flow direction.
- bcb0e5a: closer to the paper (agreed with the user on 30 Oct 2026):
  - reversed wall shear stress, sign(tau_w) in Eq. (5): where the sampled
    velocity is below the separation profile u_p U*(h u_p / nu, alpha = 0)
    under an adverse gradient, the stress points along the sampled velocity
    (before: zero stress);
  - the Stokes bound is optional, `pressure_gradient.stokes_limit`, default
    false (it is not in the paper; it came from the student's code);
  - the update kernels no longer filter (`eps` replaced by the flag);
  - all backends, unit tests with SciPy references, user guide, schema.

Tests of bcb0e5a (see `duprat_tests/README.md`):
- dp and sp builds with the CI flags; Duprat unit test 8/8 in dp and sp;
  full dp suite 59/59 (needs `LD_LIBRARY_PATH` with json-fortran).
- findent clean, flint 10.00 for all changed Fortran files.
- CUDA, HIP and Metal kernels emulated on the host, OpenCL on PoCL, all
  against the CPU kernel on 3000 random nodes: dp within 3e-8 (solver
  tolerance), sp within 2e-4, identical branch selection.
- Each branch has a single root in the scanned hill-like cases, and the
  stress passes zero continuously at the separation velocity
  (`duprat_tests/branch_scan.py`).

Coarse hill (22x8x8 elements, N = 5, statistics t = 20...40,
`analysis/wall_summary.py`, `analysis/wall_dpds_check.py`):

| run | sep | reatt | mean abs(tau_w), hill | alpha spike ratio |
|---|---|---|---|---|
| old (`feature/duprat-wall-model`) | 1.07 | 2.88 | 0.00166 | 1.7 |
| new, 50b9589 (A+B, Stokes bound) | 1.06 | 3.72 | 0.00194 | 0.6 |
| A+B without Stokes bound | 0.88 | 4.48 | 0.00330 | 1.1 |
| bcb0e5a (paper: no bound, reversed stress) | 0.79 | 5.06 | 0.00368 | 0.8 |

References: Breuer et al. 0.19 / 4.69; Duprat et al. grid A 0.38 / 4.65,
grid B 0.45 / 3.92, grid C 0.68 / 2.97. The bcb0e5a run is stable (no NaN).
Ratio of the model <|dp/ds|> to the gradient of the mean pressure:

| run | crest | bubble | recovery | windward |
|---|---|---|---|---|
| old | 0.45 | 0.21 | 0.06 | 0.33 |
| new, 50b9589 | 0.52 | 0.21 | 0.07 | 0.35 |
| A+B without Stokes bound | 0.94 | 0.56 | 0.35 | 1.04 |
| bcb0e5a | 0.88 | 0.56 | 0.63 | 1.06 |

Only one coarse run per variant with 20 time units of averaging; the
differences of about 0.1 in sep/reatt are within the statistical
uncertainty. The trend: without the bound and with the reversed stress the
bubble grows, separation moves upstream, reattachment downstream.

The Stokes bound, not the projection, was the main reason the model saw too
little pressure gradient. Earlier statements ("about a tenth because of the
projection") are to be corrected in the report and the #2888 text.

Open:
- Production hill run of bcb0e5a by the user (`periodic_hill_duprat.case`
  needs no change; `stokes_limit` defaults to false).
- (C) Driving body force: `flow_rate_force` keeps its scale only as a local
  variable in `fluid_volflow`. Options are a core change or a constant
  `body_force` in the wall-model JSON. To be decided with the user.
- Rewrite `ph_stats.py` and `ph_forces.py` for the next evaluation.
