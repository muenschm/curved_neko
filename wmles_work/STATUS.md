# Status of the work

## Upstream pull requests

- Reichardt: [ExtremeFLOW/neko#2884](https://github.com/ExtremeFLOW/neko/pull/2884),
  branch `feature/reichardt-wall-model` (head 7e9dc34).
- Duprat: [ExtremeFLOW/neko#2888](https://github.com/ExtremeFLOW/neko/pull/2888),
  branch `feature/duprat-wall-model` (head f854f9b since 16 Nov 2026: the
  reworked local mode was fast-forwarded into the PR on the user's decision,
  instead of a follow-up PR). Checked before the push on the branch merged
  with the current develop (74c61b0): no conflicts, unit suite 60/60 in dp
  and sp, coarse hill smoke run stable. New PR description and a comment
  for the colleagues in `pr_texts/`, both posted by the user on
  17 Nov 2026. Next: CI result and review comments on #2888 (the session
  cannot read ExtremeFLOW; the user forwards them).
- Review comment of timofeymukha (Gauss-Legendre tables in duprat_cpu.f90,
  "can be taken from the common library"): option A implemented in 7cd115e
  (PR head since 17 Nov 2026): the rule comes from speclib's zwgl, computed
  once in the model and passed to the kernels of all backends; new unit test
  for the rule. The commit also regenerates src/.depends (it lacked gs_ops,
  vector_math since 50b9589; the CI dependency check would have failed).
  Checks: dp/sp 59/59 on the branch, merged with develop: no conflicts,
  make depend reproduces the merged .depends, Duprat test 9/9; device
  kernels re-checked (emulation, PoCL); coarse hill smoke run stable. Reply
  text in `pr_texts/upstream_pr_duprat_reply_gauss_rule.md`, posted by the
  user on 18 Nov 2026.
- Whichever is merged second needs a merge of `develop`, and
  `WALLM_KNOWN_TYPES` then has 7 entries.
- Check on 23 Nov 2026: #2884 is NOT merged yet (develop f4519fe has 5
  entries and no Reichardt files; the Reichardt branch got a develop merge,
  f1e485d). #2888 (7cd115e) merged with develop f4519fe: no conflicts, builds
  in dp/sp, unit suite 60/60 in dp and sp, coarse Duprat hill run stable.
- Prepared and tested (locally, not pushed) the merge of #2884 into
  #2888 + develop (branch test/dup-reich, e0e7fc0). Git merges the factory
  WITHOUT conflict to `WALLM_KNOWN_TYPES(6)` with 7 names, which does not
  compile. Resolution, to be applied by whichever PR is merged second:
  1. `src/wall_models/wall_model_fctry.f90`: `WALLM_KNOWN_TYPES(7)`;
  2. user guide table, `wall_modelling.type`: `rough_log_law`, `spalding`,
     `reichardt`, `duprat` (conflict);
  3. `src/.depends` (conflict) and `.depends_device`: regenerate with
     `./regen.sh && ./configure FC=mpif90 && cd src && make depend`
     (makedepf90 3.0.1 reproduces the CI output exactly);
  4. everything else (Makefile.am, tests, configure.ac, prgm_lib, schema,
     CHANGELOG) merges correctly, each entry once.
  Result: unit suite 61/61 in dp and sp, Duprat 9/9 and Reichardt 5/5,
  coarse hill runs with both models stable.

## In progress: Duprat local-mode improvements

Branch `feature/duprat-local-dpds` in muenschm/neko, based on
`feature/duprat-wall-model`. Since 16 Nov 2026 part of #2888 (the PR branch
was fast-forwarded to f854f9b); no separate follow-up PR.

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
- f854f9b (16 Nov 2026, on the user's request): Stokes limit removed from
  the code entirely; the update kernels only project (no nu, rho, h, flag).
  Duprat test 8/8 in dp and sp, full dp suite 59/59, findent/flint clean,
  all device kernels re-checked by emulation and PoCL (unchanged results).
  Same numerics as bcb0e5a with the limit off, so the production runs
  stay valid.

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
- Production runs of bcb0e5a done (16 Nov 2026): Fritz (CPU, restart at
  t = 489.53) and Alex (GPU, CUDA, first real GPU run of the new kernels),
  both to t = 600 without NaN. Results and the comparison with the old
  Duprat version and the reference in
  `periodic_hill/results_Re10595_duprat_v2/README.md`. In short:
  reattachment 3.43 +- 0.11 / 3.64 +- 0.04 (old 3.37 / 3.21, reference
  4.2-4.3), separation unchanged 0.44; turbulence peaks closer to the
  reference at all stations; bubble C_f now at the reference level; model
  dp/ds now 0.74-1.0 of the mean gradient except the recovery region (0.24);
  hill friction and forcing about 10 % lower; windward C_f max 15 % above
  the reference.
- `periodic_hill.f90` now copies the initial velocity to the device itself
  (needed since Neko #2831; on Alex the job failed with "Initial condition
  did not set u" before).
- `analysis/` rewritten (ph_io, ph_stats, ph_forces, wall_summary,
  wall_dpds_check, duprat_model, gpu_cpu_windows) and checked against
  pySEMTools; the reader bug for files with several scalars is fixed (see
  `analysis/README.md`).
- Report updated (16 Nov 2026): new section "Duprat model closer to the
  paper" with tables and figures (new vs old vs Gloerfelt), corrected cause
  of the weak dp/ds (Stokes limit), summary and open items updated. The old
  Duprat runs were re-evaluated with the new scripts
  (`periodic_hill/results_Re10595_duprat/`); they reproduce the report.
  The PR #2888 text is rewritten (`pr_texts/upstream_pr_duprat.md`).
- Possible model output: a signed alpha (s_tau alpha) would show how often
  the reversed branch is active instantaneously.
