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

Done (commit 50b9589, not yet built or tested):
- (A) The pressure gradient is averaged across element boundaries
  (gather-scatter and `coef%mult`) before sampling.
- (B) The sampled gradient vector is filtered in time first. The filtered
  vector is then projected on the current flow direction and limited
  (update kernel with eps = 1). The device kernels are unchanged.
- User guide updated.

Open:
- Build environment after the container reset: the json-fortran and pFUnit
  source builds need the user's approval.
- Build at the CI flags in sp and dp, unit tests, a smoke run, and a short
  hill test of dp/ds against the mean pressure gradient.
- (C) Driving body force: `flow_rate_force` keeps its scale only as a local
  variable in `fluid_volflow`. Options are a core change or a constant
  `body_force` in the wall-model JSON. To be decided with the user.
- Rewrite `ph_stats.py` and `ph_forces.py` for the next evaluation.
