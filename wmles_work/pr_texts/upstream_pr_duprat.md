# Upstream PR: ExtremeFLOW/neko, base `develop`, head `muenschm:feature/duprat-wall-model`

Opened as [ExtremeFLOW/neko#2888](https://github.com/ExtremeFLOW/neko/pull/2888).

The description below follows the repository's pull request template
(`.github/pull_request_template.md`: Changes, Related Issues, Details). There
is no related issue, so that section is left out, as the template allows.

## PR title

Add the Duprat wall model with pressure-gradient effects

## PR description (paste from here)

## Changes

- Added the `duprat` wall model of Duprat et al. (2011), which accounts for the wall-tangential pressure gradient through an extended velocity scale and a pressure-gradient-dependent eddy viscosity. The pressure gradient is either prescribed as a constant or evaluated locally from the pressure field and filtered in time. Available for the CPU, CUDA, HIP, OpenCL and Metal backends.
- Updated the documentation of the `shear_stress` and `wall_model` conditions: they require `full_stress_formulation`, and they are not restricted to axis-aligned boundaries.

## Details

### Model

Duprat et al., "A wall-layer model for large-eddy simulations of turbulent flows with/out pressure gradient", Phys. Fluids 23, 015101 (2011), doi:10.1063/1.3529358.

- **Velocity scale:** the extended velocity scale u_τp = √(u_τ² + u_p²), with the pressure velocity u_p = |(ν/ρ) ∂p/∂s|^{1/3}.
- **Eddy viscosity:** a mixing-length eddy viscosity that depends on α = u_τ²/u_τp².
- **Velocity profile:** U*(y*) is obtained by integrating the thin-boundary-layer equation across the wall layer, using a composite Gauss–Legendre rule in ln(1 + y*).
- **Friction velocity:** u_τ is found with a bracketed, safeguarded Newton method. The residual is not monotone for favourable gradients.

Without a pressure gradient the model reduces to its equilibrium form.

```json
{
  "type": "wall_model",
  "zone_indices": [3, 4],
  "model": "duprat",
  "kappa": 0.41,
  "pressure_gradient": { "type": "local", "filter_time": 1.0 },
  "sampling": { "type": "gll", "value": 3 }
}
```

**Parameters:**

- `kappa` is required.
- `beta` and `A` are optional and default to the calibrated values 0.78 and 17.

**The optional `pressure_gradient` object:**

- `"type": "constant"` prescribes ∂p/∂s with `value`, which defaults to 0.
- `"type": "local"` evaluates ∂p/∂s at every wall node from the pressure field:
  1. the pressure gradient is sampled at the velocity sampling points;
  2. it is projected on the wall-parallel flow direction;
  3. it is limited so that u_p does not exceed the friction velocity of a laminar profile through the sampled velocity;
  4. it is filtered in time with an exponential moving average whose time constant is `filter_time`.

  `start_time` delays the use of the local gradient.

**Output:** ∂p/∂s and α at the wall are stored in the registry fields `duprat_dpds` and `duprat_alpha`, so they can be written or averaged, e.g. with `user_stats`.

### Files

- **Model:** `src/wall_models/duprat.f90`.
- **CPU kernels:** `src/wall_models/bcknd/cpu/duprat_cpu.f90`, for the shear stress and the pressure-gradient update.
- **Device backends:** `src/wall_models/bcknd/device/duprat_device.F90`, plus kernels for CUDA, HIP, OpenCL and Metal.
- **Registration:**
  - factory entry;
  - `Makefile.am`, `.depends`, `.depends_device`;
  - OpenCL program handle in `prgm_lib`.
- **Documentation:** case-file schema, user guide (`case-file/fluid-incompressible.md`) and changelog.
- **Tests:** a serial pFUnit suite, `tests/unit/duprat`.

### `shear_stress` and `wall_model` on curved walls

The user guide and `wall_model_bc.f90` said that both conditions only work with axis-aligned boundaries. In the current code, both require `full_stress_formulation`. The coupled vector BC projector then removes the wall-normal velocity along the resolved node normal, so the conditions apply to inclined and curved walls as well. The wall models themselves use the wall-parallel velocity at each node.

We checked this on the curved wall of the periodic hill:

- u·n stays at round-off level, ≤ 7·10⁻¹⁶ at all nodes, with the averaged node normal.
- The production runs below use the condition on the curved hill wall.

The documentation now says this and links to the page on mixed boundary-condition resolution.

### Testing

- **Unit tests:** the cases compare the CPU kernels with reference values computed independently with SciPy (adaptive quadrature of the model equation, brentq):
  - U*(y*) without, with adverse and with favourable pressure gradient;
  - the friction velocity and α on a flat wall for zero, adverse and favourable gradients;
  - custom constants;
  - an inclined wall with a wall-normal velocity component and ρ ≠ 1;
  - zero wall-parallel velocity;
  - a pressure gradient that alone accounts for the sampled velocity (vanishing stress);
  - the previous stress used as the initial guess;
  - the projection, Stokes limit and moving average of the local pressure-gradient update.

  The composite quadrature agrees with the adaptive reference to ≤ 7·10⁻¹⁰ in the tested range. Tolerances follow `NEKO_EPS`.

  Setup: pFUnit 4.15.0, gfortran 13.3, OpenMPI 4.1.6, and the CI flags `-g -O2 -pedantic -std=f2008 -w`.
  - The suite passes in single and double precision.
  - The full unit suite (`make check`) passes in both precisions, 59 of 59 tests.
  - Perturbing a reference value by a relative 1.5·10⁻⁷ makes the affected tests fail.
- **Build and lint:** `findent` with the CI options leaves the changed Fortran files unchanged. The `flint` score is 10.00 for the new files and unchanged for the modified ones.
- **Schema:** the case-file schema accepts `constant`, `local` and no `pressure_gradient`. It rejects:
  - `local` without `filter_time`;
  - `constant` with `filter_time`;
  - `local` with `value`;
  - unknown types.
- **Smoke run:** a short periodic-hill run on 4 ranks with the local pressure gradient and `user_stats` of `duprat_dpds` and `duprat_alpha` (127 steps). It runs stably and writes the statistics. The model sources are identical to those on our development branch. Against the same run there, the wall forces agree to within about 5·10⁻⁵ relative after the first steps; the small differences come from changes in `develop` since that branch was made and appear already in the first pressure solve.
- **Periodic hill, Re_h = 10595:** wall-modelled LES with 3344 elements at N = 8, local pressure gradient, averaged over t = 150 … 600. The case was run on CPU (Intel, 72 ranks) and on GPU (CUDA, one NVIDIA GPU).
  - Compared with the same case using an equilibrium law (Reichardt), the reattachment point moves downstream from x/h ≈ 2.7–2.8 to 3.2–3.4. The LES database of Gloerfelt & Cinnella (FTC 103, 2019) gives 4.3–4.5.
  - The skin-friction peak on the windward side of the hill (C_f ≈ 0.031) agrees with the reference. With Reichardt it is about 0.011.
  - CPU and GPU agree in the mean flow, the skin friction and the model's α and ∂p/∂s along the wall. There is no systematic difference in the turbulence statistics.

### Notes

- **Local pressure gradient:** in the periodic hill, the filtered ∂p/∂s that the model uses in the recirculation zone and the recovery region is roughly ten times smaller than the gradient of the mean pressure there. On the windward face it is about two thirds. Two causes contribute:
  - the projection on the instantaneous flow direction;
  - the Stokes limit, which is applied before the time filter.

  We plan to improve this in a follow-up, for example by filtering the gradient vector before projecting and limiting it. The interface would stay the same.
- **Body-force-driven flows:** in flows driven by a body force, such as a channel with a prescribed flow rate, the mean pressure gradient is not part of the pressure field. This is documented.
- **GPU status:**
  - CUDA has run the full periodic hill, and the results agree with the CPU run.
  - The OpenCL kernels were checked standalone on PoCL (a CPU OpenCL device) against reference cases, in single and double precision.
  - The HIP kernels are the CUDA kernels apart from comments. HIP and Metal have not run on hardware.
  - The model prints no device warning, like the other wall models.
- **Merge order with the Reichardt PR:** this PR and the Reichardt wall model PR both register a new model in the same lists. Whichever is merged second needs a merge of `develop`, and the length of `WALLM_KNOWN_TYPES` must then become 7.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_01Qr9jBZiogCtaH7KGcZPRSx
