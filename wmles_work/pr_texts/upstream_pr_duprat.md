# Upstream PR: ExtremeFLOW/neko, base `develop`, head `muenschm:feature/duprat-wall-model`

Opened as [ExtremeFLOW/neko#2888](https://github.com/ExtremeFLOW/neko/pull/2888).
Revised on 16 Nov 2026 after the rework of the local pressure gradient
(commits 50b9589, bcb0e5a, f854f9b). The previous version of this text is in
the git history of this file.

The description follows the repository's pull request template
(`.github/pull_request_template.md`: Changes, Related Issues, Details). There
is no related issue, so that section is left out, as the template allows.

## PR title

Add the Duprat wall model with pressure-gradient effects

## PR description (paste from here)

## Changes

- Added the `duprat` wall model of Duprat et al. (2011), which accounts for the wall-tangential pressure gradient through an extended velocity scale and a pressure-gradient-dependent eddy viscosity, including reversed near-wall flow under an adverse gradient. The pressure gradient is either prescribed as a constant or evaluated locally from the pressure field and filtered in time. Available for the CPU, CUDA, HIP, OpenCL and Metal backends.
- Updated the documentation of the `shear_stress` and `wall_model` conditions: they require `full_stress_formulation`, and they are not restricted to axis-aligned boundaries.

## Details

### Model

Duprat et al., "A wall-layer model for large-eddy simulations of turbulent flows with/out pressure gradient", Phys. Fluids 23, 015101 (2011), doi:10.1063/1.3529358.

- **Velocity scale:** the extended velocity scale u_τp = √(u_τ² + u_p²), with the pressure velocity u_p = |(ν/ρ) ∂p/∂s|^{1/3}.
- **Eddy viscosity:** a mixing-length eddy viscosity that depends on α = u_τ²/u_τp².
- **Velocity profile:** U*(y*) is obtained by integrating the thin-boundary-layer equation, dU*/dy* = (s_p (1−α)^{3/2} y* + s_τ α)/(1 + ν_t/ν), across the wall layer, using a composite Gauss–Legendre rule in ln(1 + y*) (12 sub-intervals with 5 points each; the rule is computed once with `zwgl` from `speclib` and passed to the kernels of all backends). s_p and s_τ are the signs of the pressure gradient and of the wall shear stress relative to the flow at the sampling point.
- **Reversed wall shear stress (s_τ = −1, Eq. 5 of the paper):** under an adverse gradient, the separation profile (α = 0) gives the velocity u_sep = u_p U*(h u_p/ν) at the sampling point. If the sampled velocity is smaller, the near-wall flow is reversed and the stress points along the sampled velocity. The stress passes zero continuously at u = u_sep.
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
  1. ∇p is computed element by element and averaged across element faces (gather-scatter), since the face-normal derivative is two-valued in the C⁰ basis;
  2. the gradient vector is sampled at the velocity sampling points;
  3. it is filtered in time with an exponential moving average whose time constant is `filter_time`;
  4. the filtered vector is projected on the current wall-parallel flow direction.

  `start_time` delays the use of the local gradient.

**Output:** ∂p/∂s and α at the wall are stored in the registry fields `duprat_dpds` and `duprat_alpha`, so they can be written or averaged, e.g. with `user_stats`.

### Files

- **Model:** `src/wall_models/duprat.f90`.
- **CPU kernels:** `src/wall_models/bcknd/cpu/duprat_cpu.f90`, for the shear stress and the projection of the pressure gradient.
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
  - U*(y*) without, with adverse and with favourable pressure gradient, and with reversed wall shear stress;
  - the friction velocity and α on a flat wall for zero, adverse and favourable gradients;
  - custom constants;
  - an inclined wall with a wall-normal velocity component and ρ ≠ 1;
  - zero wall-parallel velocity;
  - reversed wall shear stress: a strongly reversed case (direction and density scaling), and two cases just below and above u_sep, where the stress is small and changes sign;
  - the previous stress used as the initial guess, on both branches;
  - the projection of the filtered pressure gradient on the flow direction;
  - the Gauss–Legendre rule from `speclib` against the 5-point rule in closed form.

  The composite quadrature agrees with the adaptive reference to ≤ 7·10⁻¹⁰ in the tested range. Tolerances follow `NEKO_EPS`.

  Setup: pFUnit 4.15.0, gfortran 13.3, OpenMPI 4.1.6, and the CI flags `-g -O2 -pedantic -std=f2008 -w`, on this branch merged with the current `develop`:
  - The suite passes in single and double precision.
  - The full unit suite (`make check`) passes in both precisions.
- **Device kernels:** besides the CUDA production run below, the CUDA, HIP and Metal kernels were compiled on the host with g++ and the OpenCL kernels run on PoCL. On 3000 random wall nodes (favourable, adverse and reversed cases) they agree with the CPU kernel to 3·10⁻⁸ in double precision (the solver tolerance) and to 2·10⁻⁴ in single precision, with the same branch at every node.
- **Build and lint:** `findent` with the CI options leaves the changed Fortran files unchanged. The `flint` score is 10.00 for the new files and unchanged for the modified ones.
- **Schema:** the case-file schema accepts `constant`, `local` and no `pressure_gradient`. It rejects:
  - `local` without `filter_time`;
  - `constant` with `filter_time`;
  - `local` with `value`;
  - unknown types.
- **Smoke run:** a short periodic-hill run with the local pressure gradient on this branch merged with the current `develop`. It runs stably.
- **Periodic hill, Re_h = 10595:** wall-modelled LES with 3344 elements at N = 8, local pressure gradient, averaged over t = 150 … 600, on CPU (Intel, 72 ranks) and on GPU (CUDA, one NVIDIA GPU). Reference: the LES database of Gloerfelt & Cinnella (FTC 103, 2019), 512×256×256.
  - The model sees most of the pressure gradient: averaged over each region, the filtered ∂p/∂s is 85 % of the gradient of the mean pressure at the crest, 75–80 % in the recirculation zone, 100 % on the windward face, and 25 % in the recovery region, where the flow direction alternates.
  - The skin friction in the recirculation zone, reconstructed from the mean flow with the model, reaches the reference level (C_f ≈ −0.004).
  - The peaks of the Reynolds stresses are close to the reference; at x/h = 2 … 4 they agree within about 10 %.
  - The flow separates at x/h ≈ 0.44 and reattaches at 3.4–3.6. The reference gives −0.12 and 4.3; an equilibrium law (Reichardt) on the same case gives 0.55 and 2.7–2.8.
  - The skin-friction peak on the windward side (C_f ≈ 0.036–0.038) is about 10–15 % above the reference (0.033).
  - CPU and GPU agree in the mean flow, the separation point and the skin-friction distribution. The turbulence statistics differ by a few per cent in both directions between runs; we see no systematic difference between the code paths.

### Notes

- **Differences from the original implementation** (the code of the master's thesis this port started from):
  - The local gradient is averaged across element faces, and the gradient *vector* is filtered before the projection. Before, the projections on the instantaneous flow direction were filtered, which cancels where the flow direction changes often.
  - The original code limited ∂p/∂s so that u_p does not exceed the friction velocity of a laminar profile through the sampled velocity. This limit is not part of the paper. It turned out to cap ∂p/∂s along most of the hill wall, and the runs are stable without it, so it is not included.
  - The reversed wall shear stress (the sign(τ_w) term of Eq. 5) is included; before, the stress was set to zero where the positive branch had no solution.

  On the periodic hill, these changes bring the Reynolds stresses and the skin friction in the recirculation zone close to the reference. Compared with a version without them, the turbulence peaks are up to 17 % lower (typically 10–15 %), the reattachment moves about 0.25 downstream on average, and the separation point does not change.
- **Body-force-driven flows:** in flows driven by a body force, such as a channel with a prescribed flow rate, the mean pressure gradient is not part of the pressure field. This is documented.
- **GPU status:**
  - CUDA has run the full periodic hill, and the results agree with the CPU run.
  - The OpenCL kernels run on PoCL; the HIP and Metal kernels were checked by compiling them on the host. HIP and Metal have not run on hardware. The HIP kernels are the CUDA kernels apart from comments.
  - The model prints no device warning, like the other wall models.
- **Merge order with the Reichardt PR:** this PR and the Reichardt wall model PR both register a new model in the same lists. Whichever is merged second needs a merge of `develop`, and the length of `WALLM_KNOWN_TYPES` must then become 7.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_01Qr9jBZiogCtaH7KGcZPRSx
