# `gpu/int-D`: Group D integration (plasma, Raman, χ⁽²⁾)

Base: `gpu/14-raman` at `17c2e9b1`, which already contains the approved `gpu/13-plasma`
(`4a397683`); both descend from `gpu/12-response-traits` (`7d72431c`).
Merged in: `gpu/15-chi2` at `385c6bb9`.

After this branch every nonlinear response Luna ships has a device kernel. What is still
host-only is the *transforms*: radial, free-space and multimode.

## Commits

| commit | what |
| --- | --- |
| `3054b686` | Merge `gpu/15-chi2` into `gpu/int-D` (the merge alone, with the conflict resolutions below) |
| `43d9d63c` | Semantic sweep: dead imports, a stale docstring, the restored `nameof` assertion |
| `6cc8a1db` | Regression gate: argon plasma cases that actually ionise, and the tolerances for them |
| this file | `PR_int-D.md` |

## The merge and its conflicts

Seven conflicting regions in five files. `src/Nonlinear.jl` merged without a conflict:
`gpu/15` rewrites the χ⁽²⁾ section, `gpu/13` the plasma section and `gpu/14` the Raman
section, and in the file they are in that order with the section banners each branch added
interleaving cleanly.

| file | conflict | resolution |
| --- | --- | --- |
| `docs/src/developer/device_model.md` 135-136 | both sides rewrite adjacent rows of the kind table | take `gpu/15`'s `VectorPointwise` row (adds `Chi2Field`, `Chi2Env`) and `gpu/14`'s `Batched` row (adds `PlasmaCumtrapz`, `RamanPolarField`/`RamanPolarEnv`, `KerrFieldNoTHG`) |
| `docs/src/gpu.md` "What runs where" | one side says the χ⁽²⁾ responses are host-only, the other says plasma and Raman are | neither is true after the merge. Rewritten: the radial, free-space and multimode *transforms* are what is host-only; every response Luna ships has a kernel. `gpu/15`'s point that a χ⁽²⁾ propagation still runs on the host as a whole, because its transform does, is kept |
| `docs/src/gpu.md` refusal sentence | "(the χ⁽²⁾ responses, and anything you wrote yourself)" vs "(plasma and Raman at the moment)" | rewritten: the only response with no device kernel is one the user wrote |
| `src/Interface.jl` 373-377 (the `device` keyword docstring) | the device-capable list | rewritten: Kerr (including the no-THG form), plasma, Raman and χ⁽²⁾; false for anything user-written |
| `src/Interface.jl` 519-522 (the same list in `prop_capillary_args`) | same | same |
| `test/test_device.jl` (after the `KerrEnvTHG` testset) | both sides insert testsets at the same point | concatenated; no shared names (`plasmafield`/`ramancase` vs the χ⁽²⁾ locals) |
| `test/test_metal.jl` (helpers after `usercubic`) | `metal_plasmafield`/`metal_adkrate`/`metal_tablerate` vs `userchi2` | concatenated; independent helpers |

Checked after the merge, since three branches touched the same machinery:

- one definition each of `Nonlinear.rescale_responses` and `Nonlinear.resident_arrays_all`
  (`src/Nonlinear.jl:347, 361`), and five call sites in `src/NonlinearRHS.jl` (510, 711,
  1039, 1385, 1465), all of the form
  `rescale_responses(Tuple(resp), spec, scaling, E)` with `spec`/`scaling` keywords on the
  transform constructors;
- no duplicated helper, testset name or docstring anywhere in the merged test files;
- every merged Julia file parses (139 files under `src/`, `test/` and `examples/`).

## Semantic sweep (`43d9d63c`)

Nothing here changes what any propagation computes; the gate against the merge commit is
exactly zero on all 460 rows.

- `src/Nonlinear.jl`: removed `import FFTW`, `ldiv!` and `MArray`. `gpu/14` took the last
  FFTW user out of the file (the Raman response plans through `Utils.plan_ft`); `ldiv!` and
  `MArray` were unused before Group D. `mul!` is kept — the analytic signal and the Raman
  response use it.
- `src/Nonlinear.jl`: the `device_capable` docstring gave plasma and Raman as its examples
  of a response without a kernel. It now says a response the user wrote, which is the only
  one left.
- `test/test_interface.jl`: restored the `gpu/12-response-traits` guarantee that the
  refusal message names the response with `nameof(typeof(r))` and not with its full type.
  `gpu/13`'s review round 1 repointed that test at a user closure, whose `nameof` is a
  gensym, so the assertion had to be dropped. The fixture is now a `Columnwise` struct,
  `LongParameterResponse{A,B,C,D}`, whose instantiated type is 91 characters; the test
  asserts the name is in the message and that `NTuple` and `ComplexF64` are not.
- Consistency checks, no change needed: `Nonlinear.kind`/`device_capable` come out
  `Pointwise`/true for `Kerr_field`, `Kerr_env`, `Kerr_env_thg`, `Batched`/true for
  `Kerr_field_nothg`, `PlasmaCumtrapz` and `RamanPolarField`/`RamanPolarEnv`,
  `VectorPointwise`/true for `Chi2Field` and `Chi2Env`, `Columnwise`/false for a user
  closure. `Ionisation.device_capable` is true for `IonRateADK` and a uniformly spaced
  `IonRatePPTAccel`, false for the direct `IonRatePPT`.
- `docs/src/gpu.md` "What runs where" now lists the same set as the developer guide and as
  `Interface.jl`.
- `include("docs/make.jl")` reports the 9 unresolved `@ref`s the base branches report
  (`LinearOps.βz` ×3, `loadFFTwisdom`, `saveFFTwisdom`, `AbstractOutput`,
  `PhysData.crystal_internal_angle`, `norm_free`, `LinearOps.make_const_linop`) and none
  from Group D.

## The regression cases (`6cc8a1db`)

The four plasma cases were helium at 1 bar with 800 nJ, which does not ionise at all:
review round 1 of `gpu/13-plasma` read `stats/electrondensity` out of the baseline files
and found it exactly `0.0` in all four cases and both modes. The gate was therefore blind
to every change to the plasma response — the `0.000e+00` those branches reported on those
rows meant nothing. They are now argon at 0.1 bar, at the energies the review measured:

| case | was | is | peak ionised fraction |
| --- | --- | --- | ---: |
| `modeavg_field_plasma` | `:He`, 1 bar, 800 nJ | `:Ar`, 0.1 bar, 300 µJ | 0.78 % |
| `modeavg_field_adk` | `:He`, 1 bar, 800 nJ | `:Ar`, 0.1 bar, 300 µJ | 0.26 % |
| `modeavg_field_vector` | `:He`, 1 bar, 800 nJ | `:Ar`, 0.1 bar, 150 µJ | 0.13 % |
| `multimode_field_plasma` | `:He`, 1 bar, 800 nJ | `:Ar`, 0.1 bar, 150 µJ | 0.37 % |

Nothing else about the cases changes (core radius, length, λ0, λlims, trange, τfwhm, the
`PPT_options=NOCACHE` on the three PPT cases). They cost 1.4-16 s per mode where the
helium versions cost 0.03-0.2 s, which is most of the gate's 2 m 25 s.

New tolerances, from `test/regression/sensitivity.jl` for those four cases (100× the
one-ulp sensitivity, floor 1e-12):

| case | mode | `:Eω` | `:stats` |
| --- | --- | ---: | ---: |
| `modeavg_field_plasma` | fixed | 1.0e-12 | 2.2e-12 |
| `modeavg_field_plasma` | adaptive | 4.8e-03 | 7.3e-02 (borrowed, see below) |
| `modeavg_field_adk` | fixed | 1.0e-12 | 1.0e-12 |
| `modeavg_field_adk` | adaptive | 4.8e-10 | 1.8e-04 |
| `modeavg_field_vector` | fixed | 1.0e-12 | 8.0e-11 |
| `modeavg_field_vector` | adaptive | 6.3e-09 | 1.1e-03 |
| `multimode_field_plasma` | fixed | 6.7e-12 | 5.8e-11 |
| `multimode_field_plasma` | adaptive | 7.3e-05 | 7.3e-02 |

Ionisation feeds back into the step-size controller, so the adaptive tolerances of these
cases are loose and constrain nothing; the `:fixed` mode is what measures them (that is
what the mode is for). **`modeavg_field_plasma` in the adaptive mode has no reproducible
step sequence at all**: one ulp at the input takes it from 92 accepted steps to 98, so
`sensitivity.jl` reports a `step count` failure and `Inf` instead of a number. `Inf` is not
usable as a tolerance — `Inf <= Inf` passes, which would turn the step-count check, which
the harness means as a hard failure, into a no-op — so that one entry borrows
`multimode_field_plasma`'s 7.3e-02. The consequence is documented in `tolerances.jl` and
`test/regression/README.md`, and this branch hits it (below).

## The Group D regression record

Baselines regenerated from all three commits with the new case definitions
(`test/regression/generate.jl <commit>`), gate run on `gpu/int-D` HEAD against each.
M1 Pro, Julia 1.13.0, `-t 1`, `:estimate`, one FFTW thread, one BLAS thread, wisdom off,
`shotnoise=false`.

| baseline | result | largest difference |
| --- | --- | --- |
| `3054b686` (the `gpu/int-D` merge commit) | **460 pass, 0 fail** | `0.000e+00` on every row |
| `782f55d1` (`gpu/int-A`) | **449 pass, 1 fail** | see below |
| `fdf8dbe3` (`evanescent`) | **449 pass, 1 fail** | see below |

Every row which is not `0.000e+00` against `gpu/int-A` (`Eω` / `stats`, both modes):

| case | mode | `Eω` | `stats` | worst quantity |
| --- | --- | ---: | ---: | --- |
| `modeavg_field_plasma` | fixed | 2.511e-15 | 2.144e-14 | `stats/electrondensity` |
| `modeavg_field_plasma` | adaptive | 6.837e-05 | `Inf` | **step count: 92, baseline 93** |
| `modeavg_field_adk` | fixed | 1.245e-15 | 5.202e-15 | `stats/peak_ionisation_rate` |
| `modeavg_field_adk` | adaptive | 2.283e-13 | 8.550e-08 | `stats/peak_ionisation_rate` |
| `modeavg_field_vector` | fixed | 1.053e-15 | 6.300e-13 | `stats/transverse_integral_error_rel` |
| `modeavg_field_vector` | adaptive | 1.962e-10 | 3.294e-05 | `stats/transverse_integral_error_rel` |
| `multimode_field_plasma` | fixed | 3.785e-14 | 7.340e-13 | `stats/transverse_integral_error_abs` |
| `multimode_field_plasma` | adaptive | 7.636e-08 | 7.240e-05 | `stats/transverse_integral_error_rel` |
| `free2d_field_chi2` | fixed | 7.683e-16 | 0 | `Eω` [component 2/2, save 8/11] |
| `free2d_field_chi2` | adaptive | 4.841e-14 | 0 | `Eω` [component 2/2, save 2/11] |
| `free2d_env_chi2` | fixed | 8.348e-16 | 0 | `Eω` [component 1/2, save 11/11] |
| `free2d_env_chi2` | adaptive | 3.151e-14 | 0 | `Eω` [component 2/2, save 2/11] |

Everything else — the two Kerr field cases, the no-THG case, both Raman cases, the mixture,
both GNLSE cases, the envelope cases, the radial and 3-D free-space cases, the gradient,
the taper and the legacy-boundary case — is exactly `0.000e+00` in both classes and both
modes.

Against `evanescent` (`fdf8dbe3`) the same rows appear, plus the radial ulps that
`gpu/02-radialgrid` introduced and `gpu/int-A` recorded:

| case | mode | `Eω` | `stats` |
| --- | --- | ---: | ---: |
| `radial_field_kerr` | fixed | 1.901e-15 | 0 |
| `radial_field_kerr` | adaptive | 6.713e-11 | 0 |
| `radial_env_kerr` | fixed | 2.866e-15 | 0 |
| `radial_env_kerr` | adaptive | 4.017e-11 | 0 |

Reading the record:

- The χ⁽²⁾ rows are `gpu/15`'s FMA contraction and reproduce that branch's numbers to every
  digit (4.841e-14 and 3.151e-14 in the adaptive mode, 7.683e-16 and 8.348e-16 fixed),
  three orders inside the tolerance.
- The plasma rows are the first measurement of Group D's plasma rewrite against anything,
  because until this commit the cases carried no plasma. In the attributable `:fixed` mode
  they are 1.0e-15 to 3.8e-14 in `Eω` — the `accumulate!`-plus-correction scan against the
  old serial `Maths.cumtrapz!` loop, exactly the rounding-level difference `gpu/13`
  predicted, and two orders inside the tolerance even for the weakest of the four modes in
  `multimode_field_plasma`.
- **The one failure** is `modeavg_field_plasma` in the adaptive mode: 92 accepted steps
  against the baseline's 93. The step count is a hard failure by design. It is that
  rounding-level difference flipping one accept/reject decision in the case whose step
  sequence the sensitivity study shows is not reproducible at the one-ulp level (one ulp at
  the input moves it by six steps). The same case in the `:fixed` mode, where the step
  sequence is imposed, is 2.5e-15. Nothing else in the matrix changed its step count.

A Group E branch does not inherit this failure: it compares against `gpu/int-D`, where
every row is zero.

## Tests

All on M1 Pro / Julia 1.13.0, `-t 1` unless stated, `:estimate`, one FFTW thread, one BLAS
thread, wisdom off, through the timing wrapper.

| what | result |
| --- | --- |
| `test/test_regression.jl` ×3 baselines | 460/0, 449/1, 449/1 — the record above |
| `test/test_device.jl` (JLArrays), `-t 1` | **611 pass, 0 fail**, 31 testsets |
| `test/test_device.jl` (JLArrays), `-t 4` | **611 pass, 0 fail**, 31 testsets |
| `test/test_metal.jl` (Metal 1.11.1, hardware) | **386 pass, 0 fail**, 16 testsets |
| `test/test_interface.jl` | **340 pass, 0 fail**, 11 testsets (334 before the restored assertion) |
| `test/test_chi2.jl` | **42 pass, 0 fail**, 9 testsets |
| `test/test_raman.jl` | pass (bare `@test`s, exit 0) |
| `test/test_ionisation.jl` | **25 pass, 0 fail**, 9 testsets |
| `test/test_vectorplasma.jl` | **2 pass, 0 fail** |
| `test/test_freespace.jl` | **77 pass, 0 fail**, 55 testsets |
| `test/test_multimode.jl` | **6 pass, 0 fail**, 5 testsets |
| `test/test_output.jl` | **117 pass, 0 fail**, 8 testsets |
| `include("docs/make.jl")` | the 9 pre-existing unresolved `@ref`s, none new |
| parse check | 139 files under `src/`, `test/`, `examples/`, all parse |

The environments are the ones the test files' headers describe: a stacked environment with
`JLArrays` for `test_device.jl` and one with `Metal` for `test_metal.jl`, both with this
worktree `dev`ed in.

### Examples

Run through `Luna.run` with the plotting stack cut off (every line from the first
`Plotting.`/`plt.` call dropped, plus the `PyPlot` imports):

| example | result |
| --- | ---: |
| `low_level_interface/freespace/free2D_bbo.jl` (χ⁽²⁾, `TransFree2D`) | OK, 19.8 s |
| `low_level_interface/freespace/free2D_bbo_env.jl` (χ⁽²⁾ envelope) | OK, 30.3 s |
| `low_level_interface/Raman/Raman_modeAvg.jl` | OK, 534.2 s |
| `low_level_interface/Raman/Raman_modeAvg_env.jl` | OK, 91.7 s |
| `low_level_interface/Raman/Raman_modeAvg_noTHG.jl` | OK, 9.4 s |
| `simple_interface/RamanSSFS.jl` | OK, 21.4 s |
| `low_level_interface/plasma_ssfbs_modeAvg.jl` | OK, 25.7 s |
| `low_level_interface/basic_modeAvg.jl` (Kerr + plasma) | OK, 5.0 s |
| `low_level_interface/basic_modeAvg_field_noTHG.jl` | OK, 3.5 s |
| `low_level_interface/basic_modal.jl` (multimode, plasma) | OK, 151.7 s |

`examples/low_level_interface/polarisation/modal_nonvector_plasma.jl` does not run, and
does not run on `evanescent` either: line 37 passes `linop` to `Stats.default` two lines
before line 39 defines it (`UndefVarError: linop`). Pre-existing, unrelated to this
project, not fixed here.

## Known gaps and open questions

- **`Nonlinear.device_capable` stays derived from `kind`.** A `PlasmaCumtrapz` holding a
  rate that has no kernel (the direct `Ionisation.IonRatePPT`, or a user-written rate)
  therefore reports `device_capable == true`. `gpu/13`'s review accepted this: the
  combination is unreachable through the simple interface (`Interface.makeplasma!` only
  ever builds `IonRateADK` or a uniformly spaced cached PPT rate), and a low-level device
  run with such a rate is refused at setup by `Ionisation.device_rate` with a message
  naming the alternatives. It is a `prop_capillary`-level inaccuracy only if someone passes
  a hand-built plasma response into it, which the interface does not support.
- **The `modeavg_field_plasma` adaptive row** is expected to fail for any branch which
  moves the plasma arithmetic at rounding level and is compared against a pre-Group-D
  baseline. See above; a Group E branch compares against `gpu/int-D` and does not see it.
- **`_refuse_batched_legacy` checks the one-argument `Nonlinear.kind(r)` only**, so a
  response which is `Batched` only for `Val(2)` would slip past it. No such response
  exists, and the batched response's own shape check catches it.
- **O₂ Raman is a data gap, not a device gap.** `PhysData.raman_parameters(:O2)` has
  `# TODO τ2r` and `# TODO τ2v` (untouched since 2022), so `Raman.raman_response(t, :O2)`
  raises. `Interface.jl` does not include `:O2` in its default Raman gas list, so
  `prop_capillary` never reaches it. For `GPU_PROGRESS.md`.
- **Batched Raman and plasma buffers scale with the column count** (a Raman response holds
  `2 × (2nt × ncols)` real plus `(nt+1) × ncols` complex; the plasma response holds four
  block-sized arrays). At `nt = 8192` with 256 radial points that is ~100 MB per Raman
  response where the per-column version held ~0.4 MB. Group E should add a multi-column
  Raman gate case and decide about chunking there; no Raman example or gate case is
  multi-column today.
- **The χ⁽²⁾ responses have kernels but no device transform.** `TransFree2D`, `TransFree`
  and `TransRadial` are host-only until `gpu/20`/`gpu/21`, so a χ⁽²⁾ propagation still runs
  entirely on the host. The responses are exercised on device blocks directly in
  `test_device.jl` and `test_metal.jl`.

## How a Group E branch chooses its base

Group E branches from `gpu/int-D`, so the question their gate has to answer is "did *this
branch* change the output", i.e. the merge-base of `HEAD` with `gpu/int-D`:

```
julia --project=$PWD -t 1 test/regression/generate.jl gpu/int-D
LUNA_REGRESSION_BRANCH=gpu/int-D julia --project=$PWD -t 1 test/test_regression.jl
```

Every case must come back exactly `0.000e+00` unless the branch says which one moved and
why. The other form,

```
LUNA_REGRESSION_BASE=fdf8dbe3 julia --project=$PWD -t 1 test/test_regression.jl
```

accumulates the whole project's movement and is what a per-group record uses; on
`gpu/int-D` it is the table above, including the one step-count failure.

Both baselines were regenerated with the new (argon) case definitions, so a baseline
generated before this branch is stale for the four plasma cases and has to be regenerated;
`generate.jl` overwrites all 21 case files, so there is nothing to clean up.
