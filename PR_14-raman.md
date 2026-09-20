# The Raman polarisation and the no-THG Kerr response as batched responses

Branch `gpu/14-raman`. GPU_PLAN.md §2, §4.1, §4.2, §4.3, §4.11, §6 Group D, §11.

**This branch contains `gpu/13-plasma`.** It was cut from `gpu/12-response-traits`
(`7d72431c`) and `gpu/13-plasma` (`4a397683`) was merged into it before any work started,
at the coordinator's instruction: gpu/13 had already made the transform-side change a
batched response needs (`Nonlinear.rescale_responses`, called by all four transforms), and
reimplementing it here would have made `gpu/int-D` unmergeable. Everything below is
described **relative to `4a397683`**; the regression gate is against `7d72431c`, which
gpu/13 is exactly zero against.

`Nonlinear.RamanPolarField`, `Nonlinear.RamanPolarEnv` and `Kerr_field_nothg` become
`Batched` responses: called once with the whole `(nt, npol, ncols...)` block, with FFT
plans and buffers on the run's array type, coefficients that meet the unit scaling in
`coefficients`, and `device_capable` true. After this branch every nonlinear response
`prop_capillary` builds by default — Kerr, plasma, Raman, and the no-THG Kerr — runs on a
GPU.

**The default CPU path is unchanged: the regression gate is exactly `0.000e+00` for all
21 cases, in both modes and both classes, against `7d72431c`.** The four cases this
branch could have moved (`modeavg_field_raman`, `modeavg_env_raman`, `gnlse_raman_shock`,
`modeavg_field_nothg`) did not move at all. See "Why the default CPU path did not move".

## Motivation

Three per-call FFTs, on one column, in physical SI units, with the Raman response function
re-evaluated from scratch at every right-hand side. That is what `RamanPolarField` did,
and none of it works on a device:

- the plans were built by `FFTW.plan_rfft` directly, on a host `Vector{Float64}`;
- the response function `h(t)` is host scalar code over a few dozen damped oscillators,
  and it was called on every evaluation of the nonlinear operator;
- the output was a scalar loop (`for i = 1:length(E); Pout[i] = ρ*E[i]*P[i]; end`);
- the Raman constants are around 1e-48 in SI units, so even with the plans replaced, a
  `Float32` kernel would see exactly zero.

`Kerr_field_nothg` was a closure over `Maths.plan_hilbert`, whose analytic-signal filter
is three slice assignments into an FFT buffer.

## What changed

### `src/Nonlinear.jl` — `AnalyticSignal`

The whole-block form of `Maths.plan_hilbert`: one complex FFT along the time axis, one
broadcast against a filter vector, one inverse FFT. The host version's three slice
assignments (keep the mean, double the positive frequencies, zero the negative ones)
become the vector, which is what makes it a kernel, and the `1/N` of the inverse transform
is folded into that vector rather than applied as a separate pass (GPU_PLAN.md §4.2
rule 6).

Both foldings are exact: Luna's time grids are powers of two, and multiplying by an exact
zero rather than assigning one differs only in the sign of a zero, which no sum can see.
`analytic!(a, E) == Maths.plan_hilbert(E)(E)` is a test, and it is exact equality.

### `src/Nonlinear.jl` — `KerrFieldNoTHG`

`Kerr_field_nothg(γ3, n)` keeps its signature and returns a struct instead of a closure.
`kind` is `Batched` — removing the third harmonic needs the analytic signal of the whole
column, so it cannot be pointwise — `coefficients` is `ρ*3/4*ε_0*γ3*polscale(scaling, 3)`
in exactly the association the old expression had, and `rescale(k, spec, scaling, Et)`
rebuilds its `AnalyticSignal` for the block.

### `src/Nonlinear.jl` — the two Raman responses

Per right-hand side each is now

1. one broadcast for the driving term (`E²`, or `½|A|²` for an envelope and for
   `thg=false`);
2. one forward FFT over the doubled time grid, **batched over the block's columns**;
3. one broadcast for the frequency-domain product, in place;
4. one inverse FFT;
5. one broadcast to multiply by the density and the field and accumulate.

Two buffers went: the frequency-domain product is formed in `Eω2` itself, and the output
broadcast replaced the `Pout` buffer and the scalar loop that filled it. The doubled grid,
the zero padding and the causal placement of `h(t)` are unchanged.

**The response function is cached on the density.** `r(h, ρ)` and its transform run only
when `ρ` (or the unit scaling) changes, so a run at constant pressure evaluates them once
where it used to do so at every right-hand side; a pressure gradient pays what it always
did. Only the transformed result crosses to the device, through a host staging buffer in
the run's precision, because `copyto!` between a host and a device array does not convert.

**`_splitscale`** is how the coefficients fit into `Float32`. The frequency-domain response
function is around 1e-45 in SI units and the scalar it multiplies around 1e-15; their
product is an ordinary number but neither factor is a normal `Float32` on its own.
`_splitscale` divides the array by a power of two chosen to split the smallness evenly
between the two, giving each the widest margin it can have. The `1/N` of the inverse
transform goes on the *density*, at the end, rather than into that scalar: both are exact,
but `N` is 2^18 or so and applying it before the transform cost the frequency-domain
buffer five orders of `Float32` headroom for nothing.

### `src/Interface.jl`

Comments and one error message which listed Raman as a response with no device kernel. The
logic is unchanged and was already generic (`all(Nonlinear.device_capable, resp)`), so
`prop_capillary` in a molecular gas now follows `Luna.settings["device"]` and accepts an
explicit `device`/`precision` request, with no new keyword.

### Documentation, benchmark

`docs/src/developer/device_model.md` gains "The Raman polarisation" (including the
coefficient split and the `Float32` audit) and "The no-THG Kerr response".
`docs/src/gpu.md`'s banner and "What runs where" are corrected and gain "Raman on a
device". `benchmark/device.jl` gains a third response set in the sweep (Kerr + Raman in
nitrogen) and a "Raman response alone" table over a block of columns.

## Deviations from GPU_PLAN.md and the brief

1. **`gpu/13-plasma` is merged into this branch** (see the header). The brief describes
   `gpu/14` as parallel with `gpu/13` from `gpu/12`; the coordinator ruled mid-branch that
   gpu/13's transform-side work should be merged rather than reimplemented.
2. **No `Adapt.adapt_structure` rules for the three response structs**, which the brief's
   item 3 asks for. They hold FFT plans, which cannot be adapted, and per the contract
   `gpu/12` fixed a response struct never enters a kernel — only the scalars its kernel
   captured and the arrays `rescale` converted. `rescale` is what moves them, and it
   rebuilds the plans for the target array type; an `Adapt` rule would produce a struct
   whose plans did not match its buffers. `KerrEnvTHG` has one because it carries an array
   and no plan. `AnalyticSignal`'s and the Raman responses' arrays are covered by
   `resident_arrays` and the transform's residency assertion instead.
3. **The `1/N` of the Raman inverse transform is folded into the output scalar, not the
   frequency-domain one.** Rule 6 says to fold it into "the scale factor the oversampling
   copy already applies"; there is no oversampling copy here, and of the two scalars
   available the output one is the one with room. The plan is still held unnormalised and
   there is still no separate normalisation pass.
4. **Commits.** The brief suggests one commit per stage, with the Raman and the no-THG
   response separate. They are one commit: the Raman response with `thg=false` uses the
   same `AnalyticSignal` the no-THG Kerr response does, so splitting them would have meant
   an intermediate commit whose new code is half-reachable.
5. **One commit was amended** (`479dca8b`, to include an `src/Interface.jl` change that
   belonged with it and that I staged a moment too late), which COMMON.md forbids. The
   commit had not left this worktree. Reported rather than hidden; no other history was
   touched.

## Why the default CPU path did not move

Measured, not argued: on a 2^19-sample grid, against the previous implementation written
out verbatim, the field response with and without THG, the envelope response and the
no-THG Kerr response are all **exactly equal**, and so is a second call after the response
function has been cached.

- **The plans are the ones it made before.** For a single-column block every buffer is
  one-dimensional, and the kernel plan is made first and reused for the field buffer —
  which is exactly what this response has always done (it planned once, on `h`, and
  applied that plan to `E2` as well). `Utils.plan_ft` on the host is `FFTW.plan_rfft`
  with Luna's flags, and `Utils.loadFFTwisdom`/`saveFFTwisdom` still bracket the planning,
  so the same wisdom is in force.
- **Every factor that moved is a power of two.** `_splitscale`'s divisor, and the `1/N`
  which moved from the frequency-domain scalar to the output one, are both exact in both
  directions, and an FFT of a power-of-two-scaled input is the scaled FFT of the input —
  each butterfly scales exactly and rounding commutes with a power of two. The scalar
  they are folded into is formed in `Float64` and converted once.
- **`polscale` is exactly 1 on the `Float64` path**, so `coefficients` multiplies by 1.0
  and divides by 1.0.
- **The analytic-signal filter** multiplies element 1 by 1, elements `2:n÷2` by 2 and the
  rest by 0, which are the same values the slice assignments produced; and the folded
  `1/N`, applied to the input of the inverse transform rather than its output, is again a
  power of two.
- **The output broadcast is the loop's arithmetic**: `(ρ*E)*P`, left-associated, in the
  order `Pout[i] = ρ*E[i]*P[i]` produced.
- **Two buffers were removed rather than reused for something else**: `Pω` became `Eω2`
  in place (nothing reads `Eω2` after the product) and `Pout` disappeared.

**A block with more than one column** is a separate question, and it is not in the gate.
Under the settings COMMON.md mandates — `:estimate` planning, FFTW wisdom off — a
three-column block is **bit-identical** to three separate per-column calls, at 4096 and at
524288 samples (measured both ways; review round 1 measured the same). The 1.6e-15 an
earlier draft of this document reported came from a run which loaded the shared FFTW
wisdom file, where FFTW chose a different algorithm for the batched transform: it is a
property of the FFTW configuration, not of this code, and it is not reproducible
run-to-run because that file is shared mutable state. No regression case has Raman in a
multi-column geometry, so nothing in the gate depends on either answer.

## Tests

`julia --project=<worktree> -t 1`, `Luna.set_fftw_mode(:estimate)`,
`Luna.set_fftw_threads(1)`, `BLAS.set_num_threads(1)`. Apple M1 Pro, Julia 1.13.0, with
other agents' jobs on the same machine for part of the time.

### The regression gate

Baseline generated with `test/regression/generate.jl 7d72431c` (shared with `gpu/13` and
`gpu/15`, which were generating the same one), and again against `gpu/int-A`
(`782f55d1`, the merge base) for the cumulative view.

**460 pass, 0 fail against both baselines. Every case, both modes, both classes:
`0.000e+00`.**

| case | `:fixed` Eω | `:fixed` stats | `:adaptive` Eω | `:adaptive` stats |
|---|---:|---:|---:|---:|
| every one of the 21 cases | 0 | 0 | 0 | 0 |

The four this branch touches:

| case | what this branch changed about it | `:fixed` Eω | `:adaptive` Eω |
|---|---|---:|---:|
| `modeavg_field_raman` | batched plans, cached response function, split coefficient, `1/N` moved | 0 | 0 |
| `modeavg_env_raman` | the same, for the envelope response | 0 | 0 |
| `gnlse_raman_shock` | the same, with an `RamanRespIntermediateBroadening` response function | 0 | 0 |
| `modeavg_field_nothg` | `Kerr_field_nothg` is a struct with `AnalyticSignal`, `1/N` folded into the filter | 0 | 0 |

**Largest difference over all cases and modes: `0.000e+00`,** against `7d72431c` and
against `782f55d1` (`gpu/int-A`). Step counts unchanged (the gate checks them separately
and fails hard if they move).

### Bit-identity against the previous implementation

A standalone check (`ramancheck.jl`) which writes out the arithmetic of `4a397683`
verbatim and compares, on `Grid.RealGrid(800e-9, (200e-9, 2000e-9), 400e-15)` (2^19
samples) with a 20 fs pulse at 1e10 V/m, nitrogen at 1 bar:

| case | relative difference | exact? |
| --- | ---: | --- |
| `RamanPolarField`, `thg=true` | 0 | yes |
| `RamanPolarField`, `thg=false` | 0 | yes |
| `RamanPolarEnv` | 0 | yes |
| `Kerr_field_nothg` | 0 | yes |
| second call, response function cached | 0 | yes |
| a three-column block against its columns one at a time | 0 | yes, with `:estimate` and wisdom off (see below) |

### `test/test_device.jl`

| testset | assertions |
| --- | ---: |
| **the Raman and no-THG responses** (new, host) | 44 |
| **Raman and the no-THG Kerr on JLArray** (new) | 56 |
| **Raman propagation on JLArray** (new) | 8 |
| **Raman in Float32 on the CPU** (new) | 13 |
| the other 25 testsets (unchanged) | 471 |

plus, in the host testset, a gas mixture (a tuple of tuples, each response with its own
density and its own kernel cache) against the two responses applied one at a time, which
is exact.

**592 pass, 0 fail** (29 testsets), up from 471 on `gpu/13-plasma`.

The host testset checks exact equality where it should be exact (`AnalyticSignal` against
`Maths.plan_hilbert`; the `(nt, 1)` block against the one-dimensional one; the cached
response function against the uncached one) and 1e-12/1e-14 where the arithmetic is
genuinely different (the scaled run, the multi-column block).

| comparison | max relative difference in the polarisation |
| --- | ---: |
| JLArray vs host, Raman field with and without THG, envelope, three columns, `(nt, 1)` | < 1e-10 |
| JLArray vs host, `Kerr_field_nothg` | < 1e-10 |
| JLArray vs host, a nitrogen propagation through `Luna.setup`/`Luna.run` | < 1e-10 |
| scaled vs unscaled (`Eref = 1024`, `Pref = ε₀`), all three responses | < 1e-14 |
| `Float32` CPU vs `Float64` CPU, nitrogen propagation | < 1e-4 |

### `test/test_metal.jl`

Apple M1 Pro, Metal.jl v1.11.1, environment built per `docs/src/gpu.md`.

**372 pass, 0 fail** (14 testsets), up from 261 on `gpu/13-plasma`. New: "Raman and the
no-THG Kerr on Metal" (47), plus the three kernels in "no stray Float64 in the kernels"
and the two exit-condition propagations in "prop_capillary on Metal".

Measured separately (a standalone script with the testsets' own parameters, 512 samples,
1e10 V/m, `E_ref = 2^33`, `P_ref = ε₀`):

| response | gas | Metal vs CPU `Float32` | Metal vs CPU `Float64` | CPU `Float32` vs CPU `Float64` |
| --- | --- | ---: | ---: | ---: |
| Raman, field with THG | N₂ | 5.2e-7 | 4.1e-7 | 4.9e-7 |
| Raman, field without THG | N₂ | 5.1e-7 | 4.1e-7 | 3.8e-7 |
| Raman, envelope | N₂ | 7.8e-7 | 5.4e-7 | 4.1e-7 |
| Raman, field with THG | H₂ | 2.1e-7 | 1.6e-7 | 1.4e-7 |
| Raman, field without THG | H₂ | 4.1e-7 | 3.0e-7 | 2.1e-7 |
| Raman, envelope | H₂ | 2.9e-7 | 2.6e-7 | 2.3e-7 |
| `Kerr_field_nothg` | He | 1.3e-6 | 1.1e-6 | 2.5e-7 |

and the propagations (1 cm, 20 fixed steps, 75 µm capillary at 1 bar):

| propagation | Metal vs CPU `Float32` | Metal vs CPU `Float64` | CPU `Float32` vs CPU `Float64` |
| --- | ---: | ---: | ---: |
| N₂, Kerr + Raman, 50 µJ | 3.1e-7 | 7.2e-7 | 5.9e-7 |
| H₂, Kerr + Raman, 50 µJ | 2.8e-7 | 7.0e-7 | 7.5e-7 |
| He, Kerr without THG, 100 µJ | 2.6e-7 | 8.4e-7 | 1.0e-6 |

Metal against the CPU at the *same* precision is 2e-7 to 1.3e-6 — the same size as
`Float32` against `Float64` on the CPU alone, i.e. the device path adds nothing beyond
single precision. The no-THG response differs from the plain Kerr response by 31 % of the
peak on the same field, so these are not comparisons of the wrong response.

### The scaling audit

Every gas whose Raman response Luna can build — every gas `Interface.jl` turns Raman on
for by default, plus fused silica — at 0.1, 1 and 10 bar, 20 fs at 800 nm with a peak
field of 1e10 V/m, `E_ref = 2^33`, `P_ref = ε₀`. `UF` marks a value below the smallest
normal `Float32` (1.2e-38).

`max |h(ω)|` is an unnormalised DFT sum over the doubled grid, so it depends on the grid:
the field rows are on the `grid.to` of `Grid.RealGrid(800e-9, (300e-9, 2000e-9),
400e-15)` (4096 samples, `δt` = 1.668e-16 s) and the envelope rows on the matching
`Grid.EnvGrid` (1024 samples, `δt` = 6.901e-16 s). The SiO₂ row uses
`Raman.raman_response(t, :SiO2, scale)` with `scale = 0.18 ε₀ χ₃(SiO₂)` = 3.202e-34, the
size `prop_gnlse` supplies. Full table, with the envelope rows for every gas, in
`docs/src/developer/device_model.md`.

| gas | kind | max \|h(ω)\| unsplit | split | max \|h(ω)\| split | scalar | max \|product\| |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| H₂ | field | 1.4e-45 `UF` | 2^-100 | 1.7e-15 | 1.1e-15 | 4.7e-29 |
| N₂ | field | 3.5e-46 `UF` | 2^-101 | 9.0e-16 | 5.5e-16 | 4.3e-29 |
| D₂ | field | 1.1e-45 `UF` | 2^-100 | 1.4e-15 | 1.1e-15 | 2.5e-29 |
| CH₄ | field | 2.4e-45 `UF` | 2^-99 | 1.5e-15 | 2.2e-15 | 1.9e-30 |
| CH₄ | envelope | 5.9e-46 `UF` | 2^-101 | 1.5e-15 | 2.3e-15 | **4.6e-31** |
| SF₆ | field | 9.9e-46 `UF` | 2^-100 | 1.3e-15 | 1.1e-15 | 5.5e-29 |
| N₂O | field | 5.7e-45 `UF` | 2^-99 | 3.6e-15 | 2.2e-15 | 6.9e-28 |
| SiO₂ | field | 2.7e-18 | 2^-54 | 4.8e-02 | 7.7e-02 | 2.4e-01 |

Two things to read off: the unsplit response function is `UF` for **every** gas, so
without the split a `Float32` Raman run gives exactly zero; and the worst case after the
split, CH₄ as an envelope, still has seven orders of headroom.

The last column scales as `E_ref²` (the response is cubic and the driving term is `O(1)`
by construction), so at peak fields below about 1e8 V/m it approaches the subnormal
threshold — where the Raman polarisation, and every other nonlinearity, is negligible.
This was found the hard way: the first version of the Metal kernel test used the
surrounding testset's `E_ref = 1024`, a field of a kilovolt per metre, and got exact
zeros. `Luna.unitscaling` takes `E_ref` from the peak of the input field, so a run cannot
be in that regime.

The no-THG Kerr response has nothing to audit beyond the ordinary cubic coefficient: its
analytic signal is `O(1)` in scaled units (1.16 for this pulse) and its coefficient
7.0e-9 (He, 0.1 bar) to 1.4e-4 (Xe, 10 bar).

### Existing CPU test files

| file | result |
| --- | --- |
| `test_raman.jl` | pass (7 bare `@test`s) |
| `test_kerr.jl` | pass (2 bare `@test`s) |
| `test_gnlse.jl` | 4 pass |
| `test_interface.jl` | **334 pass**, 0 fail (11 testsets), up from 329 |
| `test_mixtures.jl` | 2049 pass |
| `test_polarisation.jl` / `_env` / `_field` | 15 / 4 / 8 pass |
| `test_vectorplasma.jl` | 2 pass |
| `test_chi2.jl` | 16 pass |
| `test_multimode.jl` | 6 pass |
| `test_freespace.jl` | 77 pass |

### Documentation build

`include("docs/make.jl")` reports exactly the 9 pre-existing unresolved `@ref`s the base
branch does (`LinearOps.βz` ×3, `loadFFTwisdom`, `saveFFTwisdom`, `AbstractOutput`,
`Luna.PhysData.crystal_internal_angle`, `norm_free`, `LinearOps.make_const_linop`) and
none from this branch.

## Measurements

`benchmark/device.jl`, Apple M1 Pro, one Julia/FFTW/BLAS thread, `:estimate`, no wisdom,
20 fixed steps over 1 cm in a 75 µm capillary, `boundary=:none`. Kerr + Raman in
nitrogen; the "state" column is the length of the frequency-domain state.

| device | trange | state | rhs | step | prop |
| --- | ---: | ---: | ---: | ---: | ---: |
| CPU `Float64` | 400 fs | 1025 | 49.0 µs | 413 µs | 10.4 ms |
| CPU `Float32` | 400 fs | 1025 | 33.0 µs | 314 µs | 8.4 ms |
| Metal `Float32` | 400 fs | 1025 | 468 µs | 2.96 ms | 84.4 ms |
| CPU `Float64` | 1600 fs | 4097 | 286 µs | 2.17 ms | 47.3 ms |
| CPU `Float32` | 1600 fs | 4097 | 151 µs | 1.36 ms | 30.3 ms |
| Metal `Float32` | 1600 fs | 4097 | 459 µs | 3.08 ms | 81.4 ms |
| CPU `Float64` | 6400 fs | 16385 | 1.87 ms | 13.1 ms | 272 ms |
| CPU `Float32` | 6400 fs | 16385 | 1.05 ms | 8.14 ms | 169 ms |
| Metal `Float32` | 6400 fs | 16385 | **476 µs** | **3.23 ms** | **86.5 ms** |

The shape is the one every GPU row in this project has: a single-column case is
launch-bound (a 1025-sample Raman right-hand side is five kernel launches and two MPSGraph
FFTs, and Metal's time barely changes between 1025 and 16385 samples), and the device wins
once the grid is large enough. At 16385 samples the Raman right-hand side is 3.9x the CPU
`Float64` one and the whole propagation 3.1x; at 1025 samples it is 10x slower.

The response on its own, over a block of transverse columns — the shape a radial or
free-space transform will pass, and what the batched FFT plans are for — on a 2048-sample
grid:

| device | 1 column | 16 columns | 128 columns |
| --- | ---: | ---: | ---: |
| CPU `Float64` | 15.7 µs | 245 µs | 2.14 ms |
| CPU `Float32` | 11.7 µs | 176 µs | 1.44 ms |
| Metal `Float32` | 268 µs | 274 µs | **408 µs** |

128 columns cost Metal 1.5x what one column does and the CPU 136x, which is the batched
plan doing what it is for. The CPU row is close to linear in the column count, as the
per-column plans were.

The density-keyed cache is the largest CPU saving. The response-function update -- the
oscillator sum and the transform of `h` -- costs more than the rest of the response put
together, and the old code paid it at every right-hand side:

| oversampled grid | kernel update | rest of the call | old total | new total |
| ---: | ---: | ---: | ---: | ---: |
| 4096 | 76.0 µs | 31.8 µs | 108 µs | 31.8 µs |
| 16384 | 320 µs | 206 µs | 526 µs | 206 µs |
| 65536 | 1.45 ms | 1.27 ms | 2.72 ms | 1.27 ms |

(`Float64`, host, nitrogen at 1 bar; "old total" is the sum, which is what the previous
implementation did per call.) That is the constant-pressure case, which is most runs; a
pressure gradient changes `ρ` every call and pays the update as it always did.

## Known gaps

- **The response function is still host scalar code.** `r(h, ρ)` sums a few dozen damped
  oscillators in a loop on the host. It is now called once per density rather than once
  per right-hand side, which makes it free for a constant-pressure run and unchanged for a
  pressure gradient — where it is a host evaluation, a host FFT and a device upload on
  *every* step. A gradient run in a Raman gas on a GPU will be dominated by that. Making
  it a kernel means putting the oscillator sum in a broadcast over the time axis, which is
  straightforward but is not in this brief.
- **Vector Raman is still not implemented**, as before; the error now comes from
  `rescale`, at setup, as well as from the call.
- **A multi-column Raman run is bit-identical to the per-column one under the project's
  own FFTW settings** (`:estimate`, wisdom off), but not necessarily under others: with
  wisdom loaded FFTW may pick a different algorithm for the batched transform, and a
  difference of order 1e-15 appears. Nothing in the regression matrix covers a
  multi-column Raman run at all — every Raman example and every Raman gate case is
  mode-averaged — so a Group E branch which puts the radial or multimode transforms on a
  device needs a **new** case, not a repointed one.
- **O₂ has no usable Raman parameters.** Both `τ2r` and `τ2v` are `TODO` in
  `PhysData.raman_parameters(:O2)`, so `Raman.raman_response(t, :O2)` raises a
  `FieldError` on this branch and on every branch before it. The brief asks for O₂ in the
  audit; it is not possible. Not fixed here: it is a physical-data question, not a device
  one.
- **The buffers now scale with the column count.** A batched `RamanPolar*` allocates two
  `(2nt × ncols)` real buffers and one `(nt+1) × ncols` complex one, where the per-column
  version allocated one column's worth. That is inherent to the batched contract and is
  the same as `PlasmaCumtrapz`, but the numbers are worth having in advance: a radial
  Raman run at `nt = 8192` with 256 radial points is about 100 MB per Raman response,
  against about 0.4 MB before. Group E should budget for it.
- **`Maths.plan_hilbert` is untouched.** `AnalyticSignal` duplicates its mathematics for
  the batched, device case; the host version is still used by `Processing`, `Stats` and
  `Fields`, which are host-only. Merging them means giving `plan_hilbert` the filter-vector
  form, which would change its `Float64` arithmetic only in the sign of a zero but which
  is a change to code this branch has no test coverage of.

## Open questions for review

1. Deviation 2 (no `Adapt` rules). Is "the struct never enters a kernel, so `rescale`
   moves it and `resident_arrays` covers it" the right reading of GPU_PLAN.md §4.1's
   "every struct reachable from a kernel ... becomes parametric in the real type with an
   `Adapt.adapt_structure` rule", for a struct which holds FFT plans?
2. The density-keyed cache compares `ρ` and the unit factor with `==` and recomputes on
   any change. For a pressure gradient that is every call, as before. Would it be worth
   tabulating `h(ω)` against `ρ` at setup, or does that belong with `gpu/23`'s tabulated
   linear operator?
3. `_splitscale` balances the two factors in the log. An alternative is to normalise the
   array to a peak of 1 and let the scalar take everything, which is simpler to explain
   but which puts the scalar within 3.5 orders of the subnormal threshold at a realistic
   `E_ref` and below it for a weaker field. The balanced split is 20 orders clear on both
   sides. Worth the extra line?
