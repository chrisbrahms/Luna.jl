# The χ⁽²⁾ responses as vector-pointwise kernels

Branch `gpu/15-chi2`, base `gpu/12-response-traits` (`7d72431c`). GPU_PLAN.md §4.1, §4.2,
§4.3 (the vector-pointwise contract), §6 Group D, §11.

This is the χ⁽²⁾ branch of Group D. `Nonlinear.Chi2Field` and `Nonlinear.Chi2Env` stop
being columnwise host loops over four heap-allocated work vectors and become
`VectorPointwise` responses: the rotation, the contracted field products and the 3×6
tensor contraction become `SVector`/`SMatrix` expressions inside the two component
broadcasts `NonlinearRHS.Et_to_Pt!` already emits for a two-component block. They run on
Metal, in `Float32`, with no buffer and no scalar indexing.

**The regression gate is 460 pass, 0 fail against both `7d72431c` and `gpu/int-A`.
`free2d_field_chi2` and `free2d_env_chi2` move at rounding level and nothing else moves
at all.** The numbers and the reason are below.

## Motivation

GPU_PLAN.md §6 Group D makes the χ⁽²⁾ examples a target for the GPU. The old
implementation could not be one: it indexed the field one time sample at a time, wrote
three `Vector{Float64}` work buffers per sample through `mul!`, and would have had to be
serialised through `Nonlinear.HostResponse` on any device. It was also the response the
vector-pointwise kind in `gpu/12-response-traits` was designed for — the Kerr vector forms
fit it, but they do not exercise it, because their two components share almost no work.

## What changed

### `src/Nonlinear.jl` — the responses

Both structs are now parametric in one real type and hold only `isbits` static matrices:

```julia
struct Chi2Field{T}
    χ2::SMatrix{3, 6, T, 18}
    toCrystal::RotMatrix3{T}
    toLab::RotMatrix3{T}
    χ2_toLab::SMatrix{3, 6, T, 18}
end
```

`Chi2Env` adds the carrier phase `C`. The four work vectors (`El`, `Ec`, `Enl`, `Pl` /
`Al`, `Ac`, `Anl`, `Pl`) are gone. `χ2` and `toLab` are unused by the response itself, as
before, and are kept because `test_chi2.jl` and
`examples/low_level_interface/freespace/chi2_benchmarking.jl` read them.

The protocol methods:

| method | what it is |
| --- | --- |
| `kind(::Chi2Field)`, `kind(::Chi2Env)` | `VectorPointwise()`, **unconditionally** |
| `coefficients(c, ρ, scaling)` | `ε_0*Luna.polscale(scaling, 2)` — quadratic in the field, so one power of `E_ref`; the density is ignored, as it was |
| `vector_kernel(c::Chi2Field, E, ρ, scaling)` | the `(ex, ey) -> SVector{2}` kernel |
| `vector_expr(c::Chi2Env, Ex, Ey, ρ, scaling)` | the same with the carrier phase as a third broadcast argument, as `KerrEnvTHG` does |
| `rescale` | converts the crystal matrices to the run's real type, moves `C`, keeps the default CPU path's pass-through |
| `resident_arrays(c::Chi2Env)` | `(c.C,)` |
| `Adapt.adapt_structure` | both |
| `device_capable` | derived from `kind`: true |

`field_products!`/`env_products!` keep their signatures and their docstrings and now write
the kernel's own expression into the caller's vector, so there is one arithmetic body per
response (GPU_PLAN.md §4.2 rule 1). The columnwise call operators are still there, still
accumulate, and go through the same kernel as the fused path.

### Unit scaling

The state is `e = E/E_ref` and the polarisation buffer holds `p = P/(P_ref E_ref)`. A
response of polynomial degree `n` carries `Luna.polscale(scaling, n) = E_ref^(n-1)/P_ref`.
χ⁽²⁾ is quadratic, so `n = 2`:

`p = ε₀χ⁽²⁾E²/(P_ref E_ref) = ε₀χ⁽²⁾(E_ref e)²/(P_ref E_ref) = (ε₀ E_ref/P_ref)·χ⁽²⁾e²`

which is `coefficients` exactly. `test_chi2.jl` checks it the way the algebra is stated —
the same physical answer comes back from a run at `E_ref = 4`, `P_ref = ε₀` after
multiplying by `P_ref E_ref` — rather than by restating the formula, so a cubic `polscale`
would fail it by a factor of `E_ref`.

### Envelope phase bookkeeping

Unchanged. `Chi2Env` still forms `cp = ½e^{+iω₀t}` and `cm = e^{-iω₀t}` from the same
carrier array and combines the sum- and difference-frequency products exactly as before;
the docstring's warning about the reference frame (an envelope linear operator which
subtracts `β0` rather than `β1·ω` puts a spurious `β0 - β1ω₀` mismatch into every χ⁽²⁾
process) is kept verbatim. `0.5*C[i]` is written `C/2`, which is the same value for every
input and avoids promoting a `ComplexF32` to `ComplexF64` on the `Float64` literal.

### Documentation

`docs/src/developer/device_model.md` gains "The χ⁽²⁾ responses" under the response
protocol: the layout, the static matrices, why `kind` is declared unconditionally, why the
kernel converts the matrices a second time, and that the responses are tested on a device
block rather than through a propagation. The existing deviation note about the two
component broadcasts is updated from "after `gpu/15`" to what is now the case.
`docs/src/gpu.md`'s "What runs where" no longer lists χ⁽²⁾ among the responses which run
on the host, and says instead that the responses have kernels but the free-space
transforms which use them do not yet.

## The price of the layout, measured in what it costs

`gpu/12`'s deviation 1 (GPU_PLAN.md §11) is that the vector form is two fused component
broadcasts, not one `reinterpret`ed `SVector{2}` write — Luna's buffers are
`(nt, npol, ncols)` and the polarisation index is the slow axis. The consequence for a
genuinely coupled response, which the Kerr vector forms are not and χ⁽²⁾ is: **the shared
intermediates are evaluated once per output component.** For these responses that is the
lab-to-crystal rotation (a 3×3 matrix-vector product) and the six contracted products.
Dead-code elimination removes the unused component of the returned `SVector`, and with it
the row of the 3×6 contraction which produced it, but not the work the two components
share.

It is still far less work than the old path, which wrote three heap vectors per time
sample; the per-call cost is in "Measurements" below. The alternative — a `Batched` kind
with a transposed scratch buffer — buys back one rotation per sample at the cost of a
field-sized buffer, a transpose and a second pass over the block, and is not obviously
faster. It is not taken here.

## Deviations from GPU_PLAN.md and the brief

1. **`kind` is declared unconditionally, not only for `Val{2}`.** The brief says to
   declare `kind(::Chi2Field, ::Val{2}) = VectorPointwise()`. Taken literally that leaves
   `kind(r)` at its `Columnwise()` default, and `gpu/12` derives
   `device_capable(r) = !(kind(r) isa Columnwise)` — so the responses would report
   themselves as having no device kernel, which is the one thing this branch changes.
   `kind(::Chi2Field) = VectorPointwise()` gives `kind(r, Val(2))` through `gpu/12`'s
   default and makes `Val(1)` an error through the check `gpu/12` added in its review round
   1 (finding 2), which is what the brief asks for ("two-component only; a scalar block
   must error"). `gpu/12`'s own comment on that check calls the unconditional declaration
   "the natural way to write a two-component-only χ⁽²⁾ response", so this is what it
   anticipated.
2. **`rescale` converts the crystal matrices *and* the kernel converts them again.**
   `gpu/12`'s structural check exempts `isbits` static matrices, so a `Chi2Field` would
   pass through `rescale` untouched and hand a `Float64` `SMatrix` to a Metal kernel. The
   `rescale` methods therefore convert them. The kernel converts again — 27 host scalars
   per right-hand side, the identity on a rescaled response — so that a low-level caller
   who assembles `Et_to_Pt!` by hand without `rescale` still gets a `Float32` kernel. It is
   belt and braces on purpose; GPU_PLAN.md §4.2 rule 3 is the rule a stray `Float64`
   breaks, and a device kernel is where it breaks loudly.
3. **The two χ⁽²⁾ regression cases are not bit-identical.** See below. The brief asks to
   prefer an association which keeps them identical "if that costs nothing"; it costs FMA.

## The regression gate

Baseline generated with `test/regression/generate.jl 7d72431c`, and again against
`gpu/int-A` (`782f55d1`, the merge base) for the cumulative view. `julia --project=. -t 1`,
`set_fftw_mode(:estimate)`, `set_fftw_threads(1)`, `BLAS.set_num_threads(1)`; Apple M1 Pro,
Julia 1.13.0, with two other agents' jobs on the same machine.

**460 pass, 0 fail against both baselines. Step counts unchanged.**

Every case is `0.000e+00` in both modes and both classes except these two:

| case | mode | `Eω` diff | `Eω` tol | stats diff | worst quantity |
|---|---|---:|---:|---:|---|
| `free2d_field_chi2` | fixed | 7.683e-16 | 1.0e-12 | 0.000e+00 | `Eω` [component 2/2, save 8/11] |
| `free2d_field_chi2` | adaptive | 4.841e-14 | 7.4e-12 | 0.000e+00 | `Eω` [component 2/2, save 2/11] |
| `free2d_env_chi2` | fixed | 8.348e-16 | 1.0e-12 | 0.000e+00 | `Eω` [component 1/2, save 11/11] |
| `free2d_env_chi2` | adaptive | 3.151e-14 | 5.0e-12 | 0.000e+00 | `Eω` [component 2/2, save 2/11] |

**Largest difference over all cases and modes: `4.841e-14`**, against `7d72431c` and
against `782f55d1` alike — the same four numbers to every digit, which is expected, since
`gpu/12` was itself `0.000e+00` against `gpu/int-A`.

### Why they moved, and why they are not made identical

The difference is one thing: **`StaticArrays` contracts its matrix-vector product with
`muladd`, and `LinearAlgebra.mul!` into a `Vector` does not.** Measured directly on the
BBO rotation, 20 000 random samples:

| comparison | elements differing |
| --- | ---: |
| `mul!(::Vector, ::RotMatrix3, ::Vector)` vs `RotMatrix3 * SVector` | 11055 / 60000 |
| `mul!(::Vector, ::SMatrix{3,6}, ::Vector)` vs `SMatrix{3,6} * SVector` | 26591 / 60000 |

and on a single sample the `SVector` product agrees bit-for-bit with an explicit `muladd`
chain, while `mul!` agrees bit-for-bit with an explicit left-associated sum of products,
with or without the leading `zero`. So the association and the order of summation are
already the old ones; only the FMA contraction is new. GPU_PLAN.md §4.2 rule 2 names "FMA
contraction" among the rounding-level differences the gate allows, and the plan describes
the χ⁽²⁾ products as "one `SMatrix` broadcast".

Making them identical would mean writing the two contractions out as explicit
non-`muladd` sums. That is not free: it gives up FMA on every backend — accuracy as well
as throughput, since an FMA product is the more accurate of the two — and it would have
to be a hand-written 3×6 expression rather than an `SMatrix` product. The measured
difference is two orders of magnitude inside the tolerance the sensitivity study derived,
and the branch keeps the `SMatrix` product.

## Tests

`julia --project=<env> -t 1`, `Luna.set_fftw_mode(:estimate)`,
`Luna.set_fftw_threads(1)`, `BLAS.set_num_threads(1)`. Apple M1 Pro, Julia 1.13.0, with
two other agents' jobs on the same machine throughout.

### `test/test_chi2.jl` (CPU)

The eight existing testsets are unchanged and still pass, including the two which read
`c.χ2`, `c.toCrystal` and `c.toLab` and the two which call `field_products!`/
`env_products!` with their own vectors. One testset is added:

| testset | assertions |
| --- | ---: |
| the eight existing ones | 16 |
| **"Chi2 response protocol"** (new) | 26 |

**42 pass, 0 fail** (9 testsets, 5 s), up from 16.

The new testset covers, for `Chi2Field` and `Chi2Env` alike: the kind and
`device_capable`; the fused block path against the columnwise call operator column by
column on an `(nt, 2, ncols)` block, as **exact equality** (`maximum(abs, P - Pref) == 0`,
which also treats `+0.0` and `-0.0` as equal); a one-component block being refused; the
density being ignored; the unit scaling round trip described above; and `rescale`
converting the crystal matrices and the carrier phase. Plus the carrier-length check.

### `test/test_device.jl` (JLArray)

Run in an environment with `JLArrays` added, which is what `Pkg.test()` gives it.

| testset | assertions |
| --- | ---: |
| **χ⁽²⁾ responses on JLArray** (new) | 12 |
| **a user χ⁽²⁾ closure through HostResponse on JLArray** (new) | 7 |
| the eighteen existing ones | 263 |

**282 pass, 0 fail** (20 testsets, ~100 s), up from 263 on the base branch.

Both new testsets work on an `(nt, 2, ncols)` block through `NonlinearRHS.Et_to_Pt!`
rather than through a propagation, because the free-space transforms which use these
responses are host-only until Group E (the brief says so, and `Luna.setup` would refuse
them on a device anyway).

| comparison | max relative difference in the block |
| --- | ---: |
| `Chi2Field` alone, JLArray vs host | < 1e-10 (tested) |
| `Chi2Field` + vector `KerrField`, fused, JLArray vs host | < 1e-10 |
| `Chi2Env` alone / with vector `KerrEnv`, JLArray vs host | < 1e-10 |
| `Chi2Field` + a user χ⁽²⁾ closure through `HostResponse`, JLArray vs host | < 1e-10 |

### `test/test_metal.jl` (hardware)

Apple M1 Pro, Metal.jl, environment built per `docs/src/gpu.md`.

| testset | assertions |
| --- | ---: |
| **χ⁽²⁾ responses on Metal** (new) | 9 |
| **a user χ⁽²⁾ closure through HostResponse on Metal** (new) | 5 |
| the twelve existing ones | 194 |

**208 pass, 0 fail** (14 testsets), up from 194 on the base branch.

Measured separately, by a standalone script with its own parameters (`nt = 1024`,
`ncols = 4`, BBO at the type-I angle, `E_ref = 1024`, `P_ref = ε₀`), not inferred from the
assertions:

| comparison | relative difference |
| --- | ---: |
| `Chi2Field`, scaled `Float64` vs unscaled `Float64` | **0** |
| `Chi2Field`, Metal vs CPU `Float32` | **0** |
| `Chi2Field`, Metal vs CPU `Float64` | 1.4e-7 |
| `Chi2Field`, CPU `Float32` vs CPU `Float64` | 1.4e-7 |
| `Chi2Env`, scaled `Float64` vs unscaled `Float64` | **0** |
| `Chi2Env`, Metal vs CPU `Float32` | **0** |
| `Chi2Env`, Metal vs CPU `Float64` | 3.2e-7 |
| `Chi2Env`, CPU `Float32` vs CPU `Float64` | 3.2e-7 |
| `Chi2Field` + user closure through `HostResponse`, Metal vs CPU `Float32` | 5.3e-8 |
| same, Metal vs CPU `Float64` | 1.1e-7 |
| same, CPU `Float32` vs CPU `Float64` | 1.1e-7 |
| the closure's own contribution (with vs without it) | 3.2e-1 |

The two zeros in the middle are the same kind as `gpu/12`'s for the vector Kerr forms and
are genuine: a vector-pointwise `Et_to_Pt!` is a pure elementwise `Float32` broadcast with
no reduction and no FFT, and both backends contract the same products with the same FMAs
in the same order. The host side is an explicit `DeviceSpec(Array, Float32)`, the device
side an `MtlArray`, and the Metal testset separately asserts the result is not all zero.
The `Float32`-against-`Float64` rows are the ordinary `Float32` rounding of the
contraction; Metal adds nothing to it. The last row says the closure is a real
contribution rather than a rounding-level one, which is what makes the `HostResponse` row
above mean something: its coefficient is of the order of `ε₀χ⁽²⁾` for BBO by construction.

The first row of each pair — scaled against unscaled `Float64`, exactly zero — is the
`polscale(scaling, 2)` algebra checked on hardware-sized data as well as in
`test_chi2.jl`.

### `test/test_freespace.jl` (CPU)

The "BBO SHG: real vs envelope" testset is the χ⁽²⁾ propagation, and the one place both
responses are run through a real free-space transform with `idcs`.

| file | result | time |
| --- | ---: | ---: |
| `test_freespace.jl` | **77 pass, 0 fail** ("BBO SHG: real vs envelope" 8/8) | 290 s |
| `test_polarisation.jl` | 15 pass, 0 fail | 0.2 s |

## Measurements

### Per-call cost of the χ⁽²⁾ block

`NonlinearRHS.Et_to_Pt!` on a `(2048, 2, 32)` block with one response and `idcs`, minimum
of 50 calls after warm-up, same script on both commits, `-t 1`, Apple M1 Pro with other
jobs on the machine:

| response | base (`7d72431c`) | this branch | allocation, base | allocation, branch |
| --- | ---: | ---: | ---: | ---: |
| `Chi2Field` | 2.897 ms | **0.588 ms** | 496 B | **0 B** |
| `Chi2Env` | 3.744 ms | **1.845 ms** | 496 B | **0 B** |

4.9× and 2.0× on the CPU, and nothing allocated per call. This is with the shared
intermediates evaluated twice, once per output component; the old path wrote three heap
vectors per *time sample*, which is what dominated it.

No propagation benchmark is reported: the free-space transforms around these responses are
unchanged, and on this machine, with two other agents running, the difference in a whole
`free2D_bbo` propagation is inside the noise.

### Documentation build

`include("docs/make.jl")` reports exactly the pre-existing unresolved `@ref`s the base
branch does and none from this branch.

## Known gaps

- **A χ⁽²⁾ propagation still runs on the host as a whole.** The responses have kernels;
  `TransFree`/`TransFree2D`/`TransRadial` do not, and `Luna.setup` refuses a device for
  them. Group E (`gpu/20-radial-device`, `gpu/21-free-device`) is what closes this, and
  only then does "χ⁽²⁾ examples on the GPU" become an end-to-end claim rather than a
  claim about the response. The brief says so explicitly and the tests are written that
  way.
- **The shared intermediates are evaluated once per output component** (the rotation and
  the six products). A `Batched` kind with a transposed scratch buffer would avoid it; it
  is not obviously faster and is not done here.
- **`Chi2Field`/`Chi2Env` still keep `χ2` and `toLab`, which the response does not use.**
  They are read by `test_chi2.jl` and by `examples/.../chi2_benchmarking.jl`, so they
  stay. They cost nothing — the struct is `isbits` and never enters a kernel whole.
- **`import LinearAlgebra: mul!` in `Nonlinear.jl` is now unused.** The χ⁽²⁾ responses
  were its last users (`ldiv!` and `MArray` on the same two import lines were already
  unused before this branch). It is left alone: removing it is a drive-by cleanup on a
  line `gpu/14-raman` is likely to need when its batched response gains an FFT, and it
  would conflict for nothing.
- **Gas χ⁽²⁾ mixtures, and a χ⁽²⁾ response with a z-dependent tensor**, are not addressed;
  neither exists in Luna.
- **`examples/low_level_interface/freespace/chi2_benchmarking.jl` is not updated.** It
  profiles the old columnwise call operator, which still works and is still what it
  measures; the numbers it produces are now the ones in the table above rather than the
  old ones. Updating it is outside this branch's scope.

## Open questions for review

1. The unconditional `kind` declaration (deviation 1). The alternative is to declare
   `kind(::Chi2Field, ::Val{2})` as the brief literally says and to override
   `device_capable` separately, which puts the two back in a position to drift apart —
   which is the thing `gpu/12` derived `device_capable` to prevent.
2. The double conversion of the crystal matrices (deviation 2). It is 27 host scalars per
   right-hand side against a class of error which only shows up on real hardware. Is the
   `rescale` conversion alone enough, given that `gpu/12`'s structural check deliberately
   exempts `isbits` static matrices and so cannot catch a missing one?
3. `field_products!`/`env_products!` are kept with their signatures and now delegate to
   the kernel's expression. They are documented, so they are public; but nothing in Luna
   calls them any more except `test_chi2.jl`. Keep, or deprecate?
4. The χ⁽²⁾ regression movement is 4.8e-14 against a 7.4e-12 tolerance and is FMA, not
   association. Is that acceptable, or should the contraction be written out to keep the
   two cases bit-identical at the cost of FMA?

## Conflicts to expect at `gpu/int-D`

`docs/src/gpu.md` lists "plasma, Raman and χ⁽²⁾" in two places; this branch removes χ⁽²⁾
from both, and `gpu/13-plasma` and `gpu/14-raman` will remove their own from the same two
sentences. Textual conflict, trivial to resolve. `src/Nonlinear.jl` is touched by all
three branches but in three disjoint regions (the χ⁽²⁾ section is between `KerrEnvTHG` and
`PlasmaCumtrapz`). `test_device.jl` and `test_metal.jl` each gain testsets at different
points in the file; `test_device.jl`'s top-of-file helper responses are the one place all
three add lines close together.

## Attribution

The commits carry `Co-Authored-By: Claude Fable 5.1`, which is what `briefs/COMMON.md`
specifies and what `PR_12-response-traits.md` (§12 of its review round 2) records as the
project's convention from that point on. The model which wrote them is Claude Opus 5
(1M context); this session's harness specifies that name instead. The session line is the
same either way.
