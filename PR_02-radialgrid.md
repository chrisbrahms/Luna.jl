# `Grid.RadialGrid`: a Luna-owned transverse grid for radially symmetric propagation

Base branch: `evanescent` (fdf8dbe3). Implements GPU_PLAN.md §4.14 on the CPU only.

## Motivation

Radial free-space propagation used a `Hankel.QDHT` directly as the transverse grid. The
`QDHT` appeared in Luna's dispatch signatures (`Luna.setup`, `LinearOps.transverse_k2`,
`NonlinearRHS.TransRadial`, `Boundaries`, `Fields`), so a user had to know Hankel.jl to set
up a radial run, and Luna called into Hankel on the per-step path. Three things follow from
that which the GPU work needs to change:

- `QDHT` is immutable with concrete `Array` fields and does not survive a `Float32`
  constructor, so it cannot describe a device grid.
- Hankel transforms along `QDHT.dim` with a `permutedims` temporary per call for
  `dim ≠ 1`. `NonlinearRHS.TransRadial` already worked around this with its own transposed
  matrices; `Boundaries.RadialCollar`, `Luna.setup` and the noise setup did not.
- Metal's accelerated GEMM needs both operands in the same element type, so the transform
  matrix has to be held per consumer in the type it multiplies, not once in the grid.

`Grid.RadialGrid` owns the sample points, the transform matrices and the integration
weights. Hankel.jl is used only inside its constructor, to get the Bessel zeros, the matrix
and the scale factors.

## What changed

### `src/Grid.jl`

New `RadialGrid <: SpaceGrid`, host-side and `Float64`:

```julia
struct RadialGrid <: SpaceGrid
    R::Float64; K::Float64; N::Int; order::Int
    r::Vector{Float64}; k::Vector{Float64}
    Tfwd::Matrix{Float64}; Tbwd::Matrix{Float64}
    wr::Vector{Float64}; wk::Vector{Float64}
end
```

- `RadialGrid(R, N; order=0)` builds a `Hankel.QDHT{order, 1}(R, N)` internally and copies
  `r`, `k`, `R`, `K`, `scaleR`, `scaleK`, forming
  `Tfwd = transpose(q.T) .* q.scaleRK` and `Tbwd = transpose(q.T) ./ q.scaleRK` — the
  matrices `NonlinearRHS.TransRadial` built for itself before (old `NonlinearRHS.jl:566-567`).
- `RadialGrid(q::Hankel.QDHT)` converts an existing transform, with `@warn maxlog=1` that
  the `QDHT` is deprecated as a Luna grid. `q.dim` is ignored.
- `Grid.HankelTransform` is an alias for `Hankel.QDHT`, so no other module imports Hankel.
- `Grid.TransverseGrid` is the union used in the free-space dispatch signatures.

Transform and integration functions, all along the **last** dimension of an array of any
rank:

- `to_kspace!(out, rg, A)`, `to_rspace!(out, rg, A)`: one `mul!` on `reshape(A, :, rg.N)`
  via `radial_matmul!`. No allocation unless `out === A`, in which case a temporary copy of
  the input is made. `Hankel.T` is symmetric, so this right multiplication is the same
  transform as Hankel's left multiplication along the radial axis, without the
  `permutedims`. Verified numerically against `Hankel.mul!`/`Hankel.ldiv!` for 1-D, 2-D and
  3-D inputs in `test/test_radialgrid.jl`.
- `to_kspace(rg, A; dim=ndims(A))`, `to_rspace(rg, A; dim=ndims(A))`: allocating forms; a
  `dim` other than the last is handled with a `permutedims` pair, for post-processing of
  `(ω, pol, r, z)` arrays.
- `integrate_r(rg, A; dim=ndims(A))`, `integrate_k(rg, A; dim=ndims(A))`: weighted sums
  over `wr`/`wk`, written as one matrix-vector product. They reproduce
  `Hankel.integrateR`/`integrateK` for the shapes Luna uses, including the 1-D and 2-D
  `(t, r)` arrays `Fields.energyfuncs` receives. Unlike `Hankel.integrateR` they drop the
  integrated dimension instead of leaving it as a singleton.
- `kperp2(rg)`, `onaxis(rg, Ak; dim)`, `symmetric(rg, A)`, `rsymmetric(rg)`: replacements
  for `q.k .^ 2`, `Hankel.onaxis`, `Hankel.symmetric`, `Hankel.Rsymmetric`. `onaxis` and
  `symmetric` are defined for 0-order grids only, as in Hankel.
- `to_dict(rg)` stores `R`, `N` and `order` only; `from_dict` rebuilds the matrices;
  `validate(::RadialGrid)` is added and is called by `from_dict`.
- Minimal `to_dict` methods for `FreeGrid` (`x`, `y`, `kx`, `ky`) and `Free2DGrid`
  (`x`, `kx`), so that `Luna.run` can describe any transverse grid in the output. Neither
  has a `from_dict` (their constructors take a window factor that `to_dict` does not
  record); `RadialGrid` does.

### Other `src/` modules

- `Luna.setup`: the two radial methods take a `Grid.RadialGrid`; `q * Eω` becomes
  `Grid.to_kspace`. A `setup(grid::Grid.TimeGrid, q::Grid.HankelTransform, args...)` method
  converts, so scripts written against the old interface still run. `doinputs_fs!` takes
  `Grid.TransverseGrid`.
- `Luna.run` writes `Grid.to_dict(Boundaries.spacegrid(transform))` to the output under the
  group `"spacegrid"` when the transform is a free-space one. Previously only the time grid
  was written.
- `LinearOps`: `transverse_k2(::Grid.RadialGrid)`, a converting method for
  `Grid.HankelTransform`, and `Grid.TransverseGrid` in the `make_const_linop`/`make_linop`
  signatures in place of the explicit `Union`.
- `NonlinearRHS`: `TransRadial.QDHT` becomes `TransRadial.rgrid`; the constructor takes the
  grid's matrices as `convert(Matrix{TT}, rgrid.Tfwd/Tbwd)` copies in the element type it
  multiplies. The noise-field `ldiv!` becomes a `radial_matmul!` with that same `Tbwd`.
  `norm_radial`/`const_norm_radial` take a `RadialGrid` and have converting methods.
  The per-step transform body is unchanged.
- `Boundaries`: `spacegrid(::TransRadial)` returns the `rgrid`; `kprofile`, `rprofile`,
  `spatialcollar` and `log_setup` dispatch on `Grid.RadialGrid`. `RadialCollar` holds the
  grid plus its own `Tfwd`/`Tbwd` in the complex spectral type and applies them with
  `Grid.radial_matmul!` in place of `Hankel.ldiv!`/`mul!`.
- `Fields`: `transform`, `prop!` and `energyfuncs` take a `RadialGrid` (with converting
  methods for a `QDHT`) and use `Grid.to_kspace`, `Grid.kperp2`, `Grid.integrate_r/k`.

### Tests and examples

- New `test/test_radialgrid.jl` (126 tests), registered in `runtests.jl`: construction
  against a `QDHT`, conversion and its warning, transform agreement with `Hankel.mul!` and
  `Hankel.ldiv!` for 1-D/2-D/3-D and for `Float64` and `ComplexF64`, in-place and aliased
  forms, transforms and integrals along a leading dimension, dimension errors, integrals
  against `Hankel.integrateR`/`integrateK` and an analytical Gaussian, Parseval,
  `onaxis`/`symmetric`/`rsymmetric` against Hankel, and `to_dict`/`from_dict`. This is the
  only file in `test/` which constructs a `QDHT`.
- `test_freespace.jl`, `test_boundaries.jl`, `test_linops.jl`, `test_fields.jl` and
  `test_modes.jl` use `Grid.RadialGrid`. Two testsets are added to `test_freespace.jl`: the
  transverse grid written to the output round-trips through `Grid.RadialGrid`, and the
  radial noise field is transformed to real space with the transform's own matrix.
- `examples/low_level_interface/freespace/radial*.jl` use `Grid.RadialGrid` and the
  `Grid.` post-processing functions. `radial_hcf.jl` also had its post-processing indices
  corrected: it indexed `Eω` as `(ω, r, z)`, but `Luna.setup` has added a polarisation axis
  for some time, so the example could not have run as written. The unused
  `import Luna: Hankel` was removed from `free2D_bbo.jl` and `free3D_bbo.jl`.

`grep -rn "Hankel" src/` shows `Hankel.` only in `src/Grid.jl`, plus `import Hankel` in
`src/Luna.jl` (kept so that `Luna.Hankel` still resolves for existing scripts) and the word
"Hankel" in comments and docstrings. No `QDHT` is constructed anywhere in `examples/` or
`docs/`.

## Numerical change and its measurement

Hankel's `permutedims` + left GEMM is replaced by a right GEMM on a reshape in three places
that were not already doing it: `Boundaries.RadialCollar`, `Luna.setup`'s input transform,
and the `TransRadial` noise setup. Folding the scalar `scaleRK` into the matrix instead of
applying it after the product also changes the order of operations. Both change rounding.
The per-step nonlinear transform is untouched.

Measured by running the radial cases of `test/test_freespace.jl` (field and envelope grid,
constant and gradient pressure, one and two polarisations, THG on and off: 12 cases, 11
saves each) on this branch and on a `evanescent` worktree, comparing
`maximum(abs, ΔEω)/maximum(abs, Eω)` at every save.

Machine: Apple M1 Pro, Julia 1.13, macOS. Settings: `julia -t 1`,
`Luna.set_fftw_mode(:estimate)`, `Luna.set_fftw_threads(1)`,
`LinearAlgebra.BLAS.set_num_threads(1)`, no shot noise, FFTW wisdom loading disabled in the
measurement script (the wisdom file is shared between concurrent runs on this machine, and
loading it would make the plans, and hence the rounding, depend on what else had run).

`maximum(abs, ΔEω)/maximum(abs, Eω)`, per case, over the 11 saves:

| case (pressure, grid, polarisations, THG) | max over saves | final save |
|---|---|---|
| const, envelope, 1 pol, no THG | 3.62e-15 | 2.90e-15 |
| const, envelope, 1 pol, THG    | 3.20e-15 | 3.20e-15 |
| const, envelope, 2 pol, no THG | 3.76e-15 | 3.36e-15 |
| const, envelope, 2 pol, THG    | 3.20e-15 | 3.20e-15 |
| const, field, 1 pol            | 3.63e-15 | 1.87e-15 |
| const, field, 2 pol            | 4.82e-15 | 4.82e-15 |
| gradient, envelope, 1 pol, no THG | 3.62e-15 | 2.90e-15 |
| gradient, envelope, 1 pol, THG    | 3.20e-15 | 3.20e-15 |
| gradient, envelope, 2 pol, no THG | 3.76e-15 | 3.36e-15 |
| gradient, envelope, 2 pol, THG    | 3.20e-15 | 3.20e-15 |
| gradient, field, 1 pol            | 3.63e-15 | 1.87e-15 |
| gradient, field, 2 pol            | 4.82e-15 | 4.82e-15 |

**Worst over all cases and all saves: 4.82e-15.** The first save (z = 0, the input field
only) is 8.7e-16 for the envelope grid and 9.7e-16 for the field grid, so the input
transform alone contributes about one ulp and the rest accumulates over 20 accepted steps.
This is well inside the 1e-12 gate and nothing needed investigating. The constant-pressure
and gradient cases agree to the last digit because the gradient case here uses a constant
refractive index; they exercise different linop code paths, not different physics.

Reproduce with the two scripts in the project scratchpad: `radial_cases.jl <out.h5>` run in
each checkout (it disables FFTW wisdom loading itself) and `compare.jl <base.h5>
<branch.h5>`. The baseline checkout used was a detached worktree of fdf8dbe3 at
`/Users/mb140/.julia/dev/Luna-gpu/baselines/base-02`.


## Tests run

All run as `julia --project=<worktree> -t 1`, with `Luna.set_fftw_mode(:estimate)`,
`Luna.set_fftw_threads(1)` and `LinearAlgebra.BLAS.set_num_threads(1)`.

| file | result | time |
|---|---|---|
| `test/test_radialgrid.jl` (new) | 126 pass | 5.8 s |
| `test/test_linops.jl` | 196 pass | 12.3 s |
| `test/test_boundaries.jl` | 182 pass | 21.9 s |
| `test/test_modes.jl` | 587 pass | 27.2 s |
| `test/test_output.jl` | 71 pass | 14.0 s |
| `test/test_fields.jl` | 179 pass | 56.1 s |
| `test/test_freespace.jl` | 77 pass | 349.6 s |

1418 tests, all passing; no failures, errors or broken tests. (Times are from a run with
several other agents' jobs on the same machine, so they are upper bounds.) `Pkg.test()` was
not run (too slow); the files above are the ones this branch touches, plus
`test_output.jl` for the new `"spacegrid"` group.


## Known gaps and deviations from GPU_PLAN.md

- The plan's §4.14 says the noise-field transform uses `to_rspace!`. It uses
  `Grid.radial_matmul!` with `TransRadial`'s own `Tbwd` instead, which is the same
  operation with the matrix in the transform's element type. Calling `to_rspace!` would
  multiply a `ComplexF64` array by the grid's `Float64` matrix for an `EnvGrid`, which is
  the mixed-type GEMM §4.14 says each consumer should avoid, and it would round differently
  from the per-step transform.
- `FreeGrid` and `Free2DGrid` get a `to_dict` but no `from_dict`: their constructors take a
  `window_factor` which is not recoverable from the fields, and the plan only asks for a
  minimal description in the output.
- The per-step `TransRadial` transform still does one `mul!` per polarisation on a `view`.
  Turning that into a single GEMM on a reshape is §4.4/§2 work for the device branch and
  would change the results, so it is not done here.
- `zmax` is untouched (`gpu/01-zmax`).
- `order ≠ 0` grids are constructed and transform correctly, but nothing else in Luna
  supports them yet; `onaxis` and `symmetric` throw a `DomainError` for them, as Hankel
  does. The field is the hook for Laguerre-Gauss modes later.
