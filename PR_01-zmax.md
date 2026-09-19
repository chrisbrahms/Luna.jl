# Remove `zmax` from the grids; pass the propagation length to `Luna.run`

Base: `evanescent` (fdf8dbe3). Branch: `gpu/01-zmax`. Implements GPU_PLAN.md §4.13.

## Motivation

`Grid.RealGrid` and `Grid.EnvGrid` carried a `zmax` field that has nothing to do with the
temporal or spectral axes they describe. It forced the propagation length to be fixed when
the grid was built, made the same number appear twice (in the grid and in the output's save
grid), and stands in the way of the device work in later branches, where a grid is uploaded
to a device and should contain only axis data. The grids now describe only the time and
frequency axes and their apodisation windows; the propagation length is an argument of
`Luna.run`.

## What changed

### `src/Grid.jl`

- `RealGrid` and `EnvGrid` lose the `zmax` field.
- New constructors `RealGrid(referenceλ, λ_lims, trange; δt=1)` and
  `EnvGrid(referenceλ, λ_lims, trange; δt=1, thg=false)`. `δt` is a keyword in both;
  it used to be the fifth positional argument of `RealGrid`.
- The old positional forms `RealGrid(zmax, referenceλ, λ_lims, trange, δt=1)` and
  `EnvGrid(zmax, referenceλ, λ_lims, trange; δt, thg)` are kept as deprecated methods that
  emit `@warn ... maxlog=1` and drop `zmax`. They are unambiguous against the new methods:
  the new ones take three positional arguments (`::Real`, `::Union{Tuple,AbstractVector}`,
  `::Real`), the deprecated ones four or five with `zmax::Real` first.
  `Base.depwarn` was not used because it is silent by default.
- The keyword constructors used by `from_dict` accept and ignore `zmax=nothing`, so output
  files written before this change still load through `Processing.makegrid`/`Grid.from_dict`
  (`to_dict` iterates `fieldnames`, so their `grid` group has a `zmax` entry).

### `src/Luna.jl`

- `run` takes `zmax` as a required keyword. It is defaulted to `nothing` internally and
  checked, so a missing `zmax` gives a message naming the keyword rather than an
  `UndefKeywordError`.
- `max_dz` defaults to `zmax/2` (previously `grid.zmax/2`).
- `zmax` is passed to `Boundaries.setup` and to `RK45.solve_precon` as `tmax`.
- If the output has a `GridCondition` save grid, `run` checks `zmax ≈ save_cond.grid[end]`
  and errors if the two disagree — the propagation length is now stated in two places and
  this stops them drifting apart.
- `run` writes `zmax` to the output as a top-level entry, so post-processing can read the
  propagation length back. The write is skipped if the entry is already there, which is the
  case when an `HDF5Output` propagation is resumed from its cache.

### `src/Output.jl`

- New `Output.hasdata(o, key)`: whether an output already holds a top-level entry. Defined
  for `MemoryOutput` and `HDF5Output`; the fallback returns `false` for outputs that are
  plain callables and cannot be queried. Used only by `run` for the `zmax` write above.

### `src/Boundaries.jl`

- `reflength`, `spectral_rate`, `temporal_rate`, `log_setup` and `setup` take `zmax` as an
  argument instead of reading `grid.zmax`. In `setup` it is a positional argument between
  `z0` and `max_dz`. The `:none`/`:legacy` branch uses `min(max_dz, zmax)` as the evanescent
  reference length, as before.

### `src/Interface.jl`

- `makegrid` drops its `flength` argument; `makeoutput` takes `flength` directly instead of
  reading it off the grid.
- `prop_capillary` and `prop_gnlse` pass `zmax=flength` to `Luna.run`.
- `prop_capillary_args`/`prop_gnlse_args` are unchanged in what they return; their
  docstrings now show that the caller has to pass `zmax=flength` to `Luna.run`.
- The simple interface's signature and behaviour are unchanged. `saveargs` already records
  `flength`, so it needed no change.

### Tests and examples

- New `test/test_grid.jl` (added to `test/runtests.jl`): the new constructors, the
  deprecated ones (warning emitted, identical grid), `to_dict`/`from_dict` round trip and
  `from_dict` with a legacy `zmax` key, `run` without `zmax` erroring, the
  `zmax`/`GridCondition` consistency check, and `output["zmax"]`.
- 24 existing test files and 49 example files updated to the new signatures. No test uses a
  deprecated constructor except the ones in `test_grid.jl` that test the deprecation.
- `docs/src` needed no edit: the grid and `Luna.run` pages are `@autodocs`/`@docs` blocks, so
  they pick up the changed docstrings. This is a deviation from item 6 of the brief only in
  the sense that there was no hand-written text to change.

## Tests

Run with `julia --project=. -t 1`, `Luna.set_fftw_mode(:estimate)`,
`Luna.set_fftw_threads(1)`, `LinearAlgebra.BLAS.set_num_threads(1)`.

PLACEHOLDER_TESTS

## Examples

All 128 `.jl` files under `examples/`, `test/` and `src/` parse. A representative subset of
examples was run up to the point where they call the plotting stack:

PLACEHOLDER_EXAMPLES

## No numerical change

Three cases were run on the base commit (`fdf8dbe3`, in a worktree at
`/Users/mb140/.julia/dev/Luna-gpu/baselines/base-02`) and on this branch, in the same
environment, with `set_fftw_mode(:estimate)`, one FFTW thread, one BLAS thread, one Julia
thread, and `rng=MersenneTwister(1234)` for the shot noise:

1. mode-averaged, field-resolved, with ADK plasma (`prop_capillary`, 125 µm, 0.1 m, Ar 1 bar);
2. radial free space (low-level, `Hankel.QDHT`, 256 points, Kerr only);
3. multimode (`prop_capillary`, `modes=2`).

PLACEHOLDER_REGRESSION

## Known gaps and open questions

- `prop_capillary`/`prop_gnlse` pass `zmax=args[2]` to `Luna.run`, i.e. they rely on
  `flength` being the second positional argument of `prop_capillary_args`/`prop_gnlse_args`.
  That is how the existing `args...` forwarding is written; splitting the signature to name
  `flength` explicitly would be a larger change to the interface than this branch should
  make.
- `Output.hasdata` is new public API introduced only so that `run` can avoid a duplicate
  write when an HDF5 propagation is resumed. An alternative would be `output("zmax", zmax;
  force=true)`, which warns on `HDF5Output` and errors on `MemoryOutput`.
- The deprecated constructors warn with `maxlog=1`, so a script that builds many grids sees
  the message once. They are meant to be removed one release after this change.
