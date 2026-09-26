# GPU and threading documentation; a test tolerance fix

Base: `gpu/31-ka-kernels` (c223fea1). Plan §4.12 and §6 Group F (`gpu/32-docs`).

## Changes

- `docs/src/gpu.md` (user page): the "work in progress, as of gpu/int-E" admonition is
  replaced by a summary of what runs on a device and what stays on the host; a note that
  the CUDA extension has not been run on CUDA hardware; a new "On the CPU: threads"
  section under Performance (the threaded broadcasts from `gpu/31`, and the finding that
  FFTW's default thread count, 4 × the Julia thread count, slows mode-averaged and small
  multimode runs by up to 3×, with a measured table and the recommended settings).
- `docs/src/developer/device_model.md`: the work-in-progress note is replaced by an
  overview of the page's order; a "Host threading" section (`Utils.threaded`, the
  thresholds and why, the nesting guard, bit-identity, the profile and speed-ups,
  KernelAbstractions measured and not used); "The cost of evaluating the closure" (the
  `gpu/30` measurements behind removing the Capillary `neff` overloads); branch names
  removed from the running text; "The modal transforms" moved from after "Testing" to
  before "The output and statistics boundary", next to the other transforms; the testing
  list names `threaded_checks.jl` and `benchmark/threaded/`.
- `docs/src/modules/Luna.md`: `Luna.set_threaded_broadcasts` added to the settings docs.
- `README.md`: a short "Running on a GPU, and threads on the CPU" section.
- `src/Luna.jl`: `set_threaded_broadcasts`'s docstring no longer `@ref`s the undocumented
  `Utils.threaded`.
- The four pre-existing unresolved `@ref`s (`loadFFTwisdom`, `saveFFTwisdom` in
  `Utils.plan_ft`, `AbstractOutput` in `Output.jl`, `PhysData.crystal_internal_angle` in
  `NonlinearRHS.jl`; none of them is in the rendered API) are plain code now. They made
  `makedocs` stop before rendering on every branch of this project.
- `test/test_device.jl` (separate commit): two testsets compared the modal transform's
  embedded quadrature error estimate, host against JLArray, relative to the estimate
  itself at 1e-10. The estimate is a difference of two quadratures and, for the
  well-resolved test case, close to rounding level, so that measured operation order. Both
  failed under `Pkg.test()` on `gpu/int-E2` (5.1 % relative difference, 3 failures in
  7072) while passing whenever the file was run on its own (checked with BLAS at 8
  threads, FFTW wisdom enabled, the six preceding test files included first, and
  `--check-bounds=yes`; the trigger inside `Pkg.test()` was not isolated). They now
  compare relative to the integral the estimate is the error of.

`CLAUDE.md` is not tracked in the repository, so it is not changed here; a corrected
version was prepared separately.

## Tests

- Docs build (`docs/make.jl`): exit 0, no errors (it failed with the four `@ref` errors
  above before). Remaining warnings: page-size thresholds, a pre-existing unescaped `$` in
  `docs/src/modules/Ionisation.md`, deployment skipped.
- `test_device.jl` on its own, after the tolerance change: 1432/1432.
- Full `Pkg.test()` on this branch: PKGTEST_RESULT.
