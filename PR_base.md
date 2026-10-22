# gpu/base: the GPU project's base with LupoLab/Luna.jl#440

`evanescent` (fdf8dbe3) + #440 ("Propagate the modified shot-noise field under the unitary
linear operator", head 6537a0f4 on jtravs/claude/luna-modified-noise-fix-40e7bb). Every GPU
branch merges this, so the GPU PRs can be based on a master which contains #440.

## Conflicts resolved

`evanescent` carries the free-space rework of #416 (polarisation dimension in the radial and
3-D transforms, the radial Hankel transform as two GEMMs, `to_time!`/`to_freq!` taking the
*forward* plan); #440 was written against master. Resolutions:

- `to_time!` calls keep `t.FT` (evanescent's convention); #440's `inv(t.FT)` would be wrong
  here, including in the auto-merged TransModeAvg noise path.
- **TransRadial**: evanescent's layout `(t, pol, r)` and GEMMs, with #440's noise: the noise
  is advanced every right-hand side, taken ω → t and k → r with the same `Tbwd` as the field
  (`Et_nl` doubles as the k-space scratch), and added in `Et_nl`.
- **TransFree**: evanescent's `(t, pol, x, y)` layout and `to_time!`; the noise buffers now
  have the field's shape (evanescent's precomputed noise buffers had no polarisation
  dimension, so `Et_nl = Eto + Et_noise` could not have broadcast).
- `Luna.setup` radial/3-D methods: `np > 1` polarisation argument kept, keyword `noise`.
- Interface: #440's `rtol`/`atol` and `save_noise` alongside evanescent's absorbing-boundary
  keywords (`boundary_kwargs`) and `linopkw` for the envelope THG option.

## Test added

`test_noise.jl` "Every transform sees E + noise(z)": for the radial, 3-D, mode-averaged and
multimode transforms, the transform with noise equals the noise-free one applied to
`Eω .+ noise_at(noise, z)`, at two z. #440's own tests cover only the capillary transforms;
this is the check that carries the free-space resolutions (and, later, every device path).

## Verification

test_noise 111/111 (incl. #440's Chen & Wise identity), test_freespace 68/68,
test_interface 301/301. (No regression gate on this branch; it arrives with gpu/00-harness.)

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_01PizrXcSBunUFRHcXRjCd6C
