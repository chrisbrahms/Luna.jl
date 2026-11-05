# gpu/base: the GPU project's base

`evanescent` without the modified shot-noise model. Every GPU branch merges this.

## History

This branch used to be `evanescent` + LupoLab/Luna.jl#440 ("Propagate the modified shot-noise
field under the unitary linear operator"), so that the GPU PRs could be based on a master
containing #440. The modified shot-noise model of #414 has since been dropped because its
physics is not correct: `revert-modified-noise` reverts #414's noise model on master (keeping
its Scans.jl queue fix, the .gitignore entries, the PPT test tolerance and the `rng` keyword,
which now seeds the input `Fields.ShotNoise`), and #440 will not be merged.

`evanescent` merged `revert-modified-noise`, and this branch then merged `evanescent` taking
its tree, so it now differs from `evanescent` only by this file. The #440 code, its fixes on
this branch (the one-mode `(nω, 1)` draw, the evanescent-component mask of the free-space draw)
and `test_noise.jl` are gone. `shotnoise=true` adds one-photon-per-mode noise to the input
field, as before #414.

Every branch downstream removes its own noise code when it merges this one.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_01WM2AxPkj1PCFvBcaaaanzV
