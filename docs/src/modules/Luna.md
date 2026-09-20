# Luna.jl

The top-level module. [`Luna.run`](@ref) is the low-level entry point which drives the
propagation; most users will go through [The simple interface](../interface.md) instead.

```@docs
Luna.run
```

## Global settings

[`Luna.settings`](@ref) holds the process-wide configuration of the FFTW planner. The
planning mode and the wisdom cache both change which plan FFTW produces, and therefore the
order in which the floating-point operations of every transform are done. Turn the wisdom
cache off when a run has to be reproducible: the wisdom file is shared by every process
using the same Julia depot, so wisdom written by an unrelated run can change the plan this
one gets.

```@docs
Luna.settings
Luna.set_fftw_mode
Luna.set_fftw_threads
Luna.set_fftw_wisdom
```
