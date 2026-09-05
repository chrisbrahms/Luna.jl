# Luna.jl
The core of the package: `setup` assembles a simulation from its composable pieces and
`run` propagates it. See [The simple interface](@ref "The simple interface") for the
high-level entry points that build these automatically.

## Setting up and running a simulation
```@docs
Luna.setup
Luna.run
```

## Global settings
```@docs
Luna.settings
Luna.set_fftw_mode
Luna.set_fftw_threads
```
