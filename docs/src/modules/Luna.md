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

## Devices and precision

Where a propagation runs and in what units. See [Running on a GPU](../gpu.md) for the
user-facing description and
[The device and precision model](../developer/device_model.md) for the internals.

```@docs
Luna.DeviceSpec
Luna.device
Luna.device_request
Luna.set_device
Luna.register_device!
Luna.DeviceHooks
Luna.device_functional
Luna.device_synchronize
Luna.device_reclaim
Luna.device_memory_status
Luna.alloc
Luna.todevice
Luna.tohost
Luna.scalar
Luna.upload_like
Luna.mask_like
Luna.assert_resident
Luna.all_resident
Luna.UnitScaling
Luna.unitscaling
Luna.polscale
Luna.GridVectors
Luna.gridvectors
Luna.HostMirror
Luna.upload!
Luna.log_device
Luna.setup
Luna.runscaling
Luna.ScaledOutput
Luna.needs_host_y
Luna.needs_host_cache
```

```@docs
Utils.Backend
Utils.CPUBackend
Utils.DeviceBackend
Utils.backend
Utils.isdevice
Utils.plan_ft
Utils.plan_ift
Utils.iplan
Utils.iscale
```
