#= CUDA glue for Luna's device backend. See `ext/LunaMetalExt.jl` for why the vendor
   operations are function-valued hooks installed from `__init__` rather than methods.

   CUDA supports `Float64`, so the registered spec is double precision: a CUDA run is
   then numerically comparable with the host run without the unit scaling. A
   single-precision CUDA run is available through
   `Luna.set_device(Luna.DeviceSpec(CUDA.CuArray, Float32))` or the `precision` keyword
   of `Luna.setup`.

   This extension is untested on hardware: there is no CUDA device on the machine this
   branch was developed on. It is written to load and register correctly; the hardware
   tests are `test/test_cuda.jl` (gated on `LUNA_TEST_CUDA=1`), which a later branch
   adds. =#
module LunaCUDAExt

import Luna
import CUDA

_functional() = CUDA.functional()

_synchronize() = CUDA.synchronize()

_reclaim() = (GC.gc(false); CUDA.reclaim())

_memory_status() = CUDA.functional() ? (CUDA.free_memory(), CUDA.total_memory()) : nothing

function __init__()
    Luna.register_device!(:cuda, Luna.DeviceSpec(CUDA.CuArray, Float64);
                          functional = _functional,
                          synchronize = _synchronize,
                          reclaim = _reclaim,
                          memory_status = _memory_status)
    return nothing
end

end
