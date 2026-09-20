#= Metal glue for Luna's device backend.

   Luna's per-step code is generic (broadcasts, planned FFTs applied with `mul!`, and
   reductions over anything following the GPUArrays interface), so this extension is
   small by design: it registers the backend and provides the four vendor operations
   which cannot be written generically.

   They are installed as function-valued hooks from `__init__` rather than defined as
   methods: a method with the same signature as a stub in Luna itself would *overwrite*
   it, which Julia rejects during precompilation.

   Metal refuses `Float64` arrays and its kernel compiler rejects any `double` which
   survives optimisation, so the registered spec is `Float32`. =#
module LunaMetalExt

import Luna
import Metal

_functional() = Metal.functional()

_synchronize() = Metal.synchronize()

# Metal's allocator keeps a pool which garbage collection alone does not return; the
# incremental collection first makes freshly dead arrays eligible.
_reclaim() = (GC.gc(false); nothing)

function _memory_status()
    Metal.functional() || return nothing
    dev = Metal.device()
    total = Int(dev.recommendedMaxWorkingSetSize)
    used = Int(dev.currentAllocatedSize)
    (total - used, total)
end

function __init__()
    Luna.register_device!(:metal, Luna.DeviceSpec(Metal.MtlArray, Float32);
                          functional = _functional,
                          synchronize = _synchronize,
                          reclaim = _reclaim,
                          memory_status = _memory_status)
    return nothing
end

end
