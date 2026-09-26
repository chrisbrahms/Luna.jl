# Wall time of the cases in cases.jl with Luna.set_threaded_broadcasts off and on, in one
# process, best of two each; the outputs are compared for bit-identity.
#   FFTWT=1 BLAST=1 julia --project=. -t 4 benchmark/threaded/speed.jl
# CASES=modeavg,radial selects cases.
import LinearAlgebra
using Luna, Printf
Luna.set_fftw_mode(:estimate); Luna.set_fftw_wisdom(false)
fftwt = parse(Int, get(ENV, "FFTWT", "1")); Luna.set_fftw_threads(fftwt)
LinearAlgebra.BLAS.set_num_threads(parse(Int, get(ENV, "BLAST", "1")))
include(joinpath(@__DIR__, "cases.jl"))
for name in Symbol.(split(get(ENV, "CASES", "modeavg,modal,radial,modeavg_long,radial_big"), ","))
    f = getfield(Main, name)
    res = Dict{Bool, Any}()
    for on in (false, true, false, true)
        Luna.set_threaded_broadcasts(on)
        t = @elapsed out = f()
        E = name in (:radial, :radial_big) ? copy(out.data["Eω"]) : copy(out["Eω"])
        res[on] = (min(t, get(res, on, (Inf,))[1]), E)
    end
    @printf("%-8s -t %d fftw %d: off %6.2f s  on %6.2f s  speed-up %.2f  identical %s\n",
            name, Threads.nthreads(), fftwt, res[false][1], res[true][1],
            res[false][1]/res[true][1], res[false][2] == res[true][2])
    flush(stdout)
end
