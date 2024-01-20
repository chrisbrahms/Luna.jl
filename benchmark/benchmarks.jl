using Luna
using BenchmarkTools
import Logging: with_logger, NullLogger

SUITE = BenchmarkGroup()

SUITE["mode_average"] = BenchmarkGroup()
SUITE["mode_average"]["with_stats"] = @benchmarkable prop_capillary(
    125e-6, 3, :He, 1;
    λ0=800e-9, τfwhm=10e-15, energy=120e-6,
    λlims=(150e-9, 4e-6), trange=0.5e-12)
SUITE["mode_average"]["no_stats"] = @benchmarkable prop_capillary(
    125e-6, 3, :He, 1;
    λ0=800e-9, τfwhm=10e-15, energy=120e-6,
    λlims=(150e-9, 4e-6), trange=0.5e-12,
    statistics=false)

SUITE["1_mode"] = BenchmarkGroup()
SUITE["1_mode"]["with_stats"] = @benchmarkable prop_capillary(
    125e-6, 3, :He, 1;
    λ0=800e-9, τfwhm=10e-15, energy=120e-6,
    λlims=(150e-9, 4e-6), trange=0.5e-12,
    modes=1)
SUITE["1_mode"]["no_stats"] = @benchmarkable prop_capillary(
    125e-6, 3, :He, 1;
    λ0=800e-9, τfwhm=10e-15, energy=120e-6,
    λlims=(150e-9, 4e-6), trange=0.5e-12,
    modes=1, statistics=false)

SUITE["4_mode"] = BenchmarkGroup()
SUITE["4_mode"]["with_stats"] = @benchmarkable prop_capillary(
    125e-6, 3, :He, 1;
    λ0=800e-9, τfwhm=10e-15, energy=120e-6,
    λlims=(150e-9, 4e-6), trange=0.5e-12,
    modes=4)
SUITE["4_mode"]["no_stats"] = @benchmarkable prop_capillary(
    125e-6, 3, :He, 1;
    λ0=800e-9, τfwhm=10e-15, energy=120e-6,
    λlims=(150e-9, 4e-6), trange=0.5e-12,
    modes=4, statistics=false)

##
# with_logger(NullLogger()) do
#     tune!(SUITE)
# end

# results = with_logger(NullLogger()) do
#     results = run(SUITE; verbose=true, seconds=60)
# end



