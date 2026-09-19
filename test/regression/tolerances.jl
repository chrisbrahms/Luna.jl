#= Per-case, per-mode tolerances for the regression gate (`test/test_regression.jl`).

   Each number is 100x the one-ulp sensitivity measured by
   `test/regression/sensitivity.jl`, with a floor of 1e-12. The sensitivity is the response
   of the case to multiplying the initial frequency-domain field by `1 + eps()`, measured
   with the metric the gate uses and maximised over `Eω`, `z` and every statistic. It is
   the floor below which no reordering of floating-point operations can be expected to
   stay.

   Regenerate with

       julia --project=<worktree> -t 1 test/regression/sensitivity.jl

   which prints the per-quantity breakdown and a `TOLERANCES` literal to paste in below.

   Reading the numbers:

   - In the `:fixed` mode the sensitivity is at rounding level (1e-15) for every case but
     two, so the tolerance is the 1e-12 floor. This is the mode that makes a difference
     attributable to the operation that produced it.
   - In the `:adaptive` mode the largest quantity is almost always `stats/dz`, followed by
     `stats/z`. The step-size controller's accept/reject decisions and its PI update
     respond to a one-ulp change far more strongly than the field does, and the response
     compounds over the steps it takes to ramp `init_dz` up to `max_dz`. `Eω` itself is
     three to six orders of magnitude tighter than the per-case number in every adaptive
     case. The adaptive tolerances are therefore loose, and for
     `modeavg_env_kerr`/`modeavg_env_raman` loose enough that the adaptive run is only a
     crash test. Read the table `test/test_regression.jl` prints, not just its pass/fail:
     the `:fixed` numbers are the gate.
   - Two `:fixed` cases are above the floor:
     - `modeavg_field_plasma`, 2.5e-14, from `stats/peak_ionisation_rate`. The PPT rate is
       exponential in the field amplitude, so it amplifies a one-ulp input change by about
       an order of magnitude.
     - `multimode_field_plasma`, 2.9e-09, from
       `stats/transverse_integral_error_rel`/`_abs`. These are `HCubature`'s own error
       estimates for the modal overlap integral. The adaptive quadrature subdivides
       differently when the integrand changes in its last bits, so the error estimate moves
       by much more than the integral does. `stats/mode_reconstruction_error` (2.7e-13) and
       everything else in that case are at rounding level.
=#
module RegressionTolerances

export TOLERANCES, tolerance

"""
    TOLERANCES

`Dict` from case name to a `Dict` from run mode (`:fixed`, `:adaptive`) to the tolerance on
the normalised maximum difference, `maximum(abs, Δ)/maximum(abs, baseline)`.

Measured on an M1 Pro, Julia 1.13.0, `-t 1`, `set_fftw_mode(:estimate)`,
`set_fftw_threads(1)`, `BLAS.set_num_threads(1)`, FFTW wisdom disabled.
"""
const TOLERANCES = Dict{String, Dict{Symbol, Float64}}(
    "modeavg_field_kerr"     => Dict(:fixed => 1.0e-12, :adaptive => 8.2e-04),
    "modeavg_field_plasma"   => Dict(:fixed => 2.5e-12, :adaptive => 1.2e-05),
    "modeavg_field_raman"    => Dict(:fixed => 1.0e-12, :adaptive => 1.6e-06),
    "modeavg_field_mixture"  => Dict(:fixed => 1.0e-12, :adaptive => 1.4e-08),
    "modeavg_env_kerr"       => Dict(:fixed => 1.0e-12, :adaptive => 7.3e-01),
    "modeavg_env_raman"      => Dict(:fixed => 1.0e-12, :adaptive => 5.9e-01),
    "gnlse_sech"             => Dict(:fixed => 1.0e-12, :adaptive => 2.8e-08),
    "multimode_field_plasma" => Dict(:fixed => 2.9e-07, :adaptive => 7.7e-06),
    "radial_field_kerr"      => Dict(:fixed => 1.0e-12, :adaptive => 1.2e-08),
    "radial_env_kerr"        => Dict(:fixed => 1.0e-12, :adaptive => 8.2e-09),
    "free3d_env_kerr"        => Dict(:fixed => 1.0e-12, :adaptive => 1.0e-12),
    "free2d_field_chi2"      => Dict(:fixed => 1.0e-12, :adaptive => 7.4e-12),
    "gradient_field_kerr"    => Dict(:fixed => 1.0e-12, :adaptive => 1.6e-04),
    "taper_field_kerr"       => Dict(:fixed => 1.0e-12, :adaptive => 3.3e-04),
    "modeavg_field_legacy"   => Dict(:fixed => 1.0e-12, :adaptive => 5.0e-06),
)

"The tolerance used for a case that is not in [`TOLERANCES`](@ref)."
const DEFAULT_TOLERANCE = 1.0e-12

"""
    tolerance(name, mode)

The regression tolerance for case `name` in run mode `mode`, or
[`DEFAULT_TOLERANCE`](@ref) if the case or the mode is not listed in
[`TOLERANCES`](@ref).
"""
function tolerance(name::AbstractString, mode::Symbol)
    haskey(TOLERANCES, name) || return DEFAULT_TOLERANCE
    get(TOLERANCES[name], mode, DEFAULT_TOLERANCE)
end

end
