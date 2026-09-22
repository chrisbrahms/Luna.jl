#= Per-case, per-mode, per-class tolerances for the regression gate
   (`test/test_regression.jl`).

   There are two tolerance classes: `:Eω`, the field itself, and `:stats`, the save grid `z`
   and every compared statistic. Each number is 100x the one-ulp sensitivity of that class,
   measured by `test/regression/sensitivity.jl`, with a floor of 1e-12. The sensitivity is
   the response of the case to multiplying the initial frequency-domain field by `1 + eps()`,
   measured with the metric the gate uses; it is the floor below which no reordering of
   floating-point operations can be expected to stay.

   Regenerate with

       julia --project=<worktree> -t 1 test/regression/sensitivity.jl

   which uses the gate's comparison, prints the per-quantity breakdown and emits a
   `TOLERANCES` literal to paste in below unchanged.

   Reading the numbers:

   - `:fixed` mode, `:Eω`: at the 1e-12 floor for every case but `multimode_field_plasma`
     (6.7e-12, in the weakest mode). That is the mode that makes a difference attributable
     to the operation that produced it, and nothing is excluded from it: the step sequence
     is imposed, so `stats/z` and `stats/dz` must match exactly.
   - `:fixed` mode, `:stats`, above the floor:
     - `modeavg_field_plasma` and `modeavg_field_adk` are both at the floor, from
       `stats/peak_ionisation_rate` (3.5e-15 and 5.7e-15). The rates are exponential in
       the field amplitude, so they amplify a one-ulp change by an order of magnitude.
     - `modeavg_field_vector` 8.0e-11 and `multimode_field_plasma` 5.8e-11, from
       `stats/transverse_integral_error_rel`/`_abs`. These are `HCubature`'s own error
       estimates for the modal overlap integral: the adaptive quadrature subdivides
       differently when the integrand changes in its last bits, so the error *estimate*
       moves by far more than the integral does. Everything else in those two cases is at or
       near the floor.
   - `:adaptive` mode: `stats/z` and `stats/dz` are excluded from the comparison
     (`RegressionCompare.STEP_STATS`), because the step-size controller responds to a one-ulp
     change far more strongly than the field does. The *number* of accepted steps is still
     checked, in both modes, as a named hard failure.
   - `:adaptive` `:stats` is set by the per-step statistics and is loose by design (up to
     7.3e-02 for `multimode_field_plasma`). `:adaptive` `:Eω` is usually two to four orders
     of magnitude tighter — which is the point of having two classes: the field is the
     quantity every later branch must not move.
   - The four ionising cases (`modeavg_field_plasma`, `modeavg_field_adk`,
     `modeavg_field_vector`, `multimode_field_plasma`, all argon at 0.1 bar since
     `gpu/int-D`) are the loose ones in the `:adaptive` mode: ionisation feeds back into
     the step-size controller, so a one-ulp change at the input moves the field by 9.2e-09
     to 7.3e-07 and the statistics by up to 7.3e-04. Their `:fixed` mode, where the step
     sequence is imposed, is at or near the floor and is what measures them.
   - The energy of `modeavg_field_plasma` is set by the same feedback. At the 300 µJ the
     `gpu/13-plasma` review recommended, the adaptive step count is not reproducible: one
     ulp at the input takes it from 92 accepted steps to 98, which is a `step count`
     failure rather than a measurable sensitivity, and a borrowed tolerance would only
     hide it (`Inf <= Inf` passes, so the step-count check, which is meant to be a hard
     failure, would become a no-op). The case runs at 175 µJ instead: the most strongly
     ionising energy whose adaptive step count is reproducible, at ±1, ±2 and ±8 ulp and
     at 1e-14. The sweep behind that choice is in `cases.jl` next to the case.
   - The two z-dependent cases, `gradient_field_kerr` and `taper_field_kerr`, used to be
     the loosest `:Eω` tolerances outside the ionising cases: 2.6e-07 and 8.7e-07 in the
     `:adaptive` mode, because the old one-point propagator rebuilt the operator at every
     stage, so where the stepper landed fed back into the field itself. `gpu/27` replaced
     that with `exp(Φ(t2) − Φ(t1))` over a table built once at setup, and the feedback is
     gone: re-measured on `gpu/int-E2` the `:adaptive` `:Eω` sensitivities are 2.5e-13
     and 7.6e-15, giving 2.5e-11 and the 1e-12 floor. The `:stats` numbers are almost
     unchanged (1.6e-04, 3.5e-05, from `stats/zdw` and `stats/peakintensity`), because
     they are still recorded once per accepted step. The `:fixed` mode was at the floor
     before and after.
   - `radial_field_raman` (added in `gpu/20-radial-device` as the matrix's only
     multi-column Raman case) is at the floor in the `:fixed` mode (2.4e-15) and 1.1e-09
     in the `:adaptive` one, in line with the two radial Kerr cases. Its baseline does not
     exist in any pre-`gpu/20` baseline directory, so the older baselines are run with
     `LUNA_REGRESSION_SKIP=radial_field_raman`.
   - `rect_modal_field` (added in `gpu/26-rectmode-fix` as the matrix's only Cartesian
     transverse domain) is the tightest multimode case: 1.8e-15 `Eω` in the `:fixed` mode
     and 2.1e-14 in the `:adaptive` one, both at the 1e-12 floor after the 100x. Unlike
     `multimode_field_plasma` there is no ionisation to feed a one-ulp change back into
     the step-size controller, and the adaptive step count (73) is reproducible. Its
     `:fixed` `:stats` 2.6e-12 is `stats/transverse_integral_error_rel`, the cubature's
     own error estimate, for the same reason as the other two adaptive-quadrature cases.
     Its baseline does not exist in any pre-`gpu/26` baseline directory, so the older
     baselines are run with `LUNA_REGRESSION_SKIP=rect_modal_field`.
=#
module RegressionTolerances

export TOLERANCES, tolerance

"""
    TOLERANCES

`Dict` from case name to a `Dict` from run mode (`:fixed`, `:adaptive`) to a `Dict` from
tolerance class (`:Eω`, `:stats`) to the tolerance on the normalised maximum difference.

Measured on an M1 Pro, Julia 1.13.0, `-t 1`, `set_fftw_mode(:estimate)`,
`set_fftw_threads(1)`, `BLAS.set_num_threads(1)`, FFTW wisdom disabled.
"""
const TOLERANCES = Dict{String, Dict{Symbol, Dict{Symbol, Float64}}}(
    "modeavg_field_kerr" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 1.0e-12, :stats => 5.9e-07)),
    "modeavg_field_nothg" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 2.3e-12, :stats => 2.1e-04)),
    "modeavg_field_plasma" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 9.2e-07, :stats => 1.2e-03)),
    "modeavg_field_raman" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 1.0e-12, :stats => 2.4e-09)),
    "modeavg_field_mixture" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 1.0e-12, :stats => 1.5e-10)),
    "modeavg_field_adk" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 4.8e-10, :stats => 1.8e-04)),
    "modeavg_field_vector" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 8.0e-11),
        :adaptive  => Dict(:Eω => 6.3e-09, :stats => 1.1e-03)),
    "modeavg_env_kerr" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 5.6e-12, :stats => 5.3e-04)),
    "modeavg_env_raman" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 2.9e-12, :stats => 4.2e-04)),
    "modeavg_env_thg" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 1.0e-12, :stats => 2.3e-07)),
    "gnlse_sech" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 3.6e-11, :stats => 1.6e-09)),
    "gnlse_raman_shock" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 1.9e-10, :stats => 1.2e-08)),
    "multimode_field_plasma" => Dict(
        :fixed     => Dict(:Eω => 6.7e-12, :stats => 5.8e-11),
        :adaptive  => Dict(:Eω => 7.3e-05, :stats => 7.3e-02)),
    "rect_modal_field" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 2.6e-12),
        :adaptive  => Dict(:Eω => 2.1e-12, :stats => 7.3e-09)),
    "radial_field_kerr" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 1.2e-08, :stats => 1.0e-12)),
    "radial_env_kerr" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 8.2e-09, :stats => 1.0e-12)),
    "radial_field_raman" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 1.1e-09, :stats => 1.0e-12)),
    "free3d_env_kerr" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 1.0e-12, :stats => 1.0e-12)),
    "free2d_field_chi2" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 7.4e-12, :stats => 1.0e-12)),
    "free2d_env_chi2" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 5.0e-12, :stats => 1.0e-12)),
    "gradient_field_kerr" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 2.5e-11, :stats => 1.6e-04)),
    "taper_field_kerr" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 1.0e-12, :stats => 3.5e-05)),
    "modeavg_field_legacy" => Dict(
        :fixed     => Dict(:Eω => 1.0e-12, :stats => 1.0e-12),
        :adaptive  => Dict(:Eω => 1.3e-11, :stats => 7.8e-08)),
)

"The tolerance used for a case, mode or class that is not in [`TOLERANCES`](@ref)."
const DEFAULT_TOLERANCE = 1.0e-12

"""
    tolerance(name, mode, class)

The regression tolerance for case `name` in run mode `mode` for tolerance class `class`
(`:Eω` or `:stats`), or [`DEFAULT_TOLERANCE`](@ref) if any of the three is not listed in
[`TOLERANCES`](@ref).
"""
function tolerance(name::AbstractString, mode::Symbol, class::Symbol)
    haskey(TOLERANCES, name) || return DEFAULT_TOLERANCE
    bymode = TOLERANCES[name]
    haskey(bymode, mode) || return DEFAULT_TOLERANCE
    get(bymode[mode], class, DEFAULT_TOLERANCE)
end

end
