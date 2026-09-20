module RK45
import Dates
import Logging
import Printf: @sprintf
import Luna.Utils: format_elapsed, isdevice

#Get Butcher tableau etc from separate file (for convenience of changing if wanted)
include("dopri.jl")

function solve(f!, y0, t, dt, tmax;
               rtol=1e-6, atol=1e-10, safety=0.9, max_dt=Inf, min_dt=0, locextrap=true,
               norm=weaknorm,
               kwargs...)
    stepper = Stepper(f!, y0, t, dt,
                      rtol=rtol, atol=atol, safety=safety, max_dt=max_dt, min_dt=min_dt, locextrap=locextrap, norm=norm)
    return solve(stepper, tmax; kwargs...)
end

function solve_precon(f!, linop, y0, t, dt, tmax;
                    rtol=1e-6, atol=1e-10, safety=0.9, max_dt=Inf, min_dt=0, locextrap=true, norm=weaknorm,
                    kwargs...)
    stepper = PreconStepper(f!, linop, y0, t, dt,
                      rtol=rtol, atol=atol, safety=safety, max_dt=max_dt, min_dt=min_dt, locextrap=locextrap, norm=norm)
    return solve(stepper, tmax; kwargs...)
end

function solve(s, tmax; stepfun=donothing!, output=false, outputN=201,
                        status_period=1, repeat_limit=10)
    if output
        #= The saved array is always a host array, whatever the solution lives on, and
           `_saveto!` copies device to host. =#
        yout = Array{eltype(s.y)}(undef, (size(s.y)..., outputN))
        tout = range(s.t, stop=tmax, length=outputN)
        saved = 1
        _saveto!(yout, 1, s.y)
    end

    steps = 0
    repeated = 0
    repeated_tot = 0

    Logging.@info "Starting propagation"
    start = Dates.now()
    tic = Dates.now()
    while s.tn <= tmax
        ok = step!(s)
        steps += 1
        if Dates.value(Dates.now()-tic) > 1000*status_period
            speed = s.tn/(Dates.value(Dates.now()-start)/1000)
            eta_in_s = (tmax-s.tn)/(speed)
            if eta_in_s > 356400
                Logging.@info @sprintf("Progress: %.2f %%, ETA: XX:XX:XX, stepsize %.2e, err %.2f, repeated %d/%d steps",
                s.tn/tmax*100, s.dt, s.err, repeated_tot, steps)
            else
                eta_in_ms = Dates.Millisecond(ceil(eta_in_s*1000))
                etad = Dates.DateTime(Dates.UTInstant(eta_in_ms))
                Logging.@info @sprintf("Progress: %.2f %%, ETA: %s, stepsize %.2e, err %.2f, repeated %d/%d steps",
                    s.tn/tmax*100, Dates.format(etad, "HH:MM:SS"), s.dt, s.err, repeated_tot, steps)
            end
            flush(stderr)
            tic = Dates.now()
        end
        if ok
            if output
                while (saved<outputN) && tout[saved+1] < s.tn
                    ti = tout[saved+1]
                    _saveto!(yout, saved+1, interpolate(s, ti))
                    saved += 1
                end
            end
            stepfun(s.yn, s.tn, s.dtn, t -> interpolate(s, t))
            repeated = 0
        else
            repeated += 1
            repeated_tot += 1
            if repeated > repeat_limit
                error("Reached limit for step repetition ($repeat_limit)")
            end
        end
    end
    totaltime = Dates.now()-start
    dtstring = format_elapsed(totaltime)
    Logging.@info @sprintf("Propagation finished in %s, %d steps",
                           dtstring, steps)

    if output
        return collect(tout), yout, steps
    else
        return nothing
    end
end

#= Write one save into the host output array. A device solution is copied down first;
   `Array(y)` is a no-op test plus a return on the host, and `solve(output=true)` is a
   low-level entry point which is not on the per-step path. =#
function _saveto!(yout, idx, y)
    selectdim(yout, ndims(yout), idx) .= y isa Array ? y : Array(y)
    nothing
end


mutable struct Stepper{T<:AbstractArray, F, nT}
    f!::F  # RHS function
    y::T  # Solution at current t
    yn::T  # Solution at t+dt
    yi::T  # Interpolant array (see interpolate())
    yerr::T  # solution error estimate (from embedded RK)
    ks::NTuple{7, T}  # k values (intermediate solutions for Runge-Kutta method)
    t::Float64  # current time (propagation variable)
    tn::Float64  # next time
    dt::Float64  # time step
    dtn::Float64  # time step for next step
    rtol::Float64  # relative tolerance on error
    atol::Float64  # absolute tolerance on error
    safety::Float64  # safety factor for stepsize control
    max_dt::Float64  # maximum value for dt (default Inf)
    min_dt::Float64  # minimum value for dt (default 0)
    locextrap::Bool  # true if using local extrapolation
    ok::Bool  # true if current step was successful. Also means "an FSAL move is owed"
              # (see step!), so nothing may write it after step! has returned.
    err::Float64  # error metric to be compared to tol
    errlast::Float64  # error of the most recent successful step
    norm::nT # function to calculate error metric, defaults to RK45.weaknorm
end

function Stepper(f!, y0, t, dt;
                 rtol=1e-6, atol=1e-10, safety=0.9, max_dt=Inf, min_dt=0,
                 locextrap=true, norm=weaknorm)
    k1 = similar(y0)
    f!(k1, y0, t)
    ks = (k1, similar(k1), similar(k1), similar(k1), similar(k1), similar(k1), similar(k1))
    yerr = similar(y0)
    return Stepper(f!, copy(y0), copy(y0), similar(y0), yerr, ks,
        float(t), float(t), float(dt), float(dt),
        float(rtol), float(atol), float(safety), float(max_dt), float(min_dt),
        locextrap, false, 0.0, 0.0, norm)
end

mutable struct PreconStepper{T<:AbstractArray, F, P, nT}
    fbar!::F  # RHS callable
    prop!::P # linear propagator callable
    y::T  # Solution at current t
    yn::T  # Solution at t+dt
    yi::T  # Interpolant array (see interpolate())
    yerr::T  # solution error estimate (from embedded RK)
    ks::NTuple{7, T}  # k values (intermediate solutions for Runge-Kutta method)
    t::Float64  # current time (propagation variable)
    tn::Float64  # next time
    dt::Float64  # time step
    dtn::Float64  # time step for next step
    rtol::Float64  # relative tolerance on error
    atol::Float64  # absolute tolerance on error
    safety::Float64  # safety factor for stepsize control
    max_dt::Float64  # maximum value for dt (default Inf)
    min_dt::Float64  # minimum value for dt (default 0)
    locextrap::Bool  # true if using local extrapolation
    ok::Bool  # true if current step was successful. Also means "an FSAL move is owed"
              # (see step!), so nothing may write it after step! has returned.
    err::Float64  # error metric to be compared to tol
    errlast::Float64  # error of the most recent successful step
    norm::nT  # function to calculate error metric, defaults to RK45.weaknorm
end

function PreconStepper(f!, linop, y0, t, dt;
                       rtol=1e-6, atol=1e-10, safety=0.9, max_dt=Inf, min_dt=0,
                       locextrap=true, norm=weaknorm)
    prop! = make_prop!(linop, y0)
    fbar! = make_fbar!(f!, prop!, y0)
    k1 = similar(y0)
    #= fbar! propagates its second argument in place, so the first RHS evaluation is made
       on a scratch copy rather than on the caller's initial condition. =#
    y0c = copy(y0)
    fbar!(k1, y0c, t, t)
    ks = (k1, similar(k1), similar(k1), similar(k1), similar(k1), similar(k1), similar(k1))
    yerr = similar(y0)

    return PreconStepper(fbar!, prop!, copy(y0), copy(y0), similar(y0), yerr, ks,
        float(t), float(t), float(dt), float(dt), float(rtol), float(atol), float(safety),
        float(max_dt), float(min_dt), locextrap, false, 0.0, 0.0, norm)
end

function step!(s)
    # FSAL: k7 of the previous step is k1 of this one. The move is deferred to here
    # (rather than done at the end of step!) because the completed step's dense output is
    # consumed between step! returning and this call, and the interpolant's k1 weight is
    # nonzero (interpC column 1; b1(1) = 35/384) -- moving it any earlier degrades every
    # interpolated save from 4th to 2nd order. s.ok is false before the first step (ks[2:7]
    # are `similar`, i.e. undef) and after a rejected step (k1 must survive for the retry).
    # For the PreconStepper this must also stay ahead of evaluate!'s rebase of ks[1] into
    # the new anchor frame, which it does by construction.
    s.ok && (s.ks[1] .= s.ks[end])
    evaluate!(s)

    # Propagate the 5th-order solution (local extrapolation, the default) or the embedded
    # 4th-order one. Both are formed explicitly rather than reusing the stage-6
    # accumulation evaluate! happens to leave in yn -- which is the 5th-order solution --
    # so this does not depend on the RHS leaving its input array alone. For locextrap the
    # two are bit-identical anyway: b5[1:6] == B[6], b5[7] == 0, same accumulation order.
    bprop = s.locextrap ? b5 : b4
    combine!(s.yn, s.y, s.ks, s.dt, bprop)

    errorestimate!(s.yerr, s.ks, s.dt)
    s.err = s.norm(s.yerr, s.y, s.yn, s.rtol, s.atol)
    s.ok = s.err <= 1
    stepcontrolPI!(s)
    if s.ok
        s.tn = s.t + s.dt
    else
        s.yn .= s.y
    end
    prop!_maybe(s) # propagate to new time to pass correct solution to stepfun
    return s.ok
end

#= Every per-step combination of the stage derivatives is one fused broadcast over the
   whole array rather than a sequence of `.+=` passes. The n-ary `+` in a broadcast is
   left-associated and the terms are written in the tableau's order, with each
   coefficient formed as `dt*b` on the host exactly as the sequential version did, so on
   the CPU in double precision the result is bit-for-bit what the loop produced. The
   scalars are converted to the state's real element type first: nothing reachable from a
   kernel argument may hold a Float64 (Metal's compiler rejects it), and for Float64 the
   conversion is the identity. =#

#= Which weights of the tableau are exactly zero. The sequential accumulation skipped
   those terms, and the fused combines below skip the same ones, so the two agree
   bitwise. Asserted rather than assumed, so that a different tableau cannot silently
   change which terms are dropped. =#
@assert all(all(!iszero, B[ii]) for ii = 1:5)
@assert iszero(B[6][2]) && !iszero(B[6][1]) && all(!iszero, B[6][3:6])
@assert iszero(b5[2]) && iszero(b5[7]) && all(!iszero, b5[[1, 3, 4, 5, 6]])
@assert iszero(b4[2]) && all(!iszero, b4[[1, 3, 4, 5, 6, 7]])
@assert iszero(errest[2]) && all(!iszero, errest[[1, 3, 4, 5, 6, 7]])

"""
    combine!(yn, y, ks, dt, b, n=7)

`yn = y + Σⱼ dt bⱼ kⱼ` over the first `n` stages, in one fused pass. `n < 7` is a
Butcher stage (`b = B[n]`); `n == 7` is the propagation, with `b` either `b5` or `b4`.
"""
function combine!(yn, y, ks, dt, b, n=7)
    R = real(eltype(yn))
    k1, k2, k3, k4, k5, k6, k7 = ks
    if n == 1
        c1 = convert(R, dt*b[1])
        @. yn = y + c1*k1
    elseif n == 2
        c1 = convert(R, dt*b[1]); c2 = convert(R, dt*b[2])
        @. yn = y + c1*k1 + c2*k2
    elseif n == 3
        c1 = convert(R, dt*b[1]); c2 = convert(R, dt*b[2]); c3 = convert(R, dt*b[3])
        @. yn = y + c1*k1 + c2*k2 + c3*k3
    elseif n == 4
        c1 = convert(R, dt*b[1]); c2 = convert(R, dt*b[2]); c3 = convert(R, dt*b[3])
        c4 = convert(R, dt*b[4])
        @. yn = y + c1*k1 + c2*k2 + c3*k3 + c4*k4
    elseif n == 5
        c1 = convert(R, dt*b[1]); c2 = convert(R, dt*b[2]); c3 = convert(R, dt*b[3])
        c4 = convert(R, dt*b[4]); c5 = convert(R, dt*b[5])
        @. yn = y + c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5
    elseif n == 6
        c1 = convert(R, dt*b[1]); c3 = convert(R, dt*b[3]); c4 = convert(R, dt*b[4])
        c5 = convert(R, dt*b[5]); c6 = convert(R, dt*b[6])
        @. yn = y + c1*k1 + c3*k3 + c4*k4 + c5*k5 + c6*k6
    else
        c1 = convert(R, dt*b[1]); c3 = convert(R, dt*b[3]); c4 = convert(R, dt*b[4])
        c5 = convert(R, dt*b[5]); c6 = convert(R, dt*b[6])
        if iszero(b[7]) # b5 is FSAL and does not use k7; b4 does
            @. yn = y + c1*k1 + c3*k3 + c4*k4 + c5*k5 + c6*k6
        else
            c7 = convert(R, dt*b[7])
            @. yn = y + c1*k1 + c3*k3 + c4*k4 + c5*k5 + c6*k6 + c7*k7
        end
    end
    yn
end

"""
Embedded error estimate `Σᵢ dt kᵢ eᵢ` in one pass. `errest[2] == 0`, which the sequential
accumulation also skipped; `dt*kᵢ*eᵢ` keeps that association so the result is unchanged.

!!! note
    The skip of the second term is written into the expression rather than branched on, as
    `combine!` branches on `iszero(b[7])`: with one tableau compiled in there is nothing to
    branch on, and the `@assert`s above hold it to DOPRI5's zero pattern. If the tableau
    ever becomes a runtime choice, this function needs the same treatment `combine!` has --
    a branch per zero pattern, or the term restored unconditionally (which costs one
    multiply-add per element and is bit-identical only when `kᵢ` is finite).
"""
function errorestimate!(yerr, ks, dt)
    R = real(eltype(yerr))
    d = convert(R, dt)
    k1, k2, k3, k4, k5, k6, k7 = ks
    e1 = convert(R, errest[1]); e3 = convert(R, errest[3]); e4 = convert(R, errest[4])
    e5 = convert(R, errest[5]); e6 = convert(R, errest[6]); e7 = convert(R, errest[7])
    @. yerr = 0 + d*k1*e1 + d*k3*e3 + d*k4*e4 + d*k5*e5 + d*k6*e6 + d*k7*e7
    yerr
end

function evaluate!(s::Stepper)
    # Set new time and stepsize values -- this happens at the beginning because
    # the interpolant still requires the old values after the step has finished
    s.dt = s.dtn
    s.t = s.tn
    s.y .= s.yn
    for ii = 1:6
        combine!(s.yn, s.y, s.ks, s.dt, B[ii], ii)
        s.f!(s.ks[ii+1], s.yn, s.t+nodes[ii]*s.dt)
    end
end

function evaluate!(s::PreconStepper)
    # Set new time and stepsize values -- this happens at the beginning because
    # the interpolant still requires the old values after the step has finished
    s.y .= s.yn
    s.prop!(s.ks[1], s.t, s.tn)
    s.dt = s.dtn
    s.t = s.tn
    for ii = 1:6
        combine!(s.yn, s.y, s.ks, s.dt, B[ii], ii)
        s.fbar!(s.ks[ii+1], s.yn, s.t, s.t+nodes[ii]*s.dt)
    end
    #= The last stage's `fbar!` propagated (clobbered) `yn` in place. That does not
       matter: `step!` rebuilds `yn` from `y` and the weight vector before anything reads
       it. =#
end

prop!_maybe(s::PreconStepper) = s.prop!(s.yn, s.t, s.tn)
prop!_maybe(s) = nothing

"""
    interpolate(s, ti)

The dense-output solution at `ti`, as one fused broadcast into the stepper's own
interpolant buffer `s.yi`.

!!! note
    The returned array is the stepper's buffer, not a fresh one: it is overwritten by the
    next call to `interpolate` and by the next step. Copy it if you need to keep it. (At
    `ti == s.t` and `ti == s.tn` the stepper's own `y`/`yn` are returned, which has always
    been the case.)
"""
function interpolate(s::Stepper, ti::Float64)
    if ti > s.tn
        error("Attempting to extrapolate!")
    end
    if ti == s.t
        return s.y
    elseif ti == s.tn
        return s.yn
    end
    interpolant!(s, ti)
end

@doc (@doc interpolate)
function interpolate(s::PreconStepper, ti::Float64)
    if ti > s.tn
        error("Attempting to extrapolate!")
    end
    if ti == s.t
        return s.y
    elseif ti == s.tn
        return s.yn
    end
    interpolant!(s, ti)
    s.prop!(s.yi, s.t, ti)
    return s.yi
end

#= `y + dt*(Σᵢ bᵢ kᵢ)` in one pass. The inner sum is left-associated over all seven
   stages in order, including the zero-weight k2 term, which is what the sequential
   fill!/.+= accumulation did, so on the CPU in double precision this is bit-identical. =#
function interpolant!(s, ti::Float64)
    σ = (ti - s.t)/s.dt
    σp = map(p -> σ^p, range(1, stop=4))
    b = sum(σp.*interpC, dims=1)
    R = real(eltype(s.y))
    d = convert(R, s.dt)
    w1 = convert(R, b[1]); w2 = convert(R, b[2]); w3 = convert(R, b[3])
    w4 = convert(R, b[4]); w5 = convert(R, b[5]); w6 = convert(R, b[6])
    w7 = convert(R, b[7])
    k1, k2, k3, k4, k5, k6, k7 = s.ks
    y = s.y
    @. s.yi = y + d*(0 + k1*w1 + k2*w2 + k3*w3 + k4*w4 + k5*w5 + k6*w6 + k7*w7)
    return s.yi
end

"Make propagator for the case of constant linear operator"
function make_prop!(linop::AbstractArray, y0)
    prop! = let linop=linop
        function prop!(y, t1, t2, bwd=false)
            dt = convert(real(eltype(y)), bwd ? (t1-t2) : (t2-t1))
            @. y *= exp(linop*dt)
        end
    end
end

"""
Make propagator for the case of non-constant linear operator.

`linop!(out, z)` is host code -- the operators in `LinearOps` are scalar loops over
`Modes.neff` -- so when the state lives on a device the operator is evaluated into a host
buffer of the state's element type and copied up, once per distinct `t2`.

This is the fallback, and the path a user-supplied `linop!` takes.
`Luna.run(...; tabulate_linop=true)` replaces the operator with a
[`LinearOps.TabulatedLinop`](@ref Luna.LinearOps.TabulatedLinop), which has its own method
of this function and does no host work per stage.
"""
function make_prop!(linop!, y0)
    linop_int = similar(y0)
    hostbuf = isdevice(y0) ? Array{eltype(y0)}(undef, size(y0)) : nothing
    lastt2 = Ref(typemin(Float64))
    prop! = let linop! = linop!, linop_int = linop_int, hostbuf = hostbuf, lastt2 = lastt2
        function prop!(y, t1, t2, bwd=false)
            #= linop is always evaluated at later time, even for backward propagation
                therefore, linop is often evaluated at the same t2 twice in a row=#
            if lastt2[] != t2
                if isnothing(hostbuf)
                    linop!(linop_int, t2)
                else
                    linop!(hostbuf, t2)
                    copyto!(linop_int, hostbuf)
                end
            end
            lastt2[] = t2
            dt = convert(real(eltype(y)), bwd ? (t1-t2) : (t2-t1))
            @. y *= exp(linop_int*dt)
        end
    end
    return prop!
end

"""
Make closure for the pre-conditioned RHS function.

!!! note
    `fbar!(out, ybar, t1, t2)` propagates `ybar` to `t2` **in place**, so the caller must
    not rely on its contents afterwards. `evaluate!(::PreconStepper)` rebuilds `yn` from
    `y` before every stage and `step!` rebuilds it afterwards, so this is safe there and
    saves one field-sized buffer.
"""
function make_fbar!(f!, prop!, y0)
    fbar! = let f! = f!, prop! = prop!
        function fbar!(out, ybar, t1, t2)
            prop!(ybar, t1, t2) # propagate to t2 (in place)
            f!(out, ybar, t2) # evaluate RHS function
            prop!(out, t1, t2, true) # propagate back to t1
        end
    end
end

#= The error norms are single reductions over the state arrays rather than scalar loops,
   so that they run on every backend. The accumulator is initialised in the state's real
   element type so that no Float64 reaches a device kernel; `rtol` and `atol` are
   converted for the same reason where they appear inside the reduction, while the final
   scalar arithmetic stays in Float64.

   A norm which indexes its arguments elementwise still works on the CPU, and a
   user-supplied one is called with the materialised error estimate exactly as before. =#

"""
One reduction over several arrays, as `op` folded over `f` applied elementwise.

The arrays are combined into a lazy `Broadcast.Broadcasted` rather than passed to
`mapreduce` directly, which matters on both backends. With an explicit `init` Base
reduces a `Broadcasted` with `mapfoldl`: a serial fold in index order which performs the
same operations in the same order as the scalar loops these replace, so the result is
bit-identical on the CPU, and which allocates only the boxed tuple result (16 bytes per
call) rather than anything field-sized. Base's *multi-array* `mapreduce`, by contrast,
materialises `map(f, As...)` first, which is a field-sized allocation per step.
`GPUArrays` has a `mapreduce` method for a `Broadcasted` of its own style, so on a device
this is its tree reduction -- whose summation order differs, as a parallel reduction's
must.
"""
@inline _zipreduce(f, op, init, arrs...) = mapreduce(
    identity, op, Broadcast.instantiate(Broadcast.broadcasted(f, arrs...)); init=init)

@inline _add3(a, b) = (a[1]+b[1], a[2]+b[2], a[3]+b[3])
@inline _max2(a, b) = (max(a[1], b[1]), max(a[2], b[2]))
@inline _maxmap(yerr, y, yn) = (abs(yerr), max(abs(y), abs(yn)))
@inline _abs2map(y, yn, yerr) = (abs2(y), abs2(yn), abs2(yerr))

"Max-ish norm (from Dane Austin's code, no idea where he got it from)."
function maxnorm(yerr, y, yn, rtol, atol)
    R = real(eltype(yerr))
    maxerr, maxy = _zipreduce(_maxmap, _max2, (zero(R), zero(R)), yerr, y, yn)
    return maxerr/(atol + rtol*maxy)
end

"Alternative form of max-ish norm."
function maxnorm_ratio(yerr, y, yn, rtol, atol)
    R = real(eltype(yerr))
    at = convert(R, atol)
    rt = convert(R, rtol)
    f = @inline function (e, a, b)
        abs(e)/(at + rt*max(abs(a), abs(b)))
    end
    return _zipreduce(f, max, zero(R), yerr, y, yn)
end

"Semi-norm as used in DifferentialEquations.jl, see Hairer, Solving Ordinary Differential
Equations: Nonstiff Problems, eq. (4.11) (p.168 of the second revised edition)."
function normnorm(yerr, y, yn, rtol, atol)
    R = real(eltype(yerr))
    at = convert(R, atol)
    rt = convert(R, rtol)
    f = @inline function (e, a, b)
        abs2(e/(at + rt*max(abs(a), abs(b))))
    end
    s = _zipreduce(f, +, zero(R), yerr, y, yn)
    sqrt(s/length(yerr))
end

"'Weak' norm as used in fnfep."
function weaknorm(yerr, y, yn, rtol, atol)
    R = real(eltype(yerr))
    sy, syn, syerr = _zipreduce(_abs2map, _add3, (zero(R), zero(R), zero(R)), y, yn, yerr)
    errwt = max(max(sqrt(sy), sqrt(syn)), atol)
    return sqrt(syerr)/rtol/errwt
end

"Simple proportional error controller, see e.g. Hairer eq. (4.13)."
function stepcontrolP!(s)
    if s.ok
        # if error is zero, there is no nonlinearity: increase step size by a lot
        s.dtn = s.err == 0 ? 1.5*s.dt : s.dt * min(5, s.safety*(s.err)^(-1/5))
    else
        if !isfinite(s.err) # check for NaN or Inf
            s.dtn = s.dt/2  # if we have one then we're in big trouble so halve the step size
        else
            s.dtn = s.dt * max(0.1, s.safety*(s.err)^(-1/5))
        end
    end
    steplims!(s)
end

"Proportional-integral error controller, aka Lund stabilisation.
See G. Söderlind and L. Wang, J. Comput. Appl. Math. 185, 225 (2006).
"
function stepcontrolPI!(s)
    β1 = 3/5 / 5
    β2 = -1/5 / 5
    ε = 0.8
    if s.ok
        s.errlast == 0 && (s.errlast = s.err) # if last error is zero, use current error instead
        if s.err == 0
            fac = 1.5 # zero error means no nonlinearity: increase step size by a lot 
        else
            fac = s.safety * (ε/s.err)^β1 * (ε/s.errlast)^β2
        end
        # (0.99 <= fac <= 1.01) && (fac = 1.0)
        s.dtn = fac * s.dt
        s.errlast = s.err
    else
        if !isfinite(s.err) # check for NaN or Inf
            s.dtn = s.dt/2  # if we have one then we're in big trouble so halve the step size
        else
            s.dtn = s.dt * max(0.1, s.safety*(s.err)^(-1/5))
        end
    end
    steplims!(s)
end

"Apply user-defined limits on step size."
function steplims!(s)
    if s.dtn > s.max_dt
        s.dtn = s.max_dt
    elseif s.dtn < s.min_dt
        s.dtn = s.min_dt
        s.ok = true
    end
end

function donothing!(y, z, dz, interpolant)
end

end