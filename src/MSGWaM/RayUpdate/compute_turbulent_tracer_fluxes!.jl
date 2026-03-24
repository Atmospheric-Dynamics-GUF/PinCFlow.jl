function compute_turbulent_tracer_fluxes! end

function compute_turbulent_tracer_fluxes!(
    state::State,
    factor::AbstractFloat,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    xr::AbstractFloat,
    yr::AbstractFloat,
    zr::AbstractFloat,
    iray::Integer,
    jray::Integer,
    kray::Integer,
)
    (; tracer_setup) = state.namelists.tracer

    compute_turbulent_tracer_fluxes!(
        state,
        tracer_setup,
        factor,
        r,
        i,
        j,
        k,
        xr,
        yr,
        zr,
        iray,
        jray,
        kray,
    )

    return
end

function compute_turbulent_tracer_fluxes!(
    state::State,
    tracer_setup::NoTracer,
    factor::AbstractFloat,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    xr::AbstractFloat,
    yr::AbstractFloat,
    zr::AbstractFloat,
    iray::Integer,
    jray::Integer,
    kray::Integer,
)
    return
end

function compute_turbulent_tracer_fluxes!(
    state::State,
    tracer_setup::TracerOn,
    factor::AbstractFloat,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    xr::AbstractFloat,
    yr::AbstractFloat,
    zr::AbstractFloat,
    iray::Integer,
    jray::Integer,
    kray::Integer,
)
    (; rays, integrals) = state.wkb
    (; coriolis_frequency) = state.namelists.atmosphere
    (; tracer_setup) = state.namelists.tracer
    (; qchi) = state.tracer.tracerwkbintegrals
    (; rhobar, thetabar) = state.atmosphere
    (; tref, kappa) = state.constants
    (; x_size, y_size) = state.namelists.domain
    (; lb) = state.turbulence.turbulenceconstants

    rhob = rhobar[iray, jray, kray]
    kr = rays.k[r, i, j, k]
    lr = rays.l[r, i, j, k]
    mr = rays.m[r, i, j, k]
    dkr = rays.dkray[r, i, j, k]
    dlr = rays.dlray[r, i, j, k]
    dmr = rays.dmray[r, i, j, k]
    n2r = interpolate_stratification(zr, state, N2())
    fc = coriolis_frequency * tref

    khr = sqrt(kr^2 + lr^2)

    omir = -sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

    wadr = rays.dens[r, i, j, k] * dmr

    if x_size > 1
        wadr *= dkr
    end
    if y_size > 1
        wadr *= dlr
    end

    n2r = interpolate_stratification(zr, state, N2())
    q10 = compute_q(state, r, i, j, k, 1.0)
    bhat0 = sqrt(
        n2r^2 * (kr^2 + lr^2) / (kr^2 + lr^2 + mr^2) * 2 * wadr / omir / rhob,
    )
    uhat0 =
        1im / mr / n2r * (omir^2 - n2r) / (omir^2 - fc^2) *
        (kr * omir + 1im * fc * lr) *
        bhat0
    vhat0 =
        1im / mr / n2r * (omir^2 - n2r) / (omir^2 - fc^2) *
        (lr * omir - 1im * fc * kr) *
        bhat0
    what0 = 1im * omir / n2r * bhat0
    pihat0 =
        1im / mr * (omir^2 - n2r^2) / n2r / thetabar[iray, jray, kray] / kappa

    dchidx = interpolate_mean_flow(xr, yr, zr, state, DChiDX())
    dchidy = interpolate_mean_flow(xr, yr, zr, state, DChiDY())
    dchidz = interpolate_mean_flow(xr, yr, zr, state, DChiDZ())

    chihat = -1im / omir * (uhat0 * dchidx + vhat0 * dchidy + what0 * dchidz)

    qchi[iray, jray, kray] += lb * imag(mr * conj(q10) * chihat) * factor

    integrals.uhat[iray, jray, kray] += uhat0 * factor
    integrals.vhat[iray, jray, kray] += vhat0 * factor
    integrals.what[iray, jray, kray] += what0 * factor
    integrals.bhat[iray, jray, kray] += bhat0 * factor
    integrals.pihat[iray, jray, kray] += pihat0 * factor

    return
end