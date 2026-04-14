function compute_gw_turbulence_integrals! end

function compute_gw_turbulence_integrals!(
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
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme compute_gw_turbulence_integrals!(
        state,
        Val(turbulence_scheme),
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

function compute_gw_turbulence_integrals!(
    state::State,
    turbulence_scheme::Val{:NoTurbulence},
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

function compute_gw_turbulence_integrals!(
    state::State,
    turbulence_scheme::Val{:TKEScheme},
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
    (; rays) = state.wkb
    (; coriolis_frequency) = state.namelists.atmosphere
    (; shear) = state.turbulence.turbulencewkbintegrals
    (; rhobar) = state.atmosphere
    (; tref) = state.constants
    (; x_size, y_size) = state.namelists.domain
    (; lv) = state.turbulence.turbulenceconstants

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
    q00 = compute_q(state, r, i, j, k, 0.0)
    bhat0 = sqrt(
        n2r^2 * (kr^2 + lr^2) / (kr^2 + lr^2 + mr^2) * 2 * wadr / omir / rhob,
    )
    uhat0 =
        1im / mr / n2r * (omir^2 - n2r) / (omir^2 - fc^2) *
        (kr * omir + 1im * fc * lr) *
        bhat0

    shear[iray, jray, kray] += lv * abs(q00) * mr^2 * abs(uhat0)^2 / 2 * factor
 
    return
end