function compute_turbulent_damping(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    xr::AbstractFloat,
    yr::AbstractFloat,
    zr::AbstractFloat,
)::AbstractFloat
    (; turbulence_scheme) = state.namelists.turbulence

    return compute_turbulent_damping(
        state,
        r,
        i,
        j,
        k,
        xr,
        yr,
        zr,
        turbulence_scheme,
    )
end

function compute_turbulent_damping(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    xr::AbstractFloat,
    yr::AbstractFloat,
    zr::AbstractFloat,
    turbulence_scheme::NoTurbulence,
)::AbstractFloat
    return 0.0
end

function compute_turbulent_damping(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    xr::AbstractFloat,
    yr::AbstractFloat,
    zr::AbstractFloat,
    turbulence_scheme::TKEScheme,
)::AbstractFloat
    (; coriolis_frequency) = state.namelists.atmosphere
    (; tref) = state.constants
    (; rays) = state.wkb
    (; lturb_ndim) = state.turbulence.turbulenceconstants

    fc = coriolis_frequency * tref

    kr = rays.k[r, i, j, k]
    lr = rays.l[r, i, j, k]
    mr = rays.m[r, i, j, k]

    khr = sqrt(kr^2 + lr^2)

    n2r = interpolate_stratification(zr, state, N2())

    delta = n2r * khr^2 / (2 * (n2r * khr^2 + fc^2 * mr^2))

    tkeloc = interpolate_tke(xr, yr, zr, state)

    gammas = mr^2 * sqrt(tkeloc) * lturb_ndim

    return gammas
end