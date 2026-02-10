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
    (; rays) = state.wkb
    (; lv, lb) = state.turbulence.turbulenceconstants
    (; coriolis_frequency) = state.namelists.atmosphere
    (; tref) = state.constants

    fc = coriolis_frequency / tref

    kr = rays.k[r, i, j, k]
    lr = rays.l[r, i, j, k]
    mr = rays.m[r, i, j, k]

    kh2 = kr^2 + lr^2

    n2r = interpolate_stratification(zr, state, N2())

    q00loc = interpolate_q(xr, yr, zr, state, Q00())
    q20loc = interpolate_q(xr, yr, zr, state, Q20())

    delta = n2r * kh2 / (2 * (n2r * kh2 + fc^2 * mr^2))

    gammas = mr^2 * q00loc * (lv * (1 - delta) + lb * delta)

    gammaw =
        mr^2 / 4 * n2r * kh2 / (n2r * kh2 + fc^2 * mr^2) *
        (lv * (1 - fc^2 / n2r) / (1 + kh2 / mr^2) - lb) *
        q20loc

    return gammas + gammaw
end