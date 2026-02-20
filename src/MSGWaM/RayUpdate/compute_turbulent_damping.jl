function compute_turbulent_damping(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    zr::AbstractFloat,
)::Tuple{<:AbstractFloat, <:AbstractFloat, <:AbstractFloat}
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

    omir = -sqrt(n2r * kh2 + fc^2 * mr^2) / sqrt(kh2 + mr^2)
    rhob = interpolate_rhobar(zr, state, Rhobar())

    q00, q10, q20 = compute_q(state, r, i, j, k)

    rays.q00[r, i, j, k] = q00
    rays.q10[r, i, j, k] = q10
    rays.q20[r, i, j, k] = q20

    delta = n2r * kh2 / (2 * (n2r * kh2 + fc^2 * mr^2))

    gammas = mr^2 * q00 * (lv * (1 - delta) + lb * delta)

    gammaw =
        mr^2 / 4 * n2r * kh2 / (n2r * kh2 + fc^2 * mr^2) *
        (lv * (1 - fc^2 / n2r) / (1 + kh2 / mr^2) - lb) *
        q20

    gammawp =
        lb * mr / omir *
        sqrt(
            n2r^2 * kh2 / (kh2 + mr^2) * rhob / 2 / omir /
            rays.dens[r, i, j, k],
        ) *
        q10

    return gammas, gammaw, gammawp
end