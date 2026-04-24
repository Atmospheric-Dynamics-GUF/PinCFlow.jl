function compute_turbulent_damping! end

function compute_turbulent_damping!(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    zr::AbstractFloat,
    dt::AbstractFloat,
)
    (; rays) = state.wkb
    (; lv, lb) = state.turbulence.turbulenceconstants
    (; coriolis_frequency) = state.namelists.atmosphere
    (; tref) = state.constants
    (; x_size, y_size) = state.namelists.domain
    (; rhobar) = state.atmosphere
    (; turbulence_damping) = state.namelists.wkb

    if !turbulence_damping
        return
    end

    fc = coriolis_frequency / tref

    kr = rays.k[r, i, j, k]
    lr = rays.l[r, i, j, k]
    mr = rays.m[r, i, j, k]

    dkr = rays.dkray[r, i, j, k]
    dlr = rays.dlray[r, i, j, k]
    dmr = rays.dmray[r, i, j, k]

    kh2 = kr^2 + lr^2

    n2r = interpolate_stratification(zr, state, N2())

    omir = -sqrt(n2r * kh2 + fc^2 * mr^2) / sqrt(kh2 + mr^2)
    rhob = rhobar[i, j, k]

    factor = dmr

    if x_size > 1
        factor *= dkr
    end
    if y_size > 1
        factor *= dlr
    end

    wadr = rays.dens[r, i, j, k] * factor

    q00, q10, q20 = compute_turbulent_velocity(state, r, i, j, k)

    delta = n2r * kh2 / (2 * (n2r * kh2 + fc^2 * mr^2))

    gammas = mr^2 * real(q00) * (lv * (1 - delta) + lb * delta)

    gammaw =
        mr^2 / 4 * n2r * kh2 / (n2r * kh2 + fc^2 * mr^2) *
        (lv * (1 - fc^2 / n2r) / (1 + kh2 / mr^2) - lb) *
        real(q20)

    gammawp =
        -lb * mr / omir / 2 *
        sqrt(n2r^2 * kh2 / (kh2 + mr^2) * rhob / 2 / abs(omir / wadr)) *
        real(1im * q10)

    wadr *= 1 - 2 * dt * (gammas + gammaw + gammawp)

    rays.dens[r, i, j, k] = wadr / factor

    return
end
