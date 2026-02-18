function compute_q end

function compute_q(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::Tuple{<:AbstractFloat, <:AbstractFloat, <:AbstractFloat}
    (; rays) = state.wkb

    if rays.dens[r, i, j, k] == 0.0
        return 0.0, 0.0, 0.0
    end

    q00 = compute_q(state, r, i, j, k, 0.0)
    q10 = compute_q(state, r, i, j, k, 1.0)
    q20 = compute_q(state, r, i, j, k, 2.0)

    return q00, q10, q20
end

function compute_q(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    beta::AbstractFloat,
)::AbstractFloat
    (; rays) = state.wkb
    (; coriolis_frequency) = state.namelists.atmosphere
    (; tref, g_ndim) = state.constants
    (; rhobar, n2) = state.atmosphere
    (; x_size, y_size) = state.namelists.domain


    (xr, yr, zr) = get_physical_position(rays, r, i, j, k)

    rhob = interpolate_rhobar(zr, state, Rhobar())
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

    int_min = 0.0
    int_max = 2 * pi
    nphi = 20
    dphi = (int_max - int_min) / nphi
    phi = 0.0
    integral = 0.0
    while phi <= int_max
        qtilde2 = compute_q(state, rhob, wadr, kr, lr, mr, n2r, fc, omir, phi)
        integral += sqrt(qtilde2) * exp(-1im * beta * phi) * dphi
        phi += dphi
    end
    if beta == 0.0
        integral /= (2 * pi)
    else
        integral /= pi
    end

    return real(integral)
end

function compute_q(
    state::State,
    rhob::AbstractFloat,
    wadr::AbstractFloat,
    kr::AbstractFloat,
    lr::AbstractFloat,
    mr::AbstractFloat,
    n2r::AbstractFloat,
    fc::AbstractFloat,
    omir::AbstractFloat,
    phi::AbstractFloat,
)::AbstractFloat
    (; ld, lv, lb) = state.turbulence.turbulenceconstants

    uhat2 =
        2 * mr^2 * wadr * (omir^2 + fc^2) / rhob / omir / (kr^2 + lr^2 + mr^2)
    u01u01 =
        -(n2r - fc^2) * (kr^2 + lr^2) * mr^2 / (kr^2 + lr^2 + mr^2)^2 *
        2 *
        wadr / omir / rhob
    bhat = sqrt(
        n2r^2 * (kr^2 + lr^2) / (kr^2 + lr^2 + mr^2) * 2 * wadr / omir / rhob,
    )

    sterm = mr^2 / 2 * (uhat2 - real(u01u01 * exp(2im * phi)))
    bterm = n2r + real(1im * mr * bhat * exp(1im * phi))

    qtilde2 = ld * (lv * sterm - lb * bterm)

    return max(0, qtilde2)
end