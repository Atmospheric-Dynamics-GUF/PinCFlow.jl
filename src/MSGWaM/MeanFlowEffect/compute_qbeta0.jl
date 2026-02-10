function compute_qbeta0 end

function compute_qbeta0(
    state::State,
    kr::AbstractFloat,
    lr::AbstractFloat,
    mr::AbstractFloat,
    omir::AbstractFloat,
    n2r::AbstractFloat,
    fc::AbstractFloat,
    wadr::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
    beta::AbstractFloat,
)::AbstractFloat
    int_min = 0.0
    int_max = 2 * pi
    nalpha = 20
    dalpha = (int_max - int_min) / nalpha
    alpha = int_min

    integral = 0.0
    while alpha <= int_max
        qalphabeta0 = compute_qalphabeta0(
            state,
            kr,
            lr,
            mr,
            omir,
            n2r,
            fc,
            wadr,
            i,
            j,
            k,
            beta,
            alpha,
        )
        integral += real(qalphabeta0 * exp(-1im * beta * alpha)) * dalpha
        alpha += dalpha
    end
    integral /= (2 * pi)

    return integral
end

function compute_qalphabeta0(
    state::State,
    kr::AbstractFloat,
    lr::AbstractFloat,
    mr::AbstractFloat,
    omir::AbstractFloat,
    n2r::AbstractFloat,
    fc::AbstractFloat,
    wadr::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
    beta::AbstractFloat,
    alpha::AbstractFloat,
)::Complex
    epsilon = 0.1
    int_min = 0.0
    int_max = 2 * pi * epsilon
    nphi = 20
    dphi = (int_max - int_min) / nphi
    phi = int_min

    integral = 0.0
    while phi <= int_max
        qtilde2 = compute_qtilde2(
            state,
            kr,
            lr,
            mr,
            omir,
            n2r,
            fc,
            wadr,
            i,
            j,
            k,
            phi,
            alpha,
            epsilon,
        )
        integral += sqrt(qtilde2) * exp(-1im * beta * phi / epsilon) * dphi
        phi += dphi
    end
    if beta == 0.0
        integral /= (2 * pi * epsilon)
    else
        integral /= (pi / beta * epsilon)
    end
    return integral
end

function compute_qtilde2(
    state::State,
    kr::AbstractFloat,
    lr::AbstractFloat,
    mr::AbstractFloat,
    omir::AbstractFloat,
    n2r::AbstractFloat,
    fc::AbstractFloat,
    wadr::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
    phi::AbstractFloat,
    alpha::AbstractFloat,
    epsilon::AbstractFloat,
)::AbstractFloat
    (; ld, lv, lb) = state.turbulence.turbulenceconstants
    (; rhobar) = state.atmosphere

    khr = sqrt(kr^2 + lr^2)

    uhat2 =
        2 * mr^2 * wadr * (omir^2 + fc^2) / rhobar[i, j, k] / omir /
        (khr^2 + mr^2)
    u10u10 =
        -(n2r - fc^2) * khr^2 * mr^2 / (khr^2 + mr^2)^2 * 2 * wadr / omir /
        rhobar[i, j, k] * exp(2im * alpha)

    b10 =
        sqrt(
            2 * wadr * mr^2 * (omir^2 + fc^2) / rhobar[i, j, k] / omir /
            (khr^2 + mr^2),
        ) * exp(1im * alpha)
    sterm = mr^2 / 2 * (uhat2 - real(u10u10 * exp(2im * phi / epsilon)))
    bterm = n2r + real(1im * mr * b10 * exp(1im * phi / epsilon))

    qsqrd = ld * (lv * sterm - lb * bterm)

    return max(0.0, qsqrd)
end