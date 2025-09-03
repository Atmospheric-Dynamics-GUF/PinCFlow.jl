function initialize_wkb_wave_packet!(
    state::State,
    omi_ini::AbstractArray{<:AbstractFloat, 4},
    wnk_ini::AbstractArray{<:AbstractFloat, 4},
    wnl_ini::AbstractArray{<:AbstractFloat, 4},
    wnm_ini::AbstractArray{<:AbstractFloat, 4},
    wad_ini::AbstractArray{<:AbstractFloat, 4},
)
    (;
        lambdax_dim,
        lambday_dim,
        lambdaz_dim,
        x0_dim,
        y0_dim,
        z0_dim,
        sigmax_dim,
        sigmay_dim,
        sigmaz_dim,
        a0,
        branch,
        wavepacketdim,
    ) = state.namelists.wavepacket
    (; lref) = state.constants
    (; k0, k1, j0, j1, i0, i1, io, jo) = state.domain
    (; bvsstrattfc, fc, rhostrattfc) = state.atmosphere
    (; x, y, ztfc) = state.grid
    (; nwm) = state.namelists.wkb

    if lambdax_dim == 0.0
        kk = 0.0
    else
        kk = 2.0 * pi / (lambdax_dim / lref)
    end

    if lambday_dim == 0.0
        ll = 0.0
    else
        ll = 2.0 * pi / (lambday_dim / lref)
    end

    mm = 2.0 * pi / (lambdaz_dim / lref)

    x0 = x0_dim / lref
    y0 = y0_dim / lref
    z0 = z0_dim / lref

    sigmax = sigmax_dim / lref
    sigmay = sigmay_dim / lref
    sigmaz = sigmaz_dim / lref

    kh = sqrt(kk^2.0 + ll^2.0)

    for iwm in 1:nwm, kz in k0:k1, jy in j0:j1, ix in i0:i1
        n2 = bvsstrattfc[ix, jy, kz]
        f = fc[jy]
        f2 = f^2.0
        omega = branch * sqrt((n2 * kh^2.0 + f2 * mm^2.0) / (kh^2.0 + mm^2.0))

        omega2 = omega^2.0

        if wavepacketdim == 1
            deltax = 0.0
            deltay = 0.0
        elseif wavepacketdim == 2
            deltax = x[io + ix] - x0
            deltay = 0.0
        elseif wavepacketdim == 3
            deltax = (x[io + ix] - x0)
            deltay = (y[jo + jy] - y0)
        end

        deltaz = ztfc[ix, jy, kz] - z0

        # Gaussian in z, Cosine in x and y 

        if sigmax == 0.0
            envel = 1.0
        elseif abs(deltax) < sigmax
            envel = 0.5 * (1.0 + cos(deltax * pi / sigmax))
        else
            envel = 0.0
        end

        if sigmay == 0.0
            envel = 1.0 * envel
        elseif abs(deltay) < sigmay
            envel = envel * 0.5 * (1.0 + cos(deltay * pi / sigmay))
        else
            envel = 0.0
        end

        envel = envel * exp(-deltaz^2.0 / (2.0 * sigmaz^2.0))

        bhat = a0 * n2 / mm * envel

        wad_ini[iwm, ix, jy, kz] =
            rhostrattfc[ix, jy, kz] * bhat^2.0 / (2.0 * n2) *
            (f2 * mm^2.0 + n2 * kh^2.0) / (omega * n2 * kh^2.0)

        omi_ini[iwm, ix, jy, kz] = omega
        wnk_ini[iwm, ix, jy, kz] = kk
        wnl_ini[iwm, ix, jy, kz] = ll
        wnm_ini[iwm, ix, jy, kz] = mm
    end

    return
end