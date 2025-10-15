function turbulence_integration! end

function turbulence_integration!(
    state::State,
    p0::Predictands,
    dt::AbstractFloat,
)
    (; turbulence_scheme) = state.namelists.turbulence

    turbulence_integration!(state, p0, dt, turbulence_scheme)

    return
end

function turbulence_integration!(
    state::State,
    p0::Predictands,
    dt::AbstractFloat,
    turbulence_scheme::NoTurbulence,
)
    return
end

function turbulence_integration!(
    state::State,
    p0::Predictands,
    dt::AbstractFloat,
    turbulence_scheme::AbstractTurbulence,
)
    set_boundaries!(state, BoundaryPredictands())

    check_tke!(state)
    turbulence_integration!(state, p0, dt * 0.5, Dissipation())
    turbulence_integration!(state, p0, dt, Advection())
    turbulence_integration!(state, p0, dt, Diffusion())
    check_tke!(state)
    turbulence_integration!(state, p0, dt * 0.5, Dissipation())

    set_boundaries!(state, BoundaryPredictands())

    return
end

function turbulence_integration!(
    state::State,
    p0::Predictands,
    dt::AbstractFloat,
    process::Dissipation,
)
    (; tke) = state.turbulence.turbulencepredictands
    (; lref) = state.constants
    (; lturb, cepsilon) = state.turbulence.turbulenceconstants
    (; i0, i1, j0, j1, k0, k1) = state.domain

    lturb_ndim = lturb / lref

    tke[i0:i1, j0:j1, k0:k1] .=
        4.0 ./
        (cepsilon / lturb_ndim * dt .+ 2 ./ sqrt.(tke[i0:i1, j0:j1, k0:k1])) .^
        2
    return
end

function turbulence_integration!(
    state::State,
    p0::Predictands,
    dt::AbstractFloat,
    process::Advection,
)
    (; nstages, stepfrac) = state.time
    (; turbulence_scheme) = state.namelists.turbulence

    for rkstage in 1:nstages
        reconstruct!(state)
        set_boundaries!(state, BoundaryReconstructions())

        compute_fluxes!(state, p0)

        set_boundaries!(state, BoundaryFluxes())

        save_backups!(state, :rho)

        update!(state, p0, dt, rkstage, turbulence_scheme)
        apply_lhs_sponge!(
            state,
            dt,
            stepfrac[rkstage] * dt,
            turbulence_scheme,
        )
        set_boundaries!(state, BoundaryPredictands())
    end

    return
end

function turbulence_integration!(
    state::State,
    p0::Predictands,
    dt::AbstractFloat,
    process::Diffusion,
)
    (; tke) = state.turbulence.turbulencepredictands
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, dz) = state.grid
    (; kek, athomas, bthomas, cthomas, fthomas, qthomas, fthomas) =
        state.turbulence.turbulenceauxiliaries

    dtdz2 = dt / (2.0 * dz^2.0)

    athomas .= 0.0
    bthomas .= 0.0
    cthomas .= 0.0
    fthomas .= 0.0
    qthomas .= 0.0

    for j in j0:j1, i in i0:i1
        for k in k0:k1
            kekd =
                (
                    jac[i, j, k - 1] * kek[i, j, k] +
                    jac[i, j, k] * kek[i, j, k - 1]
                ) / (jac[i, j, k - 1] + jac[i, j, k])
            keku =
                (
                    jac[i, j, k + 1] * kek[i, j, k] +
                    jac[i, j, k] * kek[i, j, k + 1]
                ) / (jac[i, j, k + 1] + jac[i, j, k])
            athomas[i, j, k] = -dtdz2 * kekd
            bthomas[i, j, k] = 1 + dtdz2 * keku + dtdz2 * kekd
            cthomas[i, j, k] = -dtdz2 * keku

            fthomas[i, j, k] =
                (1 - dtdz2 * keku - dtdz2 * kekd) * tke[i, j, k] +
                dtdz2 * keku * tke[i, j, k + 1] +
                dtdz2 * kekd * tke[i, j, k - 1]

        end

        athomas[i, j, k0] = 0.0
        cthomas[i, j, k1] = 0.0

        qthomas[i, j, k0] = -cthomas[i, j, k0] / bthomas[i, j, k0]
        fthomas[i, j, k0] = fthomas[i, j, k0] / bthomas[i, j, k0]

        for k in (k0 + 1):k1
            p =
                1.0 /
                (bthomas[i, j, k] + athomas[i, j, k] * qthomas[i, j, k - 1])
            qthomas[i, j, k] = -cthomas[i, j, k] * p
            fthomas[i, j, k] =
                (fthomas[i, j, k] - athomas[i, j, k] * fthomas[i, j, k - 1]) * p
        end

        for k in (k1 - 1):-1:k0
            fthomas[i, j, k] =
                fthomas[i, j, k] + qthomas[i, j, k] * fthomas[i, j, k + 1]
        end

        tke[i, j, :] .= fthomas[i, j, :]
    end

    return
end