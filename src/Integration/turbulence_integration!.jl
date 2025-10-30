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
    check_tke!(state)
    set_boundaries!(state, BoundaryPredictands())

    turbulence_integration!(state, p0, dt * 0.5, Dissipation())

    check_tke!(state)
    set_boundaries!(state, BoundaryPredictands())

    turbulence_integration!(state, p0, dt, Advection())

    check_tke!(state)
    set_boundaries!(state, BoundaryPredictands())

    turbulence_integration!(state, p0, dt, Diffusion())

    check_tke!(state)
    set_boundaries!(state, BoundaryPredictands())

    turbulence_integration!(state, p0, dt * 0.5, Dissipation())

    check_tke!(state)
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
    (; lturb_ndim) = state.turbulence.turbulenceconstants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; rhobar) = state.atmosphere

    for k in k0:k1, j in j0:j1, i in i0:i1
        tke[i, j, k] =
            4.0 /
            (
                2.0 * dt / (lturb_ndim * sqrt(rhobar[i, j, k])) +
                2.0 / sqrt(tke[i, j, k])
            )^2.0
    end

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
        reconstruct!(state, turbulence_scheme)
        set_boundaries!(state, BoundaryReconstructions())

        compute_fluxes!(state, p0, turbulence_scheme)

        set_boundaries!(state, BoundaryFluxes())

        update!(state, p0, dt, rkstage, turbulence_scheme)
        apply_lhs_sponge!(state, dt, stepfrac[rkstage] * dt, turbulence_scheme)
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
    (; nz, i0, i1, j0, j1, k0, k1) =
        state.domain
    (; nbz) = state.namelists.domain
    (; jac, dz) = state.grid
    (; ath, bth, cth, fth, qth, fth, pth, qth_bc, fth_bc) =
        state.variables.auxiliaries
    (; kek) = state.turbulence.turbulencediffusioncoefficients

    dtdz2 = dt / (2.0 * dz^2.0)

    ath .= 0.0
    bth .= 0.0
    cth .= 0.0
    fth .= 0.0
    qth .= 0.0
    pth .= 0.0
    qth_bc .= 0.0
    fth_bc .= 0.0

    for k in 1:nz
        kekd =
            (
                jac[i0:i1, j0:j1, k + nbz - 1] .* kek[i0:i1, j0:j1, k + nbz] .+
                jac[i0:i1, j0:j1, k + nbz] .* kek[i0:i1, j0:j1, k + nbz - 1]
            ) ./ (jac[i0:i1, j0:j1, k + nbz - 1] + jac[i0:i1, j0:j1, k + nbz])
        keku =
            (
                jac[i0:i1, j0:j1, k + nbz + 1] .* kek[i0:i1, j0:j1, k + nbz] .+
                jac[i0:i1, j0:j1, k + nbz] .* kek[i0:i1, j0:j1, k + nbz + 1]
            ) ./ (jac[i0:i1, j0:j1, k + 1] .+ jac[i0:i1, j0:j1, k])
        ath[:, :, k] .= .-dtdz2 .* kekd
        bth[:, :, k] .= 1 .+ dtdz2 .* keku .+ dtdz2 .* kekd
        cth[:, :, k] .= .-dtdz2 .* keku

        fth[:, :, k] =
            (1 .- dtdz2 .* keku .- dtdz2 .* kekd) .*
            tke[i0:i1, j0:j1, k + nbz] .+
            dtdz2 .* keku .* tke[i0:i1, j0:j1, k + nbz + 1] .+
            dtdz2 .* kekd .* tke[i0:i1, j0:j1, k + nbz - 1]
    end

    thomas_algorithm!(state, ath, bth, cth, fth, qth, pth, fth_bc, qth_bc)

    tke[i0:i1, j0:j1, k0:k1] .= fth

    return
end
