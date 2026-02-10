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
    turbulence_scheme::TKEScheme,
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
    (; ld) = state.turbulence.turbulenceconstants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; rhobar) = state.atmosphere

    for k in k0:k1, j in j0:j1, i in i0:i1
        tke[i, j, k] =
            (
                dt / (ld * sqrt(rhobar[i, j, k])) +
                1 / sqrt(2.0 * tke[i, j, k])
            )^(-2)
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
    (; nx, ny, nz, i0, i1, j0, j1, k0, k1) = state.domain
    (; nbz) = state.namelists.domain
    (; jac, dz) = state.grid
    (; kek) = state.turbulence.turbulencediffusioncoefficients

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)
    ath, bth, cth, fth, qth, pth, fth_bc, qth_bc =
        thomas_arrays(state, nx, ny, nz)

    ii = i0:i1
    jj = j0:j1
    kk = k0:k1

    @ivy for k in 1:nz
        knbz = k + nbz
        kekd =
            (
                jac[ii, jj, knbz - 1] .* kek[ii, jj, knbz] .+
                jac[ii, jj, knbz] .* kek[ii, jj, knbz - 1]
            ) ./ (jac[ii, jj, knbz - 1] + jac[ii, jj, knbz])
        keku =
            (
                jac[ii, jj, knbz + 1] .* kek[ii, jj, knbz] .+
                jac[ii, jj, knbz] .* kek[ii, jj, knbz + 1]
            ) ./ (jac[ii, jj, knbz + 1] .+ jac[ii, jj, knbz])
        ath[:, :, k] .= .-dtdz2 .* kekd
        bth[:, :, k] .= 1 .+ dtdz2 .* keku .+ dtdz2 .* kekd
        cth[:, :, k] .= .-dtdz2 .* keku

        fth[:, :, k] =
            (1 .- dtdz2 .* keku .- dtdz2 .* kekd) .* tke[ii, jj, knbz] .+
            dtdz2 .* keku .* tke[ii, jj, knbz + 1] .+
            dtdz2 .* kekd .* tke[ii, jj, knbz - 1]
    end

    thomas_algorithm!(state, ath, bth, cth, fth, qth, pth, fth_bc, qth_bc)

    tke[ii, jj, kk] .= fth

    return
end
