"""
```julia
turbulence_integration!(state::State, dt::AbstractFloat)
```

Integrate the turbulence energies by dispatching to the scheme-specific method.

```julia 
turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    turbulence_scheme::NoTurbulence,
)
```

Return for configurations without turbulence parameterization.

```julia 
turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    turbulence_scheme::TKEScheme,
)
```

Integrate the turbulent kinetic energy by dispatching to the specific operations.

```julia 
turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    process::Dissipation,
)
```

Integrate the dissipation contribution of the prognostic equation for the turbulent kinetic energy by dispatching to the model-specific method.

```julia 
turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    process::Dissipation,
    model::Union{PseudoIncompressible, Compressible},
)
```

Integrate the dissipation contribution of the prognostic equation for the turbulent kinetic energy for configurations in psueod-incompressible or compressible mode.

```julia 
turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    process::Dissipation,
    model::Boussinesq,
)
```

Integrate the dissipation contribution of the prognostic equation for the turbulent kinetic energy for configurations in Boussinesq mode.

```julia 
turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    process::Advection,
)
```

Integrate the advection, shear production, and buoyancy contribution terms in the prognostic equation for the turbulent kinetic energy with a Runge-Kutta time step.

```julia
turbulence_integration!(state::State, dt::AbstractFloat, process::Diffusion)
```

Integrate the turbulent diffusion term in the prognostic equation for the turbulent kinetic energy.

# Arguements:

  - `state`: Model state. 

  - `dt`: Time step.

  - `turbulence_scheme`: General turbulence parameterization configuration.

  - `process`: Terms in the prognostic equations.

  - `model`: Dynamic equations.
"""
function turbulence_integration! end

function turbulence_integration!(state::State, dt::AbstractFloat)
    (; turbulence_scheme) = state.namelists.turbulence

    turbulence_integration!(state, dt, turbulence_scheme)

    return
end

function turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    turbulence_scheme::NoTurbulence,
)
    return
end

function turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    turbulence_scheme::TKEScheme,
)
    check_tke!(state)
    set_boundaries!(state, BoundaryPredictands())

    turbulence_integration!(state, dt * 0.5, Dissipation())

    check_tke!(state)
    set_boundaries!(state, BoundaryPredictands())

    turbulence_integration!(state, dt, Advection())

    check_tke!(state)
    set_boundaries!(state, BoundaryPredictands())

    turbulence_integration!(state, dt, Diffusion())

    check_tke!(state)
    set_boundaries!(state, BoundaryPredictands())

    turbulence_integration!(state, dt * 0.5, Dissipation())

    check_tke!(state)
    set_boundaries!(state, BoundaryPredictands())

    return
end

function turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    process::Dissipation,
)
    (; model) = state.namelists.atmosphere

    turbulence_integration!(state, dt, process, model)
    return
end

function turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    process::Dissipation,
    model::Union{PseudoIncompressible, Compressible},
)
    (; tke) = state.turbulence.turbulencepredictands
    (; ld) = state.turbulence.turbulenceconstants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; rhobar) = state.atmosphere
    (; rho) = state.variables.predictands

    for k in k0:k1, j in j0:j1, i in i0:i1
        tke[i, j, k] =
            1 /
            (
                sqrt(2) * dt / (ld * sqrt(rho[i, j, k] + rhobar[i, j, k])) +
                1 / sqrt(tke[i, j, k])
            )^2.0
    end

    return
end

function turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    process::Dissipation,
    model::Boussinesq,
)
    (; tke) = state.turbulence.turbulencepredictands
    (; ld) = state.turbulence.turbulenceconstants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; rhobar) = state.atmosphere
    (; rhop) = state.variables.predictands

    for k in k0:k1, j in j0:j1, i in i0:i1
        tke[i, j, k] =
            1 /
            (
                sqrt(2) * dt / (ld * sqrt(rhop[i, j, k] + rhobar[i, j, k])) +
                1 / sqrt(tke[i, j, k])
            )^2.0
    end

    return
end

function turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    process::Advection,
)
    (; nstages, stepfrac) = state.time

    for rkstage in 1:nstages
        reconstruct!(state, TKE())

        set_boundaries!(state, BoundaryReconstructions())

        compute_fluxes!(state, TKE())

        set_boundaries!(state, BoundaryFluxes())

        update!(state, dt, rkstage, TKE())

        apply_lhs_sponge!(state, dt, stepfrac[rkstage] * dt, TKE())

        set_boundaries!(state, BoundaryPredictands())
    end

    return
end

function turbulence_integration!(state::State, dt::AbstractFloat, process::Diffusion)
    (; tke) = state.turbulence.turbulencepredictands
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; nbx, nby, nbz) = state.namelists.domain
    (; jac, dz) = state.grid
    (; kek) = state.turbulence.turbulencediffusioncoefficients
    (; ath, bth, cth, fth) = state.variables.auxiliaries

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
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

        ith = i - nbx
        jth = j - nby
        kth = k - nbz

        ath[ith, jth, kth] = -dtdz2 * kekd
        bth[ith, jth, kth] = 1 + dtdz2 * keku + dtdz2 * kekd
        cth[ith, jth, kth] = -dtdz2 * keku

        fth[ith, jth, kth] =
            (1 - dtdz2 * keku - dtdz2 * kekd) * tke[i, j, k] +
            dtdz2 * keku * tke[i, j, k + 1] +
            dtdz2 * kekd * tke[i, j, k - 1]
    end

    thomas_algorithm!(state)

    tke[i0:i1, j0:j1, k0:k1] .= fth

    return
end
