"""
```julia
turbulence_integration!(state::State, dt::AbstractFloat)
```

Integrate the turbulence energies by dispatching to the scheme-specific method.

```julia
turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    turbulence_scheme::Val{:NoTurbulence},
)
```

Return for configurations without turbulence parameterization.

```julia
turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    turbulence_scheme::Val{:TKEScheme},
)
```

Integrate the turbulent kinetic energy by dispatching to the specific operations.

```julia
turbulence_integration!(state::State, dt::AbstractFloat, process::Dissipation)
```

Integrate the dissipation contribution of the prognostic equation for the turbulent kinetic energy.

The dissipation step is given by

```math
\\left(\\rho e_\\mathrm{k}\\right) \\rightarrow \\left(\\frac{\\sqrt{2}\\Delta t}{l_d \\sqrt{\\rho}} + \\frac{1}{\\sqrt{\\rho e_\\mathrm{k}}} \\right)^{-2},
```

with turbulent mixing length ``l_d`` stored in `state.turbulence.turbulenceconstants.ld`.

```julia
turbulence_integration!(state::State, dt::AbstractFloat, process::Advection)
```

Integrate the advection, shear production, and buoyancy contribution terms in the prognostic equation for the turbulent kinetic energy with a Runge-Kutta time step.

At each Runge-Kutta stage, the mass-weighted turbulent kinetic energy is first reconstructed and its advective fluxes are calculated. Subsequently, the TKE is updated with its shear and buoyancy production terms, followed immediately by an implicit Euler step (the size of which is the fractional time step at the current Runge-Kutta stage) that accounts for the Rayleigh-damping imposed by the LHS sponge.

```julia
turbulence_integration!(state::State, dt::AbstractFloat, process::Diffusion)
```

Integrate the turbulent diffusion term in the prognostic equation for the turbulent kinetic energy using a Thomas algorithm.

The prognostic equation

```math
\\frac{\\partial e_\\mathrm{k}}{\\partial t} = \\frac{1}{J}\\frac{\\partial}{\\partial \\hat{z}}\\left(\\frac{K_\\mathrm{e_\\mathrm{k}}}{J}\\frac{\\partial e_\\mathrm{k}}{\\partial \\hat{z}}\\right)
```

is solved using the Crank-Nicolson scheme, where the system

```math
a_k \\left(e_\\mathrm{k}\\right)_{k-1}^{n+1} + b_k \\left(e_\\mathrm{k}\\right)_k^{n+1} + c_k \\left(e_\\mathrm{k}\\right)_{k+1}^{n+1} = f_k
```

is solved using a Thomas tridiagonal solver, with ``\\tilde{\\mathcal{K}}_{e_\\mathrm{k}} = \\frac{K_\\mathrm{e_\\mathrm{k}}}{J}`` and

```math
\\begin{align*}
    a_k = & -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\tilde{\\mathcal{K}}_{e_\\mathrm{k},k-1/2}}{J_k}  , \\\\
    b_k = & 1 + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\tilde{\\mathcal{K}}_{e_\\mathrm{k},k+1/2}}{J_k} + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\tilde{\\mathcal{K}}_{e_\\mathrm{k},k-1/2}}{J_k}, \\\\
    c_k = & -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\tilde{\\mathcal{K}}_{e_\\mathrm{k},k+1/2}}{J_k}  , \\\\
    f_k = & \\left( 1 - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\tilde{\\mathcal{K}}_{e_\\mathrm{k},k+1/2}}{J_k}  - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\tilde{\\mathcal{K}}_{e_\\mathrm{k},k-1/2}}{J_k}\\right) \\left(e_\\mathrm{k}\\right)_k^{n} \\\\
    & + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\tilde{\\mathcal{K}}_{e_\\mathrm{k},k+1/2}}{J_k} \\left(e_\\mathrm{k}\\right)_{k+1}^{n} + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\tilde{\\mathcal{K}}_{e_\\mathrm{k},k-1/2}}{J_k}  \\left(e_\\mathrm{k}\\right)_{k-1}^{n}.
\\end{align*}
```

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `turbulence_scheme`: General turbulence parameterization configuration.

  - `process`: Terms in the prognostic equations.

# See also

  - [`PinCFlow.Update.check_tke!`](@ref)

  - [`PinCFlow.Boundaries.set_boundaries!`](@ref)
"""
function turbulence_integration! end

function turbulence_integration!(state::State, dt::AbstractFloat)
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme turbulence_integration!(
        state,
        dt,
        Val(turbulence_scheme),
    )

    return
end

function turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    turbulence_scheme::Val{:NoTurbulence},
)
    return
end

function turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    turbulence_scheme::Val{:TKEScheme},
)
    check_tke!(state)
    set_boundaries!(state, BoundaryPredictands(), TKE())

    turbulence_integration!(state, dt * 0.5, Dissipation())

    check_tke!(state)
    set_boundaries!(state, BoundaryPredictands(), TKE())

    turbulence_integration!(state, dt, Advection())

    turbulence_integration!(state, dt, Diffusion())

    check_tke!(state)
    set_boundaries!(state, BoundaryPredictands(), TKE())

    turbulence_integration!(state, dt * 0.5, Dissipation())

    check_tke!(state)
    set_boundaries!(state, BoundaryPredictands(), TKE())

    save_turbulence_backups!(state)

    return
end

@ivy function turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    process::Dissipation,
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

@ivy function turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    process::Advection,
)
    (; nstages, stepfrac) = state.time

    for rkstage in 1:nstages
        reconstruct!(state, TKE())

        set_boundaries!(state, BoundaryReconstructions(), TKE())

        compute_fluxes!(state, TKE())

        set_boundaries!(state, BoundaryFluxes(), TKE())

        update!(state, dt, rkstage, TKE())

        apply_lhs_sponge!(state, dt, stepfrac[rkstage] * dt, TKE())

        check_tke!(state)
        set_boundaries!(state, BoundaryPredictands(), TKE())
    end

    return
end

@ivy function turbulence_integration!(
    state::State,
    dt::AbstractFloat,
    process::Diffusion,
)
    (; tke) = state.turbulence.turbulencepredictands
    (; rho) = state.variables.predictands
    (; rhobar) = state.atmosphere
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, dz) = state.grid
    (; ath, bth, cth, fth) = state.variables.auxiliaries

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)

    for k in k0:k1, j in j0:j1, i in i0:i1
        kekd =
            (
                jac[i, j, k - 1] * (
                    turbulence_diffusion_coefficient(state, i, j, k, KEK()) /
                    jac[i, j, k]
                ) +
                jac[i, j, k] * (
                    turbulence_diffusion_coefficient(
                        state,
                        i,
                        j,
                        k - 1,
                        KEK(),
                    ) / jac[i, j, k - 1]
                )
            ) / (jac[i, j, k - 1] + jac[i, j, k])
        keku =
            (
                jac[i, j, k + 1] * (
                    turbulence_diffusion_coefficient(state, i, j, k, KEK()) /
                    jac[i, j, k]
                ) +
                jac[i, j, k] * (
                    turbulence_diffusion_coefficient(
                        state,
                        i,
                        j,
                        k + 1,
                        KEK(),
                    ) / jac[i, j, k + 1]
                )
            ) / (jac[i, j, k + 1] + jac[i, j, k])

        ith = i - i0 + 1
        jth = j - j0 + 1
        kth = k - k0 + 1

        ath[ith, jth, kth] = -dtdz2 / jac[i, j, k] * kekd
        bth[ith, jth, kth] =
            1 + dtdz2 / jac[i, j, k] * keku + dtdz2 / jac[i, j, k] * kekd
        cth[ith, jth, kth] = -dtdz2 / jac[i, j, k] * keku

        fth[ith, jth, kth] =
            (1 - dtdz2 / jac[i, j, k] * keku - dtdz2 / jac[i, j, k] * kekd) *
            tke[i, j, k] / (rho[i, j, k] + rhobar[i, j, k]) +
            dtdz2 / jac[i, j, k] * keku * tke[i, j, k + 1] /
            (rho[i, j, k + 1] + rhobar[i, j, k + 1]) +
            dtdz2 / jac[i, j, k] * kekd * tke[i, j, k - 1] /
            (rho[i, j, k - 1] + rhobar[i, j, k - 1])
    end

    thomas_algorithm!(state)

    tke[i0:i1, j0:j1, k0:k1] .=
        fth .* (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])

    return
end