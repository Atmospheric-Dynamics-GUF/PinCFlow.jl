"""
```julia
correct!(
    state::State,
    dt::AbstractFloat,
    rayleigh_factor::AbstractFloat,
)
```

Correct the Exner-pressure, wind and buoyancy (density fluctuations) such that the divergence constraint is satisfied, using the Exner-pressure differences obtained from the solution to the Poisson problem.

This method calls specific methods for the correction of each individual variable.

```julia
correct!(
    state::State,
    dt::AbstractFloat,
    variable::U,
    rayleigh_factor::AbstractFloat,
)
```

Correct the zonal wind to account for the pressure differences obtained from the solution to the Poisson problem.

```julia
correct!(
    state::State,
    dt::AbstractFloat,
    variable::V,
    rayleigh_factor::AbstractFloat,
)
```

Correct the meridional wind to account for the pressure differences obtained from the solution to the Poisson problem.

```julia
correct!(
    state::State,
    dt::AbstractFloat,
    variable::W,
    rayleigh_factor::AbstractFloat,
)
```

Correct the transformed vertical wind to account for the pressure differences obtained from the solution to the Poisson problem.

```julia
correct!(
    state::State,
    dt::AbstractFloat,
    variable::RhoP,
    rayleigh_factor::AbstractFloat,
)
```

Correct the buoyancy (density fluctuations) to account for the pressure differences obtained from the solution to the Poisson problem.

```julia
correct!(state::State, variable::PiP)
```

Update the Exner-pressure fluctuations with the differences obtained from the solution to the Poisson problem.

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `variable`: Variable to correct.

  - `rayleigh_factor`: Factor by which the Rayleigh-damping coefficient is multiplied.

# See also

  - [`PinCFlow.Update.compute_pressure_gradient`](@ref)

  - [`PinCFlow.Update.compute_compressible_wind_factor`](@ref)

  - [`PinCFlow.Update.compute_buoyancy_factor`](@ref)
"""
function correct! end

function correct!(
    state::State,
    dt::AbstractFloat,
    rayleigh_factor::AbstractFloat,
)
    correct!(state, dt, U(), rayleigh_factor)
    correct!(state, dt, V(), rayleigh_factor)
    correct!(state, dt, W(), rayleigh_factor)
    correct!(state, dt, RhoP(), rayleigh_factor)
    correct!(state, PiP())
    return
end

function correct!(
    state::State,
    dt::AbstractFloat,
    variable::U,
    rayleigh_factor::AbstractFloat,
)
    (; spongelayer, sponge_uv) = state.namelists.sponge
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; kr_sp_tfc) = state.sponge
    (; corx) = state.poisson.correction
    (; dpip) = state.variables.increments
    (; u) = state.variables.predictands

    kz0 = k0
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for k in kz0:kz1, j in j0:j1, i in (i0 - 1):i1
        factor = 1.0

        if spongelayer && sponge_uv
            factor +=
                dt *
                0.5 *
                (kr_sp_tfc[i, j, k] + kr_sp_tfc[i + 1, j, k]) *
                rayleigh_factor
        end

        gradient = compute_pressure_gradient(state, dpip, (i, j, k), U())

        corx[i, j, k] = dt / factor * gradient

        jpedger = compute_compressible_wind_factor(state, (i, j, k), U())

        u[i, j, k] -= jpedger * corx[i, j, k]
    end

    return
end

function correct!(
    state::State,
    dt::AbstractFloat,
    variable::V,
    rayleigh_factor::AbstractFloat,
)
    (; spongelayer, sponge_uv) = state.namelists.sponge
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; kr_sp_tfc) = state.sponge
    (; cory) = state.poisson.correction
    (; dpip) = state.variables.increments
    (; v) = state.variables.predictands

    kz0 = k0
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for k in kz0:kz1, j in (j0 - 1):j1, i in i0:i1
        factor = 1.0

        if spongelayer && sponge_uv
            factor +=
                dt *
                0.5 *
                (kr_sp_tfc[i, j, k] + kr_sp_tfc[i, j + 1, k]) *
                rayleigh_factor
        end

        gradient = compute_pressure_gradient(state, dpip, (i, j, k), V())

        cory[i, j, k] = dt / factor * gradient

        jpedgef = compute_compressible_wind_factor(state, (i, j, k), V())

        v[i, j, k] -= jpedgef * cory[i, j, k]
    end

    return
end

function correct!(
    state::State,
    dt::AbstractFloat,
    variable::W,
    rayleigh_factor::AbstractFloat,
)
    (; spongelayer) = state.namelists.sponge
    (; zboundaries) = state.namelists.setting
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met) = state.grid
    (; bvsstrattfc) = state.atmosphere
    (; kr_sp_w_tfc) = state.sponge
    (; corx, cory) = state.poisson.correction
    (; dpip) = state.variables.increments
    (; w) = state.variables.predictands

    if zboundaries != SolidWallBoundaries()
        error("Error in correct!: Unknown zboundaries!")
    end

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 - 1 : k1

    for k in kz0:kz1, j in j0:j1, i in i0:i1
        factor = 1.0

        if spongelayer
            factor +=
                dt * (
                    jac[i, j, k + 1] * kr_sp_w_tfc[i, j, k] +
                    jac[i, j, k] * kr_sp_w_tfc[i, j, k + 1]
                ) / (jac[i, j, k] + jac[i, j, k + 1]) * rayleigh_factor
        end

        bvsstratedgeu =
            (
                jac[i, j, k + 1] * bvsstrattfc[i, j, k] +
                jac[i, j, k] * bvsstrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        gradient = compute_pressure_gradient(state, dpip, (i, j, k), W())

        jpedgeu = compute_compressible_wind_factor(state, (i, j, k), W())
        fw = compute_buoyancy_factor(state, (i, j, k), W())

        w[i, j, k] +=
            -dt / (factor + fw * bvsstratedgeu * dt^2.0) * jpedgeu * gradient -
            1.0 / (factor + fw * bvsstratedgeu * dt^2.0) *
            fw *
            bvsstratedgeu *
            dt^2.0 *
            jpedgeu *
            0.5 *
            (
                jac[i, j, k + 1] * (
                    met[i, j, k, 1, 3] * (corx[i, j, k] + corx[i - 1, j, k]) +
                    met[i, j, k, 2, 3] * (cory[i, j, k] + cory[i, j - 1, k])
                ) +
                jac[i, j, k] * (
                    met[i, j, k + 1, 1, 3] *
                    (corx[i, j, k + 1] + corx[i - 1, j, k + 1]) +
                    met[i, j, k + 1, 2, 3] *
                    (cory[i, j, k + 1] + cory[i, j - 1, k + 1])
                )
            ) / (jac[i, j, k] + jac[i, j, k + 1])
    end

    return
end

function correct!(
    state::State,
    dt::AbstractFloat,
    variable::RhoP,
    rayleigh_factor::AbstractFloat,
)
    (; nbz) = state.namelists.domain
    (; spongelayer) = state.namelists.sponge
    (; zboundaries) = state.namelists.setting
    (; g_ndim) = state.constants
    (; sizezz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met) = state.grid
    (; rhostrattfc, bvsstrattfc) = state.atmosphere
    (; kr_sp_w_tfc) = state.sponge
    (; corx, cory) = state.poisson.correction
    (; dpip) = state.variables.increments
    (; rho, rhop) = state.variables.predictands

    for k in k0:k1, j in j0:j1, i in i0:i1
        factor = 1.0

        if spongelayer
            factor += dt * kr_sp_w_tfc[i, j, k] * rayleigh_factor
        end

        lower_gradient =
            compute_pressure_gradient(state, dpip, (i, j, k - 1), W())
        upper_gradient = compute_pressure_gradient(state, dpip, (i, j, k), W())

        if ko + k == k0 && zboundaries == SolidWallBoundaries()
            lower_gradient = 0.0
        elseif ko + k == sizezz - nbz && zboundaries == SolidWallBoundaries()
            upper_gradient = 0.0
        end

        gradient = 0.5 * (lower_gradient + upper_gradient)

        fb = compute_buoyancy_factor(state, (i, j, k), RhoP())
        db =
            -1.0 / (factor + fb * bvsstrattfc[i, j, k] * dt^2.0) * (
                -fb * bvsstrattfc[i, j, k] * dt^2.0 * jac[i, j, k] * gradient +
                fb *
                bvsstrattfc[i, j, k] *
                dt *
                jac[i, j, k] *
                factor *
                0.5 *
                (
                    met[i, j, k, 1, 3] * (corx[i, j, k] + corx[i - 1, j, k]) +
                    met[i, j, k, 2, 3] * (cory[i, j, k] + cory[i, j - 1, k])
                )
            )

        rhop[i, j, k] -= (rho[i, j, k] + rhostrattfc[i, j, k]) / g_ndim * db
    end

    return
end

function correct!(state::State, variable::PiP)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; pip) = state.variables.predictands
    (; dpip) = state.variables.increments

    i = (i0 - 1):(i1 + 1)
    j = (j0 - 1):(j1 + 1)
    k = (k0 - 1):(k1 + 1)

    @views pip[i, j, k] .= pip[i, j, k] .+ dpip[i, j, k]

    return
end
