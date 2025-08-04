"""
```julia
update!(state::State, dt::AbstractFloat, m::Integer, variable::Rho)
```

Update the density, depending on whether or not the model is Boussinesq.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::Rho,
    model::Boussinesq,
)
```

Return in Boussinesq mode (the density is constant).

```julia
update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::Rho,
    model::AbstractModel,
)
```

Update the density with a Runge-Kutta step on the left-hand side of the equation (the right-hand side is zero).

```julia
update!(state::State, dt::AbstractFloat, m::Integer, variable::RhoP, side::LHS)
```

Update the density fluctuations with a Runge-Kutta step on the left-hand-side of the equation.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::RhoP,
    side::RHS,
    integration::Explicit,
)
```

Update the density fluctuations with an explicit Euler step the on right-hand side of the equation, without the Rayleigh-damping term.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::RhoP,
    side::RHS,
    integration::Implicit,
    facray::AbstractFloat,
)
```

Update the density fluctuations with an implicit Euler step on the right-hand side of the equation.

```julia
update!(state::State, dt::AbstractFloat, m::Integer, variable::U, side::LHS)
```

Update the zonal momentum with a Runge-Kutta step on the left-hand side of the equation.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::U,
    side::RHS,
    integration::Explicit,
)
```

Update the zonal wind with an explicit Euler step on the right-hand side of the equation, without the Rayleigh-damping term.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::U,
    side::RHS,
    integration::Implicit,
    facray::AbstractFloat,
)
```

Update the zonal wind with an implicit Euler step on the right-hand side of the equation.

```julia
update!(state::State, dt::AbstractFloat, m::Integer, variable::V, side::LHS)
```

Update the meridional momentum with a Runge-Kutta step on the left-hand side of the equation.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::V,
    side::RHS,
    integration::Explicit,
)
```

Update the meridional wind with an explicit Euler step on the right-hand side of the equation, without the Rayleigh-damping term.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::V,
    side::RHS,
    integration::Implicit,
    facray::AbstractFloat,
)
```

Update the meridional wind with an implicit Euler step on the right-hand side of the equation.

```julia
update!(state::State, dt::AbstractFloat, m::Integer, variable::W, side::LHS)
```

Update the transformed vertical momentum with a Runge-Kutta step on the left-hand side of the equation.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::W,
    side::RHS,
    integration::Explicit,
)
```

Update the transformed vertical wind with an explicit Euler step on the right-hand side of the equation, without the Rayleigh-damping term.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::W,
    side::RHS,
    integration::Implicit,
    facray::AbstractFloat,
)
```

Update the transformed vertical wind with an implicit Euler step on the right-hand side of the equation.

```julia
update!(state::State, dt::AbstractFloat, variable::PiP)
```

Update the Exner-pressure, depending on whether or not the model is compressible.

```julia
update!(state::State, dt::AbstractFloat, variable::PiP, model::AbstractModel)
```

Return in non-compressible modes.

```julia
update!(state::State, dt::AbstractFloat, variable::PiP, model::Compressible)
```

Update the Exner-pressure such that it is synchronized with the updated mass-weighted potential temperature.

```julia
update!(state::State, dt::AbstractFloat, m::Integer, variable::P)
```

Update the mass-weighted potential temperature, depending on whether or not the model is compressible.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::P,
    model::AbstractModel,
)
```

Return in non-compressible modes.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::P,
    model::Compressible,
)
```

Update the mass-weighted potential temperature with a Runge-Kutta step on the left-hand side of the equation (the right-hand side is zero).

```julia
update!(state::State, dt::AbstractFloat, m::Integer, tracersetup::NoTracer)
```

```julia
update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    tracersetup::AbstractTracer,
)
```

```julia
update!(state::State, dt::AbstractFloat, m::Integer, icesetup::AbstractIce)
```

```julia
update!(state::State, dt::AbstractFloat, m::Integer, icesetup::IceOn)
```

```julia
update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    turbulencesetup::NoTurbulence,
)
```

```julia
update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    turbulencesetup::AbstractTurbulence,
)
```

# Arguments

  - `state`: Model state.
  - `dt`: Time step.
  - `m`: Runge-Kutta-stage index.
  - `variable`: Variable to update.
  - `model`: Dynamic equations.
  - `side`: Side of the equation.
  - `integration`: Type of the Euler step.
  - `facray`: Factor by which the Rayleigh-damping coefficient is multiplied.

# See also

  - [`PinCFlow.Update.compute_volume_force`](@ref)
  - [`PinCFlow.Update.compute_compressible_wind_factor`](@ref)
  - [`PinCFlow.Update.compute_vertical_wind`](@ref)
  - [`PinCFlow.Update.compute_compressible_buoyancy_factor`](@ref)
  - [`PinCFlow.Update.transform`](@ref)
"""
function update! end

function update!(state::State, dt::AbstractFloat, m::Integer, variable::Rho)
    (; model) = state.namelists.setting
    update!(state, dt, m, variable, model)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::Rho,
    model::Boussinesq,
)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::Rho,
    model::AbstractModel,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; alphark, betark) = state.time
    (; drho) = state.variables.tendencies
    (; phirho) = state.variables.fluxes
    (; rho) = state.variables.predictands

    if m == 1
        drho .= 0.0
    end

    for k in k0:k1, j in j0:j1, i in i0:i1
        fl = phirho[i - 1, j, k, 1]
        fr = phirho[i, j, k, 1]
        gb = phirho[i, j - 1, k, 2]
        gf = phirho[i, j, k, 2]
        hd = phirho[i, j, k - 1, 3]
        hu = phirho[i, j, k, 3]

        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
        fluxdiff /= jac[i, j, k]

        f = -fluxdiff

        drho[i, j, k] = dt * f + alphark[m] * drho[i, j, k]
        rho[i, j, k] += betark[m] * drho[i, j, k]
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::RhoP,
    side::LHS,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; thetastrattfc) = state.atmosphere
    (; alphark, betark) = state.time
    (; drhop) = state.variables.tendencies
    (; phirhop) = state.variables.fluxes
    (; rhop) = state.variables.predictands

    if m == 1
        drhop .= 0.0
    end

    for k in k0:k1, j in j0:j1, i in i0:i1
        fl = phirhop[i - 1, j, k, 1]
        fr = phirhop[i, j, k, 1]
        gb = phirhop[i, j - 1, k, 2]
        gf = phirhop[i, j, k, 2]
        hd = phirhop[i, j, k - 1, 3]
        hu = phirhop[i, j, k, 3]

        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
        fluxdiff /= jac[i, j, k]

        heating = compute_volume_force(state, (i, j, k), P())

        f = -fluxdiff + heating / thetastrattfc[i, j, k]

        drhop[i, j, k] = dt * f + alphark[m] * drhop[i, j, k]
        rhop[i, j, k] += betark[m] * drhop[i, j, k]
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::RhoP,
    side::RHS,
    integration::Explicit,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; grid) = state
    (; g_ndim) = state.constants
    (; rhostrattfc, bvsstrattfc) = state.atmosphere
    (; predictands) = state.variables
    (; rho, rhop) = predictands

    for k in k0:k1, j in j0:j1, i in i0:i1
        jpu = compute_compressible_wind_factor(state, (i, j, k), W())
        jpd = compute_compressible_wind_factor(state, (i, j, k - 1), W())
        wvrt =
            0.5 * (
                compute_vertical_wind(i, j, k, predictands, grid) / jpu +
                compute_vertical_wind(i, j, k - 1, predictands, grid) / jpd
            )

        buoy = -g_ndim * rhop[i, j, k] / (rho[i, j, k] + rhostrattfc[i, j, k])
        fb = compute_compressible_buoyancy_factor(state, (i, j, k), RhoP())
        buoy -= dt * fb * bvsstrattfc[i, j, k] * wvrt

        rhop[i, j, k] = -buoy * (rho[i, j, k] + rhostrattfc[i, j, k]) / g_ndim
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::RhoP,
    side::RHS,
    integration::Implicit,
    facray::AbstractFloat,
)
    (; nbz) = state.namelists.domain
    (; sizezz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac, met) = state.grid
    (; zboundaries) = state.namelists.setting
    (; spongelayer) = state.namelists.sponge
    (; kr_sp_w_tfc) = state.sponge
    (; kappainv, mainv2, g_ndim) = state.constants
    (; rhostrattfc, pstrattfc, bvsstrattfc) = state.atmosphere
    (; rho, rhop, u, v, pip) = state.variables.predictands
    (; wold) = state.variables.backups

    for k in k0:k1, j in j0:j1, i in i0:i1
        rhow =
            (
                jac[i, j, k + 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhowm =
            (
                jac[i, j, k - 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k - 1]
            ) / (jac[i, j, k] + jac[i, j, k - 1])

        rhow +=
            (
                jac[i, j, k + 1] * rhostrattfc[i, j, k] +
                jac[i, j, k] * rhostrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhowm +=
            (
                jac[i, j, k - 1] * rhostrattfc[i, j, k] +
                jac[i, j, k] * rhostrattfc[i, j, k - 1]
            ) / (jac[i, j, k] + jac[i, j, k - 1])

        # Momentum is predicted before buoyancy in implicit
        # steps.
        jpu = compute_compressible_wind_factor(state, (i, j, k), W())
        jpd = compute_compressible_wind_factor(state, (i, j, k - 1), W())
        wvrt = 0.5 * (wold[i, j, k] / jpu + wold[i, j, k - 1] / jpd)

        # Compute P coefficients.
        pedgeu =
            (
                jac[i, j, k + 1] * pstrattfc[i, j, k] +
                jac[i, j, k] * pstrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        pedged =
            (
                jac[i, j, k - 1] * pstrattfc[i, j, k] +
                jac[i, j, k] * pstrattfc[i, j, k - 1]
            ) / (jac[i, j, k] + jac[i, j, k - 1])

        # Interpolate metric-tensor elements.
        met13edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 1, 3] +
                jac[i, j, k] * met[i, j, k + 1, 1, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met23edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 2, 3] +
                jac[i, j, k] * met[i, j, k + 1, 2, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met33edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 3, 3] +
                jac[i, j, k] * met[i, j, k + 1, 3, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met13edged =
            (
                jac[i, j, k - 1] * met[i, j, k, 1, 3] +
                jac[i, j, k] * met[i, j, k - 1, 1, 3]
            ) / (jac[i, j, k] + jac[i, j, k - 1])
        met23edged =
            (
                jac[i, j, k - 1] * met[i, j, k, 2, 3] +
                jac[i, j, k] * met[i, j, k - 1, 2, 3]
            ) / (jac[i, j, k] + jac[i, j, k - 1])
        met33edged =
            (
                jac[i, j, k - 1] * met[i, j, k, 3, 3] +
                jac[i, j, k] * met[i, j, k - 1, 3, 3]
            ) / (jac[i, j, k] + jac[i, j, k - 1])

        # Interpolate pressure differences.
        piredgeu =
            (
                jac[i + 1, j, k + 1] * pip[i + 1, j, k] +
                jac[i + 1, j, k] * pip[i + 1, j, k + 1]
            ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
        piledgeu =
            (
                jac[i - 1, j, k + 1] * pip[i - 1, j, k] +
                jac[i - 1, j, k] * pip[i - 1, j, k + 1]
            ) / (jac[i - 1, j, k] + jac[i - 1, j, k + 1])
        piredged =
            (
                jac[i + 1, j, k - 1] * pip[i + 1, j, k] +
                jac[i + 1, j, k] * pip[i + 1, j, k - 1]
            ) / (jac[i + 1, j, k] + jac[i + 1, j, k - 1])
        piledged =
            (
                jac[i - 1, j, k - 1] * pip[i - 1, j, k] +
                jac[i - 1, j, k] * pip[i - 1, j, k - 1]
            ) / (jac[i - 1, j, k] + jac[i - 1, j, k - 1])
        pifedgeu =
            (
                jac[i, j + 1, k + 1] * pip[i, j + 1, k] +
                jac[i, j + 1, k] * pip[i, j + 1, k + 1]
            ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
        pibedgeu =
            (
                jac[i, j - 1, k + 1] * pip[i, j - 1, k] +
                jac[i, j - 1, k] * pip[i, j - 1, k + 1]
            ) / (jac[i, j - 1, k] + jac[i, j - 1, k + 1])
        pifedged =
            (
                jac[i, j + 1, k - 1] * pip[i, j + 1, k] +
                jac[i, j + 1, k] * pip[i, j + 1, k - 1]
            ) / (jac[i, j + 1, k] + jac[i, j + 1, k - 1])
        pibedged =
            (
                jac[i, j - 1, k - 1] * pip[i, j - 1, k] +
                jac[i, j - 1, k] * pip[i, j - 1, k - 1]
            ) / (jac[i, j - 1, k] + jac[i, j - 1, k - 1])

        # Compute pressure gradients.
        pigradzedgeu =
            kappainv * mainv2 * pedgeu / rhow * (
                0.5 * met13edgeu * (piredgeu - piledgeu) / dx +
                0.5 * met23edgeu * (pifedgeu - pibedgeu) / dy +
                met33edgeu * (pip[i, j, k + 1] - pip[i, j, k]) / dz
            )
        pigradzedged =
            kappainv * mainv2 * pedged / rhowm * (
                0.5 * met13edged * (piredged - piledged) / dx +
                0.5 * met23edged * (pifedged - pibedged) / dy +
                met33edged * (pip[i, j, k] - pip[i, j, k - 1]) / dz
            )

        # Adjust at boundaries.
        if ko + k == k0 && zboundaries == SolidWallBoundaries()
            pigradzedged = 0.0
        elseif ko + k == sizezz - nbz && zboundaries == SolidWallBoundaries()
            pigradzedgeu = 0.0
        end

        # Interpolate.
        pigrad = 0.5 * (pigradzedgeu + pigradzedged)

        facw = 1.0

        if spongelayer
            facw += dt * kr_sp_w_tfc[i, j, k] * facray
        end

        # Predict buoyancy.
        buoy = -g_ndim * rhop[i, j, k] / (rho[i, j, k] + rhostrattfc[i, j, k])
        jpr = compute_compressible_wind_factor(state, (i, j, k), U())
        jpl = compute_compressible_wind_factor(state, (i - 1, j, k), U())
        jpf = compute_compressible_wind_factor(state, (i, j, k), V())
        jpb = compute_compressible_wind_factor(state, (i, j - 1, k), V())
        fb = compute_compressible_buoyancy_factor(state, (i, j, k), RhoP())
        buoy =
            1.0 / (facw + fb * bvsstrattfc[i, j, k] * dt^2.0) * (
                -fb *
                bvsstrattfc[i, j, k] *
                dt *
                jac[i, j, k] *
                (wvrt - dt * pigrad) +
                facw * buoy +
                fb *
                bvsstrattfc[i, j, k] *
                dt *
                jac[i, j, k] *
                facw *
                0.5 *
                (
                    met[i, j, k, 1, 3] *
                    (u[i, j, k] / jpr + u[i - 1, j, k] / jpl) +
                    met[i, j, k, 2, 3] *
                    (v[i, j, k] / jpf + v[i, j - 1, k] / jpb)
                )
            )

        rhop[i, j, k] = -buoy * (rho[i, j, k] + rhostrattfc[i, j, k]) / g_ndim
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::U,
    side::LHS,
)
    (; alphark, betark) = state.time
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; rhostrattfc, fc) = state.atmosphere
    (; du) = state.variables.tendencies
    (; phiu) = state.variables.fluxes
    (; rhoold, uold) = state.variables.backups
    (; rho, u, v) = state.variables.predictands

    if m == 1
        du .= 0.0
    end

    for k in k0:k1, j in j0:j1, i in (i0 - 1):i1

        # Compute zonal momentum flux divergence.
        fr = phiu[i, j, k, 1]
        fl = phiu[i - 1, j, k, 1]
        gf = phiu[i, j, k, 2]
        gb = phiu[i, j - 1, k, 2]
        hu = phiu[i, j, k, 3]
        hd = phiu[i, j, k - 1, 3]
        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz

        # Adjust zonal momentum flux divergence.
        jacedger = 0.5 * (jac[i, j, k] + jac[i + 1, j, k])
        fluxdiff /= jacedger

        # Explicit integration of Coriolis force in TFC.
        uold[i, j, k] = u[i, j, k]
        if k == k1 && ko + nzz != sizezz
            uold[i, j, k + 1] = u[i, j, k + 1]
        end
        vc = 0.5 * (v[i, j, k] + v[i, j - 1, k])
        vr = 0.5 * (v[i + 1, j, k] + v[i + 1, j - 1, k])
        volforce =
            0.5 *
            fc[j] *
            (
                (rhoold[i, j, k] + rhostrattfc[i, j, k]) * vc +
                (rhoold[i + 1, j, k] + rhostrattfc[i + 1, j, k]) * vr
            )

        # Compute force.
        force = -fluxdiff + volforce

        # Interpolate density.
        rhom_1 = 0.5 * (rhoold[i, j, k] + rhoold[i + 1, j, k])
        rhom = 0.5 * (rho[i, j, k] + rho[i + 1, j, k])
        rhostratedger = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i + 1, j, k])
        rhom_1 += rhostratedger
        rhom += rhostratedger

        # Set velocity and momentum at previous time.
        um_1 = u[i, j, k]
        momm_1 = rhom_1 * um_1

        # Compute tendency.
        du[i, j, k] = dt * force + alphark[m] * du[i, j, k]

        # Update momentum.
        momm = momm_1 + betark[m] * du[i, j, k]

        # Update wind.
        uast = momm / rhom
        u[i, j, k] = uast
    end

    # Return.
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::U,
    side::RHS,
    integration::Explicit,
)
    (; nbz) = state.namelists.domain
    (; zboundaries) = state.namelists.setting
    (; kappainv, mainv2) = state.constants
    (; sizezz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dz, met) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; rho, u, pip) = state.variables.predictands

    for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
        rhou = 0.5 * (rho[i, j, k] + rho[i + 1, j, k])
        rhostratedger = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i + 1, j, k])
        rhou += rhostratedger

        pir = pip[i + 1, j, k]
        pil = pip[i, j, k]

        # Compute values at cell edges.
        pedger = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i + 1, j, k])
        met13edger = 0.5 * (met[i, j, k, 1, 3] + met[i + 1, j, k, 1, 3])

        # Compute pressure gradient component.
        if ko + k == k0 && zboundaries == SolidWallBoundaries()
            piuuedger = 0.5 * (pip[i, j, k + 2] + pip[i + 1, j, k + 2])
            piuedger = 0.5 * (pip[i, j, k + 1] + pip[i + 1, j, k + 1])
            piedger = 0.5 * (pip[i, j, k] + pip[i + 1, j, k])
            pigrad =
                kappainv * mainv2 * pedger / rhou * (
                    (pir - pil) / dx +
                    met13edger *
                    (-piuuedger + 4.0 * piuedger - 3.0 * piedger) *
                    0.5 / dz
                )
        elseif ko + k == sizezz - nbz && zboundaries == SolidWallBoundaries()
            piddedger = 0.5 * (pip[i, j, k - 2] + pip[i + 1, j, k - 2])
            pidedger = 0.5 * (pip[i, j, k - 1] + pip[i + 1, j, k - 1])
            piedger = 0.5 * (pip[i, j, k] + pip[i + 1, j, k])
            pigrad =
                kappainv * mainv2 * pedger / rhou * (
                    (pir - pil) / dx +
                    met13edger *
                    (piddedger - 4.0 * pidedger + 3.0 * piedger) *
                    0.5 / dz
                )
        else
            piuedger = 0.5 * (pip[i, j, k + 1] + pip[i + 1, j, k + 1])
            pidedger = 0.5 * (pip[i, j, k - 1] + pip[i + 1, j, k - 1])
            pigrad =
                kappainv * mainv2 * pedger / rhou * (
                    (pir - pil) / dx +
                    met13edger * (piuedger - pidedger) * 0.5 / dz
                )
        end

        volfcx = compute_volume_force(state, (i, j, k), U())

        uhorx = u[i, j, k]

        # Update wind.
        jpr = compute_compressible_wind_factor(state, (i, j, k), U())
        uast = uhorx + dt * (-pigrad + volfcx / rhou) * jpr
        u[i, j, k] = uast
    end

    # Return.
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::U,
    side::RHS,
    integration::Implicit,
    facray::AbstractFloat,
)
    (; nbz) = state.namelists.domain
    (; zboundaries) = state.namelists.setting
    (; spongelayer, sponge_uv) = state.namelists.sponge
    (; kappainv, mainv2) = state.constants
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dz, met) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; kr_sp_tfc) = state.sponge
    (; rho, u, pip) = state.variables.predictands

    kz0 = k0
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for k in kz0:kz1, j in j0:j1, i in (i0 - 1):i1
        rhou = 0.5 * (rho[i, j, k] + rho[i + 1, j, k])
        rhostratedger = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i + 1, j, k])
        rhou += rhostratedger

        pir = pip[i + 1, j, k]
        pil = pip[i, j, k]

        # Compute values at cell edges.
        pedger = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i + 1, j, k])
        met13edger = 0.5 * (met[i, j, k, 1, 3] + met[i + 1, j, k, 1, 3])

        # Compute pressure gradient component.
        if ko + k == k0 && zboundaries == SolidWallBoundaries()
            piuuedger = 0.5 * (pip[i, j, k + 2] + pip[i + 1, j, k + 2])
            piuedger = 0.5 * (pip[i, j, k + 1] + pip[i + 1, j, k + 1])
            piedger = 0.5 * (pip[i, j, k] + pip[i + 1, j, k])
            pigradx =
                kappainv * mainv2 * pedger / rhou * (
                    (pir - pil) / dx +
                    met13edger *
                    (-piuuedger + 4.0 * piuedger - 3.0 * piedger) *
                    0.5 / dz
                )
        elseif ko + k == sizezz - nbz && zboundaries == SolidWallBoundaries()
            piddedger = 0.5 * (pip[i, j, k - 2] + pip[i + 1, j, k - 2])
            pidedger = 0.5 * (pip[i, j, k - 1] + pip[i + 1, j, k - 1])
            piedger = 0.5 * (pip[i, j, k] + pip[i + 1, j, k])
            pigradx =
                kappainv * mainv2 * pedger / rhou * (
                    (pir - pil) / dx +
                    met13edger *
                    (piddedger - 4.0 * pidedger + 3.0 * piedger) *
                    0.5 / dz
                )
        else
            piuedger = 0.5 * (pip[i, j, k + 1] + pip[i + 1, j, k + 1])
            pidedger = 0.5 * (pip[i, j, k - 1] + pip[i + 1, j, k - 1])
            pigradx =
                kappainv * mainv2 * pedger / rhou * (
                    (pir - pil) / dx +
                    met13edger * (piuedger - pidedger) * 0.5 / dz
                )
        end

        volfcx = compute_volume_force(state, (i, j, k), U())

        uhorx = u[i, j, k]

        facu = 1.0

        if spongelayer && sponge_uv
            facu +=
                dt *
                0.5 *
                (kr_sp_tfc[i, j, k] + kr_sp_tfc[i + 1, j, k]) *
                facray
        end

        # Update wind.
        jpr = compute_compressible_wind_factor(state, (i, j, k), U())
        uast = 1.0 / facu * (uhorx + dt * (-pigradx + volfcx / rhou) * jpr)
        u[i, j, k] = uast
    end

    # Return.
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::V,
    side::LHS,
)
    (; alphark, betark) = state.time
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; rhostrattfc, fc) = state.atmosphere
    (; dv) = state.variables.tendencies
    (; phiv) = state.variables.fluxes
    (; rhoold, uold, vold) = state.variables.backups
    (; rho, v) = state.variables.predictands

    if m == 1
        dv .= 0.0
    end

    for k in k0:k1, j in (j0 - 1):j1, i in i0:i1

        # Compute meridional momentum flux divergence.
        fr = phiv[i, j, k, 1]
        fl = phiv[i - 1, j, k, 1]
        gf = phiv[i, j, k, 2]
        gb = phiv[i, j - 1, k, 2]
        hu = phiv[i, j, k, 3]
        hd = phiv[i, j, k - 1, 3]
        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz

        # Adjust meridional momentum flux divergence.
        jacedgef = 0.5 * (jac[i, j, k] + jac[i, j + 1, k])
        fluxdiff /= jacedgef

        # Explicit integration of Coriolis force in TFC.
        vold[i, j, k] = v[i, j, k]
        if k == k1 && ko + nzz != sizezz
            vold[i, j, k + 1] = v[i, j, k + 1]
        end
        uc = 0.5 * (uold[i, j, k] + uold[i - 1, j, k])
        uf = 0.5 * (uold[i, j + 1, k] + uold[i - 1, j + 1, k])

        volforce =
            -0.5 * (
                fc[j] * (rhoold[i, j, k] + rhostrattfc[i, j, k]) * uc +
                fc[j + 1] *
                (rhoold[i, j + 1, k] + rhostrattfc[i, j + 1, k]) *
                uf
            )

        force = -fluxdiff + volforce

        # Interpolate density.
        rhom_1 = 0.5 * (rhoold[i, j, k] + rhoold[i, j + 1, k])
        rhom = 0.5 * (rho[i, j, k] + rho[i, j + 1, k])
        rhostratedgef = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i, j + 1, k])
        rhom_1 += rhostratedgef
        rhom += rhostratedgef

        vm_1 = v[i, j, k]
        momm_1 = rhom_1 * vm_1

        dv[i, j, k] = dt * force + alphark[m] * dv[i, j, k]

        momm = momm_1 + betark[m] * dv[i, j, k]

        # Update wind.
        vast = momm / rhom
        v[i, j, k] = vast
    end

    # Return.
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::V,
    side::RHS,
    integration::Explicit,
)
    (; nbz) = state.namelists.domain
    (; zboundaries) = state.namelists.setting
    (; kappainv, mainv2) = state.constants
    (; sizezz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dy, dz, met) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; rho, v, pip) = state.variables.predictands

    for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
        rhov = 0.5 * (rho[i, j, k] + rho[i, j + 1, k])
        rhostratedgef = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i, j + 1, k])
        rhov += rhostratedgef

        pif = pip[i, j + 1, k]
        pib = pip[i, j, k]

        # Compute values at cell edges.
        pedgef = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j + 1, k])
        met23edgef = 0.5 * (met[i, j, k, 2, 3] + met[i, j + 1, k, 2, 3])

        # Compute pressure gradient component.
        if ko + k == k0 && zboundaries == SolidWallBoundaries()
            piuuedgef = 0.5 * (pip[i, j, k + 2] + pip[i, j + 1, k + 2])
            piuedgef = 0.5 * (pip[i, j, k + 1] + pip[i, j + 1, k + 1])
            piedgef = 0.5 * (pip[i, j, k] + pip[i, j + 1, k])
            pigrad =
                kappainv * mainv2 * pedgef / rhov * (
                    (pif - pib) / dy +
                    met23edgef *
                    (-piuuedgef + 4.0 * piuedgef - 3.0 * piedgef) *
                    0.5 / dz
                )
        elseif ko + k == sizezz - nbz && zboundaries == SolidWallBoundaries()
            piddedgef = 0.5 * (pip[i, j, k - 2] + pip[i, j + 1, k - 2])
            pidedgef = 0.5 * (pip[i, j, k - 1] + pip[i, j + 1, k - 1])
            piedgef = 0.5 * (pip[i, j, k] + pip[i, j + 1, k])
            pigrad =
                kappainv * mainv2 * pedgef / rhov * (
                    (pif - pib) / dy +
                    met23edgef *
                    (piddedgef - 4.0 * pidedgef + 3.0 * piedgef) *
                    0.5 / dz
                )
        else
            piuedgef = 0.5 * (pip[i, j, k + 1] + pip[i, j + 1, k + 1])
            pidedgef = 0.5 * (pip[i, j, k - 1] + pip[i, j + 1, k - 1])
            pigrad =
                kappainv * mainv2 * pedgef / rhov * (
                    (pif - pib) / dy +
                    met23edgef * (piuedgef - pidedgef) * 0.5 / dz
                )
        end

        volfcy = compute_volume_force(state, (i, j, k), V())

        vhory = v[i, j, k]

        # Update wind.
        jpf = compute_compressible_wind_factor(state, (i, j, k), V())
        vast = vhory + dt * (-pigrad + volfcy / rhov) * jpf
        v[i, j, k] = vast
    end

    # Return.
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::V,
    side::RHS,
    integration::Implicit,
    facray::AbstractFloat,
)
    (; nbz) = state.namelists.domain
    (; zboundaries) = state.namelists.setting
    (; spongelayer, sponge_uv) = state.namelists.sponge
    (; kappainv, mainv2) = state.constants
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dy, dz, met) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; kr_sp_tfc) = state.sponge
    (; rho, v, pip) = state.variables.predictands

    kz0 = k0
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for k in kz0:kz1, j in (j0 - 1):j1, i in i0:i1
        rhov = 0.5 * (rho[i, j, k] + rho[i, j + 1, k])
        rhostratedgef = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i, j + 1, k])
        rhov += rhostratedgef

        pif = pip[i, j + 1, k]
        pib = pip[i, j, k]

        # Compute values at cell edges.
        pedgef = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j + 1, k])
        met23edgef = 0.5 * (met[i, j, k, 2, 3] + met[i, j + 1, k, 2, 3])

        # Compute pressure gradient component.
        if ko + k == k0 && zboundaries == SolidWallBoundaries()
            piuuedgef = 0.5 * (pip[i, j, k + 2] + pip[i, j + 1, k + 2])
            piuedgef = 0.5 * (pip[i, j, k + 1] + pip[i, j + 1, k + 1])
            piedgef = 0.5 * (pip[i, j, k] + pip[i, j + 1, k])
            pigrady =
                kappainv * mainv2 * pedgef / rhov * (
                    (pif - pib) / dy +
                    met23edgef *
                    (-piuuedgef + 4.0 * piuedgef - 3.0 * piedgef) *
                    0.5 / dz
                )
        elseif ko + k == sizezz - nbz && zboundaries == SolidWallBoundaries()
            piddedgef = 0.5 * (pip[i, j, k - 2] + pip[i, j + 1, k - 2])
            pidedgef = 0.5 * (pip[i, j, k - 1] + pip[i, j + 1, k - 1])
            piedgef = 0.5 * (pip[i, j, k] + pip[i, j + 1, k])
            pigrady =
                kappainv * mainv2 * pedgef / rhov * (
                    (pif - pib) / dy +
                    met23edgef *
                    (piddedgef - 4.0 * pidedgef + 3.0 * piedgef) *
                    0.5 / dz
                )
        else
            piuedgef = 0.5 * (pip[i, j, k + 1] + pip[i, j + 1, k + 1])
            pidedgef = 0.5 * (pip[i, j, k - 1] + pip[i, j + 1, k - 1])
            pigrady =
                kappainv * mainv2 * pedgef / rhov * (
                    (pif - pib) / dy +
                    met23edgef * (piuedgef - pidedgef) * 0.5 / dz
                )
        end

        volfcy = compute_volume_force(state, (i, j, k), V())

        vhory = v[i, j, k]

        facv = 1.0

        if spongelayer && sponge_uv
            facv +=
                dt *
                0.5 *
                (kr_sp_tfc[i, j, k] + kr_sp_tfc[i, j + 1, k]) *
                facray
        end

        # Update wind.
        jpf = compute_compressible_wind_factor(state, (i, j, k), V())
        vast = 1.0 / facv * (vhory + dt * (-pigrady + volfcy / rhov) * jpf)
        v[i, j, k] = vast
    end

    # Return.
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::W,
    side::LHS,
)
    (; zboundaries) = state.namelists.setting
    (; alphark, betark) = state.time
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; grid) = state
    (; dx, dy, dz, jac, met) = grid
    (; rhostrattfc, fc) = state.atmosphere
    (; dw) = state.variables.tendencies
    (; phiu, phiv, phiw) = state.variables.fluxes
    (; rhoold, uold, vold) = state.variables.backups
    (; rho, w) = state.variables.predictands

    # Initialize fields for transformation of momentum flux divergence.
    (fluxdiffu, fluxdiffv) = (zeros(2, 2) for i in 1:2)

    if m == 1
        dw .= 0.0
    end

    if zboundaries != SolidWallBoundaries()
        error("Error in update!: Unknown case zBoundary!")
    end

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 - 1 : k1

    for k in kz0:kz1, j in j0:j1, i in i0:i1
        # Compute vertical momentum flux divergence.
        fr = phiw[i, j, k, 1]
        fl = phiw[i - 1, j, k, 1]
        gf = phiw[i, j, k, 2]
        gb = phiw[i, j - 1, k, 2]
        hu = phiw[i, j, k, 3]
        hd = phiw[i, j, k - 1, 3]
        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz

        # Adjust Cartesian vertical momentum flux divergence.
        jacedgeu =
            2.0 * jac[i, j, k] * jac[i, j, k + 1] /
            (jac[i, j, k] + jac[i, j, k + 1])
        fluxdiff /= jacedgeu

        # Compute zonal momentum flux divergences.
        for ll in 0:1, mm in 0:1
            fr = phiu[i - ll, j, k + mm, 1]
            fl = phiu[i - 1 - ll, j, k + mm, 1]
            gf = phiu[i - ll, j, k + mm, 2]
            gb = phiu[i - ll, j - 1, k + mm, 2]
            hu = phiu[i - ll, j, k + mm, 3]
            hd = phiu[i - ll, j, k - 1 + mm, 3]
            fluxdiffu[ll + 1, mm + 1] =
                (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
            jacedger =
                0.5 * (jac[i - ll, j, k + mm] + jac[i + 1 - ll, j, k + mm])
            fluxdiffu[ll + 1, mm + 1] /= jacedger
        end

        # Compute meridional momentum flux divergences.
        for ll in 0:1, mm in 0:1
            fr = phiv[i, j - ll, k + mm, 1]
            fl = phiv[i - 1, j - ll, k + mm, 1]
            gf = phiv[i, j - ll, k + mm, 2]
            gb = phiv[i, j - 1 - ll, k + mm, 2]
            hu = phiv[i, j - ll, k + mm, 3]
            hd = phiv[i, j - ll, k - 1 + mm, 3]
            fluxdiffv[ll + 1, mm + 1] =
                (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
            jacedgef =
                0.5 * (jac[i, j - ll, k + mm] + jac[i, j + 1 - ll, k + mm])
            fluxdiffv[ll + 1, mm + 1] /= jacedgef
        end

        # Compute transformed vertical momentum flux divergence.
        fluxdiff = transform(
            i,
            j,
            k,
            fluxdiffu[1, 1],
            fluxdiffu[1, 2],
            fluxdiffu[2, 1],
            fluxdiffu[2, 2],
            fluxdiffv[1, 1],
            fluxdiffv[1, 2],
            fluxdiffv[2, 1],
            fluxdiffv[2, 2],
            fluxdiff,
            Transformed(),
            grid,
        )

        # Explicit integration of Coriolis force in TFC.
        vc = 0.5 * (vold[i, j, k] + vold[i, j - 1, k])
        vu = 0.5 * (vold[i, j, k + 1] + vold[i, j - 1, k + 1])
        uc = 0.5 * (uold[i, j, k] + uold[i - 1, j, k])
        uu = 0.5 * (uold[i, j, k + 1] + uold[i - 1, j, k + 1])

        volforce =
            fc[j] * (
                jac[i, j, k + 1] *
                met[i, j, k, 1, 3] *
                (rhoold[i, j, k] + rhostrattfc[i, j, k]) *
                vc +
                jac[i, j, k] *
                met[i, j, k + 1, 1, 3] *
                (rhoold[i, j, k + 1] + rhostrattfc[i, j, k + 1]) *
                vu
            ) / (jac[i, j, k] + jac[i, j, k + 1]) -
            fc[j] * (
                jac[i, j, k + 1] *
                met[i, j, k, 2, 3] *
                (rhoold[i, j, k] + rhostrattfc[i, j, k]) *
                uc +
                jac[i, j, k] *
                met[i, j, k + 1, 2, 3] *
                (rhoold[i, j, k + 1] + rhostrattfc[i, j, k + 1]) *
                uu
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        force = -fluxdiff + volforce

        # Interpolate densities.
        rhom_1 =
            (
                jac[i, j, k + 1] * rhoold[i, j, k] +
                jac[i, j, k] * rhoold[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhom =
            (
                jac[i, j, k + 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhostratedgeu =
            (
                jac[i, j, k + 1] * rhostrattfc[i, j, k] +
                jac[i, j, k] * rhostrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhom_1 += rhostratedgeu
        rhom += rhostratedgeu

        wm_1 = w[i, j, k]
        momm_1 = rhom_1 * wm_1

        dw[i, j, k] = dt * force + alphark[m] * dw[i, j, k]

        momm = momm_1 + betark[m] * dw[i, j, k]

        # Update wind.
        wast = momm / rhom
        w[i, j, k] = wast
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::W,
    side::RHS,
    integration::Explicit,
)
    (; zboundaries) = state.namelists.setting
    (; kappainv, mainv2, g_ndim) = state.constants
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac, met) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; rhopold) = state.variables.backups
    (; rho, w, pip) = state.variables.predictands

    if zboundaries != SolidWallBoundaries()
        error("Error in update!: Unknown zboundaries!")
    end

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 - 1 : k1

    for k in kz0:kz1, j in j0:j1, i in i0:i1
        rho000 = rho[i, j, k]
        rho001 = rho[i, j, k + 1]

        rhow =
            (
                jac[i, j, k + 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        rhostratedgeu =
            (
                jac[i, j, k + 1] * rhostrattfc[i, j, k] +
                jac[i, j, k] * rhostrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rho000 += rhostrattfc[i, j, k]
        rho001 += rhostrattfc[i, j, k + 1]
        rhow += rhostratedgeu

        # Compute values at cell edges.
        pedgeu =
            (
                jac[i, j, k + 1] * pstrattfc[i, j, k] +
                jac[i, j, k] * pstrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met13edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 1, 3] +
                jac[i, j, k] * met[i, j, k + 1, 1, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met23edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 2, 3] +
                jac[i, j, k] * met[i, j, k + 1, 2, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met33edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 3, 3] +
                jac[i, j, k] * met[i, j, k + 1, 3, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        piredgeu =
            (
                jac[i + 1, j, k + 1] * pip[i + 1, j, k] +
                jac[i + 1, j, k] * pip[i + 1, j, k + 1]
            ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
        piledgeu =
            (
                jac[i - 1, j, k + 1] * pip[i - 1, j, k] +
                jac[i - 1, j, k] * pip[i - 1, j, k + 1]
            ) / (jac[i - 1, j, k] + jac[i - 1, j, k + 1])
        pifedgeu =
            (
                jac[i, j + 1, k + 1] * pip[i, j + 1, k] +
                jac[i, j + 1, k] * pip[i, j + 1, k + 1]
            ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
        pibedgeu =
            (
                jac[i, j - 1, k + 1] * pip[i, j - 1, k] +
                jac[i, j - 1, k] * pip[i, j - 1, k + 1]
            ) / (jac[i, j - 1, k] + jac[i, j - 1, k + 1])

        # Compute pressure gradient component.
        pigrad =
            kappainv * mainv2 * pedgeu / rhow * (
                met13edgeu * (piredgeu - piledgeu) * 0.5 / dx +
                met23edgeu * (pifedgeu - pibedgeu) * 0.5 / dy +
                met33edgeu * (pip[i, j, k + 1] - pip[i, j, k]) / dz
            )

        volfcz = compute_volume_force(state, (i, j, k), W())

        wvert = w[i, j, k]

        buoy =
            -g_ndim * (
                jac[i, j, k + 1] * rhopold[i, j, k] / rho000 / jac[i, j, k] +
                jac[i, j, k] * rhopold[i, j, k + 1] / rho001 / jac[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        # Update wind.
        jpu = compute_compressible_wind_factor(state, (i, j, k), W())
        wast = wvert + dt * (buoy - pigrad + volfcz / rhow) * jpu
        w[i, j, k] = wast
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::W,
    side::RHS,
    integration::Implicit,
    facray::AbstractFloat,
)
    (; spongelayer) = state.namelists.sponge
    (; zboundaries) = state.namelists.setting
    (; kappainv, mainv2, g_ndim) = state.constants
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac, met) = state.grid
    (; rhostrattfc, pstrattfc, bvsstrattfc) = state.atmosphere
    (; kr_sp_w_tfc) = state.sponge
    (; rho, rhop, u, v, w, pip) = state.variables.predictands

    if zboundaries != SolidWallBoundaries()
        error("Error in update!: Unknown zboundaries!")
    end

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 - 1 : k1

    for k in kz0:kz1, j in j0:j1, i in i0:i1
        rho000 = rho[i, j, k]
        rho001 = rho[i, j, k + 1]

        rhow =
            (
                jac[i, j, k + 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        rhostratedgeu =
            (
                jac[i, j, k + 1] * rhostrattfc[i, j, k] +
                jac[i, j, k] * rhostrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rho000 += rhostrattfc[i, j, k]
        rho001 += rhostrattfc[i, j, k + 1]
        rhow += rhostratedgeu

        # Compute values at cell edges.
        pedgeu =
            (
                jac[i, j, k + 1] * pstrattfc[i, j, k] +
                jac[i, j, k] * pstrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met13edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 1, 3] +
                jac[i, j, k] * met[i, j, k + 1, 1, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met23edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 2, 3] +
                jac[i, j, k] * met[i, j, k + 1, 2, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met33edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 3, 3] +
                jac[i, j, k] * met[i, j, k + 1, 3, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        piredgeu =
            (
                jac[i + 1, j, k + 1] * pip[i + 1, j, k] +
                jac[i + 1, j, k] * pip[i + 1, j, k + 1]
            ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
        piledgeu =
            (
                jac[i - 1, j, k + 1] * pip[i - 1, j, k] +
                jac[i - 1, j, k] * pip[i - 1, j, k + 1]
            ) / (jac[i - 1, j, k] + jac[i - 1, j, k + 1])
        pifedgeu =
            (
                jac[i, j + 1, k + 1] * pip[i, j + 1, k] +
                jac[i, j + 1, k] * pip[i, j + 1, k + 1]
            ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
        pibedgeu =
            (
                jac[i, j - 1, k + 1] * pip[i, j - 1, k] +
                jac[i, j - 1, k] * pip[i, j - 1, k + 1]
            ) / (jac[i, j - 1, k] + jac[i, j - 1, k + 1])

        # Compute pressure gradient component.
        pigrad =
            kappainv * mainv2 * pedgeu / rhow * (
                met13edgeu * (piredgeu - piledgeu) * 0.5 / dx +
                met23edgeu * (pifedgeu - pibedgeu) * 0.5 / dy +
                met33edgeu * (pip[i, j, k + 1] - pip[i, j, k]) / dz
            )

        volfcz = compute_volume_force(state, (i, j, k), W())

        wvert = w[i, j, k]

        bvsstw =
            (
                jac[i, j, k + 1] * bvsstrattfc[i, j, k] +
                jac[i, j, k] * bvsstrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        facw = 1.0

        if spongelayer
            facw +=
                dt * (
                    jac[i, j, k + 1] * kr_sp_w_tfc[i, j, k] +
                    jac[i, j, k] * kr_sp_w_tfc[i, j, k + 1]
                ) / (jac[i, j, k] + jac[i, j, k + 1]) * facray
        end

        # Buoyancy is predicted after momentum in implicit steps.
        buoy =
            -g_ndim * (
                jac[i, j, k + 1] * rhop[i, j, k] / rho000 / jac[i, j, k] +
                jac[i, j, k] * rhop[i, j, k + 1] / rho001 / jac[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        jpr = compute_compressible_wind_factor(state, (i, j, k), U())
        jpl = compute_compressible_wind_factor(state, (i - 1, j, k), U())
        jpf = compute_compressible_wind_factor(state, (i, j, k), V())
        jpb = compute_compressible_wind_factor(state, (i, j - 1, k), V())
        jpur = compute_compressible_wind_factor(state, (i, j, k + 1), U())
        jpul = compute_compressible_wind_factor(state, (i - 1, j, k + 1), U())
        jpuf = compute_compressible_wind_factor(state, (i, j, k + 1), V())
        jpub = compute_compressible_wind_factor(state, (i, j - 1, k + 1), V())

        uc = 0.5 * (u[i, j, k] / jpr + u[i - 1, j, k] / jpl)
        uu = 0.5 * (u[i, j, k + 1] / jpur + u[i - 1, j, k + 1] / jpul)
        vc = 0.5 * (v[i, j, k] / jpf + v[i, j - 1, k] / jpb)
        vu = 0.5 * (v[i, j, k + 1] / jpuf + v[i, j - 1, k + 1] / jpub)

        # Update wind.
        jpu = compute_compressible_wind_factor(state, (i, j, k), W())
        fw = compute_compressible_buoyancy_factor(state, (i, j, k), W())
        wast =
            1.0 / (facw + fw * bvsstw * dt^2.0) * (
                wvert - dt * pigrad * jpu +
                dt * buoy * jpu +
                dt * volfcz / rhow * jpu +
                jpu *
                fw *
                bvsstw *
                dt^2.0 *
                (
                    jac[i, j, k + 1] *
                    (met[i, j, k, 1, 3] * uc + met[i, j, k, 2, 3] * vc) +
                    jac[i, j, k] *
                    (met[i, j, k + 1, 1, 3] * uu + met[i, j, k + 1, 2, 3] * vu)
                ) / (jac[i, j, k] + jac[i, j, k + 1])
            )
        w[i, j, k] = wast
    end

    return
end

function update!(state::State, dt::AbstractFloat, variable::PiP)
    (; model) = state.namelists.setting
    update!(state, dt, variable, model)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::PiP,
    model::AbstractModel,
)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::PiP,
    model::Compressible,
)
    (; gamma, rsp, pref) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; uold, vold, wold) = state.variables.backups
    (; pip, p) = state.variables.predictands

    for k in k0:k1, j in j0:j1, i in i0:i1
        fl = uold[i - 1, j, k]
        fr = uold[i, j, k]
        gb = vold[i, j - 1, k]
        gf = vold[i, j, k]
        hd = wold[i, j, k - 1]
        hu = wold[i, j, k]

        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
        fluxdiff /= jac[i, j, k]

        heating = compute_volume_force(state, (i, j, k), P())

        dpdpi =
            1 / (gamma - 1) * (rsp / pref)^(1 - gamma) * p[i, j, k]^(2 - gamma)

        pip[i, j, k] -= dt * (fluxdiff + heating) / dpdpi
    end

    return
end

function update!(state::State, dt::AbstractFloat, m::Integer, variable::P)
    (; model) = state.namelists.setting
    update!(state, dt, m, variable, model)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::P,
    model::AbstractModel,
)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::P,
    model::Compressible,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; alphark, betark) = state.time
    (; dp) = state.variables.tendencies
    (; phip) = state.variables.fluxes
    (; p) = state.variables.predictands

    if m == 1
        dp .= 0.0
    end

    for k in k0:k1, j in j0:j1, i in i0:i1
        fl = phip[i - 1, j, k, 1]
        fr = phip[i, j, k, 1]
        gb = phip[i, j - 1, k, 2]
        gf = phip[i, j, k, 2]
        hd = phip[i, j, k - 1, 3]
        hu = phip[i, j, k, 3]

        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
        fluxdiff /= jac[i, j, k]

        heating = compute_volume_force(state, (i, j, k), P())

        f = -fluxdiff - heating

        dp[i, j, k] = dt * f + alphark[m] * dp[i, j, k]
        p[i, j, k] += betark[m] * dp[i, j, k]
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    tracersetup::NoTracer,
)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    tracersetup::AbstractTracer,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; alphark, betark) = state.time
    (; tracertendencies, tracerpredictands, tracerfluxes) = state.tracer

    for (fd, field) in enumerate(fieldnames(TracerPredictands))
        if m == 1
            getfield(tracertendencies, fd) .= 0.0
        end

        for k in k0:k1, j in j0:j1, i in i0:i1
            fl = getfield(tracerfluxes, fd)[i - 1, j, k, 1]
            fr = getfield(tracerfluxes, fd)[i, j, k, 1]
            gb = getfield(tracerfluxes, fd)[i, j - 1, k, 2]
            gf = getfield(tracerfluxes, fd)[i, j, k, 2]
            hd = getfield(tracerfluxes, fd)[i, j, k - 1, 3]
            hu = getfield(tracerfluxes, fd)[i, j, k, 3]

            fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
            fluxdiff /= jac[i, j, k]

            f = -fluxdiff

            getfield(tracertendencies, fd)[i, j, k] =
                dt * f + alphark[m] * getfield(tracertendencies, fd)[i, j, k]
            getfield(tracerpredictands, fd)[i, j, k] +=
                betark[m] * getfield(tracertendencies, fd)[i, j, k]
        end
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    icesetup::AbstractIce,
)
    return
end

function update!(state::State, dt::AbstractFloat, m::Integer, icesetup::IceOn)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; alphark, betark) = state.time
    (; icetendencies, icepredictands, icefluxes) = state.ice

    for (fd, field) in enumerate(fieldnames(IcePredictands))
        if m == 1
            getfield(icetendencies, fd) .= 0.0
        end

        for k in k0:k1, j in j0:j1, i in i0:i1
            fl = getfield(icefluxes, fd)[i - 1, j, k, 1]
            fr = getfield(icefluxes, fd)[i, j, k, 1]
            gb = getfield(icefluxes, fd)[i, j - 1, k, 2]
            gf = getfield(icefluxes, fd)[i, j, k, 2]
            hd = getfield(icefluxes, fd)[i, j, k - 1, 3]
            hu = getfield(icefluxes, fd)[i, j, k, 3]

            fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
            fluxdiff /= jac[i, j, k]

            f = -fluxdiff

            getfield(icetendencies, fd)[i, j, k] =
                dt * f + alphark[m] * getfield(icetendencies, fd)[i, j, k]
            getfield(icepredictands, fd)[i, j, k] +=
                betark[m] * getfield(icetendencies, fd)[i, j, k]
        end
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    turbulencesetup::NoTurbulence,
)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    turbulencesetup::AbstractTurbulence,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; alphark, betark) = state.time
    (; turbulencetendencies, turbulencepredictands, turbulencefluxes) =
        state.turbulence

    for (fd, field) in enumerate(fieldnames(TurbulencePredictands))
        if m == 1
            getfield(turbulencetendencies, fd) .= 0.0
        end

        for k in k0:k1, j in j0:j1, i in i0:i1
            fl = getfield(turbulencefluxes, fd)[i - 1, j, k, 1]
            fr = getfield(turbulencefluxes, fd)[i, j, k, 1]
            gb = getfield(turbulencefluxes, fd)[i, j - 1, k, 2]
            gf = getfield(turbulencefluxes, fd)[i, j, k, 2]
            hd = getfield(turbulencefluxes, fd)[i, j, k - 1, 3]
            hu = getfield(turbulencefluxes, fd)[i, j, k, 3]

            fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
            fluxdiff /= jac[i, j, k]

            f = -fluxdiff

            getfield(turbulencetendencies, fd)[i, j, k] =
                dt * f +
                alphark[m] * getfield(turbulencetendencies, fd)[i, j, k]
            getfield(turbulencepredictands, fd)[i, j, k] +=
                betark[m] * getfield(turbulencetendencies, fd)[i, j, k]
        end
    end

    return
end
