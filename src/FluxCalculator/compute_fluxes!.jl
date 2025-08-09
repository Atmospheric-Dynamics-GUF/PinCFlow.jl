"""
```julia
compute_fluxes!(state::State, predictands::Predictands)
```

Compute fluxes by dispatching to specialized methods for each prognostic variable.

```julia
compute_fluxes!(state::State, predictands::Predictands, variable::Rho)
```

Compute the density fluxes in all three directions.

The fluxes are computed from the MUSCL reconstruction of ``\\rho / P`` and the linear interpolation of ``P \\widehat{\\boldsymbol{u}}`` at the respective cell interfaces. They are written into `state.variables.fluxes.phirho`.

```julia
compute_fluxes!(state::State, predictands::Predictands, variable::RhoP)
```

Compute the density-fluctuation fluxes in all three directions.

The computation is analogous to that of the density fluxes. The fluxes are written into `state.variables.fluxes.phirhop`.

```julia
compute_fluxes!(
    state::State,
    predictands::Predictands,
    model::AbstractModel,
    variable::P,
)
```

Return in non-compressible modes.

```julia
compute_fluxes!(
    state::State,
    predictands::Predictands,
    model::Compressible,
    variable::P,
)
```

Compute the mass-weighted potential-temperature fluxes in all three directions.

The fluxes are computed by interpolating ``P \\widehat{\\boldsymbol{u}}`` at the respective cell interfaces. They are written into `state.variables.fluxes.phip`.

```julia
compute_fluxes!(state::State, old_predictands::Predictands, variable::U)
```

Compute the sums of advective and viscous zonal-momentum fluxes in all three directions.

The advective fluxes are computed from the MUSCL reconstruction of ``\\rho u / P`` and the linear interpolation of ``P \\widehat{\\boldsymbol{u}}`` at the respective cell interfaces. The viscous fluxes are computed from linear interpolations of the corresponding elements of the viscous stress tensor. The total fluxes are written into `state.variables.fluxes.phiu`.

```julia
compute_fluxes!(state::State, old_predictands::Predictands, variable::V)
```

Compute the sums of advective and viscous meridional-momentum fluxes in all three directions.

The computation is analogous to that of the zonal-momentum fluxes. The total fluxes are written into `state.variables.fluxes.phiv`.

```julia
compute_fluxes!(state::State, old_predictands::Predictands, variable::W)
```

Compute the sums of advective and viscous vertical-momentum fluxes in all three directions.

The computation is analogous to those of the zonal-momentum and meridional-momentum fluxes. The total fluxes are written into `state.variables.fluxes.phiw`.

```julia
compute_fluxes!(state::State, predictands::Predictands, tracersetup::NoTracer)
```

Return for configurations without tracer transport.

```julia
compute_fluxes!(
    state::State,
    predictands::Predictands,
    tracersetup::AbstractTracer,
)
```

Compute the tracer fluxes in all three directions.

The computation is analogous to that of the density fluxes.

```julia
compute_fluxes!(state::State, predictands::Predictands, icesetup::NoIce)
```

Return for configurations without ice physics.

```julia
compute_fluxes!(state::State, predictands::Predictands, icesetup::AbstractIce)
```

Compute the fluxes of ice variables in all three directions.

The computation is analogous to that of the density fluxes.

```julia
compute_fluxes!(
    state::State,
    predictands::Predictands,
    turbulencesetup::NoTurbulence,
)
```

Return for configurations without turbulence physics.

```julia
compute_fluxes!(
    state::State,
    predictands::Predictands,
    turbulencesetup::AbstractTurbulence,
)
```

Compute the fluxes of turbulence variables in all three directions.

The computation is analogous to that of the density fluxes.

# Arguments

  - `state`: Model state.
  - `predictands`/`old_predictands`: The predictands that are used to compute the transporting velocities.
  - `model`: Dynamic equations.
  - `variable`: Flux variable.
  - `tracersetup`: General tracer-transport configuration.
  - `icesetup`: General ice-physics configuration.
  - `turbulencesetup`: General turbulence-physics configuration.

# See also

  - [`PinCFlow.FluxCalculator.compute_flux`](@ref)
  - [`PinCFlow.Update.compute_stress_tensor`](@ref)
"""
function compute_fluxes! end

function compute_fluxes!(state::State, predictands::Predictands)
    (; model) = state.namelists.setting

    compute_fluxes!(state, predictands, Rho())
    compute_fluxes!(state, predictands, RhoP())
    compute_fluxes!(state, predictands, U())
    compute_fluxes!(state, predictands, V())
    compute_fluxes!(state, predictands, W())

    compute_fluxes!(state, predictands, model, P())
    compute_fluxes!(state, predictands, state.namelists.tracer.tracersetup)
    compute_fluxes!(state, predictands, state.namelists.ice.icesetup)
    compute_fluxes!(
        state,
        predictands,
        state.namelists.turbulence.turbulencesetup,
    )
    return
end

function compute_fluxes!(state::State, predictands::Predictands, variable::Rho)

    # Get all necessary fields.
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid
    (; pstrattfc, rhostrattfc) = state.atmosphere
    (; rhotilde) = state.variables.reconstructions
    (; phirho) = state.variables.fluxes

    # Get old wind.
    (u0, v0, w0) = (predictands.u, predictands.v, predictands.w)

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
        rhostratedger = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i + 1, j, k])
        pedger = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i + 1, j, k])
        rhor = rhotilde[i + 1, j, k, 1, 1] + rhostratedger / pedger
        rhol = rhotilde[i, j, k, 1, 2] + rhostratedger / pedger

        pedger =
            0.5 * (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
            )
        usurf = pedger * u0[i, j, k]

        frho = compute_flux(usurf, rhol, rhor)

        phirho[i, j, k, 1] = frho
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
        rhostratedgef = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i, j + 1, k])
        pedgef = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j + 1, k])
        rhof = rhotilde[i, j + 1, k, 2, 1] + rhostratedgef / pedgef
        rhob = rhotilde[i, j, k, 2, 2] + rhostratedgef / pedgef

        pedgef =
            0.5 * (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
            )
        vsurf = pedgef * v0[i, j, k]

        grho = compute_flux(vsurf, rhob, rhof)

        phirho[i, j, k, 2] = grho
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    for k in (k0 - 1):k1, j in j0:j1, i in i0:i1
        rhostratedgeu =
            (
                jac[i, j, k + 1] * rhostrattfc[i, j, k] +
                jac[i, j, k] * rhostrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        pedgeu =
            (
                jac[i, j, k + 1] * pstrattfc[i, j, k] +
                jac[i, j, k] * pstrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhou = rhotilde[i, j, k + 1, 3, 1] + rhostratedgeu / pedgeu
        rhod = rhotilde[i, j, k, 3, 2] + rhostratedgeu / pedgeu

        pedgeu =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1])
        wsurf = pedgeu * w0[i, j, k]

        hrho = compute_flux(wsurf, rhod, rhou)

        phirho[i, j, k, 3] = hrho
    end

    # Return.
    return
end

function compute_fluxes!(state::State, predictands::Predictands, variable::RhoP)

    # Get all necessary fields.
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid
    (; pstrattfc) = state.atmosphere
    (; rhoptilde) = state.variables.reconstructions
    (; phirhop) = state.variables.fluxes

    # Get old wind.
    (u0, v0, w0) = (predictands.u, predictands.v, predictands.w)

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
        rhor = rhoptilde[i + 1, j, k, 1, 1]
        rhol = rhoptilde[i, j, k, 1, 2]

        pedger =
            0.5 * (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
            )
        usurf = pedger * u0[i, j, k]

        frhop = compute_flux(usurf, rhol, rhor)

        phirhop[i, j, k, 1] = frhop
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
        rhof = rhoptilde[i, j + 1, k, 2, 1]
        rhob = rhoptilde[i, j, k, 2, 2]

        pedgef =
            0.5 * (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
            )
        vsurf = pedgef * v0[i, j, k]

        grhop = compute_flux(vsurf, rhob, rhof)

        phirhop[i, j, k, 2] = grhop
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    for k in (k0 - 1):k1, j in j0:j1, i in i0:i1
        rhou = rhoptilde[i, j, k + 1, 3, 1]
        rhod = rhoptilde[i, j, k, 3, 2]

        pedgeu =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1])
        wsurf = pedgeu * w0[i, j, k]

        hrhop = compute_flux(wsurf, rhod, rhou)

        phirhop[i, j, k, 3] = hrhop
    end

    # Return.
    return
end

function compute_fluxes!(
    state::State,
    predictands::Predictands,
    model::AbstractModel,
    variable::P,
)
    return
end

function compute_fluxes!(
    state::State,
    predictands::Predictands,
    model::Compressible,
    variable::P,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid
    (; pstrattfc) = state.atmosphere
    (; phip) = state.variables.fluxes

    (u0, v0, w0) = (predictands.u, predictands.v, predictands.w)

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
        phip[i, j, k, 1] =
            0.5 *
            (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
            ) *
            u0[i, j, k]
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
        phip[i, j, k, 2] =
            0.5 *
            (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
            ) *
            v0[i, j, k]
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    for k in (k0 - 1):k1, j in j0:j1, i in i0:i1
        phip[i, j, k, 3] =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1]) * w0[i, j, k]
    end

    return
end

function compute_fluxes!(
    state::State,
    old_predictands::Predictands,
    variable::U,
)

    # Get all necessary fields.
    (; grid) = state
    (; re) = state.constants
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met) = grid
    (; pstrattfc, rhostrattfc) = state.atmosphere
    (; utilde) = state.variables.reconstructions
    (; phiu) = state.variables.fluxes
    (; predictands) = state.variables

    # Get old wind.
    (u0, v0, w0) = (old_predictands.u, old_predictands.v, old_predictands.w)

    kz0 = k0
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    for k in kz0:kz1, j in j0:j1, i in (i0 - 2):i1
        # The uTilde are the reconstructed specific momenta, divided by P.
        # These are to be multiplied by the linearly interpolated velocities
        # (times P) in order to obtain the desired momentum fluxes.

        ur = utilde[i + 1, j, k, 1, 1]
        ul = utilde[i, j, k, 1, 2]

        pedger =
            0.5 * (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
            )
        predger =
            0.5 * (
                jac[i + 1, j, k] * pstrattfc[i + 1, j, k] +
                jac[i + 2, j, k] * pstrattfc[i + 2, j, k]
            )
        usurf = 0.5 * (pedger * u0[i, j, k] + predger * u0[i + 1, j, k])

        frhou = compute_flux(usurf, ul, ur)

        phiu[i, j, k, 1] = frhou
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    for k in kz0:kz1, j in (j0 - 1):j1, i in (i0 - 1):i1
        # The uTilde are the reconstructed specific momenta, divided by P.
        # These are to be multiplied by the linearly interpolated velocities
        # (times P) in order to obtain the desired momentum fluxes.

        uf = utilde[i, j + 1, k, 2, 1]
        ub = utilde[i, j, k, 2, 2]

        pedgef =
            0.5 * (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
            )
        predgef =
            0.5 * (
                jac[i + 1, j, k] * pstrattfc[i + 1, j, k] +
                jac[i + 1, j + 1, k] * pstrattfc[i + 1, j + 1, k]
            )
        vsurf = 0.5 * (pedgef * v0[i, j, k] + predgef * v0[i + 1, j, k])

        grhou = compute_flux(vsurf, ub, uf)

        phiu[i, j, k, 2] = grhou
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    for k in (kz0 - 1):kz1, j in j0:j1, i in (i0 - 1):i1
        # The uTilde are the reconstructed specific momenta, divided by P.
        # These are to be multiplied by the linearly interpolated velocities
        # (times P) in order to obtain the desired momentum fluxes.

        uu = utilde[i, j, k + 1, 3, 1]
        ud = utilde[i, j, k, 3, 2]

        pedgeu =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1])
        predgeu =
            jac[i + 1, j, k] *
            jac[i + 1, j, k + 1] *
            (pstrattfc[i + 1, j, k] + pstrattfc[i + 1, j, k + 1]) /
            (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
        wsurf = 0.5 * (pedgeu * w0[i, j, k] + predgeu * w0[i + 1, j, k])

        hrhou = compute_flux(wsurf, ud, uu)

        phiu[i, j, k, 3] = hrhou
    end

    #-------------------------------------------------------------------
    #                          Viscous fluxes
    #-------------------------------------------------------------------

    if 1 / re <= eps()
        return
    end

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    for k in kz0:kz1, j in j0:j1, i in (i0 - 2):i1
        coef_v = 1 / re * rhostrattfc[i + 1, j, 1]

        frhou_visc =
            coef_v *
            jac[i + 1, j, k] *
            compute_stress_tensor(i + 1, j, k, 1, 1, predictands, grid)

        phiu[i, j, k, 1] -= frhou_visc
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    for k in kz0:kz1, j in (j0 - 1):j1, i in (i0 - 1):i1
        coef_v =
            1 / re *
            0.25 *
            (
                rhostrattfc[i, j, 1] +
                rhostrattfc[i + 1, j, 1] +
                rhostrattfc[i, j + 1, 1] +
                rhostrattfc[i + 1, j + 1, 1]
            )

        grhou_visc =
            coef_v *
            0.25 *
            (
                jac[i, j, k] *
                compute_stress_tensor(i, j, k, 1, 2, predictands, grid) +
                jac[i + 1, j, k] *
                compute_stress_tensor(i + 1, j, k, 1, 2, predictands, grid) +
                jac[i, j + 1, k] *
                compute_stress_tensor(i, j + 1, k, 1, 2, predictands, grid) +
                jac[i + 1, j + 1, k] *
                compute_stress_tensor(i + 1, j + 1, k, 1, 2, predictands, grid)
            )

        phiu[i, j, k, 2] -= grhou_visc
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    for k in (kz0 - 1):kz1, j in j0:j1, i in (i0 - 1):i1
        coef_v =
            1 / re * 0.5 * (rhostrattfc[i, j, 1] + rhostrattfc[i + 1, j, 1])

        stresstens13 =
            met[i, j, k, 1, 3] *
            compute_stress_tensor(i, j, k, 1, 1, predictands, grid) +
            met[i, j, k, 2, 3] *
            compute_stress_tensor(i, j, k, 1, 2, predictands, grid) +
            compute_stress_tensor(i, j, k, 1, 3, predictands, grid) /
            jac[i, j, k]
        stresstens13r =
            met[i + 1, j, k, 1, 3] *
            compute_stress_tensor(i + 1, j, k, 1, 1, predictands, grid) +
            met[i + 1, j, k, 2, 3] *
            compute_stress_tensor(i + 1, j, k, 1, 2, predictands, grid) +
            compute_stress_tensor(i + 1, j, k, 1, 3, predictands, grid) /
            jac[i + 1, j, k]
        stresstens13u =
            met[i, j, k + 1, 1, 3] *
            compute_stress_tensor(i, j, k + 1, 1, 1, predictands, grid) +
            met[i, j, k + 1, 2, 3] *
            compute_stress_tensor(i, j, k + 1, 1, 2, predictands, grid) +
            compute_stress_tensor(i, j, k + 1, 1, 3, predictands, grid) /
            jac[i, j, k + 1]
        stresstens13ru =
            met[i + 1, j, k + 1, 1, 3] *
            compute_stress_tensor(i + 1, j, k + 1, 1, 1, predictands, grid) +
            met[i + 1, j, k + 1, 2, 3] *
            compute_stress_tensor(i + 1, j, k + 1, 1, 2, predictands, grid) +
            compute_stress_tensor(i + 1, j, k + 1, 1, 3, predictands, grid) /
            jac[i + 1, j, k + 1]
        hrhou_visc =
            coef_v *
            0.5 *
            (
                jac[i, j, k] *
                jac[i, j, k + 1] *
                (stresstens13 + stresstens13u) /
                (jac[i, j, k] + jac[i, j, k + 1]) +
                jac[i + 1, j, k] *
                jac[i + 1, j, k + 1] *
                (stresstens13r + stresstens13ru) /
                (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
            )

        phiu[i, j, k, 3] -= hrhou_visc
    end

    # Return.
    return
end

function compute_fluxes!(
    state::State,
    old_predictands::Predictands,
    variable::V,
)

    # Get all necessary fields.
    (; grid) = state
    (; re) = state.constants
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met) = grid
    (; pstrattfc, rhostrattfc) = state.atmosphere
    (; vtilde) = state.variables.reconstructions
    (; phiv) = state.variables.fluxes
    (; predictands) = state.variables

    # Get old wind.
    (u0, v0, w0) = (old_predictands.u, old_predictands.v, old_predictands.w)

    kz0 = k0
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    for k in kz0:kz1, j in (j0 - 1):j1, i in (i0 - 1):i1
        # The vTilde are the reconstructed specific momenta, divided by P.
        # These are to be multiplied by the linearly interpolated velocities
        # (times P) in order to obtain the desired momentum fluxes.

        vr = vtilde[i + 1, j, k, 1, 1]
        vl = vtilde[i, j, k, 1, 2]

        pedger =
            0.5 * (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
            )
        pfedger =
            0.5 * (
                jac[i, j + 1, k] * pstrattfc[i, j + 1, k] +
                jac[i + 1, j + 1, k] * pstrattfc[i + 1, j + 1, k]
            )
        usurf = 0.5 * (pedger * u0[i, j, k] + pfedger * u0[i, j + 1, k])

        frhov = compute_flux(usurf, vl, vr)

        phiv[i, j, k, 1] = frhov
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    for k in kz0:kz1, j in (j0 - 2):j1, i in i0:i1
        # The vTilde are the reconstructed specific momenta, divided by P.
        # These are to be multiplied by the linearly interpolated velocities
        # (times P) in order to obtain the desired momentum fluxes.

        vf = vtilde[i, j + 1, k, 2, 1]
        vb = vtilde[i, j, k, 2, 2]

        pedgef =
            0.5 * (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
            )
        pfedgef =
            0.5 * (
                jac[i, j + 1, k] * pstrattfc[i, j + 1, k] +
                jac[i, j + 2, k] * pstrattfc[i, j + 2, k]
            )
        vsurf = 0.5 * (pedgef * v0[i, j, k] + pfedgef * v0[i, j + 1, k])

        grhov = compute_flux(vsurf, vb, vf)

        phiv[i, j, k, 2] = grhov
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    for k in (kz0 - 1):kz1, j in (j0 - 1):j1, i in i0:i1
        # The vTilde are the reconstructed specific momenta, divided by P.
        # These are to be multiplied by the linearly interpolated velocities
        # (times P) in order to obtain the desired momentum fluxes.

        vu = vtilde[i, j, k + 1, 3, 1]
        vd = vtilde[i, j, k, 3, 2]

        pedgeu =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1])
        pfedgeu =
            jac[i, j + 1, k] *
            jac[i, j + 1, k + 1] *
            (pstrattfc[i, j + 1, k] + pstrattfc[i, j + 1, k + 1]) /
            (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
        wsurf = 0.5 * (pedgeu * w0[i, j, k] + pfedgeu * w0[i, j + 1, k])

        hrhov = compute_flux(wsurf, vd, vu)

        phiv[i, j, k, 3] = hrhov
    end

    #-------------------------------------------------------------------
    #                          Viscous fluxes
    #-------------------------------------------------------------------

    if 1 / re <= eps()
        return
    end

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    for k in kz0:kz1, j in (j0 - 1):j1, i in (i0 - 1):i1
        coef_v =
            1 / re *
            0.25 *
            (
                rhostrattfc[i, j, 1] +
                rhostrattfc[i + 1, j, 1] +
                rhostrattfc[i, j + 1, 1] +
                rhostrattfc[i + 1, j + 1, 1]
            )

        frhov_visc =
            coef_v *
            0.25 *
            (
                jac[i, j, k] *
                compute_stress_tensor(i, j, k, 2, 1, predictands, grid) +
                jac[i + 1, j, k] *
                compute_stress_tensor(i + 1, j, k, 2, 1, predictands, grid) +
                jac[i, j + 1, k] *
                compute_stress_tensor(i, j + 1, k, 2, 1, predictands, grid) +
                jac[i + 1, j + 1, k] *
                compute_stress_tensor(i + 1, j + 1, k, 2, 1, predictands, grid)
            )

        phiv[i, j, k, 1] -= frhov_visc
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    for k in kz0:kz1, j in (j0 - 2):j1, i in i0:i1
        coef_v = 1 / re * rhostrattfc[i, j + 1, 1]

        grhov_visc =
            coef_v *
            jac[i, j + 1, k] *
            compute_stress_tensor(i, j + 1, k, 2, 2, predictands, grid)

        phiv[i, j, k, 2] -= grhov_visc
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    for k in (kz0 - 1):kz1, j in (j0 - 1):j1, i in i0:i1
        coef_v =
            1 / re * 0.5 * (rhostrattfc[i, j, 1] + rhostrattfc[i, j + 1, 1])

        stresstens23 =
            met[i, j, k, 1, 3] *
            compute_stress_tensor(i, j, k, 2, 1, predictands, grid) +
            met[i, j, k, 2, 3] *
            compute_stress_tensor(i, j, k, 2, 2, predictands, grid) +
            compute_stress_tensor(i, j, k, 2, 3, predictands, grid) /
            jac[i, j, k]
        stresstens23f =
            met[i, j + 1, k, 1, 3] *
            compute_stress_tensor(i, j + 1, k, 2, 1, predictands, grid) +
            met[i, j + 1, k, 2, 3] *
            compute_stress_tensor(i, j + 1, k, 2, 2, predictands, grid) +
            compute_stress_tensor(i, j + 1, k, 2, 3, predictands, grid) /
            jac[i, j + 1, k]
        stresstens23u =
            met[i, j, k + 1, 1, 3] *
            compute_stress_tensor(i, j, k + 1, 2, 1, predictands, grid) +
            met[i, j, k + 1, 2, 3] *
            compute_stress_tensor(i, j, k + 1, 2, 2, predictands, grid) +
            compute_stress_tensor(i, j, k + 1, 2, 3, predictands, grid) /
            jac[i, j, k + 1]
        stresstens23fu =
            met[i, j + 1, k + 1, 1, 3] *
            compute_stress_tensor(i, j + 1, k + 1, 2, 1, predictands, grid) +
            met[i, j + 1, k + 1, 2, 3] *
            compute_stress_tensor(i, j + 1, k + 1, 2, 2, predictands, grid) +
            compute_stress_tensor(i, j + 1, k + 1, 2, 3, predictands, grid) /
            jac[i, j + 1, k + 1]
        hrhov_visc =
            coef_v *
            0.5 *
            (
                jac[i, j, k] *
                jac[i, j, k + 1] *
                (stresstens23 + stresstens23u) /
                (jac[i, j, k] + jac[i, j, k + 1]) +
                jac[i, j + 1, k] *
                jac[i, j + 1, k + 1] *
                (stresstens23f + stresstens23fu) /
                (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
            )

        phiv[i, j, k, 3] -= hrhov_visc
    end

    # Return.
    return
end

function compute_fluxes!(
    state::State,
    old_predictands::Predictands,
    variable::W,
)

    # Get all necessary fields.
    (; grid) = state
    (; re) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met) = grid
    (; pstrattfc, rhostrattfc) = state.atmosphere
    (; wtilde) = state.variables.reconstructions
    (; phiw) = state.variables.fluxes
    (; predictands) = state.variables

    # Get old wind.
    (u0, v0, w0) = (old_predictands.u, old_predictands.v, old_predictands.w)

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    for k in (k0 - 1):k1, j in j0:j1, i in (i0 - 1):i1
        # The wTilde are the reconstructed specific momenta, divided by P.
        # These are to be multiplied by the linearly interpolated velocities
        # (times P) in order to obtain the desired momentum fluxes.

        wr = wtilde[i + 1, j, k, 1, 1]
        wl = wtilde[i, j, k, 1, 2]

        pedger =
            0.5 * (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
            )
        puedger =
            0.5 * (
                jac[i, j, k + 1] * pstrattfc[i, j, k + 1] +
                jac[i + 1, j, k + 1] * pstrattfc[i + 1, j, k + 1]
            )
        usurf =
            (
                (jac[i, j, k + 1] + jac[i + 1, j, k + 1]) *
                pedger *
                u0[i, j, k] +
                (jac[i, j, k] + jac[i + 1, j, k]) * puedger * u0[i, j, k + 1]
            ) / (
                jac[i, j, k] +
                jac[i + 1, j, k] +
                jac[i, j, k + 1] +
                jac[i + 1, j, k + 1]
            )

        frhow = compute_flux(usurf, wl, wr)

        phiw[i, j, k, 1] = frhow
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    for k in (k0 - 1):k1, j in (j0 - 1):j1, i in i0:i1
        # The wTilde are the reconstructed specific momenta, divided by P.
        # These are to be multiplied by the linearly interpolated velocities
        # (times P) in order to obtain the desired momentum fluxes.

        wf = wtilde[i, j + 1, k, 2, 1]
        wb = wtilde[i, j, k, 2, 2]

        pedgef =
            0.5 * (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
            )
        puedgef =
            0.5 * (
                jac[i, j, k + 1] * pstrattfc[i, j, k + 1] +
                jac[i, j + 1, k + 1] * pstrattfc[i, j + 1, k + 1]
            )
        vsurf =
            (
                (jac[i, j, k + 1] + jac[i, j + 1, k + 1]) *
                pedgef *
                v0[i, j, k] +
                (jac[i, j, k] + jac[i, j + 1, k]) * puedgef * v0[i, j, k + 1]
            ) / (
                jac[i, j, k] +
                jac[i, j + 1, k] +
                jac[i, j, k + 1] +
                jac[i, j + 1, k + 1]
            )

        grhow = compute_flux(vsurf, wb, wf)

        phiw[i, j, k, 2] = grhow
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    for k in (k0 - 2):k1, j in j0:j1, i in i0:i1
        # The wTilde are the reconstructed specific momenta, divided by P.
        # These are to be multiplied by the linearly interpolated velocities
        # (times P) in order to obtain the desired momentum fluxes.

        wu = wtilde[i, j, k + 1, 3, 1]
        wd = wtilde[i, j, k, 3, 2]

        pedgeu =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1])
        puedgeu =
            jac[i, j, k + 1] *
            jac[i, j, k + 2] *
            (pstrattfc[i, j, k + 1] + pstrattfc[i, j, k + 2]) /
            (jac[i, j, k + 1] + jac[i, j, k + 2])
        wsurf = 0.5 * (pedgeu * w0[i, j, k] + puedgeu * w0[i, j, k + 1])

        hrhow = compute_flux(wsurf, wd, wu)

        phiw[i, j, k, 3] = hrhow
    end

    #-------------------------------------------------------------------
    #                          Viscous fluxes
    #-------------------------------------------------------------------

    if 1 / re <= eps()
        return
    end

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    for k in (k0 - 1):k1, j in j0:j1, i in (i0 - 1):i1
        coef_v =
            1 / re * 0.5 * (rhostrattfc[i, j, 1] + rhostrattfc[i + 1, j, 1])

        frhow_visc =
            coef_v *
            0.5 *
            (
                jac[i, j, k] *
                jac[i, j, k + 1] *
                (
                    compute_stress_tensor(i, j, k, 3, 1, predictands, grid) +
                    compute_stress_tensor(i, j, k + 1, 3, 1, predictands, grid)
                ) / (jac[i, j, k] + jac[i, j, k + 1]) +
                jac[i + 1, j, k] *
                jac[i + 1, j, k + 1] *
                (
                    compute_stress_tensor(
                        i + 1,
                        j,
                        k,
                        3,
                        1,
                        predictands,
                        grid,
                    ) + compute_stress_tensor(
                        i + 1,
                        j,
                        k + 1,
                        3,
                        1,
                        predictands,
                        grid,
                    )
                ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
            )

        phiw[i, j, k, 1] -= frhow_visc
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    for k in (k0 - 1):k1, j in (j0 - 1):j1, i in i0:i1
        coef_v =
            1 / re * 0.5 * (rhostrattfc[i, j, 1] + rhostrattfc[i, j + 1, 1])

        grhow_visc =
            coef_v *
            0.5 *
            (
                jac[i, j, k] *
                jac[i, j, k + 1] *
                (
                    compute_stress_tensor(i, j, k, 3, 1, predictands, grid) +
                    compute_stress_tensor(i, j, k + 1, 3, 1, predictands, grid)
                ) / (jac[i, j, k] + jac[i, j, k + 1]) +
                jac[i, j + 1, k] *
                jac[i, j + 1, k + 1] *
                (
                    compute_stress_tensor(
                        i,
                        j + 1,
                        k,
                        3,
                        1,
                        predictands,
                        grid,
                    ) + compute_stress_tensor(
                        i,
                        j + 1,
                        k + 1,
                        3,
                        1,
                        predictands,
                        grid,
                    )
                ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
            )

        phiw[i, j, k, 2] -= grhow_visc
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    for k in (k0 - 2):k1, j in j0:j1, i in i0:i1
        coef_v = 1 / re * rhostrattfc[i, j, 1]

        hrhow_visc =
            coef_v * (
                jac[i, j, k + 1] *
                met[i, j, k + 1, 1, 3] *
                compute_stress_tensor(i, j, k + 1, 3, 1, predictands, grid) +
                jac[i, j, k + 1] *
                met[i, j, k + 1, 2, 3] *
                compute_stress_tensor(i, j, k + 1, 3, 2, predictands, grid) +
                compute_stress_tensor(i, j, k + 1, 3, 3, predictands, grid)
            )

        phiw[i, j, k, 3] -= hrhow_visc
    end

    # Return.
    return
end

function compute_fluxes!(
    state::State,
    predictands::Predictands,
    tracersetup::NoTracer,
)
    return
end

function compute_fluxes!(
    state::State,
    predictands::Predictands,
    tracersetup::AbstractTracer,
)

    # Get all necessary fields.
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid
    (; pstrattfc) = state.atmosphere
    (; tracerreconstructions, tracerfluxes) = state.tracer

    # Get old wind.
    (u0, v0, w0) = (predictands.u, predictands.v, predictands.w)

    for (fd, field) in enumerate(fieldnames(TracerPredictands))
        for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
            chir = getfield(tracerreconstructions, fd)[i + 1, j, k, 1, 1]
            chil = getfield(tracerreconstructions, fd)[i, j, k, 1, 2]

            pedger =
                0.5 * (
                    jac[i, j, k] * pstrattfc[i, j, k] +
                    jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
                )
            usurf = pedger * u0[i, j, k]

            fchi = compute_flux(usurf, chil, chir)

            getfield(tracerfluxes, fd)[i, j, k, 1] = fchi
        end

        for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
            chif = getfield(tracerreconstructions, fd)[i, j + 1, k, 2, 1]
            chib = getfield(tracerreconstructions, fd)[i, j, k, 2, 2]

            pedgef =
                0.5 * (
                    jac[i, j, k] * pstrattfc[i, j, k] +
                    jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
                )
            vsurf = pedgef * v0[i, j, k]

            gchi = compute_flux(vsurf, chib, chif)

            getfield(tracerfluxes, fd)[i, j, k, 2] = gchi
        end

        for k in (k0 - 1):k1, j in j0:j1, i in i0:i1
            chiu = getfield(tracerreconstructions, fd)[i, j, k + 1, 3, 1]
            chid = getfield(tracerreconstructions, fd)[i, j, k, 3, 2]

            pedgeu =
                jac[i, j, k] *
                jac[i, j, k + 1] *
                (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
                (jac[i, j, k] + jac[i, j, k + 1])
            wsurf = pedgeu * w0[i, j, k]

            hchi = compute_flux(wsurf, chid, chiu)

            getfield(tracerfluxes, fd)[i, j, k, 3] = hchi
        end
    end

    return
end

function compute_fluxes!(
    state::State,
    predictands::Predictands,
    icesetup::NoIce,
)
    return
end

function compute_fluxes!(
    state::State,
    predictands::Predictands,
    icesetup::AbstractIce,
)

    # Get all necessary fields.
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid
    (; pstrattfc) = state.atmosphere
    (; icereconstructions, icefluxes) = state.ice

    # Get old wind.
    (u0, v0, w0) = (predictands.u, predictands.v, predictands.w)

    for (fd, field) in enumerate(fieldnames(IcePredictands))
        for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
            chir = getfield(icereconstructions, fd)[i + 1, j, k, 1, 1]
            chil = getfield(icereconstructions, fd)[i, j, k, 1, 2]

            pedger =
                0.5 * (
                    jac[i, j, k] * pstrattfc[i, j, k] +
                    jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
                )
            usurf = pedger * u0[i, j, k]

            fchi = compute_flux(usurf, chil, chir)

            getfield(icefluxes, fd)[i, j, k, 1] = fchi
        end

        for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
            chif = getfield(icereconstructions, fd)[i, j + 1, k, 2, 1]
            chib = getfield(icereconstructions, fd)[i, j, k, 2, 2]

            pedgef =
                0.5 * (
                    jac[i, j, k] * pstrattfc[i, j, k] +
                    jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
                )
            vsurf = pedgef * v0[i, j, k]

            gchi = compute_flux(vsurf, chib, chif)

            getfield(icefluxes, fd)[i, j, k, 2] = gchi
        end

        for k in (k0 - 1):k1, j in j0:j1, i in i0:i1
            chiu = getfield(icereconstructions, fd)[i, j, k + 1, 3, 1]
            chid = getfield(icereconstructions, fd)[i, j, k, 3, 2]

            pedgeu =
                jac[i, j, k] *
                jac[i, j, k + 1] *
                (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
                (jac[i, j, k] + jac[i, j, k + 1])
            wsurf = pedgeu * w0[i, j, k]

            hchi = compute_flux(wsurf, chid, chiu)

            getfield(icefluxes, fd)[i, j, k, 3] = hchi
        end
    end

    return
end

function compute_fluxes!(
    state::State,
    predictands::Predictands,
    turbulencesetup::NoTurbulence,
)
    return
end

function compute_fluxes!(
    state::State,
    predictands::Predictands,
    turbulencesetup::AbstractTurbulence,
)

    # Get all necessary fields.
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid
    (; pstrattfc) = state.atmosphere
    (; turbulencereconstructions, turbulencefluxes) = state.turbulence

    # Get old wind.
    (u0, v0, w0) = (predictands.u, predictands.v, predictands.w)

    for (fd, field) in enumerate(fieldnames(TurbulencePredictands))
        for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
            chir = getfield(turbulencereconstructions, fd)[i + 1, j, k, 1, 1]
            chil = getfield(turbulencereconstructions, fd)[i, j, k, 1, 2]

            pedger =
                0.5 * (
                    jac[i, j, k] * pstrattfc[i, j, k] +
                    jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
                )
            usurf = pedger * u0[i, j, k]

            fchi = compute_flux(usurf, chil, chir)

            getfield(turbulencefluxes, fd)[i, j, k, 1] = fchi
        end

        for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
            chif = getfield(turbulencereconstructions, fd)[i, j + 1, k, 2, 1]
            chib = getfield(turbulencereconstructions, fd)[i, j, k, 2, 2]

            pedgef =
                0.5 * (
                    jac[i, j, k] * pstrattfc[i, j, k] +
                    jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
                )
            vsurf = pedgef * v0[i, j, k]

            gchi = compute_flux(vsurf, chib, chif)

            getfield(turbulencefluxes, fd)[i, j, k, 2] = gchi
        end

        for k in (k0 - 1):k1, j in j0:j1, i in i0:i1
            chiu = getfield(turbulencereconstructions, fd)[i, j, k + 1, 3, 1]
            chid = getfield(turbulencereconstructions, fd)[i, j, k, 3, 2]

            pedgeu =
                jac[i, j, k] *
                jac[i, j, k + 1] *
                (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
                (jac[i, j, k] + jac[i, j, k + 1])
            wsurf = pedgeu * w0[i, j, k]

            hchi = compute_flux(wsurf, chid, chiu)

            getfield(turbulencefluxes, fd)[i, j, k, 3] = hchi
        end
    end

    return
end
