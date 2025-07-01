"""
```julia
compute_fluxes!(state::State, predictands::Predictands)
```

Main function for computing fluxes for all variables. This function calls the specialized
flux calculation functions for each variable (Rho, RhoP, U, V, W, and P).

# Arguments

  - `state::State`: The current state containing grid information, variables, and other data
  - `predictands::Predictands`: The predictand variables used in flux calculations

# Returns

  - `nothing`: This function modifies the `predictands` in-place
"""
function compute_fluxes!(state::State, predictands::Predictands)
    (; model) = state.namelists.setting

    compute_fluxes!(state, predictands, Rho())
    compute_fluxes!(state, predictands, RhoP())
    compute_fluxes!(state, predictands, U())
    compute_fluxes!(state, predictands, V())
    compute_fluxes!(state, predictands, W())

    compute_fluxes!(state, predictands, model, P())

    return
end

"""
```julia
compute_fluxes!(state::State, predictands::Predictands, variable::Rho)
```

Compute fluxes for the density (Rho) variable in all three directions.

This function calculates zonal, meridional, and vertical fluxes for the density
variable, using the reconstructed density fields and the old wind fields.

# Arguments

  - `state::State`: The current state containing grid information, variables, and other data
  - `predictands::Predictands`: The predictand variables containing the wind fields (u, v, w)

# Returns

  - `nothing`: This function modifies the `phirho` flux field in the state in-place
"""
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

"""
```julia
compute_fluxes!(state::State, predictands::Predictands, variable::RhoP)
```

Compute fluxes for the density perturbations (RhoP) variable in all three directions.

This function calculates zonal, meridional, and vertical fluxes for the perturbed
variable, using the reconstructed perturbed density field and the old wind fields.

# Arguments

  - `state::State`: The current state containing grid information, variables, and other data
  - `predictands::Predictands`: The predictand variables containing the wind fields (u, v, w)

# Returns

  - `nothing`: This function modifies the `phirhop` flux field in the state in-place
"""
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

"""
```julia
compute_fluxes!(
    state::State,
    predictands::Predictands,
    model::AbstractModel,
    variable::P,
)
```

Return in non-compressible modes.
"""
function compute_fluxes!(
    state::State,
    predictands::Predictands,
    model::AbstractModel,
    variable::P,
)
    return
end

"""
```julia
compute_fluxes!(state::State, predictands::Predictands, variable::RhoP)
```

Compute fluxes for the pressure field variable in all three directions.

This function calculates zonal, meridional, and vertical fluxes for the pressure variable, using the reconstructed pressure  field and the old wind fields.

# Arguments

  - `state::State`: The current state containing grid information, variables, and other data
  - `predictands::Predictands`: The predictand variables containing the wind fields (u, v, w)

# Returns

  - `nothing`: This function modifies the `phip` flux field in the state in-place
"""
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

"""
```julia
compute_fluxes!(state::State, old_predictands::Predictands, variable::U)
compute_fluxes!(state::State, old_predictands::Predictands, variable::V)
compute_fluxes!(state::State, old_predictands::Predictands, variable::W)
```

Compute momentum fluxes for velocity components (U, V, W) in all three directions.

These methods calculate advective and viscous fluxes for momentum in the:

  - U: zonal direction (east-west)
  - V: meridional direction (north-south)
  - W: vertical direction (up-down)

The implementation follows the same pattern for all three components, with appropriate
adjustments for the direction-specific reconstructions and stress tensor components.

# Arguments

  - `state::State`: Current simulation state with grid, domain and variable data
  - `old_predictands::Predictands`: Previous timestep's wind fields
  - `variable::Union{U,V,W}`: Type parameter indicating which momentum component to compute

# Returns

  - `nothing`: Modifies the respective flux fields (`phiu`, `phiv`, `phiw`) in-place

# See also

  - [`compute_flux`](@ref)
  - [`compute_stress_tensor`](@ref)
"""
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
