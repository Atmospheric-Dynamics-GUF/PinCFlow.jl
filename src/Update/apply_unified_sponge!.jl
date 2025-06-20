"""
    apply_unified_sponge!(state::State, dt::AbstractFloat, time::AbstractFloat, variable::AbstractVariable)

Apply unified sponge layer damping to prevent spurious wave reflections.

# Mathematical Formulation

Sponge damping uses exponential relaxation: `φ_new = β*φ_old + (1-β)*φ_bg`
where `β = 1/(1 + α*dt)` and α is the spatially-varying damping coefficient.

# Implementation Strategy

  - **Variable-specific backgrounds**: Each variable relaxes toward appropriate reference state
  - **Staggered grid handling**: Damping coefficients averaged at edge locations
  - **Time-dependent perturbations**: Optional sinusoidal modulation of background states
  - **Horizontal averaging**: Option to relax toward instantaneous horizontal mean

The unified approach uses single damping coefficient field for all variables,
computed by [`compute_sponge!`](@ref) with various spatial profiles.
"""
function apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::AbstractVariable,
)
    (; model) = state.namelists.setting
    apply_unified_sponge!(state, dt, time, variable, model)
    return
end

"""
    apply_unified_sponge!(state::State, dt::AbstractFloat, time::AbstractFloat, variable::Rho, model::Boussinesq)

No-op for Boussinesq density (incompressible assumption).
"""
function apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::Rho,
    model::Boussinesq,
)
    return
end

"""
    apply_unified_sponge!(state::State, dt::AbstractFloat, time::AbstractFloat, variable::Rho, model::AbstractModel)

Apply sponge damping to density perturbations.

Relaxes toward zero background density perturbation: ρ' → 0.
"""
function apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::Rho,
    model::AbstractModel,
)
    (; spongelayer, unifiedsponge) = state.namelists.sponge
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; alphaunifiedsponge) = state.sponge
    (; rho) = state.variables.predictands

    if !spongelayer || !unifiedsponge
        return
    end

    rho_bg = 0.0
    for k in k0:k1, j in j0:j1, i in i0:i1
        alpha = alphaunifiedsponge[i, j, k]
        rho_old = rho[i, j, k]
        beta = 1.0 / (1.0 + alpha * dt)
        rho_new = (1.0 - beta) * rho_bg + beta * rho_old
        rho[i, j, k] = rho_new
    end

    return
end

"""
    apply_unified_sponge!(state::State, dt::AbstractFloat, time::AbstractFloat, variable::RhoP, model::AbstractModel)

Apply sponge damping to potential density perturbations.

Relaxes toward zero background potential density perturbation.
"""
function apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::RhoP,
    model::AbstractModel,
)
    (; spongelayer, unifiedsponge) = state.namelists.sponge
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; alphaunifiedsponge) = state.sponge
    (; rhop) = state.variables.predictands

    if !spongelayer || !unifiedsponge
        return
    end

    rho_bg = 0.0
    for k in k0:k1, j in j0:j1, i in i0:i1
        alpha = alphaunifiedsponge[i, j, k]
        rho_old = rhop[i, j, k]
        beta = 1.0 / (1.0 + alpha * dt)
        rho_new = (1.0 - beta) * rho_bg + beta * rho_old
        rhop[i, j, k] = rho_new
    end

    return
end

"""
    apply_unified_sponge!(state::State, dt::AbstractFloat, time::AbstractFloat, variable::U, model::AbstractModel)

Apply sponge damping to zonal wind with background relaxation.

# Relaxation Options

  - **Mean relaxation**: Toward horizontal mean if `relax_to_mean=true`
  - **Fixed background**: Toward `relaxation_wind[1]` with optional perturbations
  - **Edge averaging**: Uses average of adjacent cell damping coefficients

# Time Modulation

Background wind can vary sinusoidally: `u_bg = u₀(1 + A*sin(2πt/T))`
where A is amplitude and T is period.
"""
function apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::U,
    model::AbstractModel,
)
    (; sizex, sizey) = state.namelists.domain
    (;
        spongelayer,
        unifiedsponge,
        relax_to_mean,
        perturbation_period,
        perturbation_amplitude,
        relaxation_wind,
    ) = state.namelists.sponge
    (; uref, tref) = state.constants
    (; layer_comm, i0, i1, j0, j1, k0, k1) = state.domain
    (; alphaunifiedsponge, horizontal_mean) = state.sponge
    (; u) = state.variables.predictands

    if !spongelayer || !unifiedsponge
        return
    end

    horizontal_mean .= 0.0

    # Determine relaxation wind.
    if relax_to_mean
        for k in k0:k1
            @views horizontal_mean[k - k0 + 1] = sum(u[i0:i1, j0:j1, k])
        end
        MPI.Allreduce!(horizontal_mean, +, layer_comm)
        horizontal_mean ./= (sizex .* sizey)
    else
        ubg = relaxation_wind[1] / uref
        if perturbation_period > 0.0
            ubg =
                ubg * (
                    1.0 +
                    perturbation_amplitude *
                    sin(2.0 * pi * time / perturbation_period * tref)
                )
        end
    end

    # Update the zonal wind.
    for k in k0:k1
        if relax_to_mean
            ubg = horizontal_mean[k - k0 + 1]
        end
        for j in j0:j1, i in i0:i1
            alpha =
                0.5 *
                (alphaunifiedsponge[i, j, k] + alphaunifiedsponge[i + 1, j, k])
            uold = u[i, j, k]
            beta = 1.0 / (1.0 + alpha * dt)
            unew = (1.0 - beta) * ubg + beta * uold
            u[i, j, k] = unew
        end
    end

    return
end

"""
    apply_unified_sponge!(state::State, dt::AbstractFloat, time::AbstractFloat, variable::V, model::AbstractModel)

Apply sponge damping to meridional wind with background relaxation.

Similar to zonal wind but for y-component with appropriate staggering.
Uses `relaxation_wind[2]` for fixed background state.
"""
function apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::V,
    model::AbstractModel,
)
    (; sizex, sizey) = state.namelists.domain
    (;
        spongelayer,
        unifiedsponge,
        relax_to_mean,
        perturbation_period,
        perturbation_amplitude,
        relaxation_wind,
    ) = state.namelists.sponge
    (; uref, tref) = state.constants
    (; layer_comm, i0, i1, j0, j1, k0, k1) = state.domain
    (; alphaunifiedsponge, horizontal_mean) = state.sponge
    (; v) = state.variables.predictands

    if !spongelayer || !unifiedsponge
        return
    end

    horizontal_mean .= 0.0

    # Determine relaxation wind.
    if relax_to_mean
        for k in k0:k1
            @views horizontal_mean[k - k0 + 1] = sum(v[i0:i1, j0:j1, k])
        end
        MPI.Allreduce!(horizontal_mean, +, layer_comm)
        horizontal_mean ./= (sizex .* sizey)
    else
        vbg = relaxation_wind[2] / uref
        if perturbation_period > 0.0
            vbg =
                vbg * (
                    1.0 +
                    perturbation_amplitude *
                    sin(2.0 * pi * time / perturbation_period * tref)
                )
        end
    end

    # Update the meridional wind.
    for k in k0:k1
        if relax_to_mean
            vbg = horizontal_mean[k - k0 + 1]
        end
        for j in j0:j1, i in i0:i1
            alpha =
                0.5 *
                (alphaunifiedsponge[i, j, k] + alphaunifiedsponge[i, j + 1, k])
            vold = v[i, j, k]
            beta = 1.0 / (1.0 + alpha * dt)
            vnew = (1.0 - beta) * vbg + beta * vold
            v[i, j, k] = vnew
        end
    end

    return
end

"""
    apply_unified_sponge!(state::State, dt::AbstractFloat, time::AbstractFloat, variable::W, model::AbstractModel)

Apply sponge damping to vertical wind with Jacobian-weighted averaging.

Uses harmonic mean of damping coefficients weighted by Jacobian for proper
vertical staggering in terrain-following coordinates:
`α = (J[k+1]*α[k] + J[k]*α[k+1])/(J[k] + J[k+1])`

Uses `relaxation_wind[3]` for fixed background vertical velocity.
"""
function apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::W,
    model::AbstractModel,
)
    (; sizex, sizey) = state.namelists.domain
    (;
        spongelayer,
        unifiedsponge,
        relax_to_mean,
        perturbation_period,
        perturbation_amplitude,
        relaxation_wind,
    ) = state.namelists.sponge
    (; uref, tref) = state.constants
    (; layer_comm, i0, i1, j0, j1, k0, k1) = state.domain
    (; alphaunifiedsponge, horizontal_mean) = state.sponge
    (; w) = state.variables.predictands
    (; jac) = state.grid

    if !spongelayer || !unifiedsponge
        return
    end

    horizontal_mean .= 0.0

    # Determine relaxation wind.
    if relax_to_mean
        for k in k0:k1
            @views horizontal_mean[k - k0 + 1] = sum(w[i0:i1, j0:j1, k])
        end
        MPI.Allreduce!(horizontal_mean, +, layer_comm)
        horizontal_mean ./= (sizex .* sizey)
    else
        wbg = relaxation_wind[3] / uref
        if perturbation_period > 0.0
            wbg =
                wbg * (
                    1.0 +
                    perturbation_amplitude *
                    sin(2.0 * pi * time / perturbation_period * tref)
                )
        end
    end

    # Update the vertical wind.
    for k in k0:k1
        if relax_to_mean
            wbg = horizontal_mean[k - k0 + 1]
        end
        for j in j0:j1, i in i0:i1
            alpha =
                (
                    jac[i, j, k + 1] * alphaunifiedsponge[i, j, k] +
                    jac[i, j, k] * alphaunifiedsponge[i, j, k + 1]
                ) / (jac[i, j, k] + jac[i, j, k + 1])
            wold = w[i, j, k]
            beta = 1.0 / (1.0 + alpha * dt)
            wnew = (1.0 - beta) * wbg + beta * wold
            w[i, j, k] = wnew
        end
    end

    return
end

function apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::PiP,
    model::AbstractModel,
)
    return
end

function apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::PiP,
    model::Compressible,
)
    (; spongelayer, unifiedsponge) = state.namelists.sponge
    (; gamma, rsp, pref) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; alphaunifiedsponge) = state.sponge
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; rho, pip, p) = state.variables.predictands

    if !spongelayer || !unifiedsponge
        return
    end

    for k in k0:k1, j in j0:j1, i in i0:i1
        dpdpi =
            1 / (gamma - 1) * (rsp / pref)^(1 - gamma) * p[i, j, k]^(2 - gamma)
        pib =
            rhostrattfc[i, j, k] * pstrattfc[i, j, k] /
            (rho[i, j, k] + rhostrattfc[i, j, k]) / dpdpi
        alpha = alphaunifiedsponge[i, j, k]
        pipold = pip[i, j, k]
        pipnew = pipold - alpha * dt * (pstrattfc[i, j, k] / dpdpi - pib)
        pip[i, j, k] = pipnew
    end

    return
end

function apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::P,
    model::AbstractModel,
)
    return
end

"""
    apply_unified_sponge!(state::State, dt::AbstractFloat, time::AbstractFloat, variable::P, model::Compressible)

Apply sponge damping to pressure field for compressible model.

Relaxes toward hydrostatic background pressure: `p_bg = ρ₀*p₀/(ρ + ρ₀)`
"""
function apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::P,
    model::Compressible,
)
    (; spongelayer, unifiedsponge) = state.namelists.sponge
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; alphaunifiedsponge) = state.sponge
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; rho, p) = state.variables.predictands

    if !spongelayer || !unifiedsponge
        return
    end

    for k in k0:k1, j in j0:j1, i in i0:i1
        pb =
            rhostrattfc[i, j, k] * pstrattfc[i, j, k] /
            (rho[i, j, k] + rhostrattfc[i, j, k])
        alpha = alphaunifiedsponge[i, j, k]
        pold = p[i, j, k]
        beta = 1 / (1 + alpha * dt)
        pnew = (1 - beta) * pb + beta * pold
        p[i, j, k] = pnew
    end

    return
end
