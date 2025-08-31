"""
```julia
apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::AbstractVariable,
)
```

Perform an implicit substep to integrate the Rayleigh-damping term that represents the unified sponge layer in the prognostic equation for `variable` by dispatching to the appropriate model-specific method.

```julia
apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::Rho,
    model::Boussinesq,
)
```

Return in Boussinesq mode (constant density).

```julia
apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::Rho,
    model::AbstractModel,
)
```

Integrate the Rayleigh-damping term that represents the unified sponge in the continuity equation.

The update is given by

```math
\\rho \\rightarrow \\left(1 + \\alpha_\\mathrm{R} \\Delta t\\right)^{- 1} \\left(\\rho + \\alpha_\\mathrm{R} \\Delta t \\overline{\\rho}\\right),
```

where ``\\alpha_\\mathrm{R}`` is the Rayleigh-damping coefficient computed by [`PinCFlow.Update.compute_sponge!`](@ref) and ``\\Delta t`` is the time step given as input to this method.

```julia
apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::RhoP,
    model::AbstractModel,
)
```

Integrate the Rayleigh-damping term that represents the unified sponge in the auxiliary equation.

The update is given by

```math
\\rho' \\rightarrow \\left(1 + \\alpha_\\mathrm{R} \\Delta t\\right)^{- 1} \\rho'.
```

```julia
apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::U,
    model::AbstractModel,
)
```

Integrate the Rayleigh-damping term that represents the unified sponge in the zonal-momentum equation.

The update is given by

```math
u_{i + 1 / 2} \\rightarrow \\left(1 + \\alpha_{\\mathrm{R}, i + 1 / 2} \\Delta t\\right)^{- 1} \\left\\{u_{i + 1 / 2} + \\alpha_{\\mathrm{R}, i + 1 / 2} \\Delta t u_\\mathrm{r} \\left[1 + a_\\mathrm{r} \\sin \\left(\\frac{2 \\pi t}{t_\\mathrm{r}}\\right)\\right]\\right\\}.
```

If `state.namelists.sponge.relax_to_mean` is `false`, ``u_\\mathrm{r}``, ``a_\\mathrm{r}`` and ``t_\\mathrm{r}`` are given by the sponge-namelist parameters `relaxation_wind[1]`, `perturbation_amplitude` and `perturbation_period`, respectively. Otherwise, ``u_\\mathrm{r}`` is the average of ``u_{i + 1 / 2}`` across the terrain-following coordinate surface and ``a_\\mathrm{r} = 0``.

```julia
apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::V,
    model::AbstractModel,
)
```

Integrate the Rayleigh-damping term that represents the unified sponge in the meridional-momentum equation.

The update is given by

```math
v_{j + 1 / 2} \\rightarrow \\left(1 + \\alpha_{\\mathrm{R}, j + 1 / 2} \\Delta t\\right)^{- 1} \\left\\{v_{j + 1 / 2} + \\alpha_{\\mathrm{R}, j + 1 / 2} \\Delta t v_\\mathrm{r} \\left[1 + a_\\mathrm{r} \\sin \\left(\\frac{2 \\pi t}{t_\\mathrm{r}}\\right)\\right]\\right\\}.
```

The computation of the relaxation wind is analogous to that in the method for the zonal momentum, with ``v_\\mathrm{r}`` given by `state.namelists.sponge.relaxation_wind[2]`.

```julia
apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::W,
    model::AbstractModel,
)
```

Integrate the Rayleigh-damping term that represents the unified sponge in the transformed-vertical-momentum equation.

The update is given by

```math
\\widehat{w}_{k + 1 / 2} \\rightarrow \\left(1 + \\alpha_{\\mathrm{R}, k + 1 / 2} \\Delta t\\right)^{- 1} \\left\\{\\widehat{w}_{k + 1 / 2} + \\alpha_{\\mathrm{R}, k + 1 / 2} \\Delta t \\widehat{w}_\\mathrm{r} \\left[1 + a_\\mathrm{r} \\sin \\left(\\frac{2 \\pi t}{t_\\mathrm{r}}\\right)\\right]\\right\\},
```

The computation of the relaxation wind is analogous to that in the methods for the zonal and meridional momenta, with ``\\widehat{w}_\\mathrm{r}`` given by `state.namelists.sponge.relaxation_wind[3]`.

```julia
apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::PiP,
    model::AbstractModel,
)
```

Return in non-compressible modes (Exner-pressure fluctuations are only updated in the corrector step).

```julia
apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::PiP,
    model::Compressible,
)
```

Update the Exner-pressure fluctuations to account for the Rayleigh damping applied to the mass-weighted potential temperature.

The update is given by

```math
\\pi' \\rightarrow \\pi' - \\alpha_\\mathrm{R} \\Delta t P \\frac{\\partial \\pi'}{\\partial P} \\left(1 - \\frac{\\overline{\\rho}}{\\rho}\\right).
```

```julia
apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::P,
    model::AbstractModel,
)
```

Return in non-compressible modes (mass-weighted potential temperature is constant in time).

```julia
apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::P,
    model::Compressible,
)
```

Integrate the Rayleigh-damping term that represents the unified sponge in the thermodynamic-energy equation.

The update is given by

```math
P \\rightarrow \\left(1 + \\alpha_\\mathrm{R} \\Delta t\\right)^{- 1} P \\left(1 + \\alpha_\\mathrm{R} \\Delta t \\frac{\\overline{\\rho}}{\\rho}\\right).
```

```julia
apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    tracersetup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia
apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    tracersetup::AbstractTracer,
)
```

Integrate the Rayleigh-damping terms that represent the unified sponge in the tracer equations.

In each tracer equation, the update is given by

```math
\\chi \\rightarrow \\left(1 + \\alpha_\\mathrm{R} \\Delta t\\right)^{- 1} \\left(\\chi + \\alpha_\\mathrm{R} \\Delta t \\chi_0\\right),
```

where ``\\chi_0`` is the initial distribution of the tracer.

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `time`: Simulation time.

  - `variable`: Variable to apply Rayleigh damping to.

  - `model`: Dynamic equations.

  - `tracersetup`: General tracer-transport configuration.
"""
function apply_unified_sponge! end

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

function apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::Rho,
    model::Boussinesq,
)
    return
end

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

function apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    tracersetup::NoTracer,
)
    return
end

function apply_unified_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    tracersetup::AbstractTracer,
)
    (; spongelayer, unifiedsponge) = state.namelists.sponge
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; alphaunifiedsponge) = state.sponge
    (; tracerpredictands) = state.tracer
    (; initialtracer) = state.tracer.tracerauxiliaries

    if !spongelayer || !unifiedsponge
        return
    end

    for field in fieldnames(TracerPredictands)
        for k in k0:k1, j in j0:j1, i in i0:i1
            alpha = alphaunifiedsponge[i, j, k]
            chi_old = getfield(tracerpredictands, field)[i, j, k]
            beta = 1.0 / (1.0 + alpha * dt)
            chi_new = (1.0 - beta) * initialtracer[i, j, k] + beta * chi_old
            getfield(tracerpredictands, field)[i, j, k] = chi_new
        end
    end

    return
end
