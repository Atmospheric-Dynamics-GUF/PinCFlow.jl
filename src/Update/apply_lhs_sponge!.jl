"""
```julia
apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::AbstractVariable,
)
```

Perform an implicit substep to integrate the Rayleigh-damping term that represents the LHS sponge in the prognostic equation for `variable` by dispatching to the appropriate model-specific method.

```julia
apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::Rho,
    model::Boussinesq,
)
```

Return in Boussinesq mode (constant density).

```julia
apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::Rho,
    model::AbstractModel,
)
```

Integrate the Rayleigh-damping term that represents the LHS sponge in the continuity equation.

The update is given by

```math
\\rho \\rightarrow \\left(1 + \\alpha_\\mathrm{R} \\Delta t\\right)^{- 1} \\left(\\rho + \\alpha_\\mathrm{R} \\Delta t \\overline{\\rho}\\right),
```

where ``\\alpha_\\mathrm{R}`` is the Rayleigh-damping coefficient computed by [`PinCFlow.Update.compute_sponges!`](@ref) and ``\\Delta t`` is the time step given as input to this method.

```julia
apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::RhoP,
    model::AbstractModel,
)
```

Integrate the Rayleigh-damping term that represents the LHS sponge in the auxiliary equation.

The update is given by

```math
\\rho' \\rightarrow \\left(1 + \\alpha_\\mathrm{R} \\Delta t\\right)^{- 1} \\left[\\rho' + \\alpha_\\mathrm{R} \\Delta t \\overline{\\rho} \\left(1 - \\frac{P}{\\rho \\overline{\\theta}}\\right)\\right].
```

```julia
apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::U,
    model::AbstractModel,
)
```

Integrate the Rayleigh-damping term that represents the LHS sponge in the zonal-momentum equation.

The update is given by

```math
u_{i + 1 / 2} \\rightarrow \\left(1 + \\alpha_{\\mathrm{R}, i + 1 / 2} \\Delta t\\right)^{- 1} \\left\\{u_{i + 1 / 2} + \\alpha_{\\mathrm{R}, i + 1 / 2} \\Delta t u_\\mathrm{r} \\left[1 + a_\\mathrm{r} \\sin \\left(\\frac{2 \\pi t}{t_\\mathrm{r}}\\right)\\right]\\right\\}.
```

If `state.namelists.sponge.relax_to_mean` is `false`, ``u_\\mathrm{r}``, ``a_\\mathrm{r}`` and ``t_\\mathrm{r}`` are given by the sponge-namelist parameters `relaxation_wind[1]`, `perturbation_amplitude` and `perturbation_period`, respectively. Otherwise, ``u_\\mathrm{r}`` is the average of ``u_{i + 1 / 2}`` across the terrain-following coordinate surface and ``a_\\mathrm{r} = 0``.

```julia
apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::V,
    model::AbstractModel,
)
```

Integrate the Rayleigh-damping term that represents the LHS sponge in the meridional-momentum equation.

The update is given by

```math
v_{j + 1 / 2} \\rightarrow \\left(1 + \\alpha_{\\mathrm{R}, j + 1 / 2} \\Delta t\\right)^{- 1} \\left\\{v_{j + 1 / 2} + \\alpha_{\\mathrm{R}, j + 1 / 2} \\Delta t v_\\mathrm{r} \\left[1 + a_\\mathrm{r} \\sin \\left(\\frac{2 \\pi t}{t_\\mathrm{r}}\\right)\\right]\\right\\}.
```

The computation of the relaxation wind is analogous to that in the method for the zonal momentum, with ``v_\\mathrm{r}`` given by `state.namelists.sponge.relaxation_wind[2]`.

```julia
apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::W,
    model::AbstractModel,
)
```

Integrate the Rayleigh-damping term that represents the LHS sponge in the transformed-vertical-momentum equation.

The update is given by

```math
\\widehat{w}_{k + 1 / 2} \\rightarrow \\left(1 + \\alpha_{\\mathrm{R}, k + 1 / 2} \\Delta t\\right)^{- 1} \\left\\{\\widehat{w}_{k + 1 / 2} + \\alpha_{\\mathrm{R}, k + 1 / 2} \\Delta t \\widehat{w}_\\mathrm{r} \\left[1 + a_\\mathrm{r} \\sin \\left(\\frac{2 \\pi t}{t_\\mathrm{r}}\\right)\\right]\\right\\},
```

The computation of the relaxation wind is analogous to that in the methods for the zonal and meridional momenta, with ``\\widehat{w}_\\mathrm{r}`` given by `state.namelists.sponge.relaxation_wind[3]`.

```julia
apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::PiP,
    model::AbstractModel,
)
```

Return in non-compressible modes (Exner-pressure fluctuations are only updated in the corrector step).

```julia
apply_lhs_sponge!(
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
apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::P,
    model::AbstractModel,
)
```

Return in non-compressible modes (mass-weighted potential temperature is constant in time).

```julia
apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::P,
    model::Compressible,
)
```

Integrate the Rayleigh-damping term that represents the LHS sponge in the thermodynamic-energy equation.

The update is given by

```math
P \\rightarrow \\left(1 + \\alpha_\\mathrm{R} \\Delta t\\right)^{- 1} P \\left(1 + \\alpha_\\mathrm{R} \\Delta t \\frac{\\overline{\\rho}}{\\rho}\\right).
```

```julia
apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    tracer_setup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia
apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    tracer_setup::AbstractTracer,
)
```

Integrate the Rayleigh-damping terms that represent the LHS sponge in the tracer equations.

In each tracer equation, the update is given by

```math
\\left(\\rho \\chi\\right) \\rightarrow \\left(1 + \\alpha_\\mathrm{R} \\Delta t\\right)^{- 1} \\left[\\rho \\chi + \\alpha_\\mathrm{R} \\Delta t \\left(\\rho \\chi\\right)^{\\left(0\\right)}\\right].
```

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `time`: Simulation time.

  - `variable`: Variable to apply Rayleigh damping to.

  - `model`: Dynamic equations.

  - `tracer_setup`: General tracer-transport configuration.
"""
function apply_lhs_sponge! end

function apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::AbstractVariable,
)
    (; model) = state.namelists.setting
    apply_lhs_sponge!(state, dt, time, variable, model)
    return
end

function apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::Rho,
    model::Boussinesq,
)
    return
end

function apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::Rho,
    model::AbstractModel,
)
    (; use_sponge) = state.namelists.sponge
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; alphar) = state.sponge
    (; rho) = state.variables.predictands

    if !use_sponge
        return
    end

    rhobg = 0.0
    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        alpha = alphar[i, j, k]
        rhoold = rho[i, j, k]
        beta = 1.0 / (1.0 + alpha * dt)
        rhonew = (1.0 - beta) * rhobg + beta * rhoold
        rho[i, j, k] = rhonew
    end

    return
end

function apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::RhoP,
    model::Compressible,
)
    (; use_sponge) = state.namelists.sponge
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; rhostrattfc, thetastrattfc) = state.atmosphere
    (; alphar) = state.sponge
    (; rho, rhop, p) = state.variables.predictands

    if !use_sponge
        return
    end

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        rhopbg =
            rhostrattfc[i, j, k] * (
                1.0 -
                p[i, j, k] / thetastrattfc[i, j, k] /
                (rho[i, j, k] + rhostrattfc[i, j, k])
            )
        alpha = alphar[i, j, k]
        rhopold = rhop[i, j, k]
        beta = 1.0 / (1.0 + alpha * dt)
        rhopnew = (1.0 - beta) * rhopbg + beta * rhopold
        rhop[i, j, k] = rhopnew
    end

    return
end

function apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::RhoP,
    model::AbstractModel,
)
    (; use_sponge) = state.namelists.sponge
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; alphar) = state.sponge
    (; rhop) = state.variables.predictands

    if !use_sponge
        return
    end

    rhobg = 0.0
    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        alpha = alphar[i, j, k]
        rhoold = rhop[i, j, k]
        beta = 1.0 / (1.0 + alpha * dt)
        rhonew = (1.0 - beta) * rhobg + beta * rhoold
        rhop[i, j, k] = rhonew
    end

    return
end

function apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::U,
    model::AbstractModel,
)
    (; ndx, ndy) = state.namelists.domain
    (;
        use_sponge,
        relax_to_mean,
        perturbation_period,
        perturbation_amplitude,
        relaxation_wind,
    ) = state.namelists.sponge
    (; uref, tref) = state.constants
    (; layer_comm, i0, i1, j0, j1, k0, k1) = state.domain
    (; alphar, horizontal_mean) = state.sponge
    (; u) = state.variables.predictands

    if !use_sponge
        return
    end

    (ii, jj, kk) = (i0:i1, j0:j1, k0:k1)

    horizontal_mean .= 0.0

    # Determine relaxation wind.
    @ivy if relax_to_mean
        horizontal_mean .=
            sum(a -> a / ndx / ndy, u[ii, jj, kk]; dims = (1, 2))[1, 1, :]
        MPI.Allreduce!(horizontal_mean, +, layer_comm)
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
    @ivy for k in kk
        if relax_to_mean
            ubg = horizontal_mean[k - k0 + 1]
        end
        for j in jj, i in ii
            alpha = 0.5 * (alphar[i, j, k] + alphar[i + 1, j, k])
            uold = u[i, j, k]
            beta = 1.0 / (1.0 + alpha * dt)
            unew = (1.0 - beta) * ubg + beta * uold
            u[i, j, k] = unew
        end
    end

    return
end

function apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::V,
    model::AbstractModel,
)
    (; ndx, ndy) = state.namelists.domain
    (;
        use_sponge,
        relax_to_mean,
        perturbation_period,
        perturbation_amplitude,
        relaxation_wind,
    ) = state.namelists.sponge
    (; uref, tref) = state.constants
    (; layer_comm, i0, i1, j0, j1, k0, k1) = state.domain
    (; alphar, horizontal_mean) = state.sponge
    (; v) = state.variables.predictands

    if !use_sponge
        return
    end

    (ii, jj, kk) = (i0:i1, j0:j1, k0:k1)

    horizontal_mean .= 0.0

    # Determine relaxation wind.
    @ivy if relax_to_mean
        horizontal_mean .=
            sum(a -> a / ndx / ndy, v[ii, jj, kk]; dims = (1, 2))[1, 1, :]
        MPI.Allreduce!(horizontal_mean, +, layer_comm)
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
    @ivy for k in kk
        if relax_to_mean
            vbg = horizontal_mean[k - k0 + 1]
        end
        for j in jj, i in ii
            alpha = 0.5 * (alphar[i, j, k] + alphar[i, j + 1, k])
            vold = v[i, j, k]
            beta = 1.0 / (1.0 + alpha * dt)
            vnew = (1.0 - beta) * vbg + beta * vold
            v[i, j, k] = vnew
        end
    end

    return
end

function apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::W,
    model::AbstractModel,
)
    (; ndx, ndy) = state.namelists.domain
    (;
        use_sponge,
        relax_to_mean,
        perturbation_period,
        perturbation_amplitude,
        relaxation_wind,
    ) = state.namelists.sponge
    (; uref, tref) = state.constants
    (; layer_comm, i0, i1, j0, j1, k0, k1) = state.domain
    (; alphar, horizontal_mean) = state.sponge
    (; w) = state.variables.predictands
    (; jac) = state.grid

    if !use_sponge
        return
    end

    (ii, jj, kk) = (i0:i1, j0:j1, k0:k1)

    horizontal_mean .= 0.0

    # Determine relaxation wind.
    @ivy if relax_to_mean
        horizontal_mean .=
            sum(a -> a / ndx / ndy, w[ii, jj, kk]; dims = (1, 2))[1, 1, :]
        MPI.Allreduce!(horizontal_mean, +, layer_comm)
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
    @ivy for k in kk
        if relax_to_mean
            wbg = horizontal_mean[k - k0 + 1]
        end
        for j in jj, i in ii
            alpha =
                (
                    jac[i, j, k + 1] * alphar[i, j, k] +
                    jac[i, j, k] * alphar[i, j, k + 1]
                ) / (jac[i, j, k] + jac[i, j, k + 1])
            wold = w[i, j, k]
            beta = 1.0 / (1.0 + alpha * dt)
            wnew = (1.0 - beta) * wbg + beta * wold
            w[i, j, k] = wnew
        end
    end

    return
end

function apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::PiP,
    model::AbstractModel,
)
    return
end

function apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::PiP,
    model::Compressible,
)
    (; use_sponge) = state.namelists.sponge
    (; gamma, rsp, pref) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; alphar) = state.sponge
    (; rhostrattfc) = state.atmosphere
    (; rho, pip, p) = state.variables.predictands

    if !use_sponge
        return
    end

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        dpdpi =
            1 / (gamma - 1) * (rsp / pref)^(1 - gamma) * p[i, j, k]^(2 - gamma)
        pib =
            rhostrattfc[i, j, k] * p[i, j, k] /
            (rho[i, j, k] + rhostrattfc[i, j, k]) / dpdpi
        alpha = alphar[i, j, k]
        pipold = pip[i, j, k]
        pipnew = pipold - alpha * dt * (p[i, j, k] / dpdpi - pib)
        pip[i, j, k] = pipnew
    end

    return
end

function apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::P,
    model::AbstractModel,
)
    return
end

function apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::P,
    model::Compressible,
)
    (; use_sponge) = state.namelists.sponge
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; alphar) = state.sponge
    (; rhostrattfc) = state.atmosphere
    (; rho, p) = state.variables.predictands

    if !use_sponge
        return
    end

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        pb =
            rhostrattfc[i, j, k] * p[i, j, k] /
            (rho[i, j, k] + rhostrattfc[i, j, k])
        alpha = alphar[i, j, k]
        pold = p[i, j, k]
        beta = 1 / (1 + alpha * dt)
        pnew = (1 - beta) * pb + beta * pold
        p[i, j, k] = pnew
    end

    return
end

function apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    tracer_setup::NoTracer,
)
    return
end

function apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    tracer_setup::AbstractTracer,
)
    (; use_sponge) = state.namelists.sponge
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; alphar) = state.sponge
    (; tracerpredictands) = state.tracer
    (; initialtracer) = state.tracer.tracerauxiliaries

    if !use_sponge
        return
    end

    @ivy for field in fieldnames(TracerPredictands)
        for k in k0:k1, j in j0:j1, i in i0:i1
            alpha = alphar[i, j, k]
            chi_old = getfield(tracerpredictands, field)[i, j, k]
            beta = 1.0 / (1.0 + alpha * dt)
            chi_new = (1.0 - beta) * initialtracer[i, j, k] + beta * chi_old
            getfield(tracerpredictands, field)[i, j, k] = chi_new
        end
    end

    return
end
