"""
```julia
reconstruct!(state::State)
```

Reconstruct the prognostic variables at the cell interfaces of their respective grids, using the Monotonic Upstream-centered Scheme for Conservation Laws (MUSCL).

This method calls specialized methods for each prognostic variable.

```julia
reconstruct!(state::State, variable::Rho)
```

Reconstruct the density.

Since the transporting velocity is ``P \\widehat{\\boldsymbol{u}}``, the density is divided by ``P`` before reconstruction.

```julia
reconstruct!(state::State, variable::RhoP)
```

Reconstruct the density fluctuations.

Similar to the density, the density fluctuations are divided by ``P`` before reconstruction.

```julia
reconstruct!(state::State, variable::U)
```

Reconstruct the zonal momentum.

Since the transporting velocity is ``P \\widehat{\\boldsymbol{u}}``, the zonal momentum is divided by ``P`` interpolated to the respective cell interfaces before reconstruction.

```julia
reconstruct!(state::State, variable::V)
```

Reconstruct the meridional momentum.

Similar to the zonal momentum, the meridional momentum is divided by ``P`` interpolated to the respective cell interfaces before reconstruction.

```julia
reconstruct!(state::State, variable::W)
```

Reconstruct the vertical momentum.

The vertical momentum is computed with `compute_vertical_wind`, `set_zonal_boundaries_of_field!` and `set_meridional_boundaries_of_field!`. Similar to the zonal and meridional momenta, the vertical momentum is divided by ``P`` interpolated to the respective cell interfaces before reconstruction.

```julia
reconstruct!(state::State, tracer_setup::NoTracer)
```

Return for configurations without tracer transport.

```julia
reconstruct!(state::State, tracer_setup::LinearTracer)
```

Reconstruct the tracers.

Similar to the density, the tracers are divided by ``P`` before reconstruction.

# Arguments

  - `state`: Model state.

  - `variable`: The reconstructed variable.

  - `tracer_setup`: General tracer-transport configuration.

# See also

  - [`PinCFlow.FluxCalculator.apply_3d_muscl!`](@ref)

  - [`PinCFlow.Update.compute_vertical_wind`](@ref)

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)

  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)
"""
function reconstruct! end

function reconstruct!(state::State)
    (; tracer_setup) = state.namelists.tracer

    reconstruct!(state, Rho())
    reconstruct!(state, RhoP())
    reconstruct!(state, U())
    reconstruct!(state, V())
    reconstruct!(state, W())

    reconstruct!(state, tracer_setup)

    return
end

function reconstruct!(state::State, variable::Rho)
    (; limiter_type) = state.namelists.discretization
    (; k0, k1, nxx, nyy, nzz) = state.domain
    (; rho) = state.variables.predictands
    (; phi) = state.variables.auxiliaries
    (; rhotilde) = state.variables.reconstructions
    (; pbar) = state.atmosphere

    kk = (k0 - 1):(k1 + 1)

    @ivy phi[:, :, kk] .= rho[:, :, kk] ./ pbar[:, :, kk]

    apply_3d_muscl!(phi, rhotilde, nxx, nyy, nzz, limiter_type)

    return
end

function reconstruct!(state::State, variable::RhoP)
    (; limiter_type) = state.namelists.discretization
    (; k0, k1, nxx, nyy, nzz) = state.domain
    (; rhop) = state.variables.predictands
    (; phi) = state.variables.auxiliaries
    (; rhoptilde) = state.variables.reconstructions
    (; pbar) = state.atmosphere

    kk = (k0 - 1):(k1 + 1)

    @ivy phi[:, :, kk] .= rhop[:, :, kk] ./ pbar[:, :, kk]

    apply_3d_muscl!(phi, rhoptilde, nxx, nyy, nzz, limiter_type)

    return
end

function reconstruct!(state::State, variable::U)
    (; limiter_type) = state.namelists.discretization
    (; k0, k1, nxx, nyy, nzz) = state.domain
    (; rho, u) = state.variables.predictands
    (; phi) = state.variables.auxiliaries
    (; utilde) = state.variables.reconstructions
    (; rhobar, pbar) = state.atmosphere

    @ivy for k in (k0 - 1):(k1 + 1), j in 1:nyy, i in 1:(nxx - 1)
        rhoedge =
            0.5 * (
                rho[i, j, k] +
                rho[i + 1, j, k] +
                rhobar[i, j, k] +
                rhobar[i + 1, j, k]
            )
        pedge = 0.5 * (pbar[i, j, k] + pbar[i + 1, j, k])
        phi[i, j, k] = u[i, j, k] * rhoedge / pedge
    end

    apply_3d_muscl!(phi, utilde, nxx, nyy, nzz, limiter_type)

    return
end

function reconstruct!(state::State, variable::V)
    (; limiter_type) = state.namelists.discretization
    (; k0, k1, nxx, nyy, nzz) = state.domain
    (; rho, v) = state.variables.predictands
    (; phi) = state.variables.auxiliaries
    (; vtilde) = state.variables.reconstructions
    (; rhobar, pbar) = state.atmosphere

    @ivy for k in (k0 - 1):(k1 + 1), j in 1:(nyy - 1), i in 1:nxx
        rhoedge =
            0.5 * (
                rho[i, j, k] +
                rho[i, j + 1, k] +
                rhobar[i, j, k] +
                rhobar[i, j + 1, k]
            )
        pedge = 0.5 * (pbar[i, j, k] + pbar[i, j + 1, k])
        phi[i, j, k] = v[i, j, k] * rhoedge / pedge
    end

    apply_3d_muscl!(phi, vtilde, nxx, nyy, nzz, limiter_type)

    return
end

function reconstruct!(state::State, variable::W)
    (; namelists, domain, grid) = state
    (; limiter_type) = state.namelists.discretization
    (; i0, i1, j0, j1, k0, k1, nxx, nyy, nzz) = domain
    (; jac) = grid
    (; predictands) = state.variables
    (; rho, w) = predictands
    (; phi) = state.variables.auxiliaries
    (; wtilde) = state.variables.reconstructions
    (; rhobar, pbar) = state.atmosphere

    @ivy phi[:, :, (k0 - 1):(k1 + 1)] .= w[:, :, (k0 - 1):(k1 + 1)]

    @ivy for k in (k0 - 1):(k1 + 1), j in j0:j1, i in i0:i1
        phi[i, j, k] = compute_vertical_wind(i, j, k, state)
    end

    set_zonal_boundaries_of_field!(phi, namelists, domain)
    set_meridional_boundaries_of_field!(phi, namelists, domain)

    @ivy for k in (k0 - 1):(k1 + 1), j in 1:nyy, i in 1:nxx
        rhoedgeu =
            (
                jac[i, j, k + 1] * (rho[i, j, k] + rhobar[i, j, k]) +
                jac[i, j, k] * (rho[i, j, k + 1] + rhobar[i, j, k + 1])
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        pedgeu =
            (
                jac[i, j, k + 1] * pbar[i, j, k] +
                jac[i, j, k] * pbar[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        phi[i, j, k] *= rhoedgeu / pedgeu
    end

    apply_3d_muscl!(phi, wtilde, nxx, nyy, nzz, limiter_type)

    return
end

function reconstruct!(state::State, tracer_setup::NoTracer)
    return
end

function reconstruct!(state::State, tracer_setup::LinearTracer)
    (; limiter_type) = state.namelists.discretization
    (; k0, k1, nxx, nyy, nzz) = state.domain
    (; phi) = state.variables.auxiliaries
    (; pbar) = state.atmosphere
    (; tracerreconstructions, tracerpredictands) = state.tracer

    @ivy for field in 1:fieldcount(TracerPredictands)
        chi = getfield(tracerpredictands, field)[:, :, :]
        for k in (k0 - 1):(k1 + 1), j in 1:nyy, i in 1:nxx
            phi[i, j, k] = chi[i, j, k] / pbar[i, j, k]
        end
        apply_3d_muscl!(
            phi,
            getfield(tracerreconstructions, field),
            nxx,
            nyy,
            nzz,
            limiter_type,
        )
    end

    return
end
