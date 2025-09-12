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

Since the transporting velocity is ``P \\widehat{\\boldsymbol{u}}``, the density is divided by ``P`` before reconstruction. The result is written into `state.variables.reconstructions.rhotilde`.

```julia
reconstruct!(state::State, variable::RhoP)
```

Reconstruct the density fluctuations.

Similar to the density, the density fluctuations are divided by ``P`` before reconstruction. The result is written into `state.variables.reconstructions.rhoptilde`.

```julia
reconstruct!(state::State, variable::U)
```

Reconstruct the zonal momentum.

Since the transporting velocity is ``P \\widehat{\\boldsymbol{u}}``, the zonal momentum is divided by ``P`` interpolated to the respective cell interfaces before reconstruction. The result is written into `state.variables.reconstructions.utilde`.

```julia
reconstruct!(state::State, variable::V)
```

Reconstruct the meridional momentum.

Similar to the zonal momentum, the meridional momentum is divided by ``P`` interpolated to the respective cell interfaces before reconstruction. The result is written into `state.variables.reconstructions.vtilde`.

```julia
reconstruct!(state::State, variable::W)
```

Reconstruct the vertical momentum.

The vertical momentum is computed with `compute_vertical_wind`, `set_zonal_boundaries_of_field!` and `set_meridional_boundaries_of_field!`. Similar to the zonal and meridional momenta, the vertical momentum is divided by ``P`` interpolated to the respective cell interfaces before reconstruction. The result is written into `state.variables.reconstructions.wtilde`.

```julia
reconstruct!(state::State, tracersetup::NoTracer)
```

Return for configurations without tracer transport.

```julia
reconstruct!(state::State, tracersetup::AbstractTracer)
```

Reconstruct the tracers.

Similar to the density, the tracers are divided by ``P`` before reconstruction.

# Arguments

  - `state`: Model state.

  - `variable`: The reconstructed variable.

  - `tracersetup`: General tracer-transport configuration.

# See also

  - [`PinCFlow.FluxCalculator.apply_3d_muscl!`](@ref)

  - [`PinCFlow.Update.compute_vertical_wind`](@ref)

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)

  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)
"""
function reconstruct! end

function reconstruct!(state::State)
    (; tracersetup) = state.namelists.tracer

    reconstruct!(state, Rho())
    reconstruct!(state, RhoP())
    reconstruct!(state, U())
    reconstruct!(state, V())
    reconstruct!(state, W())

    reconstruct!(state, tracersetup)

    return
end

function reconstruct!(state::State, variable::Rho)
    (; limitertype) = state.namelists.discretization
    (; k0, k1, nxx, nyy, nzz) = state.domain
    (; rho) = state.variables.predictands
    (; phi) = state.variables.auxiliaries
    (; rhotilde) = state.variables.reconstructions
    (; pstrattfc) = state.atmosphere

    @ivy for kz in (k0 - 1):(k1 + 1), jy in 1:nyy, ix in 1:nxx
        phi[ix, jy, kz] = rho[ix, jy, kz] / pstrattfc[ix, jy, kz]
    end
    apply_3d_muscl!(phi, rhotilde, nxx, nyy, nzz, limitertype)

    return
end

function reconstruct!(state::State, variable::RhoP)
    (; limitertype) = state.namelists.discretization
    (; k0, k1, nxx, nyy, nzz) = state.domain
    (; rhop) = state.variables.predictands
    (; phi) = state.variables.auxiliaries
    (; rhoptilde) = state.variables.reconstructions
    (; pstrattfc) = state.atmosphere

    @ivy for kz in (k0 - 1):(k1 + 1), jy in 1:nyy, ix in 1:nxx
        phi[ix, jy, kz] = rhop[ix, jy, kz] / pstrattfc[ix, jy, kz]
    end
    apply_3d_muscl!(phi, rhoptilde, nxx, nyy, nzz, limitertype)

    return
end

function reconstruct!(state::State, variable::U)
    (; limitertype) = state.namelists.discretization
    (; k0, k1, nxx, nyy, nzz) = state.domain
    (; rho, u) = state.variables.predictands
    (; phi) = state.variables.auxiliaries
    (; utilde) = state.variables.reconstructions
    (; rhostrattfc, pstrattfc) = state.atmosphere

    @ivy for kz in (k0 - 1):(k1 + 1), jy in 1:nyy, ix in 1:(nxx - 1)
        rhoedge =
            0.5 * (
                rho[ix, jy, kz] +
                rho[ix + 1, jy, kz] +
                rhostrattfc[ix, jy, kz] +
                rhostrattfc[ix + 1, jy, kz]
            )
        pedge = 0.5 * (pstrattfc[ix, jy, kz] + pstrattfc[ix + 1, jy, kz])
        phi[ix, jy, kz] = u[ix, jy, kz] * rhoedge / pedge
    end

    apply_3d_muscl!(phi, utilde, nxx, nyy, nzz, limitertype)

    return
end

function reconstruct!(state::State, variable::V)
    (; limitertype) = state.namelists.discretization
    (; k0, k1, nxx, nyy, nzz) = state.domain
    (; rho, v) = state.variables.predictands
    (; phi) = state.variables.auxiliaries
    (; vtilde) = state.variables.reconstructions
    (; rhostrattfc, pstrattfc) = state.atmosphere

    @ivy for kz in (k0 - 1):(k1 + 1), jy in 1:(nyy - 1), ix in 1:nxx
        rhoedge =
            0.5 * (
                rho[ix, jy, kz] +
                rho[ix, jy + 1, kz] +
                rhostrattfc[ix, jy, kz] +
                rhostrattfc[ix, jy + 1, kz]
            )
        pedge = 0.5 * (pstrattfc[ix, jy, kz] + pstrattfc[ix, jy + 1, kz])
        phi[ix, jy, kz] = v[ix, jy, kz] * rhoedge / pedge
    end

    apply_3d_muscl!(phi, vtilde, nxx, nyy, nzz, limitertype)

    return
end

function reconstruct!(state::State, variable::W)
    (; namelists, domain, grid) = state
    (; limitertype) = state.namelists.discretization
    (; i0, i1, j0, j1, k0, k1, nxx, nyy, nzz) = domain
    (; jac) = grid
    (; predictands) = state.variables
    (; rho, w) = predictands
    (; phi) = state.variables.auxiliaries
    (; wtilde) = state.variables.reconstructions
    (; rhostrattfc, pstrattfc) = state.atmosphere

    @. @ivy phi[:, :, (k0 - 1):(k1 + 1)] = w[:, :, (k0 - 1):(k1 + 1)]

    @ivy for kz in (k0 - 1):(k1 + 1), jy in j0:j1, ix in i0:i1
        phi[ix, jy, kz] = compute_vertical_wind(ix, jy, kz, state)
    end

    set_zonal_boundaries_of_field!(phi, namelists, domain)
    set_meridional_boundaries_of_field!(phi, namelists, domain)

    @ivy for kz in (k0 - 1):(k1 + 1), jy in 1:nyy, ix in 1:nxx
        rhoedgeu =
            (
                jac[ix, jy, kz + 1] *
                (rho[ix, jy, kz] + rhostrattfc[ix, jy, kz]) +
                jac[ix, jy, kz] *
                (rho[ix, jy, kz + 1] + rhostrattfc[ix, jy, kz + 1])
            ) / (jac[ix, jy, kz] + jac[ix, jy, kz + 1])
        pedgeu =
            (
                jac[ix, jy, kz + 1] * pstrattfc[ix, jy, kz] +
                jac[ix, jy, kz] * pstrattfc[ix, jy, kz + 1]
            ) / (jac[ix, jy, kz] + jac[ix, jy, kz + 1])
        phi[ix, jy, kz] *= rhoedgeu / pedgeu
    end

    apply_3d_muscl!(phi, wtilde, nxx, nyy, nzz, limitertype)

    return
end

function reconstruct!(state::State, tracersetup::NoTracer)
    return
end

function reconstruct!(state::State, tracersetup::AbstractTracer)
    (; limitertype) = state.namelists.discretization
    (; k0, k1, nxx, nyy, nzz) = state.domain
    (; phi) = state.variables.auxiliaries
    (; pstrattfc) = state.atmosphere
    (; tracerreconstructions, tracerpredictands) = state.tracer

    @ivy for (fd, field) in enumerate(fieldnames(TracerPredictands))
        for kz in (k0 - 1):(k1 + 1), jy in 1:nyy, ix in 1:nxx
            phi[ix, jy, kz] =
                getfield(tracerpredictands, fd)[ix, jy, kz] /
                pstrattfc[ix, jy, kz]
        end
        apply_3d_muscl!(
            phi,
            getfield(tracerreconstructions, fd),
            nxx,
            nyy,
            nzz,
            limitertype,
        )
    end

    return
end
