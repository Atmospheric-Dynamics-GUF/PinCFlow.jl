"""
    reconstruct!(state::State)
    reconstruct!(state::State, variable::Union{Rho,RhoP})
    reconstruct!(state::State, variable::Union{U,V,W})

Perform MUSCL reconstruction for state variables in a computational domain.

These functions handle the reconstruction of various flow quantities using the MUSCL scheme.
The reconstruction is performed separately for each variable.

# Arguments

  - `state::State`: The complete state object containing all simulation variables
  - `variable`: Type parameter indicating which variable to reconstruct (Rho, RhoP, U, V, or W)

# Details

  - For density (Rho) and density perturbation (RhoP): Scales by reference pressure
  - For velocities (U,V,W): Includes density averaging at appropriate cell faces and pressure scaling
  - W-component additionally requires special vertical wind computation and boundary handling

The functions modify the corresponding 'tilde' variables in `state.variables.reconstructions`.
"""
function reconstruct!(state::State)
    (; tracersetup) = state.namelists.tracer
    (; icesetup) = state.namelists.ice
    (; turbulencesetup) = state.namelists.turbulence

    reconstruct!(state, Rho())
    reconstruct!(state, RhoP())
    reconstruct!(state, U())
    reconstruct!(state, V())
    reconstruct!(state, W())

    reconstruct!(state, tracersetup)
    reconstruct!(state, icesetup)
    reconstruct!(state, turbulencesetup)

    return
end

# Density reconstruction methods share implementation pattern
"""
    reconstruct!(state::State, variable::Rho)
    reconstruct!(state::State, variable::RhoP)

Reconstruct density or density perturbation fields.
Scales the field by the reference pressure before applying MUSCL reconstruction.
"""
function reconstruct!(state::State, variable::Rho)
    (; limitertype) = state.namelists.discretization
    (; k0, k1, nxx, nyy, nzz) = state.domain
    (; rho) = state.variables.predictands
    (; phi) = state.variables.auxiliaries
    (; rhotilde) = state.variables.reconstructions
    (; pstrattfc) = state.atmosphere

    for kz in (k0 - 1):(k1 + 1), jy in 1:nyy, ix in 1:nxx
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

    for kz in (k0 - 1):(k1 + 1), jy in 1:nyy, ix in 1:nxx
        phi[ix, jy, kz] = rhop[ix, jy, kz] / pstrattfc[ix, jy, kz]
    end
    apply_3d_muscl!(phi, rhoptilde, nxx, nyy, nzz, limitertype)

    return
end

"""
    reconstruct!(state::State, variable::U)
    reconstruct!(state::State, variable::V)
    reconstruct!(state::State, variable::W)

Reconstruct velocity components (u,v,w).
Includes:

  - Density averaging at appropriate cell faces
  - Pressure scaling
  - For W: Additional vertical wind computation and boundary condition handling
"""
function reconstruct!(state::State, variable::U)
    (; limitertype) = state.namelists.discretization
    (; k0, k1, nxx, nyy, nzz) = state.domain
    (; rho, u) = state.variables.predictands
    (; phi) = state.variables.auxiliaries
    (; utilde) = state.variables.reconstructions
    (; rhostrattfc, pstrattfc) = state.atmosphere

    for kz in (k0 - 1):(k1 + 1), jy in 1:nyy, ix in 1:(nxx - 1)
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

    for kz in (k0 - 1):(k1 + 1), jy in 1:(nyy - 1), ix in 1:nxx
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

    @views phi[:, :, (k0 - 1):(k1 + 1)] .= w[:, :, (k0 - 1):(k1 + 1)]
    for kz in (k0 - 1):(k1 + 1), jy in j0:j1, ix in i0:i1
        phi[ix, jy, kz] = compute_vertical_wind(ix, jy, kz, predictands, grid)
    end
    set_zonal_boundaries_of_field!(phi, namelists, domain)
    set_meridional_boundaries_of_field!(phi, namelists, domain)
    for kz in (k0 - 1):(k1 + 1), jy in 1:nyy, ix in 1:nxx
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

    for (fd, field) in enumerate(fieldnames(TracerPredictands))
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

function reconstruct!(state::State, icesetup::NoIce)
    return
end

function reconstruct!(state::State, icesetup::AbstractIce)
    (; limitertype) = state.namelists.discretization
    (; k0, k1, nxx, nyy, nzz) = state.domain
    (; phi) = state.variables.auxiliaries
    (; pstrattfc) = state.atmosphere
    (; icereconstructions, icepredictands) = state.ice

    for (fd, field) in enumerate(fieldnames(IcePredictands))
        for kz in (k0 - 1):(k1 + 1), jy in 1:nyy, ix in 1:nxx
            phi[ix, jy, kz] =
                getfield(icepredictands, fd)[ix, jy, kz] / pstrattfc[ix, jy, kz]
        end
        apply_3d_muscl!(
            phi,
            getfield(icereconstructions, fd),
            nxx,
            nyy,
            nzz,
            limitertype,
        )
    end

    return
end

function reconstruct!(state::State, turbulencesetup::NoTurbulence)
    return
end

function reconstruct!(state::State, turbulencesetup::AbstractTurbulence)
    (; limitertype) = state.namelists.discretization
    (; k0, k1, nxx, nyy, nzz) = state.domain
    (; phi) = state.variables.auxiliaries
    (; pstrattfc) = state.atmosphere
    (; turbulencereconstructions, turbulencepredictands) = state.turbulence

    for (fd, field) in enumerate(fieldnames(TurbulencePredictands))
        for kz in (k0 - 1):(k1 + 1), jy in 1:nyy, ix in 1:nxx
            phi[ix, jy, kz] =
                getfield(turbulencepredictands, fd)[ix, jy, kz] /
                pstrattfc[ix, jy, kz]
        end
        apply_3d_muscl!(
            phi,
            getfield(turbulencereconstructions, fd),
            nxx,
            nyy,
            nzz,
            limitertype,
        )
    end

    return
end
