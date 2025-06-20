"""
    compute_rhs!(state, b, model::Boussinesq)

Compute the right-hand side of the Poisson equation for Boussinesq approximation.

For the Boussinesq model, the RHS is computed from the velocity divergence scaled by
stratification parameters. The function also computes convergence tolerance based on
the L2 norm of the divergence field.

# Arguments

  - `state::State`: Simulation state containing velocity fields and grid information
  - `b::AbstractArray{<:AbstractFloat, 3}`: Output array for RHS values
  - `model::Boussinesq`: Boussinesq model type dispatch

# Returns

  - `tolref::AbstractFloat`: Reference tolerance for iterative solver convergence
"""
function compute_rhs!(
    state::State,
    b::AbstractArray{<:AbstractFloat, 3},
    model::Boussinesq,
)
    (; sizex, sizey, sizez) = state.namelists.domain
    (; ma, kappa) = state.constants
    (; comm, nx, ny, nz, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; thetastrattfc) = state.atmosphere
    (; u, v, w) = state.variables.predictands

    # Initialize summation fields.
    divmax = 0.0
    divsum = 0.0
    divsum_local = 0.0
    divl2 = 0.0
    divl2_local = 0.0
    divl2_norm = 0.0
    divl2_norm_local = 0.0

    # Calculate RHS for TFC.
    for k in k0:k1, j in j0:j1, i in i0:i1
        # Store velocities at cell edges.
        ur = (jac[i, j, k] + jac[i + 1, j, k]) / 2 * u[i, j, k]
        ul = (jac[i, j, k] + jac[i - 1, j, k]) / 2 * u[i - 1, j, k]
        vf = (jac[i, j, k] + jac[i, j + 1, k]) / 2 * v[i, j, k]
        vb = (jac[i, j, k] + jac[i, j - 1, k]) / 2 * v[i, j - 1, k]
        wu =
            2 * jac[i, j, k] * jac[i, j, k + 1] /
            (jac[i, j, k] + jac[i, j, k + 1]) * w[i, j, k]
        wd =
            2 * jac[i, j, k] * jac[i, j, k - 1] /
            (jac[i, j, k] + jac[i, j, k - 1]) * w[i, j, k - 1]
        # Determine indices for RHS.
        ib = i - i0 + 1
        jb = j - j0 + 1
        kb = k - k0 + 1
        # Compute RHS.
        bu =
            (ur - ul) / dx / jac[i, j, k] * ma^2 * kappa /
            thetastrattfc[i, j, k]
        bv =
            (vf - vb) / dy / jac[i, j, k] * ma^2 * kappa /
            thetastrattfc[i, j, k]
        bw =
            (wu - wd) / dz / jac[i, j, k] * ma^2 * kappa /
            thetastrattfc[i, j, k]
        b[ib, jb, kb] = bu + bv + bw
        bl2loc = bu^2 + bv^2 + bw^2
        # Compute check sum for solvability criterion.
        divsum_local += b[ib, jb, kb]
        divl2_local += b[ib, jb, kb]^2
        divl2_norm_local += bl2loc
        if abs(b[ib, jb, kb]) > divmax
            divmax = abs(b[ib, jb, kb])
        end
    end

    # MPI: sum divSum_local over all procs
    divsum = MPI.Allreduce(divsum_local, +, comm)

    # MPI: sum divL2_local over all procs
    divl2 = MPI.Allreduce(divl2_local, +, comm)

    # MPI: sum divL2_norm_local over all procs
    divl2_norm = MPI.Allreduce(divl2_norm_local, +, comm)

    # scale div
    divl2_local = sqrt(divl2_local / nx / ny / nz)
    divl2 = sqrt(divl2 / sizex / sizey / sizez)

    divl2_norm_local = sqrt(divl2_norm_local / nx / ny / nz)
    divl2_norm = sqrt(divl2_norm / sizex / sizey / sizez)

    if divl2_norm != 0.0
        tolref = divl2 / divl2_norm
    else
        if divl2 == 0.0
            tolref = 1.0
        else
            error("Error in compute_rhs: divl2_norm = 0 while divl2 != 0!")
        end
    end

    return tolref
end

"""
    compute_rhs!(state, b, model::PseudoIncompressible)

Compute the right-hand side of the Poisson equation for the pseudo-incompressible approximation.

This version accounts for density variations by including pressure-weighted velocity
divergence terms.

# Arguments

  - `state::State`: Simulation state
  - `b::AbstractArray{<:AbstractFloat, 3}`: Output array for RHS values
  - `model::PseudoIncompressible`: Model type for dispatch

# Returns

  - `tolref::AbstractFloat`: Reference tolerance for iterative solver convergence
"""
function compute_rhs!(
    state::State,
    b::AbstractArray{<:AbstractFloat, 3},
    model::PseudoIncompressible,
)
    (; sizex, sizey, sizez) = state.namelists.domain
    (; ma, kappa) = state.constants
    (; comm, nx, ny, nz, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; u, v, w) = state.variables.predictands

    # Initialize summation fields.
    divmax = 0.0
    divsum = 0.0
    divsum_local = 0.0
    divl2 = 0.0
    divl2_local = 0.0
    divl2_norm = 0.0
    divl2_norm_local = 0.0

    # Calculate RHS for TFC.
    for k in k0:k1, j in j0:j1, i in i0:i1
        # Calculate scaling factor.
        fcscal = sqrt(pstrattfc[i, j, k]^2.0 / rhostrattfc[i, j, k])
        # Store velocities at cell edges.
        ur = u[i, j, k]
        ul = u[i - 1, j, k]
        vf = v[i, j, k]
        vb = v[i, j - 1, k]
        wu = w[i, j, k]
        wd = w[i, j, k - 1]
        # Calculate P at cell edges.
        pedger =
            0.5 * (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
            )
        pedgel =
            0.5 * (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i - 1, j, k] * pstrattfc[i - 1, j, k]
            )
        pedgef =
            0.5 * (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
            )
        pedgeb =
            0.5 * (
                jac[i, j, k] * pstrattfc[i, j, k] +
                jac[i, j - 1, k] * pstrattfc[i, j - 1, k]
            )
        pedgeu =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1])
        pedged =
            jac[i, j, k] *
            jac[i, j, k - 1] *
            (pstrattfc[i, j, k] + pstrattfc[i, j, k - 1]) /
            (jac[i, j, k] + jac[i, j, k - 1])
        # Determine indices for RHS.
        ib = i - i0 + 1
        jb = j - j0 + 1
        kb = k - k0 + 1
        # Compute RHS.
        bu = (pedger * ur - pedgel * ul) / dx / jac[i, j, k] * ma^2.0 * kappa
        bv = (pedgef * vf - pedgeb * vb) / dy / jac[i, j, k] * ma^2.0 * kappa
        bw = (pedgeu * wu - pedged * wd) / dz / jac[i, j, k] * ma^2.0 * kappa
        divsum_local += bu + bv + bw
        bu /= fcscal
        bv /= fcscal
        bw /= fcscal
        b[ib, jb, kb] = bu + bv + bw
        # Compute check sum for solvability criterion.
        divl2_local += b[ib, jb, kb]^2.0
        bl2loc = bu^2.0 + bv^2.0 + bw^2.0
        divl2_norm_local += bl2loc
        if abs(b[ib, jb, kb]) > divmax
            divmax = abs(b[ib, jb, kb])
        end
    end

    # MPI: sum divSum_local over all procs
    divsum = MPI.Allreduce(divsum_local, +, comm)

    # MPI: sum divL2_local over all procs
    divl2 = MPI.Allreduce(divl2_local, +, comm)

    # MPI: sum divL2_norm_local over all procs
    divl2_norm = MPI.Allreduce(divl2_norm_local, +, comm)

    # scale div
    divl2_local = sqrt(divl2_local / nx / ny / nz)
    divl2 = sqrt(divl2 / sizex / sizey / sizez)

    divl2_norm_local = sqrt(divl2_norm_local / nx / ny / nz)
    divl2_norm = sqrt(divl2_norm / sizex / sizey / sizez)

    if divl2_norm != 0.0
        tolref = divl2 / divl2_norm
    else
        if divl2 == 0.0
            tolref = 1.0
        else
            error("Error in compute_rhs: divl2_norm = 0 while divl2 != 0!")
        end
    end

    return tolref
end

"""
    compute_rhs!(state, b, model::Compressible)

Compute the RHS for the Poisson equation for the fully compressible equations including heating terms.

The compressible formulation includes additional source terms from diabatic heating
computed through the volume force calculation.

# Arguments

  - `state::State`: Simulation state
  - `b::AbstractArray{<:AbstractFloat, 3}`: Output RHS array
  - `model::Compressible`: Compressible model type

# Returns

  - `tolref::AbstractFloat`: Reference tolerance for convergence
"""
function compute_rhs!(
    state::State,
    b::AbstractArray{<:AbstractFloat, 3},
    model::Compressible,
)
    (; sizex, sizey, sizez) = state.namelists.domain
    (; ma, kappa) = state.constants
    (; comm, nx, ny, nz, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; u, v, w) = state.variables.predictands

    # Initialize summation fields.
    divmax = 0.0
    divsum = 0.0
    divsum_local = 0.0
    divl2 = 0.0
    divl2_local = 0.0
    divl2_norm = 0.0
    divl2_norm_local = 0.0

    # Calculate RHS for TFC.
    for k in k0:k1, j in j0:j1, i in i0:i1
        # Calculate scaling factor.
        fcscal = sqrt(pstrattfc[i, j, k]^2.0 / rhostrattfc[i, j, k])
        # Store velocities at cell edges.
        ur = u[i, j, k]
        ul = u[i - 1, j, k]
        vf = v[i, j, k]
        vb = v[i, j - 1, k]
        wu = w[i, j, k]
        wd = w[i, j, k - 1]
        # Determine indices for RHS.
        ib = i - i0 + 1
        jb = j - j0 + 1
        kb = k - k0 + 1
        # Compute the heating.
        heating = compute_volume_force(state, (i, j, k), P()) * ma^2.0 * kappa
        # Compute RHS.
        bu = (ur - ul) / dx / jac[i, j, k] * ma^2.0 * kappa
        bv = (vf - vb) / dy / jac[i, j, k] * ma^2.0 * kappa
        bw = (wu - wd) / dz / jac[i, j, k] * ma^2.0 * kappa
        divsum_local += bu + bv + bw + heating
        bu /= fcscal
        bv /= fcscal
        bw /= fcscal
        heating /= fcscal
        b[ib, jb, kb] = bu + bv + bw + heating
        # Compute check sum for solvability criterion.
        divl2_local += b[ib, jb, kb]^2.0
        bl2loc = bu^2.0 + bv^2.0 + bw^2.0 + heating^2.0
        divl2_norm_local += bl2loc
        if abs(b[ib, jb, kb]) > divmax
            divmax = abs(b[ib, jb, kb])
        end
    end

    # MPI: sum divSum_local over all procs
    divsum = MPI.Allreduce(divsum_local, +, comm)

    # MPI: sum divL2_local over all procs
    divl2 = MPI.Allreduce(divl2_local, +, comm)

    # MPI: sum divL2_norm_local over all procs
    divl2_norm = MPI.Allreduce(divl2_norm_local, +, comm)

    # scale div
    divl2_local = sqrt(divl2_local / nx / ny / nz)
    divl2 = sqrt(divl2 / sizex / sizey / sizez)

    divl2_norm_local = sqrt(divl2_norm_local / nx / ny / nz)
    divl2_norm = sqrt(divl2_norm / sizex / sizey / sizez)

    if divl2_norm != 0.0
        tolref = divl2 / divl2_norm
    else
        if divl2 == 0.0
            tolref = 1.0
        else
            error("Error in compute_rhs: divl2_norm = 0 while divl2 != 0!")
        end
    end

    return tolref
end
