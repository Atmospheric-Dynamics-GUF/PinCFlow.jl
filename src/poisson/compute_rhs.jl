function compute_rhs!(
  state::State,
  b::AbstractArray{<:AbstractFloat, 3},
  dt::AbstractFloat,
  model::PseudoIncompressible,
)
  (; sizex, sizey, sizez) = state.namelists.domain
  (; ma, kappa) = state.constants
  (; master, comm, nx, ny, nz) = state.domain
  (; dx, dy, dz, jac) = state.grid
  (; rhostrattfc, pstrattfc) = state.atmosphere
  (; u, v, w) = state.variables.predictands

  # Initialize summation fields.
  divsum = 0.0
  divl2 = 0.0
  divmax = 0.0
  divsum_local = 0.0
  divl2_local = 0.0
  divl2_norm = 0.0
  divl2_norm_local = 0.0

  # Calculate RHS for TFC.
  for k in 1:nz
    for j in 1:ny
      for i in 1:nx
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
        # Compute RHS.
        bu = (pedger * ur - pedgel * ul) / dx / jac[i, j, k] * ma^2.0 * kappa
        bv = (pedgef * vf - pedgeb * vb) / dy / jac[i, j, k] * ma^2.0 * kappa
        bw = (pedgeu * wu - pedged * wd) / dz / jac[i, j, k] * ma^2.0 * kappa
        divsum_local = divsum_local + bu + bv + bw
        bu = bu / fcscal
        bv = bv / fcscal
        bw = bw / fcscal
        b[i, j, k] = bu + bv + bw
        # Compute check sum for solvability criterion.
        divl2_local = divl2_local + b[i, j, k]^2.0
        bl2loc = bu^2.0 + bv^2.0 + bw^2.0
        divl2_norm_local = divl2_norm_local + bl2loc
        if abs(b[i, j, k]) > divmax
          divmax = abs(b[i, j, k])
        end
      end
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

  b_norm = divl2

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
