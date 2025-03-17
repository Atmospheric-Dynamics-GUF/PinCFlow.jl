function apply_unified_sponge!(
  state::State,
  dt::AbstractFloat,
  time::AbstractFloat,
  variable::Rho,
  model::PseudoIncompressible,
)
  (; spongelayer, unifiedsponge) = state.namelists.sponge
  (; nx, ny, nz) = state.domain
  (; alphaunifiedsponge) = state.sponge
  (; rho) = state.variables.predictands

  if !spongelayer || !unifiedsponge
    return
  end

  rho_bg = 0.0
  for k in 1:nz
    for j in 1:ny
      for i in 1:nx
        alpha = alphaunifiedsponge[i, j, k]
        rho_old = rho[i, j, k]
        beta = 1.0 / (1.0 + alpha * dt)
        rho_new = (1.0 - beta) * rho_bg + beta * rho_old
        rho[i, j, k] = rho_new
      end
    end
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
  (; nx, ny, nz) = state.domain
  (; alphaunifiedsponge) = state.sponge
  (; rhop) = state.variables.predictands

  if !spongelayer || !unifiedsponge
    return
  end

  rho_bg = 0.0
  for k in 1:nz
    for j in 1:ny
      for i in 1:nx
        alpha = alphaunifiedsponge[i, j, k]
        rho_old = rhop[i, j, k]
        beta = 1.0 / (1.0 + alpha * dt)
        rho_new = (1.0 - beta) * rho_bg + beta * rho_old
        rhop[i, j, k] = rho_new
      end
    end
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
  (; backgroundflow_dim) = state.namelists.atmosphere
  (;
    spongelayer,
    unifiedsponge,
    relax_to_mean,
    relaxation_period,
    relaxation_amplitude,
  ) = state.namelists.sponge
  (; uref, tref) = state.constants
  (; comm, nx, ny, nz) = state.domain
  (; alphaunifiedsponge) = state.sponge
  (; u) = state.variables.predictands

  if !spongelayer || !unifiedsponge
    return
  end

  (sum_local, sum_global) = (zeros(nz) for i in 1:2)

  # Determine relaxation wind.
  if relax_to_mean
    for k in 1:nz
      sum_local[k] = sum(u[1:nx, 1:ny, k])
    end
    sum_global[:] = MPI.Allreduce(sum_local, sum, comm)
    sum_global[:] = sum_global ./ (sizex .* sizey)
  else
    ubg = backgroundflow_dim[1] / uref
    if relaxation_period > 0.0
      ubg =
        ubg * (
          1.0 +
          relaxation_amplitude *
          sin(2.0 * pi * time / relaxation_period * tref)
        )
    end
  end

  # Update the zonal wind.
  for k in 1:nz
    if relax_to_mean
      ubg = sum_global[k]
    end
    for j in 1:ny
      for i in 1:nx
        alpha =
          0.5 * (alphaunifiedsponge[i, j, k] + alphaunifiedsponge[i + 1, j, k])
        uold = u[i, j, k]
        beta = 1.0 / (1.0 + alpha * dt)
        unew = (1.0 - beta) * ubg + beta * uold
        u[i, j, k] = unew
      end
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
  (; backgroundflow_dim) = state.namelists.atmosphere
  (;
    spongelayer,
    unifiedsponge,
    relax_to_mean,
    relaxation_period,
    relaxation_amplitude,
  ) = state.namelists.sponge
  (; uref, tref) = state.constants
  (; comm, nx, ny, nz) = state.domain
  (; alphaunifiedsponge) = state.sponge
  (; v) = state.variables.predictands

  if !spongelayer || !unifiedsponge
    return
  end

  (sum_local, sum_global) = (zeros(nz) for i in 1:2)

  # Determine relaxation wind.
  if relax_to_mean
    for k in 1:nz
      sum_local[k] = sum(v[1:nx, 1:ny, k])
    end
    sum_global[:] = MPI.Allreduce(sum_local, sum, comm)
    sum_global[:] = sum_global ./ (sizex .* sizey)
  else
    vbg = backgroundflow_dim[2] / uref
    if relaxation_period > 0.0
      vbg =
        vbg * (
          1.0 +
          relaxation_amplitude *
          sin(2.0 * pi * time / relaxation_period * tref)
        )
    end
  end

  # Update the meridional wind.
  for k in 1:nz
    if relax_to_mean
      vbg = sum_global[k]
    end
    for j in 1:ny
      for i in 1:nx
        alpha =
          0.5 * (alphaunifiedsponge[i, j, k] + alphaunifiedsponge[i, j + 1, k])
        vold = v[i, j, k]
        beta = 1.0 / (1.0 + alpha * dt)
        vnew = (1.0 - beta) * vbg + beta * vold
        v[i, j, k] = vnew
      end
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
  (; backgroundflow_dim) = state.namelists.atmosphere
  (;
    spongelayer,
    unifiedsponge,
    relax_to_mean,
    relaxation_period,
    relaxation_amplitude,
  ) = state.namelists.sponge
  (; uref, tref) = state.constants
  (; comm, nx, ny, nz) = state.domain
  (; alphaunifiedsponge) = state.sponge
  (; w) = state.variables.predictands
  (; jac) = state.grid

  if !spongelayer || !unifiedsponge
    return
  end

  (sum_local, sum_global) = (zeros(nz) for i in 1:2)

  # Determine relaxation wind.
  if relax_to_mean
    for k in 1:nz
      sum_local[k] = sum(w[1:nx, 1:ny, k])
    end
    sum_global[:] = MPI.Allreduce(sum_local, sum, comm)
    sum_global[:] = sum_global ./ (sizex .* sizey)
  else
    wbg = backgroundflow_dim[3] / uref
    if relaxation_period > 0.0
      wbg =
        wbg * (
          1.0 +
          relaxation_amplitude *
          sin(2.0 * pi * time / relaxation_period * tref)
        )
    end
  end

  # Update the vertical wind.
  for k in 1:nz
    if relax_to_mean
      wbg = sum_global[k]
    end
    for j in 1:ny
      for i in 1:nx
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
  end

  return
end
