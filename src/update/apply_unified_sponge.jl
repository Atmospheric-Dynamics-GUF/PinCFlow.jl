function apply_unified_sponge!(
  state::State,
  dt::AbstractFloat,
  time::AbstractFloat,
  variable::Rho,
  model::PseudoIncompressible,
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
  (; backgroundflow_dim) = state.namelists.atmosphere
  (;
    spongelayer,
    unifiedsponge,
    relax_to_mean,
    relaxation_period,
    relaxation_amplitude,
  ) = state.namelists.sponge
  (; uref, tref) = state.constants
  (; comm, i0, i1, j0, j1, k0, k1, local_sum, global_sum) = state.domain
  (; alphaunifiedsponge) = state.sponge
  (; u) = state.variables.predictands

  if !spongelayer || !unifiedsponge
    return
  end

  # Determine relaxation wind.
  if relax_to_mean
    for k in k0:k1
      @views local_sum[k - k0 + 1] = sum(u[i0:i1, j0:j1, k])
    end
    MPI.Allreduce!(local_sum, global_sum, .+, comm)
    global_sum ./= (sizex .* sizey)
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
  for k in k0:k1
    if relax_to_mean
      ubg = global_sum[k - k0 + 1]
    end
    for j in j0:j1, i in i0:i1
      alpha =
        0.5 * (alphaunifiedsponge[i, j, k] + alphaunifiedsponge[i + 1, j, k])
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
  (; backgroundflow_dim) = state.namelists.atmosphere
  (;
    spongelayer,
    unifiedsponge,
    relax_to_mean,
    relaxation_period,
    relaxation_amplitude,
  ) = state.namelists.sponge
  (; uref, tref) = state.constants
  (; comm, i0, i1, j0, j1, k0, k1, local_sum, global_sum) = state.domain
  (; alphaunifiedsponge) = state.sponge
  (; v) = state.variables.predictands

  if !spongelayer || !unifiedsponge
    return
  end

  # Determine relaxation wind.
  if relax_to_mean
    for k in k0:k1
      @views local_sum[k - k0 + 1] = sum(v[i0:i1, j0:j1, k])
    end
    MPI.Allreduce!(local_sum, global_sum, .+, comm)
    global_sum ./= (sizex .* sizey)
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
  for k in k0:k1
    if relax_to_mean
      vbg = global_sum[k - k0 + 1]
    end
    for j in j0:j1, i in i0:i1
      alpha =
        0.5 * (alphaunifiedsponge[i, j, k] + alphaunifiedsponge[i, j + 1, k])
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
  (; backgroundflow_dim) = state.namelists.atmosphere
  (;
    spongelayer,
    unifiedsponge,
    relax_to_mean,
    relaxation_period,
    relaxation_amplitude,
  ) = state.namelists.sponge
  (; uref, tref) = state.constants
  (; comm, i0, i1, j0, j1, k0, k1, local_sum, global_sum) = state.domain
  (; alphaunifiedsponge) = state.sponge
  (; w) = state.variables.predictands
  (; jac) = state.grid

  if !spongelayer || !unifiedsponge
    return
  end

  # Determine relaxation wind.
  if relax_to_mean
    for k in k0:k1
      @views local_sum[k - k0 + 1] = sum(w[i0:i1, j0:j1, k])
    end
    MPI.Allreduce!(local_sum, global_sum, .+, comm)
    global_sum ./= (sizex .* sizey)
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
  for k in k0:k1
    if relax_to_mean
      wbg = global_sum[k - k0 + 1]
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
