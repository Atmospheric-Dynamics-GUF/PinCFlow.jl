abstract type AbstractSpongeVariable end
struct SpongeDensity <: AbstractSpongeVariable end
struct SpongeDensityFluctuations <: AbstractSpongeVariable end
struct SpongeWind <: AbstractSpongeVariable end

function compute_sponge!(state::State, dt::AbstractFloat)
  (; nx, ny, nz) = state.domain
  (; ztfc, lz) = state.grid
  (; kr_sp_tfc, kr_sp_w_tfc, zsponge) = state.sponge
  (; unifiedsponge, spongetype, spongealphaz_fac) = state.namelists.sponge

  if unifiedsponge
    compute_unified_sponge!(state, dt, spongetype)
  else
    alpspg = spongealphaz_fac / dt

    for k in 1:nz
      for j in 1:ny
        for i in 1:nx
          if ztfc[i, j, k] >= zsponge
            kr_sp_tfc[i, j, k] =
              alpspg *
              sin(0.5 * pi * (ztfc[i, j, k] - zsponge) / (lz[1] - zsponge))^2.0
            kr_sp_w_tfc[i, j, k] = kr_sp_tfc[i, j, k] / jac[i, j, k]
          end
        end
      end
    end

    kr_sp_tfc[:, :, 0] = kr_sp_tfc[:, :, 1]
    kr_sp_tfc[:, :, d.nz + 1] = kr_sp_tfc[:, :, d.nz]
    kr_sp_w_tfc[:, :, 0] = kr_sp_w_tfc[:, :, 1]
    kr_sp_w_tfc[:, :, d.nz + 1] = kr_sp_w_tfc[:, :, d.nz]
  end

  return
end

function compute_unified_sponge!(
  state::State,
  dt::AbstractFloat,
  spongetype::ExponentialSponge,
)
  (; sizex, sizey, nbx, nby) = state.namelists.domain
  (; is, js, nx, ny, nz) = state.domain
  (; x, y, ztfc, lx, ly, lz) = state.grid
  (; tref) = state.constants
  (; lateralsponge, spongealphaz_dim) = state.namelists.sponge
  (; alphaunifiedsponge, dxsponge, dysponge, dzsponge) = state.sponge

  spongealphaz = spongealphaz_dim * tref

  if lateralsponge
    i00 = is + nbx - 1
    j00 = js + nby - 1
    if sizex > 1 && sizey > 1
      spongealphaz = spongealphaz / 3.0
      spongealphax = spongealphaz
      spongealphay = spongealphaz
    elseif sizex > 1
      spongealphaz = spongealphaz / 2.0
      spongealphax = spongealphaz
      spongealphay = 0.0
    elseif sizey > 1
      spongealphaz = spongealphaz / 2.0
      spongealphax = 0.0
      spongealphay = spongealphaz
    end
  end

  alphaunifiedsponge .= 0.0

  for k in 1:nz
    for j in 0:(ny + 1)
      for i in 0:(nx + 1)
        height = ztfc[i, j, k]

        if sizez > 1
          alphaunifiedsponge[i, j, k] =
            alphaunifiedsponge[i, j, k] +
            spongealphaz * exp((height - lz[1]) / dzsponge)
        end
        if lateralsponge
          if sizex > 1
            if x[i00 + i] <= 0.5 * (lx[0] + lx[1])
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                spongealphax * exp((lx[0] - x[i00 + i]) / dxsponge)
            else
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                spongealphax * exp((x[i00 + i] - lx[1]) / dxsponge)
            end
          end
          if sizey > 1
            if y[j00 + j] <= 0.5 * (ly[0] + ly[1])
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                spongealphay * exp((ly[0] - y[j00 + j]) / dysponge)
            else
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                spongealphay * exp((y[j00 + j] - ly[1]) / dysponge)
            end
          end
        end
      end
    end
  end

  alphaunifiedsponge[:, :, 0] = alphaunifiedsponge[:, :, 1]
  alphaunifiedsponge[:, :, d.nz + 1] = alphaunifiedsponge[:, :, d.nz]

  return
end

function compute_unified_sponge!(
  state::State,
  dt::AbstractFloat,
  spongetype::COSMOSponge,
)
  (; sizex, sizey, nbx, nby) = state.namelists.domain
  (; is, js, nx, ny, nz) = state.domain
  (; x, y, ztfc, lx, ly, lz) = state.grid
  (; tref) = state.constants
  (; lateralsponge, cosmosteps) = state.namelists.sponge
  (; alphaunifiedsponge, dxsponge, dysponge, dzsponge) = state.sponge

  spongealphaz = spongealphaz_dim * tref

  if lateralsponge
    i00 = is + nbx - 1
    j00 = js + nby - 1
    if sizex > 1 && sizey > 1
      spongealphaz = spongealphaz / 3.0
      spongealphax = spongealphaz
      spongealphay = spongealphaz
    elseif sizex > 1
      spongealphaz = spongealphaz / 2.0
      spongealphax = spongealphaz
      spongealphay = 0.0
    elseif sizey > 1
      spongealphaz = spongealphaz / 2.0
      spongealphax = 0.0
      spongealphay = spongealphaz
    end
  end

  alphaunifiedsponge .= 0.0

  for k in 1:nz
    for j in 0:(ny + 1)
      for i in 0:(nx + 1)
        height = ztfc[i, j, k]

        if sizez > 1
          if height >= zsponge
            alphaunifiedsponge[i, j, k] =
              alphaunifiedsponge[i, j, k] +
              0.5 / cosmosteps / dt *
              (1.0 - cos(pi * (height - zsponge) / dzsponge))
          end
        end
        if lateralsponge
          if sizex > 1
            if x[i00 + i] <= xsponge0
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                0.5 / cosmosteps / dt *
                (1.0 - cos(pi * (xsponge0 - x[i00 + i]) / dxsponge))
            elseif x[i00 + i] >= xsponge1
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                0.5 / cosmosteps / dt *
                (1.0 - cos(pi * (x[i00 + i] - xsponge1) / dxsponge))
            end
          end
          if sizey > 1
            if y[j00 + j] <= ysponge0
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                0.5 / cosmosteps / dt *
                (1.0 - cos(pi * (ysponge0 - y[j00 + j]) / dysponge))
            elseif y[j00 + j] >= ysponge1
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                0.5 / cosmosteps / dt *
                (1.0 - cos(pi * (y[j00 + j] - ysponge1) / dysponge))
            end
          end
        end
      end
    end
  end

  alphaunifiedsponge[:, :, 0] = alphaunifiedsponge[:, :, 1]
  alphaunifiedsponge[:, :, d.nz + 1] = alphaunifiedsponge[:, :, d.nz]

  return
end

function compute_unified_sponge!(
  state::State,
  dt::AbstractFloat,
  spongetype::PolynomialSponge,
)
  (; sizex, sizey, nbx, nby) = state.namelists.domain
  (; is, js, nx, ny, nz) = state.domain
  (; x, y, ztfc, lx, ly, lz) = state.grid
  (; tref) = state.constants
  (; lateralsponge, spongealphaz, spongeorder) = state.namelists.sponge
  (; alphaunifiedsponge, dxsponge, dysponge, dzsponge) = state.sponge

  spongealphaz = spongealphaz_dim * tref

  if lateralsponge
    i00 = is + nbx - 1
    j00 = js + nby - 1
    if sizex > 1 && sizey > 1
      spongealphaz = spongealphaz / 3.0
      spongealphax = spongealphaz
      spongealphay = spongealphaz
    elseif sizex > 1
      spongealphaz = spongealphaz / 2.0
      spongealphax = spongealphaz
      spongealphay = 0.0
    elseif sizey > 1
      spongealphaz = spongealphaz / 2.0
      spongealphax = 0.0
      spongealphay = spongealphaz
    end
  end

  alphaunifiedsponge .= 0.0

  for k in 1:nz
    for j in 0:(ny + 1)
      for i in 0:(nx + 1)
        height = ztfc[i, j, k]

        if sizez > 1
          if height >= zsponge
            alphaunifiedsponge[i, j, k] =
              alphaunifiedsponge[i, j, k] +
              spongealphaz * ((height - zsponge) / dzsponge)^spongeorder
          end
        end
        if lateralsponge
          if sizex > 1
            if x[i00 + i] <= xsponge0
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                spongealphax * ((xsponge0 - x[i00 + i]) / dxsponge)^spongeorder
            elseif x[i00 + i] >= xsponge1
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                spongealphax * ((x[i00 + i] - xsponge1) / dxsponge)^spongeorder
            end
          end
          if sizey > 1
            if y[j00 + j] <= ysponge0
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                spongealphay * ((ysponge0 - y[j00 + j]) / dysponge)^spongeorder
            elseif y[j00 + j] >= ysponge1
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                spongealphay * ((y[j00 + j] - ysponge1) / dysponge)^spongeorder
            end
          end
        end
      end
    end
  end

  alphaunifiedsponge[:, :, 0] = alphaunifiedsponge[:, :, 1]
  alphaunifiedsponge[:, :, d.nz + 1] = alphaunifiedsponge[:, :, d.nz]

  return
end

function compute_unified_sponge!(
  state::State,
  dt::AbstractFloat,
  spongetype::SinusoidalSponge,
)
  (; sizex, sizey, nbx, nby) = state.namelists.domain
  (; is, js, nx, ny, nz) = state.domain
  (; x, y, ztfc, lx, ly, lz) = state.grid
  (; tref) = state.constants
  (; lateralsponge, spongealphaz) = state.namelists.sponge
  (; alphaunifiedsponge, dxsponge, dysponge, dzsponge) = state.sponge

  spongealphaz = spongealphaz_dim * tref

  if lateralsponge
    i00 = is + nbx - 1
    j00 = js + nby - 1
    if sizex > 1 && sizey > 1
      spongealphaz = spongealphaz / 3.0
      spongealphax = spongealphaz
      spongealphay = spongealphaz
    elseif sizex > 1
      spongealphaz = spongealphaz / 2.0
      spongealphax = spongealphaz
      spongealphay = 0.0
    elseif sizey > 1
      spongealphaz = spongealphaz / 2.0
      spongealphax = 0.0
      spongealphay = spongealphaz
    end
  end

  alphaunifiedsponge .= 0.0

  for k in 1:nz
    for j in 0:(ny + 1)
      for i in 0:(nx + 1)
        height = ztfc[i, j, k]

        if sizez > 1
          if height >= zsponge
            alphaunifiedsponge[i, j, k] =
              alphaunifiedsponge[i, j, k] +
              spongealphaz * sin(0.5 * pi * (height - zsponge) / dzsponge)^2.0
          end
        end
        if lateralsponge
          if sizex > 1
            if x[i00 + i] <= xsponge0
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                spongealphax *
                sin(0.5 * pi * (xsponge0 - x[i00 + i]) / dxsponge)^2.0
            elseif x[i00 + i] >= xsponge1
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                spongealphax *
                sin(0.5 * pi * (x[i00 + i] - xsponge1) / dxsponge)^2.0
            end
          end
          if sizey > 1
            if y[j00 + j] <= ysponge0
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                spongealphay *
                sin(0.5 * pi * (ysponge0 - y[j00 + j]) / dysponge)^2.0
            elseif y[j00 + j] >= ysponge1
              alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                spongealphay *
                sin(0.5 * pi * (y[j00 + j] - ysponge1) / dysponge)^2.0
            end
          end
        end
      end
    end
  end

  alphaunifiedsponge[:, :, 0] = alphaunifiedsponge[:, :, 1]
  alphaunifiedsponge[:, :, d.nz + 1] = alphaunifiedsponge[:, :, d.nz]

  return
end

function apply_unified_sponge!(
  state::State,
  dt::AbstractFloat,
  time::AbstractFloat,
  variable::SpongeDensity,
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
  variable::SpongeDensityFluctuations,
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
  variable::SpongeWind,
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
  (; u, v, w) = state.variables.predictands
  (; jac) = state.grid

  (sum_local, sum_global) = (zeros(nz) for i in 1:2)

  # Apply the sponge to the zonal wind.

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

  # Apply the sponge to the meridional wind.

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

  # Apply the sponge to the vertical wind.

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
