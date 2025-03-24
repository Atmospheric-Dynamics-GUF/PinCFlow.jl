function compute_sponge!(state::State, dt::AbstractFloat)
  (; nx, ny, nz) = state.domain
  (; ztfc, lz, jac) = state.grid
  (; kr_sp_tfc, kr_sp_w_tfc, zsponge) = state.sponge
  (; unifiedsponge, spongetype, spongealphaz_fac) = state.namelists.sponge

  if unifiedsponge
    compute_sponge!(state, dt, spongetype)
  else
    alpspg = spongealphaz_fac / dt

    for k in 1:nz, j in 1:ny, i in 1:nx
      if ztfc[i, j, k] >= zsponge
        kr_sp_tfc[i, j, k] =
          alpspg *
          sin(0.5 * pi * (ztfc[i, j, k] - zsponge) / (lz[1] - zsponge))^2.0
        kr_sp_w_tfc[i, j, k] = kr_sp_tfc[i, j, k] / jac[i, j, k]
      end
    end

    @views kr_sp_tfc[:, :, 0] .= kr_sp_tfc[:, :, 1]
    @views kr_sp_tfc[:, :, nz + 1] .= kr_sp_tfc[:, :, nz]
    @views kr_sp_w_tfc[:, :, 0] .= kr_sp_w_tfc[:, :, 1]
    @views kr_sp_w_tfc[:, :, nz + 1] .= kr_sp_w_tfc[:, :, nz]
  end

  return
end

function compute_sponge!(
  state::State,
  dt::AbstractFloat,
  spongetype::ExponentialSponge,
)
  (; sizex, sizey, sizez, nbx, nby) = state.namelists.domain
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

  for k in 1:nz, j in 0:(ny + 1), i in 0:(nx + 1)
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

  @views alphaunifiedsponge[:, :, 0] .= alphaunifiedsponge[:, :, 1]
  @views alphaunifiedsponge[:, :, nz + 1] .= alphaunifiedsponge[:, :, nz]

  return
end

function compute_sponge!(
  state::State,
  dt::AbstractFloat,
  spongetype::COSMOSponge,
)
  (; sizex, sizey, sizez, nbx, nby) = state.namelists.domain
  (; is, js, nx, ny, nz) = state.domain
  (; x, y, ztfc, lx, ly, lz) = state.grid
  (; tref) = state.constants
  (; lateralsponge, cosmosteps) = state.namelists.sponge
  (;
    alphaunifiedsponge,
    dxsponge,
    dysponge,
    dzsponge,
    xsponge0,
    xsponge1,
    ysponge0,
    ysponge1,
    zsponge,
  ) = state.sponge

  if lateralsponge
    i00 = is + nbx - 1
    j00 = js + nby - 1
  end

  alphaunifiedsponge .= 0.0

  for k in 1:nz, j in 0:(ny + 1), i in 0:(nx + 1)
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

  @views alphaunifiedsponge[:, :, 0] .= alphaunifiedsponge[:, :, 1]
  @views alphaunifiedsponge[:, :, nz + 1] .= alphaunifiedsponge[:, :, nz]

  return
end

function compute_sponge!(
  state::State,
  dt::AbstractFloat,
  spongetype::PolynomialSponge,
)
  (; sizex, sizey, sizez, nbx, nby) = state.namelists.domain
  (; is, js, nx, ny, nz) = state.domain
  (; x, y, ztfc, lx, ly, lz) = state.grid
  (; tref) = state.constants
  (; lateralsponge, spongealphaz_dim, spongeorder) = state.namelists.sponge
  (;
    alphaunifiedsponge,
    dxsponge,
    dysponge,
    dzsponge,
    xsponge0,
    xsponge1,
    ysponge0,
    ysponge1,
    zsponge,
  ) = state.sponge

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

  for k in 1:nz, j in 0:(ny + 1), i in 0:(nx + 1)
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

  @views alphaunifiedsponge[:, :, 0] .= alphaunifiedsponge[:, :, 1]
  @views alphaunifiedsponge[:, :, nz + 1] .= alphaunifiedsponge[:, :, nz]

  return
end

function compute_sponge!(
  state::State,
  dt::AbstractFloat,
  spongetype::SinusoidalSponge,
)
  (; sizex, sizey, sizez, nbx, nby) = state.namelists.domain
  (; is, js, nx, ny, nz) = state.domain
  (; x, y, ztfc, lx, ly, lz) = state.grid
  (; tref) = state.constants
  (; lateralsponge, spongealphaz_dim) = state.namelists.sponge
  (;
    alphaunifiedsponge,
    dxsponge,
    dysponge,
    dzsponge,
    xsponge0,
    xsponge1,
    ysponge0,
    ysponge1,
    zsponge,
  ) = state.sponge

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

  for k in 1:nz, j in 0:(ny + 1), i in 0:(nx + 1)
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

  @views alphaunifiedsponge[:, :, 0] .= alphaunifiedsponge[:, :, 1]
  @views alphaunifiedsponge[:, :, nz + 1] .= alphaunifiedsponge[:, :, nz]

  return
end
