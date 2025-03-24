function write_output(state::State, time::AbstractFloat, iout::Integer)

  # Get all necessary fields.
  (; namelists, domain, grid) = state
  (; sizex, sizey, sizez) = namelists.domain
  (; prepare_restart, atmvarout) = namelists.output
  (; comm, master, nx, ny, nz, local_array, global_array) = domain
  (; tref, lref, rhoref, thetaref, uref) = state.constants
  (; x, y, ztfc) = grid
  (; rhostrattfc, thetastrattfc, bvsstrattfc, pstrattfc) = state.atmosphere
  (; predictands) = state.variables
  (; rho, u, v, w, pip) = predictands

  # Advance output counter.
  iout += 1

  # Open the dataset.
  if master
    dataset = NCDataset("pincflow_output.nc", "a")
  end

  # Write the time.
  if master
    dataset["t"][iout] = time * tref
  end

  # Write the horizontal grid.
  if master
    @views dataset["x"][:] = x[1:sizex] .* lref
    @views dataset["y"][:] = y[1:sizey] .* lref
  end

  # Write the vertical grid.
  @views local_array[:, :, :] = ztfc[1:nx, 1:ny, 1:nz] .* lref
  compute_global_array!(namelists, domain)
  if master
    dataset["z"][:, :, :] = global_array
  end

  # Write the background density.
  @views local_array[:, :, :] = rhostrattfc[1:nx, 1:ny, 1:nz] .* rhoref
  compute_global_array!(namelists, domain)
  if master
    dataset["rhobar"][:, :, :] = global_array
  end

  # Write the background potential temperature.
  @views local_array[:, :, :] = thetastrattfc[1:nx, 1:ny, 1:nz] .* thetaref
  compute_global_array!(namelists, domain)
  if master
    dataset["thetabar"][:, :, :] = global_array
  end

  # Write the buoyancy frequency.
  @views local_array[:, :, :] = bvsstrattfc[1:nx, 1:ny, 1:nz] ./ tref .^ 2
  compute_global_array!(namelists, domain)
  if master
    dataset["n2"][:, :, :] = global_array
  end

  # Write the mass-weighted potential temperature.
  @views local_array[:, :, :] = pstrattfc[1:nx, 1:ny, 1:nz] ./ tref .^ 2
  compute_global_array!(namelists, domain)
  if master
    dataset["p"][:, :, :] = global_array
  end

  # Write the density fluctuations.
  if prepare_restart || RhoP() in atmvarout
    @views local_array[:, :, :] = rho[1:nx, 1:ny, 1:nz] .* rhoref
    compute_global_array!(namelists, domain)
    if master
      dataset["rhop"][:, :, :, iout] = global_array
    end
  end

  # Write the zonal winds.
  if U() in atmvarout
    @views local_array[:, :, :] =
      (u[1:nx, 1:ny, 1:nz] .+ u[0:(nx - 1), 1:ny, 1:nz]) ./ 2 .* uref
    compute_global_array!(namelists, domain)
    if master
      dataset["u"][:, :, :, iout] = global_array
    end
  end

  # Write the staggered zonal winds.
  if prepare_restart || US() in atmvarout
    @views local_array[:, :, :] = u[1:nx, 1:ny, 1:nz] .* uref
    compute_global_array!(namelists, domain)
    if master
      dataset["us"][:, :, :, iout] = global_array
    end
  end

  # Write the meridional winds.
  if V() in atmvarout
    @views local_array[:, :, :] =
      (v[1:nx, 1:ny, 1:nz] .+ v[1:nx, 0:(ny - 1), 1:nz]) ./ 2 .* uref
    compute_global_array!(namelists, domain)
    if master
      dataset["v"][:, :, :, iout] = global_array
    end
  end

  # Write the staggered meridional winds.
  if prepare_restart || VS() in atmvarout
    @views local_array[:, :, :] = v[1:nx, 1:ny, 1:nz] .* uref
    compute_global_array!(namelists, domain)
    if master
      dataset["vs"][:, :, :, iout] = global_array
    end
  end

  # Write the vertical winds.
  if W() in atmvarout
    for k in 1:nz, j in 1:ny, i in 1:nx
      local_array[i, j, k] =
        (
          compute_vertical_wind(i, j, k, predictands, grid) +
          compute_vertical_wind(i, j, k - 1, predictands, grid)
        ) / 2 * uref
    end
    compute_global_array!(namelists, domain)
    if master
      dataset["w"][:, :, :, iout] = global_array
    end
  end

  # Write the staggered vertical winds.
  if WS() in atmvarout
    for k in 1:nz, j in 1:ny, i in 1:nx
      local_array[i, j, k] =
        compute_vertical_wind(i, j, k, predictands, grid) * uref
    end
    compute_global_array!(namelists, domain)
    if master
      dataset["ws"][:, :, :, iout] = global_array
    end
  end

  # Write the transformed vertical winds.
  if WTFC() in atmvarout
    @views local_array[:, :, :] =
      (w[1:nx, 1:ny, 1:nz] .+ w[1:nx, 1:ny, 0:(nz - 1)]) ./ 2 .* uref
    compute_global_array!(namelists, domain)
    if master
      dataset["wtfc"][:, :, :, iout] = global_array
    end
  end

  # Write the staggered transformed vertical winds.
  if prepare_restart || WSTFC() in atmvarout
    @views local_array[:, :, :] = w[1:nx, 1:ny, 1:nz] .* uref
    compute_global_array!(namelists, domain)
    if master
      dataset["wstfc"][:, :, :, iout] = global_array
    end
  end

  # Write the potential-temperature fluctuations.
  if ThetaP() in atmvarout
    @views local_array[:, :, :] =
      (
        pstrattfc[1:nx, 1:ny, 1:nz] ./
        (rhostrattfc[1:nx, 1:ny, 1:nz] .+ rho[1:nx, 1:ny, 1:nz]) .-
        thetastrattfc[1:nx, 1:ny, 1:nz]
      ) .* thetaref
    if master
      dataset["thetap"][:, :, :, iout] = global_array
    end
  end

  # Write the Exner-pressure fluctuations.
  if prepare_restart || PiP() in atmvarout
    @views local_array[:, :, :] = pip[1:nx, 1:ny, 1:nz]
    compute_global_array!(namelists, domain)
    if master
      dataset["pip"][:, :, :, iout] = global_array
    end
  end

  # Close the dataset.
  if master
    close(dataset)
  end

  # Return.
  return iout
end
