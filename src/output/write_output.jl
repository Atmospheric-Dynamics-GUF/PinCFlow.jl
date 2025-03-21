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
    dataset["x"][:] = view(x, 1:sizex) .* lref
    dataset["y"][:] = view(y, 1:sizey) .* lref
  end

  # Write the vertical grid.
  local_array[:, :, :] = view(ztfc, 1:nx, 1:ny, 1:nz) .* lref
  compute_global_array!(namelists, domain)
  if master
    dataset["z"][:, :, :] = global_array
  end

  # Write the background density.
  local_array[:, :, :] = view(rhostrattfc, 1:nx, 1:ny, 1:nz) .* rhoref
  compute_global_array!(namelists, domain)
  if master
    dataset["rhobar"][:, :, :] = global_array
  end

  # Write the background potential temperature.
  local_array[:, :, :] = view(thetastrattfc, 1:nx, 1:ny, 1:nz) .* thetaref
  compute_global_array!(namelists, domain)
  if master
    dataset["thetabar"][:, :, :] = global_array
  end

  # Write the buoyancy frequency.
  local_array[:, :, :] = view(bvsstrattfc, 1:nx, 1:ny, 1:nz) ./ tref .^ 2
  compute_global_array!(namelists, domain)
  if master
    dataset["n2"][:, :, :] = global_array
  end

  # Write the mass-weighted potential temperature.
  local_array[:, :, :] = view(pstrattfc, 1:nx, 1:ny, 1:nz) ./ tref .^ 2
  compute_global_array!(namelists, domain)
  if master
    dataset["p"][:, :, :] = global_array
  end

  # Write the density fluctuations.
  if prepare_restart || RhoP() in atmvarout
    local_array[:, :, :] = view(rho, 1:nx, 1:ny, 1:nz) .* rhoref
    compute_global_array!(namelists, domain)
    if master
      dataset["rhop"][:, :, :, iout] = global_array
    end
  end

  # Write the zonal winds.
  if U() in atmvarout
    local_array[:, :, :] =
      (view(u, 1:nx, 1:ny, 1:nz) .+ view(u, 0:(nx - 1), 1:ny, 1:nz)) ./ 2 .*
      uref
    compute_global_array!(namelists, domain)
    if master
      dataset["u"][:, :, :, iout] = global_array
    end
  end

  # Write the staggered zonal winds.
  if prepare_restart || US() in atmvarout
    local_array[:, :, :] = view(u, 1:nx, 1:ny, 1:nz) .* uref
    compute_global_array!(namelists, domain)
    if master
      dataset["us"][:, :, :, iout] = global_array
    end
  end

  # Write the meridional winds.
  if V() in atmvarout
    local_array[:, :, :] =
      (view(v, 1:nx, 1:ny, 1:nz) .+ view(v, 1:nx, 0:(ny - 1), 1:nz)) ./ 2 .*
      uref
    compute_global_array!(namelists, domain)
    if master
      dataset["v"][:, :, :, iout] = global_array
    end
  end

  # Write the staggered meridional winds.
  if prepare_restart || VS() in atmvarout
    local_array[:, :, :] = view(v, 1:nx, 1:ny, 1:nz) .* uref
    compute_global_array!(namelists, domain)
    if master
      dataset["vs"][:, :, :, iout] = global_array
    end
  end

  # Write the vertical winds.
  if W() in atmvarout
    for k in 1:nz
      for j in 1:ny
        for i in 1:nx
          local_array[i, j, k] =
            (
              compute_vertical_wind(i, j, k, predictands, grid) +
              compute_vertical_wind(i, j, k - 1, predictands, grid)
            ) / 2 * uref
        end
      end
    end
    compute_global_array!(namelists, domain)
    if master
      dataset["w"][:, :, :, iout] = global_array
    end
  end

  # Write the staggered vertical winds.
  if WS() in atmvarout
    for k in 1:nz
      for j in 1:ny
        for i in 1:nx
          local_array[i, j, k] =
            compute_vertical_wind(i, j, k, predictands, grid) * uref
        end
      end
    end
    compute_global_array!(namelists, domain)
    if master
      dataset["ws"][:, :, :, iout] = global_array
    end
  end

  # Write the transformed vertical winds.
  if WTFC() in atmvarout
    local_array[:, :, :] =
      (view(w, 1:nx, 1:ny, 1:nz) .+ view(w, 1:nx, 1:ny, 0:(nz - 1))) ./ 2 .*
      uref
    compute_global_array!(namelists, domain)
    if master
      dataset["wtfc"][:, :, :, iout] = global_array
    end
  end

  # Write the staggered transformed vertical winds.
  if prepare_restart || WSTFC() in atmvarout
    local_array[:, :, :] = view(w, 1:nx, 1:ny, 1:nz) .* uref
    compute_global_array!(namelists, domain)
    if master
      dataset["wstfc"][:, :, :, iout] = global_array
    end
  end

  # Write the potential-temperature fluctuations.
  if ThetaP() in atmvarout
    local_array =
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
    local_array[:, :, :] = view(pip, 1:nx, 1:ny, 1:nz)
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
