function write_output(
  state::State,
  time::AbstractFloat,
  iout::Integer,
  cpu_start_time::DateTime,
)

  # Get all necessary fields.
  (; namelists, domain, grid) = state
  (; sizex, sizey) = namelists.domain
  (; prepare_restart, atmvarout) = namelists.output
  (; master, nx, ny, nz, i0, i1, j0, j1, k0, k1, local_array, global_array) =
    domain
  (; tref, lref, rhoref, thetaref, uref) = state.constants
  (; x, y, ztfc) = grid
  (; rhostrattfc, thetastrattfc, bvsstrattfc, pstrattfc) = state.atmosphere
  (; predictands) = state.variables
  (; rho, u, v, w, pip) = predictands

  # Print information.
  if master
    println("")
    println(repeat("-", 80))
    println("Output into file pincflow_output.nc")
    println("Physical time: ", time * tref, " s")
    println("CPU time: ", canonicalize(now() - cpu_start_time))
    println(repeat("-", 80))
  end

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
  if master && iout == 1
    @views dataset["x"][:] .= x[i0:(i0 + sizex - 1)] .* lref
    @views dataset["y"][:] .= y[j0:(j0 + sizey - 1)] .* lref
  end

  # Write the vertical grid.
  if iout == 1
    @views local_array[:, :, :] .= ztfc[i0:i1, j0:j1, k0:k1] .* lref
    compute_global_array!(namelists, domain)
    if master
      dataset["z"][:, :, :] .= global_array
    end
  end

  # Write the background density.
  if iout == 1
    @views local_array[:, :, :] .= rhostrattfc[i0:i1, j0:j1, k0:k1] .* rhoref
    compute_global_array!(namelists, domain)
    if master
      dataset["rhobar"][:, :, :] .= global_array
    end
  end

  # Write the background potential temperature.
  if iout == 1
    @views local_array[:, :, :] .=
      thetastrattfc[i0:i1, j0:j1, k0:k1] .* thetaref
    compute_global_array!(namelists, domain)
    if master
      dataset["thetabar"][:, :, :] .= global_array
    end
  end

  # Write the buoyancy frequency.
  if iout == 1
    @views local_array[:, :, :] .= bvsstrattfc[i0:i1, j0:j1, k0:k1] ./ tref .^ 2
    compute_global_array!(namelists, domain)
    if master
      dataset["n2"][:, :, :] .= global_array
    end
  end

  # Write the mass-weighted potential temperature.
  if iout == 1
    @views local_array[:, :, :] .= pstrattfc[i0:i1, j0:j1, k0:k1] ./ tref .^ 2
    compute_global_array!(namelists, domain)
    if master
      dataset["p"][:, :, :] .= global_array
    end
  end

  # Write the density fluctuations.
  if prepare_restart || RhoP() in atmvarout
    @views local_array[:, :, :] .= rho[i0:i1, j0:j1, k0:k1] .* rhoref
    compute_global_array!(namelists, domain)
    if master
      dataset["rhop"][:, :, :, iout] .= global_array
    end
  end

  # Write the zonal winds.
  if U() in atmvarout
    @views local_array[:, :, :] .=
      (u[i0:i1, j0:j1, k0:k1] .+ u[(i0 - 1):(i1 - 1), j0:j1, k0:k1]) ./ 2 .*
      uref
    compute_global_array!(namelists, domain)
    if master
      dataset["u"][:, :, :, iout] .= global_array
    end
  end

  # Write the staggered zonal winds.
  if prepare_restart || US() in atmvarout
    @views local_array[:, :, :] .= u[i0:i1, j0:j1, k0:k1] .* uref
    compute_global_array!(namelists, domain)
    if master
      dataset["us"][:, :, :, iout] .= global_array
    end
  end

  # Write the meridional winds.
  if V() in atmvarout
    @views local_array[:, :, :] .=
      (v[i0:i1, j0:j1, k0:k1] .+ v[i0:i1, (j0 - 1):(j1 - 1), k0:k1]) ./ 2 .*
      uref
    compute_global_array!(namelists, domain)
    if master
      dataset["v"][:, :, :, iout] .= global_array
    end
  end

  # Write the staggered meridional winds.
  if prepare_restart || VS() in atmvarout
    @views local_array[:, :, :] .= v[i0:i1, j0:j1, k0:k1] .* uref
    compute_global_array!(namelists, domain)
    if master
      dataset["vs"][:, :, :, iout] .= global_array
    end
  end

  # Write the vertical winds.
  if W() in atmvarout
    for k in 1:nz, j in 1:ny, i in 1:nx
      local_array[i, j, k] =
        (
          compute_vertical_wind(
            i + i0 - 1,
            j + j0 - 1,
            k + k0 - 1,
            predictands,
            grid,
          ) + compute_vertical_wind(
            i + i0 - 1,
            j + j0 - 1,
            k + k0 - 2,
            predictands,
            grid,
          )
        ) / 2 * uref
    end
    compute_global_array!(namelists, domain)
    if master
      dataset["w"][:, :, :, iout] .= global_array
    end
  end

  # Write the staggered vertical winds.
  if WS() in atmvarout
    for k in 1:nz, j in 1:ny, i in 1:nx
      local_array[i, j, k] =
        compute_vertical_wind(
          i + i0 - 1,
          j + j0 - 1,
          k + k0 - 1,
          predictands,
          grid,
        ) * uref
    end
    compute_global_array!(namelists, domain)
    if master
      dataset["ws"][:, :, :, iout] .= global_array
    end
  end

  # Write the transformed vertical winds.
  if WTFC() in atmvarout
    @views local_array[:, :, :] .=
      (w[i0:i1, j0:j1, k0:k1] .+ w[i0:i1, j0:j1, (k0 - 1):(k1 - 1)]) ./ 2 .*
      uref
    compute_global_array!(namelists, domain)
    if master
      dataset["wtfc"][:, :, :, iout] .= global_array
    end
  end

  # Write the staggered transformed vertical winds.
  if prepare_restart || WSTFC() in atmvarout
    @views local_array[:, :, :] .= w[i0:i1, j0:j1, k0:k1] .* uref
    compute_global_array!(namelists, domain)
    if master
      dataset["wstfc"][:, :, :, iout] .= global_array
    end
  end

  # Write the potential-temperature fluctuations.
  if ThetaP() in atmvarout
    @views local_array[:, :, :] .=
      (
        pstrattfc[i0:i1, j0:j1, k0:k1] ./
        (rhostrattfc[i0:i1, j0:j1, k0:k1] .+ rho[i0:i1, j0:j1, k0:k1]) .-
        thetastrattfc[i0:i1, j0:j1, k0:k1]
      ) .* thetaref
    if master
      dataset["thetap"][:, :, :, iout] .= global_array
    end
  end

  # Write the Exner-pressure fluctuations.
  if prepare_restart || PiP() in atmvarout
    @views local_array[:, :, :] .= pip[i0:i1, j0:j1, k0:k1]
    compute_global_array!(namelists, domain)
    if master
      dataset["pip"][:, :, :, iout] .= global_array
    end
  end

  # Close the dataset.
  if master
    close(dataset)
  end

  # Return.
  return iout
end
