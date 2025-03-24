function read_input!(state::State)

  # Get all necessary fields.
  (; namelists, domain) = state
  (; iin) = namelists.output
  (; comm, local_array, global_array) = domain
  (; tref, lref, rhoref, thetaref, uref) = state.constants
  (; rho, rhop, u, v, w, pip) = predictands

  # Open the dataset.
  if master
    dataset = NCDataset("pincflow_input.nc", "r")
  end

  # Read the time.
  if master
    time = dataset["t"][iin] / tref
  end
  time = MPI.bcast(time, comm)

  # Read the density fluctuations.
  if master
    @views global_array[:, :, :] .= dataset["rhop"][:, :, :, iin] ./ rhoref
  end
  compute_local_array!(namelists, domain)
  rho[1:nx, 1:ny, 1:nz] .= local_array
  rhop[1:nx, 1:ny, 1:nz] .= local_array

  # Read the staggered zonal winds.
  if master
    @views global_array[:, :, :] .= dataset["us"][:, :, :, iin] ./ uref
  end
  compute_local_array!(namelists, domain)
  u[1:nx, 1:ny, 1:nz] .= local_array

  # Read the staggered meridional winds.
  if master
    @views global_array[:, :, :] .= dataset["vs"][:, :, :, iin] ./ uref
  end
  compute_local_array!(namelists, domain)
  v[1:nx, 1:ny, 1:nz] .= local_array

  # Read the staggered transformed vertical winds.
  if master
    @views global_array[:, :, :] .= dataset["wstfc"][:, :, :, iin] ./ uref
  end
  compute_local_array!(namelists, domain)
  w[1:nx, 1:ny, 1:nz] .= local_array

  # Read the Exner-pressure fluctuations.
  if master
    @views global_array[:, :, :] .= dataset["pip"][:, :, :, iin]
  end
  compute_local_array!(namelists, domain)
  pip[1:nx, 1:ny, 1:nz] .= local_array

  # Close the dataset.
  if master
    close(dataset)
  end

  # Return.
  return time
end
