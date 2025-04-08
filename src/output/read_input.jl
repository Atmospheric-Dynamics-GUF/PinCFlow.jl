function read_input!(state::State)

  # Get all necessary fields.
  (; iin, folder) = state.namelists.output
  (; comm, nx, ny, nz, io, jo, i0, i1, j0, j1, k0, k1) = state.domain
  (; tref, rhoref, uref) = state.constants
  (; rho, rhop, u, v, w, pip) = state.variables.predictands

  # Open the file. Note: Fused in-place assignments cannot be used here!
  time = h5open(folder * "/pincflow_input.h5", "r", comm) do file

    # Read the time.
    time = file["t"][iin] / tref

    # Read the density fluctuations.
    @views rhop[i0:i1, j0:j1, k0:k1] =
      file["rhop"][(io + 1):(io + nx), (jo + ny):(jo + j1), 1:nz, iin] ./ rhoref
    @views rho[i0:i1, j0:j1, k0:k1] .= rhop[i0:i1, j0:j1, k0:k1]

    # Read the staggered zonal wind.
    @views u[i0:i1, j0:j1, k0:k1] =
      file["us"][(io + 1):(io + nx), (jo + ny):(jo + j1), 1:nz, iin] ./ uref

    # Read the staggered meridional wind.
    @views v[i0:i1, j0:j1, k0:k1] =
      file["vs"][(io + 1):(io + nx), (jo + ny):(jo + j1), 1:nz, iin] ./ uref

    # Read the staggered transformed vertical wind.
    @views w[i0:i1, j0:j1, k0:k1] =
      file["wstfc"][(io + 1):(io + nx), (jo + ny):(jo + j1), 1:nz, iin] ./ uref

    # Read the Exner-pressure fluctuations.
    @views pip[i0:i1, j0:j1, k0:k1] =
      file["pip"][(io + 1):(io + nx), (jo + ny):(jo + j1), 1:nz, iin]

    # Return.
    return time
  end

  # Return.
  return time
end
