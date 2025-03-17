function compute_time_step(state::State)
  (; grid) = state
  (; cfl, dtmin_dim, dtmax_dim, adaptive_time_step) =
    state.namelists.discretization
  (; tref, small, re) = state.constants
  (; master, comm, root, nx, ny, nz) = state.domain
  (; dx, dy, dz, jac) = grid
  (; predictands) = state.variables
  (; u, v, w) = predictands

  #-------------------------------------------
  #              Fixed time step
  #-------------------------------------------

  if !adaptive_time_step
    dt = dtmax_dim / tref

    if master
      println("dt = dtfix = ", dt * tref, " seconds")
    end
  else

    #-------------------------------------------
    #           Variable time step
    #-------------------------------------------

    #----------------------
    #     CFL condition
    #----------------------

    umax = maximum(abs.(u[1:nx, 1:ny, 1:nz])) + small
    vmax = maximum(abs.(v[1:nx, 1:ny, 1:nz])) + small
    wmax = maximum(abs.(w[1:nx, 1:ny, 1:nz])) + small

    dtconv_loc = cfl * min(dx / umax, dy / vmax, dz / wmax)

    for k in 1:nz
      for j in 1:ny
        for i in 1:nx
          dtconv_loc = min(
            dtconv_loc,
            cfl * jac[i, j, k] * dz / (
              abs(
                0.5 * (
                  compute_vertical_wind(i, j, k, predictands, grid) +
                  compute_vertical_wind(i, j, k - 1, predictands, grid)
                ),
              ) + small
            ),
          )
        end
      end
    end
    dtconv = MPI.Reduce(dtconv_loc, min, comm; root = root)
    MPI.bcast(dtconv, comm; root = root)

    #---------------------------
    #   von Neumann condition
    #----------------------------

    dtvisc = 0.5 * min(dx^2, dy^2, dz^2) * re

    dtvisc_loc = dtvisc
    for k in 1:(nz)
      for j in 1:(ny)
        for i in 1:(nx)
          dtvisc_loc = min(dtvisc_loc, 0.5 * (jac[i, j, k] * dz)^2.0 * re)
        end
      end
    end
    dtvisc = MPI.Reduce(dtvisc_loc, min, comm; root = root)
    MPI.bcast(dtvisc, comm; root = root)

    #----------------------------
    #    Maximal time step
    #----------------------------

    dtmax = dtmax_dim / tref

    #-------------------------------
    #        Make your choice
    #-------------------------------

    dt = min(dtvisc, dtconv, dtmax)

    #-----------------------------------------
    #     Inform on time step restrictions
    #-----------------------------------------

    if master
      println("dtvisc = ", dtvisc * tref, " seconds")
      println("dtconv = ", dtconv * tref, " seconds")
      println("dtmax = ", dtmax * tref, " seconds")
      println("")

      if dt == dtmax
        println("--> dt = dtmax = ", dt * tref, " seconds")
      elseif dt == dtconv
        println("--> dt = dtconv = ", dt * tref, " seconds")
      elseif dt == dtvisc
        println("--> dt = dtvisc = ", dt * tref, " seconds")
      else
        println("--> dt = ????? = ", dt * tref, " seconds")
      end
      println("")
    end
  end

  if dt * tref < dtmin_dim
    println(
      "Error in compute_time_step: dt = ",
      dt * tref,
      " < ",
      dtmin_dim,
      " = dtmin",
    )
    exit()
  end

  return dt
end
