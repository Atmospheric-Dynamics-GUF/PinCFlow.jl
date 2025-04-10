function compute_time_step(state::State)
    (; grid) = state
    (; cfl, dtmin_dim, dtmax_dim, adaptive_time_step) =
        state.namelists.discretization
    (; tref, re) = state.constants
    (; master, comm, i0, i1, j0, j1, k0, k1) = state.domain
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

        @views umax = maximum(abs.(u[i0:i1, j0:j1, k0:k1])) + eps()
        @views vmax = maximum(abs.(v[i0:i1, j0:j1, k0:k1])) + eps()
        @views wmax = maximum(abs.(w[i0:i1, j0:j1, k0:k1])) + eps()

        dtconv_loc = cfl * min(dx / umax, dy / vmax, dz / wmax)

        for k in k0:k1, j in j0:j1, i in i0:i1
            dtconv_loc = min(
                dtconv_loc,
                cfl * jac[i, j, k] * dz / (
                    abs(
                        0.5 * (
                            compute_vertical_wind(i, j, k, predictands, grid) + compute_vertical_wind(
                                i,
                                j,
                                k - 1,
                                predictands,
                                grid,
                            )
                        ),
                    ) + eps()
                ),
            )
        end
        dtconv = MPI.Allreduce(dtconv_loc, min, comm)

        #---------------------------
        #   von Neumann condition
        #----------------------------

        dtvisc = 0.5 * min(dx^2, dy^2, dz^2) * re

        dtvisc_loc = dtvisc
        for k in k0:k1, j in j0:j1, i in i0:i1
            dtvisc_loc = min(dtvisc_loc, 0.5 * (jac[i, j, k] * dz)^2.0 * re)
        end
        dtvisc = MPI.Allreduce(dtvisc_loc, min, comm)

        #----------------------------
        #    Maximal time step
        #----------------------------

        dtmax = dtmax_dim / tref

        #---------------------------------
        #    Gravity wave time period
        #---------------------------------

        # dtwave = 1. / ()

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
