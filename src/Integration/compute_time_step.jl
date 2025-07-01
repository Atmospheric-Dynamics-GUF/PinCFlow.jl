"""
```julia
compute_time_step(state::State) -> Float64
```

Compute adaptive time step based on stability constraints.

Calculates the maximum allowable time step by evaluating multiple stability
criteria including CFL condition, viscous stability, and wave propagation
constraints. Returns the most restrictive time step.

# Arguments

  - `state::State`: Complete simulation state containing grid, flow fields, and parameters

# Returns

  - `Float64`: Time step in non-dimensional units

# Stability Criteria

## CFL Condition (Advective Stability)

  - Horizontal: `dt ≤ cfl * min(dx/|u|, dy/|v|)`
  - Vertical: `dt ≤ cfl * jac * dz / |w_physical|`
  - Uses terrain-following coordinate transformation for vertical velocity

## von Neumann Condition (Viscous Stability)

  - `dt ≤ 0.5 * Re * min(dx², dy², (jac*dz)²)`
  - Ensures stability of viscous diffusion terms

## WKB-CFL Condition (for gravity wave test cases)

  - `dt ≤ cfl_wave * min(dx/cgx, dy/cgy, jac*dz/cgz)`
  - Based on group velocity components from WKB ray tracing
  - Only applied for `AbstractWKBTestCase` types

## Maximum Time Step

  - User-specified upper limit: `dt ≤ dtmax`

# Adaptive vs Fixed Time Step

  - If `adaptive_time_step = false`: returns fixed `dtmax`
  - If `adaptive_time_step = true`: returns minimum of all constraints

# MPI Communication

  - Uses `MPI.Allreduce` with `min` operation to find global minimum
  - Ensures all processes use the same time step

# Output

  - Master process prints all constraint values and selected time step
  - Identifies which constraint is most restrictive

# Error Conditions

  - Throws error if computed time step falls below `dtmin`
  - Prevents simulation from becoming unstable due to overly small time steps
"""
function compute_time_step(state::State)
    (; grid) = state
    (; cfl, cfl_wave, dtmin_dim, dtmax_dim, adaptive_time_step) =
        state.namelists.discretization
    (; tref, re) = state.constants
    (; master, comm, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = grid
    (; predictands) = state.variables
    (; u, v, w) = predictands
    (; testcase) = state.namelists.setting
    (; sizex, sizey) = state.namelists.domain
    (; cgx_max, cgy_max, cgz_max) = state.wkb

    #-------------------------------------------
    #              Fixed time step
    #-------------------------------------------

    if !adaptive_time_step
        dt = dtmax_dim / tref

        if master
            println("dt = dtfix = ", dt * tref, " seconds")
            println("")
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

        #----------------------------------
        #         WKB-CFL criterion
        #----------------------------------

        if typeof(testcase) <: AbstractWKBTestCase
            dtwkb_loc = jac[i0, j0, k0] * dz / (cgz_max[i0, j0, k0] + eps())

            kz0 = ko == 0 ? k0 - 1 : k0
            kz1 = k1

            for kz in kz0:kz1, jy in j0:j1, ix in i0:i1
                @views dtwkb_loc = min(
                    dtwkb_loc,
                    minimum(
                        jac[
                            (ix - 1):(ix + 1),
                            (jy - 1):(jy + 1),
                            (kz - 1):(kz + 1),
                        ],
                    ) * dz / (cgz_max[ix, jy, kz] + eps()),
                )
            end

            if sizex > 1
                dtwkb_loc = min(dtwkb_loc, dx / (cgx_max[] + eps()))
            end
            if sizey > 1
                dtwkb_loc = min(dtwkb_loc, dy / (cgy_max[] + eps()))
            end

            dtwkb_loc *= cfl_wave

            # find global minimum

            dtwkb = MPI.Allreduce(dtwkb_loc, min, comm)
        end
        #-------------------------------
        #        Make your choice
        #-------------------------------

        if typeof(testcase) <: AbstractWKBTestCase
            dt = min(dtvisc, dtconv, dtmax, dtwkb)
        else
            dt = min(dtvisc, dtconv, dtmax)
        end

        #-----------------------------------------
        #     Inform on time step restrictions
        #-----------------------------------------

        if master
            println("dtvisc = ", dtvisc * tref, " seconds")
            println("dtconv = ", dtconv * tref, " seconds")
            println("dtmax = ", dtmax * tref, " seconds")
            if typeof(testcase) <: AbstractWKBTestCase
                println("dtwkb = ", dtwkb * tref, " seconds")
            end
            println("")

            if dt == dtmax
                println("=> dt = dtmax = ", dt * tref, " seconds")
            elseif dt == dtconv
                println("=> dt = dtconv = ", dt * tref, " seconds")
            elseif dt == dtvisc
                println("=> dt = dtvisc = ", dt * tref, " seconds")
            elseif dt == dtwkb
                println("=> dt = dtwkb = ", dt * tref, " seconds")
            else
                println("=> dt = ??? = ", dt * tref, " seconds")
            end
            println("")
        end
    end

    if dt * tref < dtmin_dim
        error(
            "Error in compute_time_step: dt = ",
            dt * tref,
            " < ",
            dtmin_dim,
            " = dtmin",
        )
    end

    return dt
end
