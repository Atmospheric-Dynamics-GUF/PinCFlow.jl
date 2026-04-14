"""
```julia
compute_time_step(state::State)::AbstractFloat
```

Compute and return an adaptive time step based on several stability criteria.

If `state.namelists.discretization.adaptive_time_step` is set to `true`, the returned time step is given by

```math
\\Delta t = \\min \\left(\\Delta t_\\mathrm{CFL}, \\Delta t_\\mathrm{WKB}, \\Delta t_\\mathrm{viscous}, \\Delta t_{\\max}\\right),
```

where ``\\Delta t_\\mathrm{CFL}`` and ``\\Delta t_\\mathrm{WKB}`` are computed from CFL conditions with respect to the resolved flow and unresolved gravity waves, respectively, ``\\Delta t_\\mathrm{viscous}`` is determined from a von Neumann condition that takes the viscosity into account and ``\\Delta t_{\\max}`` is an upper limit specified in `state.namelists.discretization`. Otherwise, the returned time step is equal to ``\\Delta t_{\\max}``. If ``\\Delta t`` is smaller than ``\\Delta t_{\\min}`` (also specified in `state.namelists.discretization`), an error is thrown.

The individual stability criteria are as follows.

  - CFL condition with respect to the resolved flow (where ``u_{\\max}``, ``v_{\\max}``, and ``\\hat{w}_{\\max}`` are the maximum absolute velocities in zonal, meridional and transformed vertical direction, respectively):

    ```math
    \\Delta t_\\mathrm{CFL} = \\mu_\\mathrm{CFL} \\min \\left(\\frac{\\Delta \\hat{x}}{u_{\\max}}, \\frac{\\Delta \\hat{y}}{v_{\\max}}, \\frac{\\Delta \\hat{z}}{\\hat{w}_{\\max}}\\right)
    ```

  - CFL condition with respect to the group velocities of unresolved gravity waves (where ``c_{\\mathrm{g}, x, \\max}``, ``c_{\\mathrm{g}, y, \\max}``, and ``c_{\\mathrm{g}, z, \\max}`` are the maximum absolute group velocities in zonal, meridional, and vertical direction, respectively, and ``\\Delta z_{\\min}`` is the minimum layer depth):

    ```math
    \\Delta t_\\mathrm{WKB} = \\mu_\\mathrm{WKB} \\min \\left(\\frac{\\Delta \\hat{x}}{c_{\\mathrm{g}, x, \\max}}, \\frac{\\Delta \\hat{y}}{c_{\\mathrm{g}, y, \\max}}, \\frac{\\Delta z_{\\min}}{c_{\\mathrm{g}, z, \\max}}\\right)
    ```

  - Von Neumann condition (with ``\\mathrm{Re}`` being the Reynolds number):

    ```math
    \\Delta t_\\mathrm{viscous} = \\frac{\\mathrm{Re}}{2} \\min \\left[\\left(\\Delta \\hat{x}\\right)^2, \\left(\\Delta \\hat{y}\\right)^2, \\left(\\Delta z_{\\min}\\right)^2\\right]
    ```

# Arguments

  - `state`: Model state.

# See also

  - [`PinCFlow.Update.compute_vertical_wind`](@ref)
"""
function compute_time_step end

function compute_time_step(state::State)::AbstractFloat
    (; grid) = state
    (; cfl_number, wkb_cfl_number, dtmin, dtmax, adaptive_time_step) =
        state.namelists.discretization
    (; tref, re) = state.constants
    (; master, comm, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, dzcmin) = grid
    (; predictands) = state.variables
    (; u, v, w) = predictands
    (; x_size, y_size) = state.namelists.domain
    (; wkb_mode) = state.namelists.wkb
    (; cgx_max, cgy_max, cgz_max) = state.wkb

    @ivy if !adaptive_time_step
        dt = dtmax / tref

        if master
            println("dt = dtfix = ", dt * tref, " seconds")
            println("")
        end
    else
        #----------------------
        #     CFL condition
        #----------------------

        umax = maximum(abs, u[i0:i1, j0:j1, k0:k1]) + eps()
        vmax = maximum(abs, v[i0:i1, j0:j1, k0:k1]) + eps()
        wmax = maximum(abs, w[i0:i1, j0:j1, k0:k1]) + eps()

        dtconv = cfl_number * min(dx / umax, dy / vmax, dz / wmax)

        dtconv = MPI.Allreduce(dtconv, min, comm)

        #---------------------------
        #   von Neumann condition
        #----------------------------

        dtvisc = 0.5 * min(dx^2, dy^2, dzcmin^2) * re

        #----------------------------------
        #         WKB-CFL criterion
        #----------------------------------

        if wkb_mode in (:SingleColumn, :MultiColumn)
            dtwkb = dzcmin / (cgz_max[] + eps())

            if x_size > 1
                dtwkb = min(dtwkb, dx / (cgx_max[] + eps()))
            end

            if y_size > 1
                dtwkb = min(dtwkb, dy / (cgy_max[] + eps()))
            end

            dtwkb *= wkb_cfl_number

            dtwkb = MPI.Allreduce(dtwkb, min, comm)
        end

        #-------------------------------
        #        Make your choice
        #-------------------------------

        if wkb_mode in (:SingleColumn, :MultiColumn)
            dt = min(dtvisc, dtconv, dtmax / tref, dtwkb)
        else
            dt = min(dtvisc, dtconv, dtmax / tref)
        end

        #-----------------------------------------
        #     Inform on time step restrictions
        #-----------------------------------------

        if master
            println("dtvisc = ", dtvisc * tref, " seconds")
            println("dtconv = ", dtconv * tref, " seconds")
            println("dtmax = ", dtmax, " seconds")
            if wkb_mode in (:SingleColumn, :MultiColumn)
                println("dtwkb = ", dtwkb * tref, " seconds")
            end
            println("")

            if dt == dtmax / tref
                println("=> dt = dtmax = ", dt * tref, " seconds")
            elseif dt == dtconv
                println("=> dt = dtconv = ", dt * tref, " seconds")
            elseif dt == dtvisc
                println("=> dt = dtvisc = ", dt * tref, " seconds")
            elseif wkb_mode in (:SingleColumn, :MultiColumn) && dt == dtwkb
                println("=> dt = dtwkb = ", dt * tref, " seconds")
            else
                println("=> dt = ??? = ", dt * tref, " seconds")
            end
            println("")
        end
    end

    if dt * tref < dtmin
        error(
            "Error in compute_time_step: dt = ",
            dt * tref,
            " < ",
            dtmin,
            " = dtmin",
        )
    end

    return dt
end
