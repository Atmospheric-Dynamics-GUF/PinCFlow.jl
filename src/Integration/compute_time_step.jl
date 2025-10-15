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

  - CFL condition with respect to the resolved flow (where ``w`` is computed with `compute_vertical_wind`):

    ```math
    \\Delta t_\\mathrm{CFL} = \\mu_\\mathrm{CFL} \\min\\limits_\\mathrm{global} \\left[\\frac{\\Delta \\widehat{x}}{u_{\\max}}, \\frac{\\Delta \\widehat{y}}{v_{\\max}}, \\min \\left(\\frac{J \\Delta \\widehat{z}}{w}\\right)\\right]
    ```

  - CFL condition with respect to the group velocities of unresolved gravity waves (where ``J_{\\min}`` is the minimum Jacobian in a one-grid-cell radius and ``c_{\\mathrm{g} z}`` is the maximum vertical group velocity within a grid cell):

    ```math
    \\Delta t_\\mathrm{WKB} = \\mu_\\mathrm{WKB} \\min\\limits_\\mathrm{global} \\left[\\frac{\\Delta \\widehat{x}}{c_{\\mathrm{g} x, \\max}}, \\frac{\\Delta \\widehat{y}}{c_{\\mathrm{g} y, \\max}}, \\min \\left(\\frac{J_{\\min} \\Delta \\widehat{z}}{c_{\\mathrm{g} z}}\\right)\\right]
    ```

  - Von Neumann condition (with ``\\mathrm{Re}`` being the Reynolds number):

    ```math
    \\Delta t_\\mathrm{viscous} = \\frac{\\mathrm{Re}}{2} \\min\\limits_\\mathrm{global} \\left[\\left(\\Delta \\widehat{x}\\right)^2, \\left(\\Delta \\widehat{y}\\right)^2, \\left(J \\Delta \\widehat{z}\\right)^2\\right]
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
    (; dx, dy, dz, jac) = grid
    (; predictands) = state.variables
    (; u, v, w, rho) = predictands
    (; test_case) = state.namelists.setting
    (; x_size, y_size) = state.namelists.domain
    (; cgx_max, cgy_max, cgz_max) = state.wkb
    (; tke) = state.turbulence.turbulencepredictands
    (; rhobar) = state.atmosphere
    (; turbulence_scheme) = state.namelists.turbulence

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

        for k in k0:k1, j in j0:j1, i in i0:i1
            dtconv = min(
                dtconv,
                cfl_number * jac[i, j, k] * dz / (
                    abs(
                        0.5 * (
                            compute_vertical_wind(i, j, k, state) +
                            compute_vertical_wind(i, j, k - 1, state)
                        ),
                    ) + eps()
                ),
            )
        end
        dtconv = MPI.Allreduce(dtconv, min, comm)

        #---------------------------
        #   von Neumann condition
        #----------------------------

        dtvisc = 0.5 * min(dx^2, dy^2, dz^2) * re

        for k in k0:k1, j in j0:j1, i in i0:i1
            dtvisc = min(dtvisc, 0.5 * (jac[i, j, k] * dz)^2.0 * re)
        end
        dtvisc = MPI.Allreduce(dtvisc, min, comm)

        #----------------------------------
        #         WKB-CFL criterion
        #----------------------------------

        if typeof(test_case) <: AbstractWKBTestCase
            dtwkb = jac[i0, j0, k0] * dz / (cgz_max[i0, j0, k0] + eps())

            kmin = ko == 0 ? k0 - 1 : k0
            kmax = k1

            for k in kmin:kmax, j in j0:j1, i in i0:i1
                dtwkb = min(
                    dtwkb,
                    minimum(
                        jac[(i - 1):(i + 1), (j - 1):(j + 1), (k - 1):(k + 1)],
                    ) * dz / (cgz_max[i, j, k] + eps()),
                )
            end

            if x_size > 1
                dtwkb = min(dtwkb, dx / (cgx_max[] + eps()))
            end
            if y_size > 1
                dtwkb = min(dtwkb, dy / (cgy_max[] + eps()))
            end

            dtwkb *= wkb_cfl_number

            # find global minimum

            dtwkb = MPI.Allreduce(dtwkb, min, comm)
        end

        #-------------------------------
        #     Turbulence criterion 
        #---------------------

        if turbulence_scheme != NoTurbulence()
            uturb =
                maximum(
                    abs,
                    sqrt.(
                        tke[i0:i1, j0:j1, k0:k1] ./ (
                            rho[i0:i1, j0:j1, k0:k1] .+
                            rhobar[i0:i1, j0:j1, k0:k1]
                        )
                    ),
                ) + eps()

            dtturb = cfl_number * min(dx / uturb, dy / uturb, dz / uturb)

            dtturb = MPI.Allreduce(dtturb, min, comm)
        end
        #-------------------------------
        #        Make your choice
        #-------------------------------

        if typeof(test_case) <: AbstractWKBTestCase
            dt = min(dtvisc, dtconv, dtmax / tref, dtwkb)
        else
            dt = min(dtvisc, dtconv, dtmax / tref)
        end

        if turbulence_scheme != NoTurbulence()
            dt = min(dt, dtturb)
        end

        #-----------------------------------------
        #     Inform on time step restrictions
        #-----------------------------------------

        if master
            println("dtvisc = ", dtvisc * tref, " seconds")
            println("dtconv = ", dtconv * tref, " seconds")
            println("dtmax = ", dtmax, " seconds")
            if typeof(test_case) <: AbstractWKBTestCase
                println("dtwkb = ", dtwkb * tref, " seconds")
            end
            if turbulence_scheme != NoTurbulence()
                println("dtturb = ", dtturb * tref, " seconds")
            end
            println("")

            if dt == dtmax / tref
                println("=> dt = dtmax = ", dt * tref, " seconds")
            elseif dt == dtconv
                println("=> dt = dtconv = ", dt * tref, " seconds")
            elseif dt == dtvisc
                println("=> dt = dtvisc = ", dt * tref, " seconds")
            elseif dt == dtwkb
                println("=> dt = dtwkb = ", dt * tref, " seconds")
            elseif dt == dtturb
                println("=> dt = dtturb = ", dt * tref, " seconds")
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
