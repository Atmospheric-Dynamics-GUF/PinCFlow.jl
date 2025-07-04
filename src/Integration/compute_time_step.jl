"""
```julia
compute_time_step(state::State)
```

Compute adaptive time step based on several stability criteria.

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

  - CFL condition with respect to the (grid-cell) maxima of the group velocities of unresolved gravity waves (with ``J_{\\min}`` being the minimum Jacobian in a one-grid-cell radius):

```math
\\Delta t_\\mathrm{WKB} = \\mu_\\mathrm{WKB} \\min\\limits_\\mathrm{global} \\left[\\frac{\\Delta \\widehat{x}}{c_{\\mathrm{g} x, \\max}}, \\frac{\\Delta \\widehat{y}}{c_{\\mathrm{g} y, \\max}}, \\min \\left(\\frac{J_{\\min} \\Delta \\widehat{z}}{c_{\\mathrm{g} z}}\\right)\\right]
```

  - Von Neumann condition (with ``\\mathrm{Re}`` being the Reynolds number):

```math
\\Delta t_\\mathrm{viscous} = \\frac{\\mathrm{Re}}{2} \\min\\limits_\\mathrm{global} \\left[\\left(\\Delta \\widehat{x}\\right)^2, \\left(\\Delta \\widehat{y}\\right)^2, \\left(J \\Delta \\widehat{z}\\right)^2\\right]
```

# Arguments

  - `state`: Model state.

# Returns

  - `::Float64`: Time step.

# See also

  - [`PinCFlow.Update.compute_vertical_wind`](@ref)
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
