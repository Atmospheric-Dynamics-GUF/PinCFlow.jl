
# ============================================================
# Prescribed lower-boundary predictands
#
# The incident fields are provided by the generic boundary
# source callback:
#
#     source(x, y, z, t) -> (u, w, b, p)
#
# An optional v component is also supported.
#
# All source fields are dimensional:
#   u, v, w : m/s
#   b       : m/s^2
#   p       : m^2/s^2 (kinematic pressure)
#
# No source-specific physical parameters are used here.
# ============================================================

function set_prescribed_lower_predictands!(state::State)
    (; domain, grid, constants, atmosphere) = state
    (; ko, i0, i1, j0, j1, k0) = domain

    # Only the bottom MPI subdomains have a physical lower boundary.
    ko == 0 || return nothing

    (; nbx, nby, nbz) = state.namelists.domain
    (; lref, uref, thetaref, g, kappa, rsp) = constants
    (; x, y, dx, dy, zc, zctilde, hb) = grid
    (; rhobar, thetabar) = atmosphere
    (; rhop, u, v, w, pip) = state.variables.predictands

    # The physical vertical velocity equals transformed w only
    # for a flat lower boundary. General terrain is not handled yet.
    all(iszero, hb) ||
        throw(ArgumentError(
            "Prescribed lower-boundary forcing currently requires flat topography."
        ))

    ii = (i0 - nbx):(i1 + nbx)
    jj = (j0 - nby):(j1 + nby)

    @inbounds for j in jj, i in ii
        xc = x[i] * lref
        yc = y[j] * lref

        # Staggered horizontal-velocity positions.
        xu = (x[i] + 0.5 * dx) * lref
        yv = (y[j] + 0.5 * dy) * lref

        # ----------------------------------------------------
        # 1. Prescribe vertical velocity on the lower face.
        # ----------------------------------------------------

        zface = zctilde[i, j, k0 - 1] * lref

        incident_face = evaluate_boundary_source(
            state, xc, yc, zface
        )

        w[i, j, k0 - 1] = incident_face.w / uref

        # ----------------------------------------------------
        # 2. Cell-centered ghost values.
        #
        # For a variable q:
        # qghost = qinc_ghost + s * (qinside - qinc_inside)
        #
        # s = +1 for u, v, pip
        # s = -1 for rhop
        # ----------------------------------------------------

        for n in 1:nbz
            kg = k0 - n
            km = k0 + n - 1

            zg = zc[i, j, kg] * lref
            zm = zc[i, j, km] * lref

            # Cell-centered incident fields.
            sg = evaluate_boundary_source(state, xc, yc, zg)
            sm = evaluate_boundary_source(state, xc, yc, zm)

            # Incident u at the zonal velocity faces.
            sug = evaluate_boundary_source(state, xu, yc, zg)
            sum = evaluate_boundary_source(state, xu, yc, zm)

            # Incident v at the meridional velocity faces.
            svg = evaluate_boundary_source(state, xc, yv, zg)
            svm = evaluate_boundary_source(state, xc, yv, zm)

            vg = hasproperty(svg, :v) ? svg.v : 0.0
            vm = hasproperty(svm, :v) ? svm.v : 0.0

            # Even residual reflection: horizontal velocities.
            u[i, j, kg] =
                u[i, j, km] + (sug.u - sum.u) / uref

            v[i, j, kg] =
                v[i, j, km] + (vg - vm) / uref

            # Incident density fluctuations from buoyancy.
            # rhop is nondimensionalized by rhoref.
            rhop_inc_g = -rhobar[i, j, kg] * sg.b / g
            rhop_inc_m = -rhobar[i, j, km] * sm.b / g

            # Odd residual reflection: density fluctuations.
            rhop[i, j, kg] =
                -rhop[i, j, km] + rhop_inc_g + rhop_inc_m

            # Convert kinematic pressure to Exner pressure.
            theta_g = thetabar[i, j, kg] * thetaref
            theta_m = thetabar[i, j, km] * thetaref

            pip_inc_g = kappa * sg.p / (rsp * theta_g)
            pip_inc_m = kappa * sm.p / (rsp * theta_m)

            # Even residual reflection: pressure.
            pip[i, j, kg] =
                pip[i, j, km] + pip_inc_g - pip_inc_m
        end

        # ----------------------------------------------------
        # 3. Ghost values of staggered vertical velocity.
        #
        # The first lower index is the prescribed boundary
        # face. Ghost values begin one index below it.
        # ----------------------------------------------------

        for n in 2:nbz
            kg = k0 - n
            km = k0 + n - 2

            zg = zctilde[i, j, kg] * lref
            zm = zctilde[i, j, km] * lref

            sg = evaluate_boundary_source(state, xc, yc, zg)
            sm = evaluate_boundary_source(state, xc, yc, zm)

            # Odd reflection of the residual vertical velocity.
            w[i, j, kg] =
                -w[i, j, km] + (sg.w + sm.w) / uref
        end
    end

    return nothing
end
