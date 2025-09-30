"""
```julia
compute_gw_tendencies!(state::State)
```

Compute the gravity-wave impact on the momentum and mass-weighted potential temperature.

Calculates the tendencies that are to be added to the equations for momentum and mass-weighted potential temperature. These are given by

```math
\\begin{align*}
    \\left(\\frac{\\partial \\rho_\\mathrm{b} u_\\mathrm{b}}{\\partial t}\\right)_\\mathrm{w} & = - \\frac{\\rho_\\mathrm{b}}{\\overline{\\rho}} \\left[\\frac{\\left(\\overline{\\rho} \\left\\langle u' u' \\right\\rangle\\right)_{i + 1} - \\left(\\overline{\\rho} \\left\\langle u' u' \\right\\rangle\\right)_{i - 1}}{2 \\Delta \\widehat{x}} + G^{13} \\frac{\\left(\\overline{\\rho} \\left\\langle u' u' \\right\\rangle\\right)_{k + 1} - \\left(\\overline{\\rho} \\left\\langle u' u' \\right\\rangle\\right)_{k - 1}}{2 \\Delta \\widehat{z}}\\right.\\\\
    & \\qquad \\qquad + \\frac{\\left(\\overline{\\rho} \\left\\langle u' v' \\right\\rangle\\right)_{j + 1} - \\left(\\overline{\\rho} \\left\\langle u' v' \\right\\rangle\\right)_{j - 1}}{2 \\Delta \\widehat{y}} + G^{23} \\frac{\\left(\\overline{\\rho} \\left\\langle u' v' \\right\\rangle\\right)_{k + 1} - \\left(\\overline{\\rho} \\left\\langle u' v' \\right\\rangle\\right)_{k - 1}}{2 \\Delta \\widehat{z}}\\\\
    & \\qquad \\qquad + \\left.\\frac{\\left(\\overline{\\rho} \\left\\langle u' w' \\right\\rangle\\right)_{k + 1} - \\left(\\overline{\\rho} \\left\\langle u' w' \\right\\rangle\\right)_{k - 1}}{2 J \\Delta \\widehat{z}}\\right] - \\rho_\\mathrm{b} \\frac{f}{\\overline{\\theta}} \\left\\langle \\theta' v' \\right\\rangle,\\\\
    \\left(\\frac{\\partial \\rho_\\mathrm{b} v_\\mathrm{b}}{\\partial t}\\right)_\\mathrm{w} & = - \\frac{\\rho_\\mathrm{b}}{\\overline{\\rho}} \\left[\\frac{\\left(\\overline{\\rho} \\left\\langle v' u' \\right\\rangle\\right)_{i + 1} - \\left(\\overline{\\rho} \\left\\langle v' u' \\right\\rangle\\right)_{i - 1}}{2 \\Delta \\widehat{x}} + G^{13} \\frac{\\left(\\overline{\\rho} \\left\\langle v' u' \\right\\rangle\\right)_{k + 1} - \\left(\\overline{\\rho} \\left\\langle v' u' \\right\\rangle\\right)_{k - 1}}{2 \\Delta \\widehat{z}}\\right.\\\\
    & \\qquad \\qquad + \\frac{\\left(\\overline{\\rho} \\left\\langle v' v' \\right\\rangle\\right)_{j + 1} - \\left(\\overline{\\rho} \\left\\langle v' v' \\right\\rangle\\right)_{j - 1}}{2 \\Delta \\widehat{y}} + G^{23} \\frac{\\left(\\overline{\\rho} \\left\\langle v' v' \\right\\rangle\\right)_{k + 1} - \\left(\\overline{\\rho} \\left\\langle v' v' \\right\\rangle\\right)_{k - 1}}{2 \\Delta \\widehat{z}}\\\\
    & \\qquad \\qquad + \\left.\\frac{\\left(\\overline{\\rho} \\left\\langle v' w' \\right\\rangle\\right)_{k + 1} - \\left(\\overline{\\rho} \\left\\langle v' w' \\right\\rangle\\right)_{k - 1}}{2 J \\Delta \\widehat{z}}\\right] + \\rho_\\mathrm{b} \\frac{f}{\\overline{\\theta}} \\left\\langle \\theta' u' \\right\\rangle,\\\\
    \\left(\\frac{\\partial \\rho_\\mathrm{b} \\widehat{w}_\\mathrm{b}}{\\partial t}\\right)_\\mathrm{w} & = G^{13} \\left(\\frac{\\partial \\rho_\\mathrm{b} u_\\mathrm{b}}{\\partial t}\\right)_\\mathrm{w} + G^{23} \\left(\\frac{\\partial \\rho_\\mathrm{b} v_\\mathrm{b}}{\\partial t}\\right)_\\mathrm{w},\\\\
    \\left(\\frac{\\partial P_\\mathrm{b}}{\\partial t}\\right)_\\mathrm{w} & = - \\rho_\\mathrm{b} \\left[\\frac{\\left(\\overline{\\rho} \\left\\langle \\theta' u' \\right\\rangle\\right)_{i + 1} - \\left(\\overline{\\rho} \\left\\langle \\theta' u' \\right\\rangle\\right)_{i - 1}}{2 \\Delta \\widehat{x}} + G^{13} \\frac{\\left(\\overline{\\rho} \\left\\langle \\theta' u' \\right\\rangle\\right)_{k + 1} - \\left(\\overline{\\rho} \\left\\langle \\theta' u' \\right\\rangle\\right)_{k - 1}}{2 \\Delta \\widehat{z}}\\right.\\\\
    & \\qquad \\qquad + \\left.\\frac{\\left(\\overline{\\rho} \\left\\langle \\theta' v' \\right\\rangle\\right)_{j + 1} - \\left(\\overline{\\rho} \\left\\langle \\theta' v' \\right\\rangle\\right)_{j - 1}}{2 \\Delta \\widehat{y}} + G^{23} \\frac{\\left(\\overline{\\rho} \\left\\langle \\theta' v' \\right\\rangle\\right)_{k + 1} - \\left(\\overline{\\rho} \\left\\langle \\theta' v' \\right\\rangle\\right)_{k - 1}}{2 \\Delta \\widehat{z}}\\right],
\\end{align*}
```

where ``\\left(u_\\mathrm{b}, v_\\mathrm{b}, \\widehat{w}_\\mathrm{b}\\right)`` are the components of the transformed (i.e. terrain-following) resolved wind, ``\\rho_\\mathrm{b}`` is the resolved density (including the reference part ``\\overline{\\rho}``) and ``P_\\mathrm{b}`` is the resolved mass-weighted potential temperature. For a documentation of the fluxes, see [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_integrals!`](@ref). Below `state.namelists.wkb.zmin_wkb_dim`, all tendencies are set to zero.

# Arguments

  - `state::State`: Model state.
"""
function compute_gw_tendencies! end

function compute_gw_tendencies!(state::State)
    (; ndx, ndy) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; zmin_wkb_dim) = state.namelists.wkb
    (; tref, lref) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, ztfc, jac, met) = state.grid
    (; rhostrattfc, thetastrattfc) = state.atmosphere
    (; rho) = state.variables.predictands
    (; integrals, tendencies) = state.wkb

    # Set the Coriolis parameter.
    fc = coriolis_frequency * tref

    for field in fieldnames(WKBTendencies)
        getfield(tendencies, field) .= 0.0
    end

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        if ztfc[i, j, k] < zmin_wkb_dim / lref
            continue
        end

        rhotot = rho[i, j, k] + rhostrattfc[i, j, k]

        # Compute the drag on the zonal wind.

        tendencies.dudt[i, j, k] =
            -rhotot / rhostrattfc[i, j, k] / jac[i, j, k] *
            (integrals.uw[i, j, k + 1] - integrals.uw[i, j, k - 1]) / (2.0 * dz)

        if ndx > 1
            tendencies.dudt[i, j, k] -=
                rhotot / rhostrattfc[i, j, k] * (
                    (integrals.uu[i + 1, j, k] - integrals.uu[i - 1, j, k]) /
                    (2.0 * dx) +
                    met[i, j, k, 1, 3] *
                    (integrals.uu[i, j, k + 1] - integrals.uu[i, j, k - 1]) /
                    (2.0 * dz)
                )
        end

        if ndy > 1
            tendencies.dudt[i, j, k] -=
                rhotot / rhostrattfc[i, j, k] * (
                    (integrals.uv[i, j + 1, k] - integrals.uv[i, j - 1, k]) /
                    (2.0 * dy) +
                    met[i, j, k, 2, 3] *
                    (integrals.uv[i, j, k + 1] - integrals.uv[i, j, k - 1]) /
                    (2.0 * dz)
                )
        end

        tendencies.dudt[i, j, k] -=
            rhotot * fc / thetastrattfc[i, j, k] * integrals.vtheta[i, j, k]

        # Compute the drag on the meridional wind.

        tendencies.dvdt[i, j, k] =
            -rhotot / rhostrattfc[i, j, k] / jac[i, j, k] *
            (integrals.vw[i, j, k + 1] - integrals.vw[i, j, k - 1]) / (2.0 * dz)

        if ndx > 1
            tendencies.dvdt[i, j, k] -=
                rhotot / rhostrattfc[i, j, k] * (
                    (integrals.uv[i + 1, j, k] - integrals.uv[i - 1, j, k]) /
                    (2.0 * dx) +
                    met[i, j, k, 1, 3] *
                    (integrals.uv[i, j, k + 1] - integrals.uv[i, j, k - 1]) /
                    (2.0 * dz)
                )
        end

        if ndy > 1
            tendencies.dvdt[i, j, k] -=
                rhotot / rhostrattfc[i, j, k] * (
                    (integrals.vv[i, j + 1, k] - integrals.vv[i, j - 1, k]) /
                    (2.0 * dy) +
                    met[i, j, k, 2, 3] *
                    (integrals.vv[i, j, k + 1] - integrals.vv[i, j, k - 1]) /
                    (2.0 * dz)
                )
        end

        tendencies.dvdt[i, j, k] +=
            rhotot * fc / thetastrattfc[i, j, k] * integrals.utheta[i, j, k]

        # Compute the heating.

        if fc != 0.0 && (ndx > 1 || ndy > 1)
            if ndx > 1
                tendencies.dthetadt[i, j, k] -=
                    rhotot * (
                        (
                            integrals.utheta[i + 1, j, k] -
                            integrals.utheta[i - 1, j, k]
                        ) / (2.0 * dx) +
                        met[i, j, k, 1, 3] * (
                            integrals.utheta[i, j, k + 1] -
                            integrals.utheta[i, j, k - 1]
                        ) / (2.0 * dz)
                    )
            end

            if ndy > 1
                tendencies.dthetadt[i, j, k] -=
                    rhotot * (
                        (
                            integrals.vtheta[i, j + 1, k] -
                            integrals.vtheta[i, j - 1, k]
                        ) / (2.0 * dy) +
                        met[i, j, k, 2, 3] * (
                            integrals.vtheta[i, j, k + 1] -
                            integrals.vtheta[i, j, k - 1]
                        ) / (2.0 * dz)
                    )
            end
        end

        compute_leading_order_tracer_forcing!(
            state,
            i,
            j,
            k,
            state.namelists.tracer.tracersetup,
        )
    end
end
