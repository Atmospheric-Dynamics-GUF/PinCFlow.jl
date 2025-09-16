"""
```julia
compute_gw_tendencies!(state::State)
```

Compute the gravity-wave impact on the momentum and mass-weighted potential temperature.

Calculates the tendencies that are to be added to the equations for momentum and mass-weighted potential temperature. These are given by

```math
\\begin{align*}
    \\left(\\frac{\\partial \\rho_\\mathrm{b} u_\\mathrm{b}}{\\partial t}\\right)_\\mathrm{w} & = - \\frac{\\rho_\\mathrm{b}}{\\overline{\\rho}} \\left(\\frac{M_{u u, i + 1} - M_{u u, i - 1}}{2 \\Delta \\widehat{x}} + G^{13} \\frac{M_{u u, k + 1} - M_{u u, k - 1}}{2 \\Delta \\widehat{z}} + \\frac{M_{u v, j + 1} - M_{u v, j - 1}}{2 \\Delta \\widehat{y}} + G^{23} \\frac{M_{u v, k + 1} - M_{u v, k - 1}}{2 \\Delta \\widehat{z}} + \\frac{M_{u w, k + 1} - M_{u w, k - 1}}{2 J \\Delta \\widehat{z}}\\right) + E_u,\\\\
    \\left(\\frac{\\partial \\rho_\\mathrm{b} v_\\mathrm{b}}{\\partial t}\\right)_\\mathrm{w} & = - \\frac{\\rho_\\mathrm{b}}{\\overline{\\rho}} \\left(\\frac{M_{v u, i + 1} - M_{v u, i - 1}}{2 \\Delta \\widehat{x}} + G^{13} \\frac{M_{v u, k + 1} - M_{v u, k - 1}}{2 \\Delta \\widehat{z}} + \\frac{M_{v v, j + 1} - M_{v v, j - 1}}{2 \\Delta \\widehat{y}} + G^{23} \\frac{M_{v v, k + 1} - M_{v v, k - 1}}{2 \\Delta \\widehat{z}} + \\frac{M_{v w, k + 1} - M_{v w, k - 1}}{2 J \\Delta \\widehat{z}}\\right) + E_v,\\\\
    \\left(\\frac{\\partial \\rho_\\mathrm{b} \\widehat{w}_\\mathrm{b}}{\\partial t}\\right)_\\mathrm{w} & = G^{13} \\left(\\frac{\\partial \\rho_\\mathrm{b} u_\\mathrm{b}}{\\partial t}\\right)_\\mathrm{w} + G^{23} \\left(\\frac{\\partial \\rho_\\mathrm{b} v_\\mathrm{b}}{\\partial t}\\right)_\\mathrm{w},\\\\
    \\left(\\frac{\\partial P_\\mathrm{b}}{\\partial t}\\right)_\\mathrm{w} & = - \\rho_\\mathrm{b} \\left(\\frac{T_{u, i + 1} - T_{u, i - 1}}{2 \\Delta \\widehat{x}} + G^{13} \\frac{T_{u, k + 1} - T_{u, k - 1}}{2 \\Delta \\widehat{z}} + \\frac{T_{v, j + 1} - T_{u, j - 1}}{2 \\Delta \\widehat{y}} + G^{23} \\frac{T_{v, k + 1} - T_{v, k - 1}}{2 \\Delta \\widehat{z}}\\right),
\\end{align*}
```

where ``\\left(u_\\mathrm{b}, v_\\mathrm{b}, \\widehat{w}_\\mathrm{b}\\right)`` are the components of the transformed (i.e. terrain-following) resolved wind, ``\\rho_\\mathrm{b}`` is the resolved density (including the reference part ``\\overline{\\rho}``) and ``P_\\mathrm{b}`` is the resolved mass-weighted potential temperature. For a documentation of the fluxes ``\\left(M_{u u}, M_{u v}, M_{u w}, M_{v v}, M_{v w}\\right)``, ``\\left(E_u, E_v\\right)`` and ``\\left(T_u, T_v\\right)``, see [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_integrals!`](@ref). Below `state.namelists.wkb.zmin_wkb_dim`, all tendencies are set to zero.

# Arguments

  - `state::State`: Model state.
"""
function compute_gw_tendencies! end

function compute_gw_tendencies!(state::State)
    (; sizex, sizey) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; zmin_wkb_dim) = state.namelists.wkb
    (; tref, lref) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, ztfc, jac, met) = state.grid
    (; rhostrattfc) = state.atmosphere
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

        if sizex > 1
            tendencies.dudt[i, j, k] -=
                rhotot / rhostrattfc[i, j, k] * (
                    (integrals.uu[i + 1, j, k] - integrals.uu[i - 1, j, k]) /
                    (2.0 * dx) +
                    met[i, j, k, 1, 3] *
                    (integrals.uu[i, j, k + 1] - integrals.uu[i, j, k - 1]) /
                    (2.0 * dz)
                )
        end

        if sizey > 1
            tendencies.dudt[i, j, k] -=
                rhotot / rhostrattfc[i, j, k] * (
                    (integrals.uv[i, j + 1, k] - integrals.uv[i, j - 1, k]) /
                    (2.0 * dy) +
                    met[i, j, k, 2, 3] *
                    (integrals.uv[i, j, k + 1] - integrals.uv[i, j, k - 1]) /
                    (2.0 * dz)
                )
        end

        tendencies.dudt[i, j, k] += rhotot * integrals.etx[i, j, k]

        # Compute the drag on the meridional wind.

        tendencies.dvdt[i, j, k] =
            -rhotot / rhostrattfc[i, j, k] / jac[i, j, k] *
            (integrals.vw[i, j, k + 1] - integrals.vw[i, j, k - 1]) / (2.0 * dz)

        if sizex > 1
            tendencies.dvdt[i, j, k] -=
                rhotot / rhostrattfc[i, j, k] * (
                    (integrals.uv[i + 1, j, k] - integrals.uv[i - 1, j, k]) /
                    (2.0 * dx) +
                    met[i, j, k, 1, 3] *
                    (integrals.uv[i, j, k + 1] - integrals.uv[i, j, k - 1]) /
                    (2.0 * dz)
                )
        end

        if sizey > 1
            tendencies.dvdt[i, j, k] -=
                rhotot / rhostrattfc[i, j, k] * (
                    (integrals.vv[i, j + 1, k] - integrals.vv[i, j - 1, k]) /
                    (2.0 * dy) +
                    met[i, j, k, 2, 3] *
                    (integrals.vv[i, j, k + 1] - integrals.vv[i, j, k - 1]) /
                    (2.0 * dz)
                )
        end

        tendencies.dvdt[i, j, k] += rhotot * integrals.ety[i, j, k]

        # Compute the heating.

        if fc != 0.0 && (sizex > 1 || sizey > 1)
            if sizex > 1
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

            if sizey > 1
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
    end
end
