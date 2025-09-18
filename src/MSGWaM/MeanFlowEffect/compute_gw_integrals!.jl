"""
```julia
compute_gw_integrals!(state::State)
```

Compute the gravity-wave integrals needed for the computation of the mean-flow impact by dispatching to a WKB-mode-specific method.

```julia
compute_gw_integrals!(state::State, wkb_mode::MultiColumn)
```

Compute the gravity-wave integrals needed for the computation of the mean-flow impact in multi-column mode.

This method computes the sums

```math
\\begin{align*}
    M_{u u} & = \\overline{\\rho} \\sum_{\\lambda, \\mu, \\nu, \\alpha} \\left(F u_\\mathrm{w} u_\\mathrm{w}^*\\right)_{i + \\lambda, j + \\mu, k + \\nu, \\alpha},\\\\
    M_{u v} & = \\overline{\\rho} \\sum_{\\lambda, \\mu, \\nu, \\alpha} \\left(F u_\\mathrm{w} v_\\mathrm{w}^*\\right)_{i + \\lambda, j + \\mu, k + \\nu, \\alpha},\\\\
    M_{u w} & = \\overline{\\rho} \\sum_{\\lambda, \\mu, \\nu, \\alpha} \\left(F u_\\mathrm{w} w_\\mathrm{w}^*\\right)_{i + \\lambda, j + \\mu, k + \\nu, \\alpha},\\\\
    M_{v v} & = \\overline{\\rho} \\sum_{\\lambda, \\mu, \\nu, \\alpha} \\left(F v_\\mathrm{w} v_\\mathrm{w}^*\\right)_{i + \\lambda, j + \\mu, k + \\nu, \\alpha},\\\\
    M_{v w} & = \\overline{\\rho} \\sum_{\\lambda, \\mu, \\nu, \\alpha} \\left(F v_\\mathrm{w} w_\\mathrm{w}^*\\right)_{i + \\lambda, j + \\mu, k + \\nu, \\alpha},\\\\
    T_u & = \\sum_{\\lambda, \\mu, \\nu, \\alpha} \\left(F u_\\mathrm{w} \\theta_\\mathrm{w}^*\\right)_{i + \\lambda, j + \\mu, k + \\nu, \\alpha},\\\\
    T_v & = \\sum_{\\lambda, \\mu, \\nu, \\alpha} \\left(F v_\\mathrm{w} \\theta_\\mathrm{w}^*\\right)_{i + \\lambda, j + \\mu, k + \\nu, \\alpha},\\\\
    E_u & = \\frac{f}{\\overline{\\theta}} \\sum_{\\lambda, \\mu, \\nu, \\alpha} \\left(F v_\\mathrm{w} \\theta_\\mathrm{w}^*\\right)_{i + \\lambda, j + \\mu, k + \\nu, \\alpha},\\\\
    E_v & = \\frac{f}{\\overline{\\theta}} \\sum_{\\lambda, \\mu, \\nu, \\alpha} \\left(F u_\\mathrm{w} \\theta_\\mathrm{w}^*\\right)_{i + \\lambda, j + \\mu, k + \\nu, \\alpha},\\\\
    E & = \\sum_{\\lambda, \\mu, \\nu, \\alpha} \\left(F \\mathcal{E}\\right)_{i + \\lambda, j + \\mu, k + \\nu, \\alpha}.
\\end{align*}
```

Therein, ``\\left(\\lambda, \\mu, \\nu\\right)`` are index shifts to ray volumes that are at least partially within the grid cell at ``\\left(i, j, k\\right)``, ``F_{i + \\lambda, j + \\mu, k + \\nu, \\alpha}`` are the corresponding ray volume fractions and ``\\left(u_\\mathrm{w}, v_\\mathrm{w}, w_\\mathrm{w}, \\theta_\\mathrm{w}, \\mathcal{E}\\right)_{i + \\lambda, j + \\mu, k + \\nu, \\alpha}`` are the wave amplitudes of the wind, those of the potential temperature and the wave-energy densities. The background density, Coriolis parameter and background potential temperature are denoted by ``\\overline{\\rho}``, ``f`` and ``\\overline{\\theta}``, respectively. The computation is based on the relations

```math
\\begin{align*}
    \\overline{\\rho} u_{\\mathrm{w}, \\alpha} u_{\\mathrm{w}, \\alpha}^* & = \\left(k_\\alpha \\widehat{c}_{\\mathrm{g} x, \\alpha} - \\mathrm{sgn} \\left(\\left|f\\right|\\right) \\frac{k_\\alpha \\widehat{c}_{\\mathrm{g} x, \\alpha} + l_\\alpha \\widehat{c}_{\\mathrm{g} y, \\alpha}}{1 - \\left(\\widehat{\\omega}_\\alpha / f\\right)^2}\\right) \\mathcal{A}_\\alpha,\\\\
    \\overline{\\rho} u_{\\mathrm{w}, \\alpha} v_{\\mathrm{w}, \\alpha}^* & = l_\\alpha \\widehat{c}_{\\mathrm{g} x, \\alpha} \\mathcal{A}_\\alpha,\\\\
    \\overline{\\rho} u_{\\mathrm{w}, \\alpha} w_{\\mathrm{w}, \\alpha}^* & = \\frac{k_\\alpha \\widehat{c}_{\\mathrm{g} z, \\alpha}}{1 - \\left(f / \\widehat{\\omega}_\\alpha\\right)^2} \\mathcal{A}_\\alpha,\\\\
    \\overline{\\rho} v_{\\mathrm{w}, \\alpha} v_{\\mathrm{w}, \\alpha}^* & = \\left(l_\\alpha \\widehat{c}_{\\mathrm{g} y, \\alpha} - \\mathrm{sgn} \\left(\\left|f\\right|\\right) \\frac{k_\\alpha \\widehat{c}_{\\mathrm{g} x, \\alpha} + l_\\alpha \\widehat{c}_{\\mathrm{g} y, \\alpha}}{1 - \\left(\\widehat{\\omega}_\\alpha / f\\right)^2}\\right) \\mathcal{A}_\\alpha,\\\\
    \\overline{\\rho} v_{\\mathrm{w}, \\alpha} w_{\\mathrm{w}, \\alpha}^* & = \\frac{l_\\alpha \\widehat{c}_{\\mathrm{g} z, \\alpha}}{1 - \\left(f / \\widehat{\\omega}_\\alpha\\right)^2} \\mathcal{A}_\\alpha,\\\\
    \\overline{\\rho} u_{\\mathrm{w}, \\alpha} \\theta_{\\mathrm{w}, \\alpha}^* & = \\frac{f}{g \\overline{\\theta}} \\frac{l_\\alpha m_\\alpha N_\\alpha^2}{\\widehat{\\omega}_\\alpha \\left|\\boldsymbol{k}_\\alpha\\right|^2},\\\\
    \\overline{\\rho} v_{\\mathrm{w}, \\alpha} \\theta_{\\mathrm{w}, \\alpha}^* & = - \\frac{f}{g \\overline{\\theta}} \\frac{k_\\alpha m_\\alpha N_\\alpha^2}{\\widehat{\\omega}_\\alpha \\left|\\boldsymbol{k}_\\alpha\\right|^2},\\\\
    \\mathcal{E}_\\alpha & = \\widehat{\\omega}_\\alpha \\mathcal{A}_\\alpha,
\\end{align*}
```

where ``\\boldsymbol{k}_\\alpha = \\left(k_\\alpha, l_\\alpha, m_\\alpha\\right)^\\mathrm{T}``, ``\\widehat{\\boldsymbol{c}}_{\\mathrm{g}, \\alpha} = \\left(\\widehat{c}_{\\mathrm{g} x, \\alpha}, \\widehat{c}_{\\mathrm{g} y, \\alpha}, \\widehat{c}_{\\mathrm{g} z, \\alpha}\\right)^\\mathrm{T}``, ``\\widehat{\\omega}_\\alpha``, ``\\mathcal{A}_\\alpha``, ``g`` and ``N_\\alpha^2`` are the wavevector, intrinsic group velocity, intrinsic frequency, wave-action density, gravitational acceleration and squared buoyancy frequency interpolated to the ray-volume position, respectively. The components of the intrinsic group velocity are given by

```math
\\begin{align*}
    \\widehat{c}_{\\mathrm{g} x, \\alpha} & = \\frac{k_\\alpha \\left(N_\\alpha^2 - \\widehat{\\omega}_\\alpha^2\\right)}{\\widehat{\\omega}_\\alpha \\left|\\boldsymbol{k}_\\alpha\\right|^2},\\\\
    \\widehat{c}_{\\mathrm{g} y, \\alpha} & = \\frac{l_\\alpha \\left(N_\\alpha^2 - \\widehat{\\omega}_\\alpha^2\\right)}{\\widehat{\\omega}_\\alpha \\left|\\boldsymbol{k}_\\alpha\\right|^2},\\\\
    \\widehat{c}_{\\mathrm{g} z, \\alpha} & = - \\frac{m_\\alpha \\left(\\widehat{\\omega}_\\alpha^2 - f^2\\right)}{\\widehat{\\omega}_\\alpha \\left|\\boldsymbol{k}_\\alpha\\right|^2}.
\\end{align*}
```

```julia
compute_gw_integrals!(state::State, wkb_mode::SingleColumn)
```

Compute the gravity-wave integrals needed for the computation of the mean-flow impact in single-column mode.

This method computes the sums ``M_{u w}``, ``M_{v w}``, ``T_u``, ``T_v``, ``E_u``, ``E_v`` and ``E`` (see above for details).

```julia
compute_gw_integrals!(state::State, wkb_mode::SteadyState)
```

Compute the gravity-wave integrals needed for the computation of the mean-flow impact in steady-state mode.

This method computes the sums ``M_{u w}`` and ``M_{v w}`` (see above for details). In contrast to the multi-column and single-column modes, the steady-state mode uses the pseudo-momentum approximation

```math
\\begin{align*}
    \\overline{\\rho} u_{\\mathrm{w}, \\alpha} w_{\\mathrm{w}, \\alpha}^* & = k_\\alpha \\widehat{c}_{\\mathrm{g} z, \\alpha} \\mathcal{A}_\\alpha,\\\\
    \\overline{\\rho} v_{\\mathrm{w}, \\alpha} w_{\\mathrm{w}, \\alpha}^* & = l_\\alpha \\widehat{c}_{\\mathrm{g} z, \\alpha} \\mathcal{A}_\\alpha.
\\end{align*}
```

# Arguments

  - `state::State`: Model state.

  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_stratification`](@ref)

  - [`PinCFlow.MSGWaM.MeanFlowEffect.compute_horizontal_cell_indices`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.get_next_half_level`](@ref)
"""
function compute_gw_integrals! end

function compute_gw_integrals!(state::State)
    (; wkb_mode) = state.namelists.wkb
    compute_gw_integrals!(state, wkb_mode)
    return
end

function compute_gw_integrals!(state::State, wkb_mode::MultiColumn)
    (; domain, grid) = state
    (; sizex, sizey) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; branchr) = state.namelists.wkb
    (; tref, g_ndim) = state.constants
    (; i0, i1, j0, j1, k0, k1, io, jo) = domain
    (; dx, dy, dz, x, y, ztildetfc, jac) = grid
    (; rhostrattfc, thetastrattfc) = state.atmosphere
    (; nray, rays, integrals) = state.wkb

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    for field in fieldnames(WKBIntegrals)
        getfield(integrals, field) .= 0.0
    end

    set_tracer_field_zero!(state)

    @ivy for k in (k0 - 1):(k1 + 1),
        j in (j0 - 1):(j1 + 1),
        i in (i0 - 1):(i1 + 1)

        for r in 1:nray[i, j, k]
            if rays.dens[r, i, j, k] == 0
                continue
            end

            xr = rays.x[r, i, j, k]
            yr = rays.y[r, i, j, k]
            zr = rays.z[r, i, j, k]

            dxr = rays.dxray[r, i, j, k]
            dyr = rays.dyray[r, i, j, k]
            dzr = rays.dzray[r, i, j, k]

            kr = rays.k[r, i, j, k]
            lr = rays.l[r, i, j, k]
            mr = rays.m[r, i, j, k]

            dkr = rays.dkray[r, i, j, k]
            dlr = rays.dlray[r, i, j, k]
            dmr = rays.dmray[r, i, j, k]

            khr = sqrt(kr^2 + lr^2)

            n2r = interpolate_stratification(zr, state, N2())

            omir =
                branchr * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            cgirx = kr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
            cgiry = lr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
            cgirz = -mr * (omir^2 - fc^2) / (omir * (khr^2 + mr^2))

            (imin, imax, jmin, jmax) =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for iray in imin:imax
                if sizex > 1
                    dxi = (
                        min(xr + dxr / 2, x[io + iray] + dx / 2) -
                        max(xr - dxr / 2, x[io + iray] - dx / 2)
                    )

                    fcpspx = dkr * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jray in jmin:jmax
                    if sizey > 1
                        dyi = (
                            min(yr + dyr / 2, y[jo + jray] + dy / 2) -
                            max(yr - dyr / 2, y[jo + jray] - dy / 2)
                        )

                        fcpspy = dlr * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kmin = get_next_half_level(iray, jray, zr - dzr / 2, state)
                    kmax = get_next_half_level(iray, jray, zr + dzr / 2, state)

                    for kray in kmin:kmax
                        dzi =
                            min((zr + dzr / 2), ztildetfc[iray, jray, kray]) -
                            max((zr - dzr / 2), ztildetfc[iray, jray, kray - 1])

                        fcpspz = dmr * dzi / jac[iray, jray, kray] / dz

                        wadr = fcpspx * fcpspy * fcpspz * rays.dens[r, i, j, k]

                        if sizex > 1
                            if fc != 0
                                integrals.uu[iray, jray, kray] +=
                                    wadr * (
                                        kr * cgirx -
                                        (kr * cgirx + lr * cgiry) /
                                        (1 - (omir / fc)^2)
                                    )
                            else
                                integrals.uu[iray, jray, kray] +=
                                    wadr * kr * cgirx
                            end
                        end

                        if sizex > 1 || sizey > 1
                            integrals.uv[iray, jray, kray] += wadr * cgirx * lr
                        end

                        integrals.uw[iray, jray, kray] +=
                            wadr * kr * cgirz / (1 - (fc / omir)^2)

                        if sizey > 1
                            if fc != 0
                                integrals.vv[iray, jray, kray] +=
                                    wadr * (
                                        lr * cgiry -
                                        (kr * cgirx + lr * cgiry) /
                                        (1 - (omir / fc)^2)
                                    )
                            else
                                integrals.vv[iray, jray, kray] +=
                                    wadr * lr * cgiry
                            end
                        end

                        integrals.vw[iray, jray, kray] +=
                            wadr * lr * cgirz / (1 - (fc / omir)^2)

                        if fc != 0
                            integrals.etx[iray, jray, kray] +=
                                wadr * fc^2 * n2r * kr * mr / (
                                    rhostrattfc[iray, jray, kray] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )

                            integrals.ety[iray, jray, kray] +=
                                wadr * fc^2 * n2r * lr * mr / (
                                    rhostrattfc[iray, jray, kray] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )
                        end

                        integrals.e[iray, jray, kray] += wadr * omir

                        compute_leading_order_tracer_fluxes!(
                            state,
                            state.namelists.tracer.tracersetup,
                            fc,
                            omir,
                            kr,
                            lr,
                            mr,
                            wadr,
                            xr,
                            yr,
                            zr,
                            iray,
                            jray,
                            kray,
                        )
                    end
                end
            end
        end
    end

    @ivy if fc != 0
        for k in k0:k1, j in j0:j1, i in i0:i1
            integrals.utheta[i, j, k] =
                thetastrattfc[i, j, k] / fc * integrals.ety[i, j, k]
            integrals.vtheta[i, j, k] =
                -thetastrattfc[i, j, k] / fc * integrals.etx[i, j, k]
        end
    end
end

function compute_gw_integrals!(state::State, wkb_mode::SingleColumn)
    (; domain, grid) = state
    (; sizex, sizey) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; branchr) = state.namelists.wkb
    (; g_ndim, tref) = state.constants
    (; i0, i1, j0, j1, k0, k1, io, jo) = domain
    (; dx, dy, dz, x, y, ztildetfc, jac) = grid
    (; rhostrattfc, thetastrattfc) = state.atmosphere
    (; nray, rays, integrals) = state.wkb

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    for field in fieldnames(WKBIntegrals)
        getfield(integrals, field) .= 0.0
    end

    @ivy for k in (k0 - 1):(k1 + 1),
        j in (j0 - 1):(j1 + 1),
        i in (i0 - 1):(i1 + 1)

        for r in 1:nray[i, j, k]
            if rays.dens[r, i, j, k] == 0
                continue
            end

            xr = rays.x[r, i, j, k]
            yr = rays.y[r, i, j, k]
            zr = rays.z[r, i, j, k]

            dxr = rays.dxray[r, i, j, k]
            dyr = rays.dyray[r, i, j, k]
            dzr = rays.dzray[r, i, j, k]

            kr = rays.k[r, i, j, k]
            lr = rays.l[r, i, j, k]
            mr = rays.m[r, i, j, k]

            dkr = rays.dkray[r, i, j, k]
            dlr = rays.dlray[r, i, j, k]
            dmr = rays.dmray[r, i, j, k]

            khr = sqrt(kr^2 + lr^2)

            n2r = interpolate_stratification(zr, state, N2())

            omir =
                branchr * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            cgirz = -mr * (omir^2 - fc^2) / (omir * (khr^2 + mr^2))

            imin, imax, jmin, jmax =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for iray in imin:imax
                if sizex > 1
                    dxi = (
                        min(xr + dxr / 2, x[io + iray] + dx / 2) -
                        max(xr - dxr / 2, x[io + iray] - dx / 2)
                    )

                    fcpspx = dkr * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jray in jmin:jmax
                    if sizey > 1
                        dyi = (
                            min(yr + dyr / 2, y[jo + jray] + dy / 2) -
                            max(yr - dyr / 2, y[jo + jray] - dy / 2)
                        )

                        fcpspy = dlr * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kmin = get_next_half_level(iray, jray, zr - dzr / 2, state)
                    kmax = get_next_half_level(iray, jray, zr + dzr / 2, state)

                    for kray in kmin:kmax
                        dzi =
                            min((zr + dzr / 2), ztildetfc[iray, jray, kray]) -
                            max((zr - dzr / 2), ztildetfc[iray, jray, kray - 1])

                        fcpspz = dmr * dzi / jac[iray, jray, kray] / dz

                        wadr = fcpspx * fcpspy * fcpspz * rays.dens[r, i, j, k]

                        integrals.uw[iray, jray, kray] +=
                            wadr * kr * cgirz / (1 - (fc / omir)^2)

                        integrals.vw[iray, jray, kray] +=
                            wadr * lr * cgirz / (1 - (fc / omir)^2)

                        if fc != 0
                            integrals.etx[iray, jray, kray] +=
                                wadr * fc^2 * n2r * kr * mr / (
                                    rhostrattfc[iray, jray, kray] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )

                            integrals.ety[iray, jray, kray] +=
                                wadr * fc^2 * n2r * lr * mr / (
                                    rhostrattfc[iray, jray, kray] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )
                        end

                        integrals.e[iray, jray, kray] += wadr * omir

                        compute_leading_order_tracer_fluxes!(
                            state,
                            state.namelists.tracer.tracersetup,
                            fc,
                            omir,
                            0.0,
                            0.0,
                            wnrm,
                            wadr,
                            xr,
                            yr,
                            zr,
                            iray,
                            jray,
                            kray,
                        )
                    end
                end
            end
        end
    end

    @ivy if fc != 0
        for k in k0:k1, j in j0:j1, i in i0:i1
            integrals.utheta[i, j, k] =
                thetastrattfc[i, j, k] / fc * integrals.ety[i, j, k]
            integrals.vtheta[i, j, k] =
                -thetastrattfc[i, j, k] / fc * integrals.etx[i, j, k]
        end
    end
end

function compute_gw_integrals!(state::State, wkb_mode::SteadyState)
    (; domain, grid) = state
    (; coriolis_frequency) = state.namelists.atmosphere
    (; tref) = state.constants
    (; i0, i1, j0, j1, k0, k1, io, jo) = state.domain
    (; dx, dy, dz, x, y, ztildetfc, jac) = state.grid
    (; sizex, sizey) = state.namelists.domain
    (; branchr) = state.namelists.wkb
    (; nray, rays, integrals) = state.wkb

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    for field in fieldnames(WKBIntegrals)
        getfield(integrals, field) .= 0.0
    end

    @ivy for k in (k0 - 1):(k1 + 1),
        j in (j0 - 1):(j1 + 1),
        i in (i0 - 1):(i1 + 1)

        for r in 1:nray[i, j, k]
            if rays.dens[r, i, j, k] == 0
                continue
            end

            xr = rays.x[r, i, j, k]
            yr = rays.y[r, i, j, k]
            zr = rays.z[r, i, j, k]

            dxr = rays.dxray[r, i, j, k]
            dyr = rays.dyray[r, i, j, k]
            dzr = rays.dzray[r, i, j, k]

            kr = rays.k[r, i, j, k]
            lr = rays.l[r, i, j, k]
            mr = rays.m[r, i, j, k]

            dkr = rays.dkray[r, i, j, k]
            dlr = rays.dlray[r, i, j, k]
            dmr = rays.dmray[r, i, j, k]

            khr = sqrt(kr^2 + lr^2)

            n2r = interpolate_stratification(zr, state, N2())

            omir =
                branchr * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            cgirz = -mr * (omir^2 - fc^2) / (omir * (khr^2 + mr^2))

            imin, imax, jmin, jmax =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for iray in imin:imax
                if sizex > 1
                    dxi = (
                        min(xr + dxr / 2, x[io + iray] + dx / 2) -
                        max(xr - dxr / 2, x[io + iray] - dx / 2)
                    )

                    fcpspx = dkr * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jray in jmin:jmax
                    if sizey > 1
                        dyi = (
                            min(yr + dyr / 2, y[jo + jray] + dy / 2) -
                            max(yr - dyr / 2, y[jo + jray] - dy / 2)
                        )

                        fcpspy = dlr * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kmin = get_next_half_level(iray, jray, zr - dzr / 2, state)
                    kmax = get_next_half_level(iray, jray, zr + dzr / 2, state)

                    for kray in kmin:kmax
                        dzi =
                            min((zr + dzr / 2), ztildetfc[iray, jray, kray]) -
                            max((zr - dzr / 2), ztildetfc[iray, jray, kray - 1])

                        fcpspz = dmr * dzi / jac[iray, jray, kray] / dz

                        wadr = fcpspx * fcpspy * fcpspz * rays.dens[r, i, j, k]

                        integrals.uw[iray, jray, kray] += wadr * kr * cgirz

                        integrals.vw[iray, jray, kray] += wadr * lr * cgirz

                        compute_leading_order_tracer_fluxes!(
                            state,
                            state.namelists.tracer.tracersetup,
                            fc,
                            omir,
                            wnrk,
                            wnrl,
                            wnrm,
                            wadr,
                            xr,
                            yr,
                            zr,
                            iray,
                            jray,
                            kray,
                        )
                    end
                end
            end
        end
    end
end
