"""
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

function compute_gw_integrals!(state::State, wkb_mode::MultiColumn)
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

    for kzrv in (k0 - 1):(k1 + 1),
        jyrv in (j0 - 1):(j1 + 1),
        ixrv in (i0 - 1):(i1 + 1)

        for iray in 1:nray[ixrv, jyrv, kzrv]
            if rays.dens[iray, ixrv, jyrv, kzrv] == 0
                continue
            end

            xr = rays.x[iray, ixrv, jyrv, kzrv]
            yr = rays.y[iray, ixrv, jyrv, kzrv]
            zr = rays.z[iray, ixrv, jyrv, kzrv]

            dxr = rays.dxray[iray, ixrv, jyrv, kzrv]
            dyr = rays.dyray[iray, ixrv, jyrv, kzrv]
            dzr = rays.dzray[iray, ixrv, jyrv, kzrv]

            kr = rays.k[iray, ixrv, jyrv, kzrv]
            lr = rays.l[iray, ixrv, jyrv, kzrv]
            mr = rays.m[iray, ixrv, jyrv, kzrv]

            dkr = rays.dkray[iray, ixrv, jyrv, kzrv]
            dlr = rays.dlray[iray, ixrv, jyrv, kzrv]
            dmr = rays.dmray[iray, ixrv, jyrv, kzrv]

            khr = sqrt(kr^2 + lr^2)

            n2r = interpolate_stratification(zr, state, N2())

            omir =
                branchr * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            cgirx = kr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
            cgiry = lr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
            cgirz = -mr * (omir^2 - fc^2) / (omir * (khr^2 + mr^2))

            (ixmin, ixmax, jymin, jymax) =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for ix in ixmin:ixmax
                if sizex > 1
                    dxi = (
                        min(xr + dxr / 2, x[io + ix] + dx / 2) -
                        max(xr - dxr / 2, x[io + ix] - dx / 2)
                    )

                    fcpspx = dkr * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jy in jymin:jymax
                    if sizey > 1
                        dyi = (
                            min(yr + dyr / 2, y[jo + jy] + dy / 2) -
                            max(yr - dyr / 2, y[jo + jy] - dy / 2)
                        )

                        fcpspy = dlr * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kzmin =
                        get_next_half_level(ix, jy, zr - dzr / 2, domain, grid; dkd = 1, dku = 1)
                    kzmax =
                        get_next_half_level(ix, jy, zr + dzr / 2, domain, grid; dkd = 1, dku = 1)

                    for kz in kzmin:kzmax
                        dzi =
                            min((zr + dzr / 2), ztildetfc[ix, jy, kz]) -
                            max((zr - dzr / 2), ztildetfc[ix, jy, kz - 1])

                        fcpspz = dmr * dzi / jac[ix, jy, kz] / dz

                        wadr =
                            fcpspx *
                            fcpspy *
                            fcpspz *
                            rays.dens[iray, ixrv, jyrv, kzrv]

                        if sizex > 1
                            if fc != 0
                                integrals.uu[ix, jy, kz] +=
                                    wadr * (
                                        kr * cgirx -
                                        (kr * cgirx + lr * cgiry) /
                                        (1 - (omir / fc)^2)
                                    )
                            else
                                integrals.uu[ix, jy, kz] += wadr * kr * cgirx
                            end
                        end

                        if sizex > 1 || sizey > 1
                            integrals.uv[ix, jy, kz] += wadr * cgirx * lr
                        end

                        integrals.uw[ix, jy, kz] +=
                            wadr * kr * cgirz / (1 - (fc / omir)^2)

                        if sizey > 1
                            if fc != 0
                                integrals.vv[ix, jy, kz] +=
                                    wadr * (
                                        lr * cgiry -
                                        (kr * cgirx + lr * cgiry) /
                                        (1 - (omir / fc)^2)
                                    )
                            else
                                integrals.vv[ix, jy, kz] += wadr * lr * cgiry
                            end
                        end

                        integrals.vw[ix, jy, kz] +=
                            wadr * lr * cgirz / (1 - (fc / omir)^2)

                        if fc != 0
                            integrals.etx[ix, jy, kz] +=
                                wadr * fc^2 * n2r * kr * mr / (
                                    rhostrattfc[ix, jy, kz] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )

                            integrals.ety[ix, jy, kz] +=
                                wadr * fc^2 * n2r * lr * mr / (
                                    rhostrattfc[ix, jy, kz] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )
                        end

                        integrals.e[ix, jy, kz] += wadr * omir
                    end
                end
            end
        end
    end

    if fc != 0
        for kz in k0:k1, jy in j0:j1, ix in i0:i1
            integrals.utheta[ix, jy, kz] =
                thetastrattfc[ix, jy, kz] / fc * integrals.ety[ix, jy, kz]
            integrals.vtheta[ix, jy, kz] =
                -thetastrattfc[ix, jy, kz] / fc * integrals.etx[ix, jy, kz]
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

    for kzrv in (k0 - 1):(k1 + 1),
        jyrv in (j0 - 1):(j1 + 1),
        ixrv in (i0 - 1):(i1 + 1)

        for iray in 1:nray[ixrv, jyrv, kzrv]
            if rays.dens[iray, ixrv, jyrv, kzrv] == 0
                continue
            end

            xr = rays.x[iray, ixrv, jyrv, kzrv]
            yr = rays.y[iray, ixrv, jyrv, kzrv]
            zr = rays.z[iray, ixrv, jyrv, kzrv]

            dxr = rays.dxray[iray, ixrv, jyrv, kzrv]
            dyr = rays.dyray[iray, ixrv, jyrv, kzrv]
            dzr = rays.dzray[iray, ixrv, jyrv, kzrv]

            kr = rays.k[iray, ixrv, jyrv, kzrv]
            lr = rays.l[iray, ixrv, jyrv, kzrv]
            mr = rays.m[iray, ixrv, jyrv, kzrv]

            dkr = rays.dkray[iray, ixrv, jyrv, kzrv]
            dlr = rays.dlray[iray, ixrv, jyrv, kzrv]
            dmr = rays.dmray[iray, ixrv, jyrv, kzrv]

            khr = sqrt(kr^2 + lr^2)

            n2r = interpolate_stratification(zr, state, N2())

            omir =
                branchr * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            cgirz = -mr * (omir^2 - fc^2) / (omir * (khr^2 + mr^2))

            ixmin, ixmax, jymin, jymax =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for ix in ixmin:ixmax
                if sizex > 1
                    dxi = (
                        min(xr + dxr / 2, x[io + ix] + dx / 2) -
                        max(xr - dxr / 2, x[io + ix] - dx / 2)
                    )

                    fcpspx = dkr * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jy in jymin:jymax
                    if sizey > 1
                        dyi = (
                            min(yr + dyr / 2, y[jo + jy] + dy / 2) -
                            max(yr - dyr / 2, y[jo + jy] - dy / 2)
                        )

                        fcpspy = dlr * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kzmin =
                        get_next_half_level(ix, jy, zr - dzr / 2, domain, grid)
                    kzmax =
                        get_next_half_level(ix, jy, zr + dzr / 2, domain, grid)

                    for kz in kzmin:kzmax
                        dzi =
                            min((zr + dzr / 2), ztildetfc[ix, jy, kz]) -
                            max((zr - dzr / 2), ztildetfc[ix, jy, kz - 1])

                        fcpspz = dmr * dzi / jac[ix, jy, kz] / dz

                        wadr =
                            fcpspx *
                            fcpspy *
                            fcpspz *
                            rays.dens[iray, ixrv, jyrv, kzrv]

                        integrals.uw[ix, jy, kz] +=
                            wadr * kr * cgirz / (1 - (fc / omir)^2)

                        integrals.vw[ix, jy, kz] +=
                            wadr * lr * cgirz / (1 - (fc / omir)^2)

                        if fc != 0
                            integrals.etx[ix, jy, kz] +=
                                wadr * fc^2 * n2r * kr * mr / (
                                    rhostrattfc[ix, jy, kz] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )

                            integrals.ety[ix, jy, kz] +=
                                wadr * fc^2 * n2r * lr * mr / (
                                    rhostrattfc[ix, jy, kz] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )
                        end

                        integrals.e[ix, jy, kz] += wadr * omir
                    end
                end
            end
        end
    end

    if fc != 0
        for kz in k0:k1, jy in j0:j1, ix in i0:i1
            integrals.utheta[ix, jy, kz] =
                thetastrattfc[ix, jy, kz] / fc * integrals.ety[ix, jy, kz]
            integrals.vtheta[ix, jy, kz] =
                -thetastrattfc[ix, jy, kz] / fc * integrals.etx[ix, jy, kz]
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

    for kzrv in (k0 - 1):(k1 + 1),
        jyrv in (j0 - 1):(j1 + 1),
        ixrv in (i0 - 1):(i1 + 1)

        for iray in 1:nray[ixrv, jyrv, kzrv]
            if rays.dens[iray, ixrv, jyrv, kzrv] == 0
                continue
            end

            xr = rays.x[iray, ixrv, jyrv, kzrv]
            yr = rays.y[iray, ixrv, jyrv, kzrv]
            zr = rays.z[iray, ixrv, jyrv, kzrv]

            dxr = rays.dxray[iray, ixrv, jyrv, kzrv]
            dyr = rays.dyray[iray, ixrv, jyrv, kzrv]
            dzr = rays.dzray[iray, ixrv, jyrv, kzrv]

            kr = rays.k[iray, ixrv, jyrv, kzrv]
            lr = rays.l[iray, ixrv, jyrv, kzrv]
            mr = rays.m[iray, ixrv, jyrv, kzrv]

            dkr = rays.dkray[iray, ixrv, jyrv, kzrv]
            dlr = rays.dlray[iray, ixrv, jyrv, kzrv]
            dmr = rays.dmray[iray, ixrv, jyrv, kzrv]

            khr = sqrt(kr^2 + lr^2)

            n2r = interpolate_stratification(zr, state, N2())

            omir =
                branchr * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            cgirz = -mr * (omir^2 - fc^2) / (omir * (khr^2 + mr^2))

            ixmin, ixmax, jymin, jymax =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for ix in ixmin:ixmax
                if sizex > 1
                    dxi = (
                        min(xr + dxr / 2, x[io + ix] + dx / 2) -
                        max(xr - dxr / 2, x[io + ix] - dx / 2)
                    )

                    fcpspx = dkr * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jy in jymin:jymax
                    if sizey > 1
                        dyi = (
                            min(yr + dyr / 2, y[jo + jy] + dy / 2) -
                            max(yr - dyr / 2, y[jo + jy] - dy / 2)
                        )

                        fcpspy = dlr * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kzmin =
                        get_next_half_level(ix, jy, zr - dzr / 2, domain, grid)
                    kzmax =
                        get_next_half_level(ix, jy, zr + dzr / 2, domain, grid)

                    for kz in kzmin:kzmax
                        dzi =
                            min((zr + dzr / 2), ztildetfc[ix, jy, kz]) -
                            max((zr - dzr / 2), ztildetfc[ix, jy, kz - 1])

                        fcpspz = dmr * dzi / jac[ix, jy, kz] / dz

                        wadr =
                            fcpspx *
                            fcpspy *
                            fcpspz *
                            rays.dens[iray, ixrv, jyrv, kzrv]

                        integrals.uw[ix, jy, kz] += wadr * kr * cgirz

                        integrals.vw[ix, jy, kz] += wadr * lr * cgirz
                    end
                end
            end
        end
    end
end
