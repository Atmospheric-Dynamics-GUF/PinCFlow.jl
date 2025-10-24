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
    \\overline{\\rho} \\left\\langle u' u' \\right\\rangle & = \\overline{\\rho} \\sum_{r, \\lambda, \\mu, \\nu} \\left(F u_\\mathrm{w} u_\\mathrm{w}^*\\right)_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    \\overline{\\rho} \\left\\langle u' v' \\right\\rangle & = \\overline{\\rho} \\sum_{r, \\lambda, \\mu, \\nu} \\left(F u_\\mathrm{w} v_\\mathrm{w}^*\\right)_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    \\overline{\\rho} \\left\\langle u' w' \\right\\rangle & = \\overline{\\rho} \\sum_{r, \\lambda, \\mu, \\nu} \\left(F u_\\mathrm{w} w_\\mathrm{w}^*\\right)_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    \\overline{\\rho} \\left\\langle v' v' \\right\\rangle & = \\overline{\\rho} \\sum_{r, \\lambda, \\mu, \\nu} \\left(F v_\\mathrm{w} v_\\mathrm{w}^*\\right)_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    \\overline{\\rho} \\left\\langle v' w' \\right\\rangle & = \\overline{\\rho} \\sum_{r, \\lambda, \\mu, \\nu} \\left(F v_\\mathrm{w} w_\\mathrm{w}^*\\right)_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    \\left\\langle \\theta' u' \\right\\rangle & = \\sum_{r, \\lambda, \\mu, \\nu} \\left(F \\theta_\\mathrm{w} u_\\mathrm{w}^*\\right)_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    \\left\\langle \\theta' v' \\right\\rangle & = \\sum_{r, \\lambda, \\mu, \\nu} \\left(F \\theta_\\mathrm{w} v_\\mathrm{w}^*\\right)_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    \\mathcal{E} & = \\sum_{r, \\lambda, \\mu, \\nu} \\left(F \\mathcal{A} \\widehat{\\omega}\\right)_{r, i + \\lambda, j + \\mu, k + \\nu}.
\\end{align*}
```

Therein, ``\\left(\\lambda, \\mu, \\nu\\right)`` are index shifts to ray volumes that are at least partially within the grid cell at ``\\left(i, j, k\\right)``, ``F_{r, i + \\lambda, j + \\mu, k + \\nu}`` are the corresponding ray volume fractions and ``\\left(u_\\mathrm{w}, v_\\mathrm{w}, w_\\mathrm{w}, \\theta_\\mathrm{w}\\right)_{r, i + \\lambda, j + \\mu, k + \\nu}`` are the wave amplitudes of the wind and the potential temperature. The computation is based on the relations

```math
\\begin{align*}
    \\overline{\\rho} u_{\\mathrm{w}, r} u_{\\mathrm{w}, r}^* & = \\left(k_r \\widehat{c}_{\\mathrm{g} x, r} - \\mathrm{sgn} \\left(\\left|f\\right|\\right) \\frac{k_r \\widehat{c}_{\\mathrm{g} x, r} + l_r \\widehat{c}_{\\mathrm{g} y, r}}{1 - \\left(\\widehat{\\omega}_r / f\\right)^2}\\right) \\mathcal{A}_r,\\\\
    \\overline{\\rho} u_{\\mathrm{w}, r} v_{\\mathrm{w}, r}^* & = l_r \\widehat{c}_{\\mathrm{g} x, r} \\mathcal{A}_r,\\\\
    \\overline{\\rho} u_{\\mathrm{w}, r} w_{\\mathrm{w}, r}^* & = \\frac{k_r \\widehat{c}_{\\mathrm{g} z, r}}{1 - \\left(f / \\widehat{\\omega}_r\\right)^2} \\mathcal{A}_r,\\\\
    \\overline{\\rho} v_{\\mathrm{w}, r} v_{\\mathrm{w}, r}^* & = \\left(l_r \\widehat{c}_{\\mathrm{g} y, r} - \\mathrm{sgn} \\left(\\left|f\\right|\\right) \\frac{k_r \\widehat{c}_{\\mathrm{g} x, r} + l_r \\widehat{c}_{\\mathrm{g} y, r}}{1 - \\left(\\widehat{\\omega}_r / f\\right)^2}\\right) \\mathcal{A}_r,\\\\
    \\overline{\\rho} v_{\\mathrm{w}, r} w_{\\mathrm{w}, r}^* & = \\frac{l_r \\widehat{c}_{\\mathrm{g} z, r}}{1 - \\left(f / \\widehat{\\omega}_r\\right)^2} \\mathcal{A}_r,\\\\
    \\theta_{\\mathrm{w}, r} u_{\\mathrm{w}, r}^* & = \\frac{f \\overline{\\theta}}{g \\overline{\\rho}} \\frac{l_r m_r N_r^2}{\\widehat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2} \\mathcal{A}_r,\\\\
    \\theta_{\\mathrm{w}, r} v_{\\mathrm{w}, r}^* & = - \\frac{f \\overline{\\theta}}{g \\overline{\\rho}} \\frac{k_r m_r N_r^2}{\\widehat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2} \\mathcal{A}_r,
\\end{align*}
```

where ``N_r^2`` is the squared buoyancy frequency interpolated to the ray-volume position. The components of the intrinsic group velocity are given by

```math
\\begin{align*}
    \\widehat{c}_{\\mathrm{g} x, r} & = \\frac{k_r \\left(N_r^2 - \\widehat{\\omega}_r^2\\right)}{\\widehat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2},\\\\
    \\widehat{c}_{\\mathrm{g} y, r} & = \\frac{l_r \\left(N_r^2 - \\widehat{\\omega}_r^2\\right)}{\\widehat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2},\\\\
    \\widehat{c}_{\\mathrm{g} z, r} & = - \\frac{m_r \\left(\\widehat{\\omega}_r^2 - f^2\\right)}{\\widehat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2}.
\\end{align*}
```

```julia
compute_gw_integrals!(state::State, wkb_mode::SingleColumn)
```

Compute the gravity-wave integrals needed for the computation of the mean-flow impact in single-column mode.

This method computes ``\\overline{\\rho} \\left\\langle u' w' \\right\\rangle``, ``\\overline{\\rho} \\left\\langle v' w' \\right\\rangle``, ``\\left\\langle \\theta' u' \\right\\rangle``, ``\\left\\langle \\theta' v' \\right\\rangle`` and ``\\mathcal{E}`` (see above for details).

```julia
compute_gw_integrals!(state::State, wkb_mode::SteadyState)
```

Compute the gravity-wave integrals needed for the computation of the mean-flow impact in steady-state mode.

This method computes the sums ``\\overline{\\rho} \\left\\langle u' w' \\right\\rangle`` and ``\\overline{\\rho} \\left\\langle v' w' \\right\\rangle`` (see above for details). In contrast to the multi-column and single-column modes, the steady-state mode uses the pseudo-momentum approximation

```math
\\begin{align*}
    \\overline{\\rho} u_{\\mathrm{w}, r} w_{\\mathrm{w}, r}^* & = k_r \\widehat{c}_{\\mathrm{g} z, r} \\mathcal{A}_r,\\\\
    \\overline{\\rho} v_{\\mathrm{w}, r} w_{\\mathrm{w}, r}^* & = l_r \\widehat{c}_{\\mathrm{g} z, r} \\mathcal{A}_r.
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
    (; x_size, y_size) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; branch) = state.namelists.wkb
    (; tref, g_ndim) = state.constants
    (; i0, i1, j0, j1, k0, k1) = domain
    (; dx, dy, dz, x, y, zctilde, jac) = grid
    (; rhobar, thetabar) = state.atmosphere
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

            omir = branch * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            cgirx = kr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
            cgiry = lr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
            cgirz = -mr * (omir^2 - fc^2) / (omir * (khr^2 + mr^2))

            (imin, imax, jmin, jmax) =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for iray in imin:imax
                if x_size > 1
                    dxi = (
                        min(xr + dxr / 2, x[iray] + dx / 2) -
                        max(xr - dxr / 2, x[iray] - dx / 2)
                    )

                    fcpspx = dkr * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jray in jmin:jmax
                    if y_size > 1
                        dyi = (
                            min(yr + dyr / 2, y[jray] + dy / 2) -
                            max(yr - dyr / 2, y[jray] - dy / 2)
                        )

                        fcpspy = dlr * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kmin = get_next_half_level(
                        iray,
                        jray,
                        zr - dzr / 2,
                        state;
                        dkd = 1,
                    )
                    kmax = get_next_half_level(
                        iray,
                        jray,
                        zr + dzr / 2,
                        state;
                        dkd = 1,
                    )

                    for kray in kmin:kmax
                        dzi =
                            min((zr + dzr / 2), zctilde[iray, jray, kray]) -
                            max((zr - dzr / 2), zctilde[iray, jray, kray - 1])

                        fcpspz = dmr * dzi / jac[iray, jray, kray] / dz

                        wadr = fcpspx * fcpspy * fcpspz * rays.dens[r, i, j, k]

                        if x_size > 1
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

                        if x_size > 1 || y_size > 1
                            integrals.uv[iray, jray, kray] += wadr * cgirx * lr
                        end

                        integrals.uw[iray, jray, kray] +=
                            wadr * kr * cgirz / (1 - (fc / omir)^2)

                        if y_size > 1
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
                            integrals.vtheta[iray, jray, kray] -=
                                wadr *
                                fc *
                                n2r *
                                kr *
                                mr *
                                thetabar[iray, jray, kray] / (
                                    rhobar[iray, jray, kray] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )

                            integrals.utheta[iray, jray, kray] +=
                                wadr *
                                fc *
                                n2r *
                                lr *
                                mr *
                                thetabar[iray, jray, kray] / (
                                    rhobar[iray, jray, kray] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )
                        end

                        integrals.e[iray, jray, kray] += wadr * omir

                        compute_leading_order_tracer_fluxes!(
                            state,
                            state.namelists.tracer.tracer_setup,
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
end

function compute_gw_integrals!(state::State, wkb_mode::SingleColumn)
    (; domain, grid) = state
    (; x_size, y_size) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; branch) = state.namelists.wkb
    (; g_ndim, tref) = state.constants
    (; i0, i1, j0, j1, k0, k1) = domain
    (; dx, dy, dz, x, y, zctilde, jac) = grid
    (; rhobar, thetabar) = state.atmosphere
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

            omir = branch * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            cgirz = -mr * (omir^2 - fc^2) / (omir * (khr^2 + mr^2))

            imin, imax, jmin, jmax =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for iray in imin:imax
                if x_size > 1
                    dxi = (
                        min(xr + dxr / 2, x[iray] + dx / 2) -
                        max(xr - dxr / 2, x[iray] - dx / 2)
                    )

                    fcpspx = dkr * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jray in jmin:jmax
                    if y_size > 1
                        dyi = (
                            min(yr + dyr / 2, y[jray] + dy / 2) -
                            max(yr - dyr / 2, y[jray] - dy / 2)
                        )

                        fcpspy = dlr * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kmin = get_next_half_level(
                        iray,
                        jray,
                        zr - dzr / 2,
                        state;
                        dkd = 1,
                    )
                    kmax = get_next_half_level(
                        iray,
                        jray,
                        zr + dzr / 2,
                        state;
                        dkd = 1,
                    )

                    for kray in kmin:kmax
                        dzi =
                            min((zr + dzr / 2), zctilde[iray, jray, kray]) -
                            max((zr - dzr / 2), zctilde[iray, jray, kray - 1])

                        fcpspz = dmr * dzi / jac[iray, jray, kray] / dz

                        wadr = fcpspx * fcpspy * fcpspz * rays.dens[r, i, j, k]

                        integrals.uw[iray, jray, kray] +=
                            wadr * kr * cgirz / (1 - (fc / omir)^2)

                        integrals.vw[iray, jray, kray] +=
                            wadr * lr * cgirz / (1 - (fc / omir)^2)

                        if fc != 0
                            integrals.vtheta[iray, jray, kray] -=
                                wadr *
                                fc *
                                n2r *
                                kr *
                                mr *
                                thetabar[iray, jray, kray] / (
                                    rhobar[iray, jray, kray] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )

                            integrals.utheta[iray, jray, kray] +=
                                wadr *
                                fc *
                                n2r *
                                lr *
                                mr *
                                thetabar[iray, jray, kray] / (
                                    rhobar[iray, jray, kray] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )
                        end

                        integrals.e[iray, jray, kray] += wadr * omir

                        compute_leading_order_tracer_fluxes!(
                            state,
                            state.namelists.tracer.tracer_setup,
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
end

function compute_gw_integrals!(state::State, wkb_mode::SteadyState)
    (; domain, grid) = state
    (; coriolis_frequency) = state.namelists.atmosphere
    (; tref) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, x, y, zctilde, jac) = state.grid
    (; x_size, y_size) = state.namelists.domain
    (; branch) = state.namelists.wkb
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

            omir = branch * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            cgirz = -mr * (omir^2 - fc^2) / (omir * (khr^2 + mr^2))

            imin, imax, jmin, jmax =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for iray in imin:imax
                if x_size > 1
                    dxi = (
                        min(xr + dxr / 2, x[iray] + dx / 2) -
                        max(xr - dxr / 2, x[iray] - dx / 2)
                    )

                    fcpspx = dkr * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jray in jmin:jmax
                    if y_size > 1
                        dyi = (
                            min(yr + dyr / 2, y[jray] + dy / 2) -
                            max(yr - dyr / 2, y[jray] - dy / 2)
                        )

                        fcpspy = dlr * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kmin = get_next_half_level(
                        iray,
                        jray,
                        zr - dzr / 2,
                        state;
                        dkd = 1,
                    )
                    kmax = get_next_half_level(
                        iray,
                        jray,
                        zr + dzr / 2,
                        state;
                        dkd = 1,
                    )

                    for kray in kmin:kmax
                        dzi =
                            min((zr + dzr / 2), zctilde[iray, jray, kray]) -
                            max((zr - dzr / 2), zctilde[iray, jray, kray - 1])

                        fcpspz = dmr * dzi / jac[iray, jray, kray] / dz

                        wadr = fcpspx * fcpspy * fcpspz * rays.dens[r, i, j, k]

                        integrals.uw[iray, jray, kray] += wadr * kr * cgirz

                        integrals.vw[iray, jray, kray] += wadr * lr * cgirz

                        compute_leading_order_tracer_fluxes!(
                            state,
                            state.namelists.tracer.tracer_setup,
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
