function compute_next_order_tracer_fluxes! end

function compute_next_order_tracer_fluxes!(state::State, dt::AbstractFloat)
    (; tracer_setup) = state.namelists.tracer

    compute_next_order_tracer_fluxes!(state, dt, tracer_setup)

    return
end

function compute_next_order_tracer_fluxes!(
    state::State,
    dt::AbstractFloat,
    tracer_setup::NoTracer,
)
    return
end

function compute_next_order_tracer_fluxes!(
    state::State,
    dt::AbstractFloat,
    tracer_setup::TracerOn,
)
    (; domain, grid) = state
    (; x_size, y_size) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; branch) = state.namelists.wkb
    (; tref, g_ndim) = state.constants
    (; i0, i1, j0, j1, k0, k1) = domain
    (; dx, dy, dz, x, y, zctilde, jac) = grid
    (; rhobar, thetabar) = state.atmosphere
    (; nray, rays) = state.wkb
    (; uhat, vhat, what, bhat, pihat) = state.wkb.integrals
    (; uhatold, vhatold, whatold, bhatold, pihatold, mat) =
        state.wkb.wkbauxiliaries

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

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
                    factor = dxi / dx
                else
                    fcpspx = 1.0
                    factor = 1.0
                end

                for jray in jmin:jmax
                    if y_size > 1
                        dyi = (
                            min(yr + dyr / 2, y[jray] + dy / 2) -
                            max(yr - dyr / 2, y[jray] - dy / 2)
                        )

                        fcpspy = dlr * dyi / dy
                        factor *= dyi / dy
                    else
                        fcpspy = 1.0
                        factor *= 1.0
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
                        factor *= dzi / jac[iray, jray, kray] / dz

                        wadr = fcpspx * fcpspy * fcpspz * rays.dens[r, i, j, k]

                        mat .= [
                            -1im*omir -fc 0 0 1im*kr
                            fc -1im*omir 0 0 1im*lr
                            0 0 -1im*omir -sqrt(n2r) 1im
                            0 0 sqrt(n2r) -1im*omir 0
                            1im*kr 1im*lr 1im*mr 0 0
                        ]

                        uhatr = interpolate_scalar(state, xr, yr, zr, uhat)
                        uhatoldr =
                            interpolate_scalar(state, xr, yr, zr, uhatold)
                        vhatr = interpolate_scalar(state, xr, yr, zr, vhat)
                        vhatoldr =
                            interpolate_scalar(state, xr, yr, zr, vhatold)
                        whatr = interpolate_scalar(state, xr, yr, zr, what)
                        whatoldr =
                            interpolate_scalar(state, xr, yr, zr, whatold)
                        bhatr = interpolate_scalar(state, xr, yr, zr, bhat)
                        bhatoldr =
                            interpolate_scalar(state, xr, yr, zr, bhatold)
                        pihatr = interpolate_scalar(state, xr, yr, zr, pihat)
                        pihatoldr =
                            interpolate_scalar(state, xr, yr, zr, pihatold)

                        duhatdt = (uhatr - uhatoldr) / dt
                    end
                end
            end
        end
    end

    return
end
