function compute_next_order_tracer_fluxes! end

function compute_next_order_tracer_fluxes!(state::State, dt::AbstractFloat)
    (; tracer_setup) = state.namelists.tracer

    @dispatch_tracer_setup compute_next_order_tracer_fluxes!(state, dt, Val(tracer_setup))

    return
end

function compute_next_order_tracer_fluxes!(
    state::State,
    dt::AbstractFloat,
    tracer_setup::Val{:NoTracer},
)
    return
end

function compute_next_order_tracer_fluxes!(
    state::State,
    dt::AbstractFloat,
    tracer_setup::Val{:TracerOn},
)
    (; domain, grid) = state
    (; x_size, y_size) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; branch) = state.namelists.wkb
    (; tref, g_ndim, kappa) = state.constants
    (; i0, i1, j0, j1, k0, k1) = domain
    (; dx, dy, dz, x, y, zctilde, jac) = grid
    (; rhobar, thetabar, pbar) = state.atmosphere
    (; nray, rays) = state.wkb
    (; uhat, vhat, what, bhat, pihat, chihat) = state.wkb.integrals
    (; uhatold, vhatold, whatold, bhatold, chihatold, mat) =
        state.wkb.wkbauxiliaries
    (; rho, pip) = state.variables.predictands
    (; uchi1, vchi1, wchi1) = state.tracer.tracerwkbintegrals

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

            khr = sqrt(kr^2 + lr^2)

            n2r = interpolate_stratification(zr, state, N2())

            omir = branch * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            (imin, imax, jmin, jmax) =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for iray in imin:imax
                if x_size > 1
                    dxi = (
                        min(xr + dxr / 2, x[iray] + dx / 2) -
                        max(xr - dxr / 2, x[iray] - dx / 2)
                    )

                    factor = dxi / dx
                else
                    factor = 1.0
                end

                for jray in jmin:jmax
                    if y_size > 1
                        dyi = (
                            min(yr + dyr / 2, y[jray] + dy / 2) -
                            max(yr - dyr / 2, y[jray] - dy / 2)
                        )

                        factor *= dyi / dy
                    else
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

                        factor *= dzi / jac[iray, jray, kray] / dz

                        uhatr =
                            interpolate_scalar(xr, yr, zr, state, uhat, None())
                        vhatr =
                            interpolate_scalar(xr, yr, zr, state, vhat, None())
                        whatr =
                            interpolate_scalar(xr, yr, zr, state, what, None())
                        pihatr =
                            interpolate_scalar(xr, yr, zr, state, pihat, None())
                        bhatr =
                            interpolate_scalar(xr, yr, zr, state, bhat, None())
                        chihatr = interpolate_scalar(
                            xr,
                            yr,
                            zr,
                            state,
                            chihat,
                            None(),
                        )

                        uhatoldr = interpolate_scalar(
                            xr,
                            yr,
                            zr,
                            state,
                            uhatold,
                            None(),
                        )
                        vhatoldr = interpolate_scalar(
                            xr,
                            yr,
                            zr,
                            state,
                            vhatold,
                            None(),
                        )
                        whatoldr = interpolate_scalar(
                            xr,
                            yr,
                            zr,
                            state,
                            whatold,
                            None(),
                        )
                        bhatoldr = interpolate_scalar(
                            xr,
                            yr,
                            zr,
                            state,
                            bhatold,
                            None(),
                        )
                        chihatoldr = interpolate_scalar(
                            xr,
                            yr,
                            zr,
                            state,
                            chihatold,
                            None(),
                        )

                        ur = interpolate_mean_flow(xr, yr, zr, state, U())
                        vr = interpolate_mean_flow(xr, yr, zr, state, V())
                        rhor =
                            interpolate_scalar(xr, yr, zr, state, rho, None())

                        rhobarr = interpolate_scalar(
                            xr,
                            yr,
                            zr,
                            state,
                            rhobar,
                            None(),
                        )
                        pbarr =
                            interpolate_scalar(xr, yr, zr, state, pbar, None())
                        thetabarr = interpolate_scalar(
                            xr,
                            yr,
                            zr,
                            state,
                            thetabar,
                            None(),
                        )

                        thetar = pbarr / (rhor + rhobarr) - thetabarr

                        dthetadxr = interpolate_scalar(
                            xr,
                            yr,
                            zr,
                            state,
                            pbar ./ (rho .+ rhobar) .- thetabar,
                            DX(),
                        )
                        dthetadyr = interpolate_scalar(
                            xr,
                            yr,
                            zr,
                            state,
                            pbar ./ (rho .+ rhobar) .- thetabar,
                            DY(),
                        )
                        dthetadzr = interpolate_scalar(
                            xr,
                            yr,
                            zr,
                            state,
                            pbar ./ (rho .+ rhobar) .- thetabar,
                            DZ(),
                        )

                        thetahatr = bhatr * thetabarr / g_ndim

                        dupdxr = interpolate_scalar(
                            xr,
                            yr,
                            zr,
                            state,
                            rhobar .* thetabar .* uhat,
                            DX(),
                        )
                        dvpdyr = interpolate_scalar(
                            xr,
                            yr,
                            zr,
                            state,
                            rhobar .* thetabar .* vhat,
                            DY(),
                        )
                        dwpdzr = interpolate_scalar(
                            xr,
                            yr,
                            zr,
                            state,
                            rhobar .* thetabar .* what,
                            DZ(),
                        )

                        dudxr = interpolate_mean_flow(xr, yr, zr, state, DUDX())
                        dudyr = interpolate_mean_flow(xr, yr, zr, state, DUDY())
                        dudzr = interpolate_mean_flow(xr, yr, zr, state, DUDZ())
                        dvdxr = interpolate_mean_flow(xr, yr, zr, state, DVDX())
                        dvdyr = interpolate_mean_flow(xr, yr, zr, state, DVDY())
                        dvdzr = interpolate_mean_flow(xr, yr, zr, state, DVDZ())
                        dpidxr =
                            interpolate_scalar(xr, yr, zr, state, pip, DX())
                        dpidyr =
                            interpolate_scalar(xr, yr, zr, state, pip, DY())
                        dpidzr =
                            interpolate_scalar(xr, yr, zr, state, pip, DZ())
                        dchidxr =
                            interpolate_mean_flow(xr, yr, zr, state, DChiDX())
                        dchidyr =
                            interpolate_mean_flow(xr, yr, zr, state, DChiDY())
                        dchidzr =
                            interpolate_mean_flow(xr, yr, zr, state, DChiDZ())

                        duhatdxr =
                            interpolate_scalar(xr, yr, zr, state, uhat, DX())
                        duhatdyr =
                            interpolate_scalar(xr, yr, zr, state, uhat, DY())
                        dvhatdxr =
                            interpolate_scalar(xr, yr, zr, state, vhat, DX())
                        dvhatdyr =
                            interpolate_scalar(xr, yr, zr, state, vhat, DY())
                        dwhatdxr =
                            interpolate_scalar(xr, yr, zr, state, what, DX())
                        dwhatdyr =
                            interpolate_scalar(xr, yr, zr, state, what, DY())
                        dpihatdxr =
                            interpolate_scalar(xr, yr, zr, state, pihat, DX())
                        dpihatdyr =
                            interpolate_scalar(xr, yr, zr, state, pihat, DY())
                        dpihatdzr =
                            interpolate_scalar(xr, yr, zr, state, pihat, DZ())
                        dbhatdxr =
                            interpolate_scalar(xr, yr, zr, state, bhat, DX())
                        dbhatdyr =
                            interpolate_scalar(xr, yr, zr, state, bhat, DY())
                        dchihatdxr =
                            interpolate_scalar(xr, yr, zr, state, chihat, DX())
                        dchihatdyr =
                            interpolate_scalar(xr, yr, zr, state, chihat, DY())

                        duhatdt = (uhatr - uhatoldr) / dt
                        dvhatdt = (vhatr - vhatoldr) / dt
                        dwhatdt = (whatr - whatoldr) / dt
                        dbhatdt = (bhatr - bhatoldr) / dt
                        dchihatdt = (chihatr - chihatoldr) / dt

                        ru =
                            -(duhatdt + ur * duhatdxr + vr * duhatdyr) -
                            (uhatr * dudxr + vhatr * dudyr + whatr * dudzr) -
                            1 / kappa * (
                                thetabarr * dpihatdxr +
                                1im * kr * thetar * pihatr +
                                thetahatr * dpidxr
                            )
                        rv =
                            -(dvhatdt + ur * dvhatdxr + vr * dvhatdyr) -
                            (uhatr * dvdxr + vhatr * dvdyr + whatr * dvdzr) -
                            1 / kappa * (
                                thetabarr * dpihatdyr +
                                1im * lr * thetar * pihatr +
                                thetahatr * dpidyr
                            )
                        rw =
                            -(dwhatdt + ur * dwhatdxr + vr * dwhatdyr) -
                            1 / kappa * (
                                thetabarr * dpihatdzr +
                                1im * mr * thetar * pihatr +
                                thetahatr * dpidzr
                            )
                        rb =
                            -(dbhatdt + ur * dbhatdxr + vr * dbhatdyr) - (
                                uhatr / thetabarr * dthetadxr +
                                vhatr / thetabarr * dthetadyr +
                                whatr / thetabarr * dthetadzr
                            )

                        rpi =
                            -1 / (rhobarr * thetabarr) *
                            (dupdxr + dvpdyr + dwpdzr)

                        mat .= [
                            -1im*omir -fc 0 0 1im*kr
                            fc -1im*omir 0 0 1im*lr
                            0 0 -1im*omir -sqrt(n2r) 1im*mr
                            0 0 sqrt(n2r) -1im*omir 0
                            1im*kr 1im*lr 1im*mr 0 0
                        ]

                        matinv = pinv(mat; atol = 1e-12)

                        next_order_wave_amplitudes =
                            matinv * [ru, rv, rw, rb / sqrt(n2r), rpi]

                        uhatr2 = next_order_wave_amplitudes[1]
                        vhatr2 = next_order_wave_amplitudes[2]
                        whatr2 = next_order_wave_amplitudes[3]

                        chihatr2 =
                            -1im / omir * (
                                (
                                    dchihatdt +
                                    ur * dchihatdxr +
                                    vr * dchihatdyr
                                ) + (
                                    uhatr2 * dchidxr +
                                    vhatr2 * dchidyr +
                                    whatr2 * dchidzr
                                )
                            )

                        uchi1[iray, jray, kray] +=
                            real(
                                uhatr * conj(chihatr2) + uhatr2 * conj(chihatr),
                            ) * factor
                        vchi1[iray, jray, kray] +=
                            real(
                                vhatr * conj(chihatr2) + vhatr2 * conj(chihatr),
                            ) * factor
                        wchi1[iray, jray, kray] +=
                            real(
                                whatr * conj(chihatr2) + whatr2 * conj(chihatr),
                            ) * factor
                    end
                end
            end
        end
    end

    return
end
