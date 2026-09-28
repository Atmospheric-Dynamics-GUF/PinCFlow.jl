"""
```julia
initialize_rays!(state::State)
```

Complete the initialization of MS-GWaM by dispatching to a WKB-mode-specific method.

```julia
initialize_rays!(state::State, wkb_mode::NoWKB)
```

Return for non-WKB configurations.

```julia
initialize_rays!(
    state::State,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)
```

Complete the initialization of MS-GWaM.

In each grid cell, `wave_modes` wave modes are computed, using `state.namelists.wkb.initial_wave_field`, as well as `activate_orographic_source!` for mountain waves. For each of these modes, `nrx * nry * nrz * nrk * nrl * nrm` ray volumes are then defined such that they evenly divide the volume one would get for `nrx = nry = nrz = nrk = nrl = nrm = 1` (the parameters are taken from `state.namelists.wkb`). Finally, the maximum group velocities are determined for the corresponding CFL condition that is used in the computation of the time step.

# Arguments

  - `state`: Model state.

  - `wkb_mode`: Approximations used by MS-GWaM.

# See also

  - [`PinCFlow.MSGWaM.RaySources.activate_orographic_source!`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_stratification`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_mean_flow`](@ref)
"""
function initialize_rays! end

function initialize_rays!(state::State)
    (; wkb_mode) = state.namelists.wkb
    initialize_rays!(state, wkb_mode)
    return
end

function initialize_rays!(state::State, wkb_mode::NoWKB)
    return
end

function initialize_rays!(
    state::State,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)
    (; source_mode) = state.namelists.wkb
    (; triad_mode, ray_volume_ini) = state.namelists.triad

    initialize_rays!(state, wkb_mode, triad_mode, ray_volume_ini, source_mode)

    return
end

function initialize_rays!(
    state::State,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
    triad_mode::NoTriad,
    ray_volume_ini::GaussianDist,
    source_mode::Union{NoRaySource, OrographicSource},
)
    error("GaussianDist ray-volume initialization is only supported with Triad2D.")
end

function initialize_rays!(
    state::State,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
    triad_mode::Union{NoTriad, Triad2D},
    ray_volume_ini::UniformDist,
    source_mode::Union{NoRaySource, OrographicSource},
)
    (; x_size, y_size) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (;
        nrx,
        nry,
        nrz,
        nrk,
        nrl,
        nrm,
        wave_modes,
        dkr_factor,
        dlr_factor,
        dmr_factor,
        initial_wave_field,
    ) = state.namelists.wkb
    (; lref, tref, rhoref, uref) = state.constants
    (; comm, master, nxx, nyy, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, x, y, zc, jac) = state.grid
    (;
        nray_max,
        nray_wrk,
        n_sfc,
        nray,
        rays,
        surface_indices,
        cgx_max,
        cgy_max,
        cgz_max,
    ) = state.wkb

    if x_size == 1 && nrk != 1
        error(
            "Error in initialize_rays!: nrk must be 1 when x_size == 1. ",
            "Otherwise identical zero-width k ray volumes are initialized.",
        )
    end

    if y_size == 1 && nrl != 1
        error(
            "Error in initialize_rays!: nrl must be 1 when y_size == 1. ",
            "Otherwise identical zero-width l ray volumes are initialized.",
        )
    end
    # Set Coriolis parameter.
    fc = coriolis_frequency * tref
    # Initialize local arrays.
    omi_ini = zeros(wave_modes, nxx, nyy, nzz)
    wnk_ini = zeros(wave_modes, nxx, nyy, nzz)
    wnl_ini = zeros(wave_modes, nxx, nyy, nzz)
    wnm_ini = zeros(wave_modes, nxx, nyy, nzz)
    wad_ini = zeros(wave_modes, nxx, nyy, nzz)

    # Compute initial wavenumbers, intrinsic frequencies and wave-action
    # densities with initial_wave_field.
    if wkb_mode != SteadyState()
        for k in k0:k1, j in j0:j1, i in i0:i1, alpha in 1:wave_modes
            (kdim, ldim, mdim, omegadim, adim) = initial_wave_field(
                alpha,
                x[i] * lref,
                y[j] * lref,
                zc[i, j, k] * lref,
            )
            wnk_ini[alpha, i, j, k] = kdim * lref
            wnl_ini[alpha, i, j, k] = ldim * lref
            wnm_ini[alpha, i, j, k] = mdim * lref
            omi_ini[alpha, i, j, k] = omegadim * tref
            wad_ini[alpha, i, j, k] = adim / rhoref / uref^2 / tref
        end
    else
        println(
            "Warning: MS-GWaM's steady-state mode currently ignores non-orographic initializations!",
        )
        println("")
    end

    # Add orographic wave modes.
    if source_mode isa OrographicSource
        activate_orographic_source!(
            state,
            omi_ini,
            wnk_ini,
            wnl_ini,
            wnm_ini,
            wad_ini,
        )
    end

    # Set initial spectral extents (these will be overwritten in the loop).
    dk_ini_nd = 0.0
    dl_ini_nd = 0.0
    dm_ini_nd = 0.0

    # Set vertical index bounds.
    kmin = ko == 0 ? k0 - 1 : k0
    kmax = k1

    # Loop over all grid cells with ray volumes.
    @ivy for k in kmin:kmax, j in j0:j1, i in i0:i1
        r = 0
        s = 0

        # Loop over all ray volumes within a spatial cell.
        for ix in 1:nrx,
            ik in 1:nrk,
            jy in 1:nry,
            jl in 1:nrl,
            kz in 1:nrz,
            km in 1:nrm,
            alpha in 1:wave_modes

            # Set ray-volume indices.
            if ko == 0 && k == k0 - 1
                s += 1

                # Set surface indices.
                surface_indices.ixs[s] = ix
                surface_indices.jys[s] = jy
                surface_indices.kzs[s] = kz
                surface_indices.iks[s] = ik
                surface_indices.jls[s] = jl
                surface_indices.kms[s] = km
                surface_indices.alphas[s] = alpha

                # Set surface ray-volume index.
                if wad_ini[alpha, i, j, k] == 0.0
                    surface_indices.rs[s, i, j] = -1
                    continue
                else
                    r += 1
                    surface_indices.rs[s, i, j] = r
                end
            else
                if wad_ini[alpha, i, j, k] == 0.0
                    continue
                end
                r += 1
            end

            # Set ray-volume positions.
            rays.x[r, i, j, k] = (x[i] - 0.5 * dx + (ix - 0.5) * dx / nrx)
            rays.y[r, i, j, k] = (y[j] - 0.5 * dy + (jy - 0.5) * dy / nry)
            rays.z[r, i, j, k] = (
                zc[i, j, k] - 0.5 * jac[i, j, k] * dz +
                (kz - 0.5) * jac[i, j, k] * dz / nrz
            )

            xr = rays.x[r, i, j, k]
            yr = rays.y[r, i, j, k]
            zr = rays.z[r, i, j, k]

            # Check if ray volume is too low.
            if zr < -dz
                error(
                    "Error in initialize_rays!: Ray volume",
                    r,
                    "at",
                    i,
                    j,
                    k,
                    "is too low!",
                )
            end

            # Compute local stratification.
            n2r = interpolate_stratification(zr, state, N2())

            # Set spatial extents.
            rays.dxray[r, i, j, k] = dx / nrx
            rays.dyray[r, i, j, k] = dy / nry
            rays.dzray[r, i, j, k] = jac[i, j, k] * dz / nrz

            wnk0 = wnk_ini[alpha, i, j, k]
            wnl0 = wnl_ini[alpha, i, j, k]
            wnm0 = wnm_ini[alpha, i, j, k]

            # Ensure correct wavenumber extents.
            wnh0 = sqrt(wnk0^2 + wnl0^2)

            if x_size == 1
                dk_ini_nd = 0.0
            else
                dk_ini_nd = dkr_factor[alpha] * wnh0

                if dk_ini_nd <= 0.0
                    error(
                        "Error in initialize_rays!: dk_ini_nd <= 0 for mode ",
                        alpha,
                        " with x_size > 1.",
                    )
                end
            end

            if y_size == 1
                dl_ini_nd = 0.0
            else
                dl_ini_nd = dlr_factor[alpha] * wnh0

                if dl_ini_nd <= 0.0
                    error(
                        "Error in initialize_rays!: dl_ini_nd <= 0 for mode ",
                        alpha,
                        " with y_size > 1.",
                    )
                end
            end

            if wnm0 == 0.0
                error(
                    "Error in initialize_rays!: wnm0 = 0 for mode ",
                    alpha,
                    ".",
                )
            end

            dm_ini_nd = dmr_factor[alpha] * abs(wnm0)

            if dm_ini_nd <= 0.0
                error(
                    "Error in initialize_rays!: dm_ini_nd <= 0 for mode ",
                    alpha,
                    ".",
                )
            end
            # Set ray-volume wavenumbers.
            rays.k[r, i, j, k] =
                (wnk0 - 0.5 * dk_ini_nd + (ik - 0.5) * dk_ini_nd / nrk)
            rays.l[r, i, j, k] =
                (wnl0 - 0.5 * dl_ini_nd + (jl - 0.5) * dl_ini_nd / nrl)
            rays.m[r, i, j, k] =
                (wnm0 - 0.5 * dm_ini_nd + (km - 0.5) * dm_ini_nd / nrm)

            # Set spectral extents.
            rays.dkray[r, i, j, k] = dk_ini_nd / nrk
            rays.dlray[r, i, j, k] = dl_ini_nd / nrl
            rays.dmray[r, i, j, k] = dm_ini_nd / nrm

            # Set spectral volume.
            pspvol = dm_ini_nd
            if x_size > 1
                pspvol = pspvol * dk_ini_nd
            end
            if y_size > 1
                pspvol = pspvol * dl_ini_nd
            end

            # Set phase-space wave-action density.
            rays.dens[r, i, j, k] = wad_ini[alpha, i, j, k] / pspvol

            # Interpolate winds to ray-volume position.
            uxr = interpolate_mean_flow(xr, yr, zr, state, U())
            vyr = interpolate_mean_flow(xr, yr, zr, state, V())
            wzr = interpolate_mean_flow(xr, yr, zr, state, W())

            wnrk = rays.k[r, i, j, k]
            wnrl = rays.l[r, i, j, k]
            wnrm = rays.m[r, i, j, k]
            wnrh = sqrt(wnrk^2 + wnrl^2)
            omir = omi_ini[alpha, i, j, k]

            # Compute maximum group velocities.
            cgirx = wnrk * (n2r - omir^2) / (omir * (wnrh^2 + wnrm^2))
            if abs(uxr + cgirx) > abs(cgx_max[])
                cgx_max[] = abs(uxr + cgirx)
            end
            cgiry = wnrl * (n2r - omir^2) / (omir * (wnrh^2 + wnrm^2))
            if abs(vyr + cgiry) > abs(cgy_max[])
                cgy_max[] = abs(vyr + cgiry)
            end
            cgirz = -wnrm * (omir^2 - fc^2) / (omir * (wnrh^2 + wnrm^2))
            if abs(wzr + cgirz) > abs(cgz_max[i, j, k])
                cgz_max[i, j, k] = max(cgz_max[i, j, k], abs(wzr + cgirz))
            end
        end

        # Set ray-volume count.
        nray[i, j, k] = r
        if r > nray_wrk
            error(
                "Error in initialize_rays!: nray",
                [i, j, k],
                " > nray_wrk =",
                nray_wrk,
            )
        end

        # Check if surface ray-volume count is correct.
        if ko == 0 && k == k0 - 1
            if s != n_sfc
                error(
                    "Error in initialize_rays!: s =",
                    s,
                    "/= n_sfc =",
                    n_sfc,
                    "at (i, j, k) = ",
                    (i, j, k),
                )
            end
        end
    end

    # Compute global ray-volume count.
    @ivy local_sum = sum(nray[i0:i1, j0:j1, kmin:kmax])
    global_sum = MPI.Allreduce(local_sum, +, comm)

    # Print information.
    if master
        println("MS-GWaM:")
        println("Global ray-volume count: ", global_sum)
        println("Maximum number of ray volumes per cell: ", nray_max)
        println("")
    end

    return
end


function initialize_rays!(
    state::State,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
    triad_mode::Union{Triad2D},
    ray_volume_ini::GaussianDist,
    source_mode::Union{NoRaySource, OrographicSource},
)
    (; x_size, y_size) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere

    (;
        nrx,
        nry,
        nrz,
        nrk,
        nrl,
        nrm,
        wave_modes,
        dkr_factor,
        dlr_factor,
        dmr_factor,
        initial_wave_field,
    ) = state.namelists.wkb

    (; lref, tref, rhoref, uref) = state.constants

    (;
        comm,
        master,
        nxx,
        nyy,
        nzz,
        ko,
        i0,
        i1,
        j0,
        j1,
        k0,
        k1,
    ) = state.domain

    (; dx, dy, dz, x, y, zc, jac) = state.grid

    (;
        nray_max,
        nray_wrk,
        n_sfc,
        nray,
        rays,
        surface_indices,
        cgx_max,
        cgy_max,
        cgz_max,
    ) = state.wkb

    (; m_sigma_cutoff) = state.namelists.triad

    if x_size == 1 && nrk != 1
        error(
            "Error in initialize_rays!: nrk must be 1 when x_size == 1. ",
            "Otherwise identical zero-width k ray volumes are initialized.",
        )
    end

    if y_size == 1 && nrl != 1
        error(
            "Error in initialize_rays!: nrl must be 1 when y_size == 1. ",
            "Otherwise identical zero-width l ray volumes are initialized.",
        )
    end

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    # Initialize arrays for the initial wave properties.
    omi_ini = zeros(wave_modes, nxx, nyy, nzz)
    wnk_ini = zeros(wave_modes, nxx, nyy, nzz)
    wnl_ini = zeros(wave_modes, nxx, nyy, nzz)
    wnm_ini = zeros(wave_modes, nxx, nyy, nzz)
    wad_ini = zeros(wave_modes, nxx, nyy, nzz)

    # ------------------------------------------------------------------
    # Initialize the carrier-wave properties at Eulerian cell centres.
    #
    # The carrier wavenumbers and frequencies continue to be initialized
    # in the original way. The wave-action density of the interior rays
    # will subsequently be re-evaluated at the actual physical sub-ray
    # positions.
    # ------------------------------------------------------------------

    if wkb_mode != SteadyState()
        for k in k0:k1, j in j0:j1, i in i0:i1, alpha in 1:wave_modes
            (kdim, ldim, mdim, omegadim, adim) =
                initial_wave_field(
                    alpha,
                    x[i] * lref,
                    y[j] * lref,
                    zc[i, j, k] * lref,
                )

            wnk_ini[alpha, i, j, k] = kdim * lref
            wnl_ini[alpha, i, j, k] = ldim * lref
            wnm_ini[alpha, i, j, k] = mdim * lref
            omi_ini[alpha, i, j, k] = omegadim * tref
            wad_ini[alpha, i, j, k] = adim / rhoref / uref^2 / tref
        end
    else
        if master
            println(
                "Warning: MS-GWaM's steady-state mode currently ignores non-orographic initializations!",
            )
            println("")
        end
    end

    # Set the initial properties of unresolved orographic gravity waves
    # in the surface launch layer k0 - 1.
    if source_mode isa OrographicSource
        activate_orographic_source!(
            state,
            omi_ini,
            wnk_ini,
            wnl_ini,
            wnm_ini,
            wad_ini,
        )
    end

    # Gaussian weights for the spectral m sub-rays.
    #
    # The complete spectral interval represents
    #
    #     m0 - m_sigma_cutoff * sigma_m
    #     <= m <=
    #     m0 + m_sigma_cutoff * sigma_m.
    #
    # The weights are normalized such that
    #
    #     sum(m_weights) / nrm = 1,
    #
    # so that the total initialized wave action is unchanged.

    if m_sigma_cutoff <= 0.0
        error("m_sigma_cutoff must be positive.")
    end

    m_weights = zeros(nrm)

    for km in 1:nrm
        xi = -m_sigma_cutoff + (km - 0.5) * 2 * m_sigma_cutoff / nrm
        m_weights[km] = exp(-0.5 * xi^2)
    end

    m_weights .*= nrm / sum(m_weights)

    # Initialize spectral ray-volume extents.
    dk_ini_nd = 0.0
    dl_ini_nd = 0.0
    dm_ini_nd = 0.0

    # Include the surface launch layer on the lowest MPI subdomain.
    kmin = ko == 0 ? k0 - 1 : k0
    kmax = k1

    # ------------------------------------------------------------------
    # Initialize ray volumes.
    # ------------------------------------------------------------------

    @ivy for k in kmin:kmax, j in j0:j1, i in i0:i1

        r = 0
        s = 0

        for ix in 1:nrx,
            ik in 1:nrk,
            jy in 1:nry,
            jl in 1:nrl,
            kz in 1:nrz,
            km in 1:nrm,
            alpha in 1:wave_modes

            # ----------------------------------------------------------
            # Candidate physical sub-ray centre.
            #
            # This is computed before deciding whether the ray exists,
            # because the initial wave-action density is now evaluated
            # at the actual sub-ray centre rather than the Eulerian
            # cell centre.
            # ----------------------------------------------------------

            xr = x[i] - 0.5 * dx + (ix - 0.5) * dx / nrx
            yr = y[j] - 0.5 * dy + (jy - 0.5) * dy / nry
            zr =
                zc[i, j, k] -
                0.5 * jac[i, j, k] * dz +
                (kz - 0.5) * jac[i, j, k] * dz / nrz

            surface_ray = ko == 0 && k == k0 - 1

            # ----------------------------------------------------------
            # Determine whether this ray volume should be initialized.
            # ----------------------------------------------------------

            if surface_ray

                # Surface rays are treated by the existing orographic
                # source mechanism.
                s += 1

                surface_indices.ixs[s] = ix
                surface_indices.jys[s] = jy
                surface_indices.kzs[s] = kz
                surface_indices.iks[s] = ik
                surface_indices.jls[s] = jl
                surface_indices.kms[s] = km
                surface_indices.alphas[s] = alpha

                wad_r = wad_ini[alpha, i, j, k]

                if wad_r == 0.0
                    surface_indices.rs[s, i, j] = -1
                    continue
                end

                r += 1
                surface_indices.rs[s, i, j] = r

            else

                # Steady-state mode does not contain an independently
                # initialized interior wave packet.
                if wkb_mode == SteadyState()
                    continue
                end

                # ------------------------------------------------------
                # Evaluate the physical-space wave-action density at
                # the actual sub-ray centre.
                #
                # Only the action is re-evaluated here. The nominal
                # carrier wave numbers remain those initialized at the
                # Eulerian cell centre.
                # ------------------------------------------------------

                (_, _, _, _, adim_r) =
                    initial_wave_field(
                        alpha,
                        xr * lref,
                        yr * lref,
                        zr * lref,
                    )

                wad_r = adim_r / rhoref / uref^2 / tref

                # This test must use the ray-centre value rather than
                # wad_ini at the Eulerian cell centre. This is important
                # for explicitly truncated Gaussian packets.
                if wad_r == 0.0
                    continue
                end

                r += 1
            end

            if r > nray_wrk
                error(
                    "Error in initialize_rays!: Number of ray volumes exceeds nray_wrk = ",
                    nray_wrk,
                    " at grid cell ",
                    (i, j, k),
                )
            end

            # ----------------------------------------------------------
            # Set physical ray-volume position.
            # ----------------------------------------------------------

            rays.x[r, i, j, k] = xr
            rays.y[r, i, j, k] = yr
            rays.z[r, i, j, k] = zr

            if zr < -dz
                error(
                    "Error in initialize_rays!: Ray volume ",
                    r,
                    " at ",
                    (i, j, k),
                    " is too low!",
                )
            end

            # Interpolate stratification at the ray centre.
            n2r = interpolate_stratification(zr, state, N2())

            # ----------------------------------------------------------
            # Set physical ray-volume extent.
            # ----------------------------------------------------------

            rays.dxray[r, i, j, k] = dx / nrx
            rays.dyray[r, i, j, k] = dy / nry
            rays.dzray[r, i, j, k] = jac[i, j, k] * dz / nrz

            # ----------------------------------------------------------
            # Get nominal carrier wave numbers.
            # ----------------------------------------------------------

            wnk0 = wnk_ini[alpha, i, j, k]
            wnl0 = wnl_ini[alpha, i, j, k]
            wnm0 = wnm_ini[alpha, i, j, k]

            # ----------------------------------------------------------
            # Compute spectral ray-volume extents.
            #
            # dkr_factor, dlr_factor and dmr_factor are mode-dependent
            # vectors of length wave_modes.
            # ----------------------------------------------------------

            dkr_factor_alpha = dkr_factor[alpha]
            dlr_factor_alpha = dlr_factor[alpha]
            dmr_factor_alpha = dmr_factor[alpha]

            wnh0 = sqrt(wnk0^2 + wnl0^2)

            if x_size == 1
                dk_ini_nd = 0.0
            else
                dk_ini_nd = dkr_factor_alpha * wnh0

                if dk_ini_nd <= 0.0
                    error(
                        "Error in initialize_rays!: dk_ini_nd <= 0 for mode ",
                        alpha,
                        " with x_size > 1.",
                    )
                end
            end

            if y_size == 1
                dl_ini_nd = 0.0
            else
                dl_ini_nd = dlr_factor_alpha * wnh0

                if dl_ini_nd <= 0.0
                    error(
                        "Error in initialize_rays!: dl_ini_nd <= 0 for mode ",
                        alpha,
                        " with y_size > 1.",
                    )
                end
            end

            if wnm0 == 0.0
                error(
                    "Error in initialize_rays!: wnm0 = 0 for mode ",
                    alpha,
                    ".",
                )
            end

            dm_ini_nd = dmr_factor_alpha * abs(wnm0)

            if dm_ini_nd <= 0.0
                error(
                    "Error in initialize_rays!: dm_ini_nd <= 0 for mode ",
                    alpha,
                    ".",
                )
            end
            # ----------------------------------------------------------
            # Set spectral ray-volume position.
            #
            # For m, the nrm equal-width sub-rays span the complete
            # interval m0 ± 2.5 sigma_m.
            # ----------------------------------------------------------

            rays.k[r, i, j, k] =
                wnk0 -
                0.5 * dk_ini_nd +
                (ik - 0.5) * dk_ini_nd / nrk

            rays.l[r, i, j, k] =
                wnl0 -
                0.5 * dl_ini_nd +
                (jl - 0.5) * dl_ini_nd / nrl

            rays.m[r, i, j, k] =
                wnm0 -
                0.5 * dm_ini_nd +
                (km - 0.5) * dm_ini_nd / nrm

            # ----------------------------------------------------------
            # Set spectral ray-volume extent.
            # ----------------------------------------------------------

            rays.dkray[r, i, j, k] = dk_ini_nd / nrk
            rays.dlray[r, i, j, k] = dl_ini_nd / nrl
            rays.dmray[r, i, j, k] = dm_ini_nd / nrm

            # ----------------------------------------------------------
            # Compute the complete spectral volume represented by the
            # initialized packet at this physical position.
            # ----------------------------------------------------------

            pspvol = dm_ini_nd

            if x_size > 1
                pspvol *= dk_ini_nd
            end

            if y_size > 1
                pspvol *= dl_ini_nd
            end

            # ----------------------------------------------------------
            # Set phase-space wave-action density.
            #
            # Interior analytic packets:
            #   Gaussian distribution in m.
            #
            # Orographic surface source:
            #   retain the original uniform distribution.
            # ----------------------------------------------------------

            if surface_ray
                rays.dens[r, i, j, k] = wad_r / pspvol
            else
                rays.dens[r, i, j, k] =
                    wad_r / pspvol * m_weights[km]
            end

            # ----------------------------------------------------------
            # Compute maximum group velocities for WKB CFL estimate.
            # ----------------------------------------------------------

            uxr = interpolate_mean_flow(xr, yr, zr, state, U())
            vyr = interpolate_mean_flow(xr, yr, zr, state, V())
            wzr = interpolate_mean_flow(xr, yr, zr, state, W())

            wnrk = rays.k[r, i, j, k]
            wnrl = rays.l[r, i, j, k]
            wnrm = rays.m[r, i, j, k]

            wnrh = sqrt(wnrk^2 + wnrl^2)
            omir = omi_ini[alpha, i, j, k]

            # Compute maximum zonal group velocity.
            cgirx =
                wnrk *
                (n2r - omir^2) /
                (omir * (wnrh^2 + wnrm^2))

            if abs(uxr + cgirx) > abs(cgx_max[])
                cgx_max[] = abs(uxr + cgirx)
            end

            # Compute maximum meridional group velocity.
            cgiry =
                wnrl *
                (n2r - omir^2) /
                (omir * (wnrh^2 + wnrm^2))

            if abs(vyr + cgiry) > abs(cgy_max[])
                cgy_max[] = abs(vyr + cgiry)
            end

            # Compute maximum vertical group velocity.
            cgirz =
                -wnrm *
                (omir^2 - fc^2) /
                (omir * (wnrh^2 + wnrm^2))

            if abs(wzr + cgirz) > abs(cgz_max[i, j, k])
                cgz_max[i, j, k] =
                    max(cgz_max[i, j, k], abs(wzr + cgirz))
            end
        end

        # Set number of ray volumes in this Eulerian grid cell.
        nray[i, j, k] = r

        if r > nray_wrk
            error(
                "Error in initialize_rays!: nray = ",
                r,
                " > nray_wrk = ",
                nray_wrk,
            )
        end

        # Verify number of surface ray-volume slots.
        if ko == 0 && k == k0 - 1
            if s != n_sfc
                error(
                    "Error in initialize_rays!: Number of surface ray volumes ",
                    s,
                    " != n_sfc = ",
                    n_sfc,
                )
            end
        end
    end

    # ------------------------------------------------------------------
    # Global ray-volume count.
    # ------------------------------------------------------------------

    @ivy local_sum = sum(nray[i0:i1, j0:j1, kmin:kmax])
    global_sum = MPI.Allreduce(local_sum, +, comm)

    if master
        println("MS-GWaM:")
        println("Global ray-volume count: ", global_sum)
        println("Maximum number of ray volumes per cell: ", nray_max)
        println("")
    end

    return
end


# =============================================================================
# SpatialGaussianDist
#
# Physical-space initialization:
#   Wave-action density is evaluated independently at every physical ray-volume
#   centre.
#
# Spectral initialization:
#   Uniform distribution over the spectral ray volumes.
#
# Supported triad modes:
#   NoTriad, Triad2D
# =============================================================================

function initialize_rays!(
    state::State,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
    triad_mode::Union{NoTriad, Triad2D},
    ray_volume_ini::SpatialGaussianDist,
    source_mode::Union{NoRaySource, OrographicSource},
)
    (; x_size, y_size) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere

    (;
        nrx,
        nry,
        nrz,
        nrk,
        nrl,
        nrm,
        wave_modes,
        dkr_factor,
        dlr_factor,
        dmr_factor,
        initial_wave_field,
    ) = state.namelists.wkb

    (; lref, tref, rhoref, uref) = state.constants

    (;
        comm,
        master,
        nxx,
        nyy,
        nzz,
        ko,
        i0,
        i1,
        j0,
        j1,
        k0,
        k1,
    ) = state.domain

    (; dx, dy, dz, x, y, zc, jac) = state.grid

    (;
        nray_max,
        nray_wrk,
        n_sfc,
        nray,
        rays,
        surface_indices,
        cgx_max,
        cgy_max,
        cgz_max,
    ) = state.wkb

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    # Initialize arrays for the initial wave properties.
    omi_ini = zeros(wave_modes, nxx, nyy, nzz)
    wnk_ini = zeros(wave_modes, nxx, nyy, nzz)
    wnl_ini = zeros(wave_modes, nxx, nyy, nzz)
    wnm_ini = zeros(wave_modes, nxx, nyy, nzz)
    wad_ini = zeros(wave_modes, nxx, nyy, nzz)

    # Carrier-wave properties are evaluated at the Eulerian cell centres.
    # The interior wave-action density itself is subsequently re-evaluated
    # at each physical ray-volume centre.
    if wkb_mode != SteadyState()
        for k in k0:k1, j in j0:j1, i in i0:i1, alpha in 1:wave_modes
            (kdim, ldim, mdim, omegadim, adim) = initial_wave_field(
                alpha,
                x[i] * lref,
                y[j] * lref,
                zc[i, j, k] * lref,
            )

            wnk_ini[alpha, i, j, k] = kdim * lref
            wnl_ini[alpha, i, j, k] = ldim * lref
            wnm_ini[alpha, i, j, k] = mdim * lref
            omi_ini[alpha, i, j, k] = omegadim * tref
            wad_ini[alpha, i, j, k] = adim / rhoref / uref^2 / tref
        end
    else
        if master
            println(
                "Warning: MS-GWaM's steady-state mode currently ignores non-orographic initializations!",
            )
            println("")
        end
    end

    # Add orographic wave modes.
    if source_mode isa OrographicSource
        activate_orographic_source!(
            state,
            omi_ini,
            wnk_ini,
            wnl_ini,
            wnm_ini,
            wad_ini,
        )
    end

    dk_ini_nd = 0.0
    dl_ini_nd = 0.0
    dm_ini_nd = 0.0

    kmin = ko == 0 ? k0 - 1 : k0
    kmax = k1

    @ivy for k in kmin:kmax, j in j0:j1, i in i0:i1
        r = 0
        s = 0

        for ix in 1:nrx,
            ik in 1:nrk,
            jy in 1:nry,
            jl in 1:nrl,
            kz in 1:nrz,
            km in 1:nrm,
            alpha in 1:wave_modes

            # ----------------------------------------------------------
            # Candidate physical ray-volume centre.
            #
            # This must be computed before the existence test because
            # the physical packet is sampled at the actual ray centre.
            # ----------------------------------------------------------

            xr = x[i] - 0.5 * dx + (ix - 0.5) * dx / nrx
            yr = y[j] - 0.5 * dy + (jy - 0.5) * dy / nry
            zr =
                zc[i, j, k] -
                0.5 * jac[i, j, k] * dz +
                (kz - 0.5) * jac[i, j, k] * dz / nrz

            surface_ray = ko == 0 && k == k0 - 1

            # ----------------------------------------------------------
            # Determine whether the candidate ray exists.
            # ----------------------------------------------------------

            if surface_ray
                s += 1

                surface_indices.ixs[s] = ix
                surface_indices.jys[s] = jy
                surface_indices.kzs[s] = kz
                surface_indices.iks[s] = ik
                surface_indices.jls[s] = jl
                surface_indices.kms[s] = km
                surface_indices.alphas[s] = alpha

                # Orographic launch retains its original cell-based value.
                wad_r = wad_ini[alpha, i, j, k]

                if wad_r == 0.0
                    surface_indices.rs[s, i, j] = -1
                    continue
                end

                r += 1
                surface_indices.rs[s, i, j] = r

            else
                # Steady-state mode has no independently initialized
                # interior wave packet.
                if wkb_mode == SteadyState()
                    continue
                end

                # Evaluate action directly at the physical ray centre.
                (_, _, _, _, adim_r) = initial_wave_field(
                    alpha,
                    xr * lref,
                    yr * lref,
                    zr * lref,
                )

                wad_r = adim_r / rhoref / uref^2 / tref

                # The existence test must use the ray-centre value.
                if wad_r == 0.0
                    continue
                end

                r += 1
            end

            if r > nray_wrk
                error(
                    "Error in initialize_rays!: Number of ray volumes exceeds nray_wrk = ",
                    nray_wrk,
                    " at grid cell ",
                    (i, j, k),
                )
            end

            # Set physical ray-volume position.
            rays.x[r, i, j, k] = xr
            rays.y[r, i, j, k] = yr
            rays.z[r, i, j, k] = zr

            if zr < -dz
                error(
                    "Error in initialize_rays!: Ray volume ",
                    r,
                    " at ",
                    (i, j, k),
                    " is too low!",
                )
            end

            n2r = interpolate_stratification(zr, state, N2())

            # Physical extents.
            rays.dxray[r, i, j, k] = dx / nrx
            rays.dyray[r, i, j, k] = dy / nry
            rays.dzray[r, i, j, k] = jac[i, j, k] * dz / nrz

            # Carrier wave numbers.
            wnk0 = wnk_ini[alpha, i, j, k]
            wnl0 = wnl_ini[alpha, i, j, k]
            wnm0 = wnm_ini[alpha, i, j, k]

            wnh0 = sqrt(wnk0^2 + wnl0^2)

            # Spectral packet widths.
            if x_size == 1
                dk_ini_nd = 0.0
            else
                dk_ini_nd = dkr_factor[alpha] * wnh0

                if dk_ini_nd <= 0.0
                    error(
                        "Error in initialize_rays!: dk_ini_nd <= 0 for mode ",
                        alpha,
                        " with x_size > 1.",
                    )
                end
            end

            if y_size == 1
                dl_ini_nd = 0.0
            else
                dl_ini_nd = dlr_factor[alpha] * wnh0

                if dl_ini_nd <= 0.0
                    error(
                        "Error in initialize_rays!: dl_ini_nd <= 0 for mode ",
                        alpha,
                        " with y_size > 1.",
                    )
                end
            end

            if wnm0 == 0.0
                error(
                    "Error in initialize_rays!: wnm0 = 0 for mode ",
                    alpha,
                    ".",
                )
            end

            dm_ini_nd = dmr_factor[alpha] * abs(wnm0)

            if dm_ini_nd <= 0.0
                error(
                    "Error in initialize_rays!: dm_ini_nd <= 0 for mode ",
                    alpha,
                    ".",
                )
            end

            # Spectral ray-volume positions.
            rays.k[r, i, j, k] =
                wnk0 -
                0.5 * dk_ini_nd +
                (ik - 0.5) * dk_ini_nd / nrk

            rays.l[r, i, j, k] =
                wnl0 -
                0.5 * dl_ini_nd +
                (jl - 0.5) * dl_ini_nd / nrl

            rays.m[r, i, j, k] =
                wnm0 -
                0.5 * dm_ini_nd +
                (km - 0.5) * dm_ini_nd / nrm

            # Spectral ray-volume extents.
            rays.dkray[r, i, j, k] = dk_ini_nd / nrk
            rays.dlray[r, i, j, k] = dl_ini_nd / nrl
            rays.dmray[r, i, j, k] = dm_ini_nd / nrm

            # Complete spectral volume.
            pspvol = dm_ini_nd

            if x_size > 1
                pspvol *= dk_ini_nd
            end

            if y_size > 1
                pspvol *= dl_ini_nd
            end

            # Uniform spectral distribution, but with physical-space
            # action evaluated at the actual ray-volume centre.
            rays.dens[r, i, j, k] = wad_r / pspvol

            # ----------------------------------------------------------
            # Group velocities for WKB CFL estimate.
            # ----------------------------------------------------------

            uxr = interpolate_mean_flow(xr, yr, zr, state, U())
            vyr = interpolate_mean_flow(xr, yr, zr, state, V())
            wzr = interpolate_mean_flow(xr, yr, zr, state, W())

            wnrk = rays.k[r, i, j, k]
            wnrl = rays.l[r, i, j, k]
            wnrm = rays.m[r, i, j, k]
            wnrh = sqrt(wnrk^2 + wnrl^2)
            omir = omi_ini[alpha, i, j, k]

            cgirx =
                wnrk *
                (n2r - omir^2) /
                (omir * (wnrh^2 + wnrm^2))

            if abs(uxr + cgirx) > abs(cgx_max[])
                cgx_max[] = abs(uxr + cgirx)
            end

            cgiry =
                wnrl *
                (n2r - omir^2) /
                (omir * (wnrh^2 + wnrm^2))

            if abs(vyr + cgiry) > abs(cgy_max[])
                cgy_max[] = abs(vyr + cgiry)
            end

            cgirz =
                -wnrm *
                (omir^2 - fc^2) /
                (omir * (wnrh^2 + wnrm^2))

            if abs(wzr + cgirz) > abs(cgz_max[i, j, k])
                cgz_max[i, j, k] =
                    max(cgz_max[i, j, k], abs(wzr + cgirz))
            end
        end

        nray[i, j, k] = r

        if r > nray_wrk
            error(
                "Error in initialize_rays!: nray = ",
                r,
                " > nray_wrk = ",
                nray_wrk,
            )
        end

        if ko == 0 && k == k0 - 1
            if s != n_sfc
                error(
                    "Error in initialize_rays!: Number of surface ray volumes ",
                    s,
                    " != n_sfc = ",
                    n_sfc,
                    " at grid cell ",
                    (i, j, k),
                )
            end
        end
    end

    @ivy local_sum = sum(nray[i0:i1, j0:j1, kmin:kmax])
    global_sum = MPI.Allreduce(local_sum, +, comm)

    if master
        println("MS-GWaM:")
        println("Global ray-volume count: ", global_sum)
        println("Maximum number of ray volumes per cell: ", nray_max)
        println("")
    end

    return
end

# =============================================================================
# ContinuousSpectralSource
#
# Initialize a continuously forced lower-boundary spectral source.
#
# Only the artificial launch layer k0 - 1 is populated initially.
# The physical domain k0:k1 is initially empty.
#
# The source spectrum is obtained from initial_wave_field at the lower physical
# boundary and distributed in m using the same Gaussian discretization as
# GaussianDist:
#
#   - dmr_factor controls the complete spectral width,
#   - m_sigma_cutoff specifies the represented Gaussian interval,
#   - nrm gives the number of spectral sub-rays in m,
#   - nrz gives the number of vertically staggered source ray volumes.
#
# The Gaussian weights are normalized such that the total wave action
# represented by the source layer is independent of nrm. Likewise, subdivision
# into nrz physical ray volumes changes only the vertical resolution of the
# source reservoir, not its phase-space wave-action density.
#
# The source slots are stored in surface_indices and will subsequently be
# handled by activate_continuous_spectral_source!.
# =============================================================================

function initialize_rays!(
    state::State,
    wkb_mode::Union{SingleColumn, MultiColumn},
    triad_mode::Union{NoTriad, Triad2D},
    ray_volume_ini::GaussianDist,
    source_mode::ContinuousSpectralSource,
)
    (; x_size, y_size) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere

    (;
        nrx,
        nry,
        nrz,
        nrk,
        nrl,
        nrm,
        wave_modes,
        dkr_factor,
        dlr_factor,
        dmr_factor,
        branch,
        initial_wave_field,
    ) = state.namelists.wkb

    (; m_sigma_cutoff) = state.namelists.triad
    (; lref, tref, rhoref, uref) = state.constants

    (;
        comm,
        master,
        nxx,
        nyy,
        ko,
        i0,
        i1,
        j0,
        j1,
        k0,
        k1,
    ) = state.domain

    (; dx, dy, dz, x, y, zc, zctilde, jac) = state.grid

    (;
        nray_wrk,
        n_sfc,
        nray,
        rays,
        surface_indices,
        cgx_max,
        cgy_max,
        cgz_max,
    ) = state.wkb

    # ------------------------------------------------------------------
    # Validate source configuration.
    # ------------------------------------------------------------------

    if nrz < 1
        error(
            "Error in initialize_rays!: nrz must be >= 1 for ",
            "ContinuousSpectralSource.",
        )
    end

    if nrm < 1
        error(
            "Error in initialize_rays!: nrm must be >= 1 for ",
            "ContinuousSpectralSource.",
        )
    end

    if x_size == 1 && nrk != 1
        error(
            "Error in initialize_rays!: nrk must be 1 when x_size == 1. ",
            "Otherwise identical zero-width k ray volumes are initialized.",
        )
    end

    if y_size == 1 && nrl != 1
        error(
            "Error in initialize_rays!: nrl must be 1 when y_size == 1. ",
            "Otherwise identical zero-width l ray volumes are initialized.",
        )
    end

    if m_sigma_cutoff <= 0.0
        error(
            "Error in initialize_rays!: m_sigma_cutoff must be positive ",
            "for ContinuousSpectralSource.",
        )
    end

    # Set nondimensional Coriolis parameter.
    fc = coriolis_frequency * tref

    # ------------------------------------------------------------------
    # Construct Gaussian weights in m.
    #
    # xi spans the interval
    #
    #     -m_sigma_cutoff <= xi <= m_sigma_cutoff
    #
    # at the centres of the nrm spectral ray volumes.
    #
    # The normalization
    #
    #     sum(m_weights) / nrm = 1
    #
    # ensures that the complete Gaussian packet represents the wave action
    # supplied by initial_wave_field.
    # ------------------------------------------------------------------

    m_weights = zeros(nrm)

    for km in 1:nrm
        xi =
            -m_sigma_cutoff +
            (km - 0.5) * 2.0 * m_sigma_cutoff / nrm

        m_weights[km] = exp(-0.5 * xi^2)
    end

    m_weights .*= nrm / sum(m_weights)

    # ------------------------------------------------------------------
    # Source carrier-wave properties.
    #
    # initial_wave_field defines the spectrum at the lower physical
    # boundary. The actual source ray volumes are stored below this
    # boundary, in k0 - 1.
    #
    # These arrays are needed only on the lowest vertical MPI subdomain.
    # ------------------------------------------------------------------

    wnk_src = zeros(wave_modes, nxx, nyy)
    wnl_src = zeros(wave_modes, nxx, nyy)
    wnm_src = zeros(wave_modes, nxx, nyy)
    wad_src = zeros(wave_modes, nxx, nyy)

    if ko == 0
        for j in j0:j1, i in i0:i1, alpha in 1:wave_modes

            # Interface between the artificial source layer k0 - 1
            # and the first physical layer k0.
            z_source = zctilde[i, j, k0 - 1]

            (kdim, ldim, mdim, _, adim) = initial_wave_field(
                alpha,
                x[i] * lref,
                y[j] * lref,
                z_source * lref,
            )

            wnk_src[alpha, i, j] = kdim * lref
            wnl_src[alpha, i, j] = ldim * lref
            wnm_src[alpha, i, j] = mdim * lref

            wad_src[alpha, i, j] =
                adim / rhoref / uref^2 / tref
        end
    end

    # ------------------------------------------------------------------
    # Initialize only the source layer k0 - 1.
    #
    # All physical cells k0:k1 remain empty at t = 0.
    # ------------------------------------------------------------------

    if ko == 0

        k = k0 - 1

        @ivy for j in j0:j1, i in i0:i1

            r = 0
            s = 0

            for ix in 1:nrx,
                ik in 1:nrk,
                jy in 1:nry,
                jl in 1:nrl,
                kz in 1:nrz,
                km in 1:nrm,
                alpha in 1:wave_modes

                # ------------------------------------------------------
                # Register persistent source slot.
                #
                # The complete tuple
                #
                #   (ix, jy, kz, ik, jl, km, alpha)
                #
                # identifies the source ray. In particular, kz is retained
                # so that nrz > 1 gives vertically staggered source rays.
                # ------------------------------------------------------

                s += 1

                surface_indices.ixs[s] = ix
                surface_indices.jys[s] = jy
                surface_indices.kzs[s] = kz
                surface_indices.iks[s] = ik
                surface_indices.jls[s] = jl
                surface_indices.kms[s] = km
                surface_indices.alphas[s] = alpha

                wad0 = wad_src[alpha, i, j]

                if wad0 == 0.0
                    surface_indices.rs[s, i, j] = -1
                    continue
                end

                r += 1

                if r > nray_wrk
                    error(
                        "Error in initialize_rays!: Number of source ray ",
                        "volumes exceeds nray_wrk = ",
                        nray_wrk,
                        " at source cell ",
                        (i, j, k),
                    )
                end

                surface_indices.rs[s, i, j] = r

                # ------------------------------------------------------
                # Physical position.
                #
                # nrz equal-width ray volumes uniformly subdivide the
                # artificial source layer.
                # ------------------------------------------------------

                rays.x[r, i, j, k] =
                    x[i] -
                    0.5 * dx +
                    (ix - 0.5) * dx / nrx

                rays.y[r, i, j, k] =
                    y[j] -
                    0.5 * dy +
                    (jy - 0.5) * dy / nry

                rays.z[r, i, j, k] =
                    zc[i, j, k] -
                    0.5 * jac[i, j, k] * dz +
                    (kz - 0.5) * jac[i, j, k] * dz / nrz

                zr = rays.z[r, i, j, k]

                # ------------------------------------------------------
                # Physical ray-volume extents.
                #
                # No additional 1/nrz factor enters rays.dens. The smaller
                # dzray automatically reduces the action carried by each
                # individual vertical subvolume.
                # ------------------------------------------------------

                rays.dxray[r, i, j, k] = dx / nrx
                rays.dyray[r, i, j, k] = dy / nry
                rays.dzray[r, i, j, k] =
                    jac[i, j, k] * dz / nrz

                # ------------------------------------------------------
                # Carrier wavenumbers.
                # ------------------------------------------------------

                wnk0 = wnk_src[alpha, i, j]
                wnl0 = wnl_src[alpha, i, j]
                wnm0 = wnm_src[alpha, i, j]

                wnh0 = sqrt(wnk0^2 + wnl0^2)

                if wnh0 <= 0.0
                    error(
                        "Error in initialize_rays!: Horizontal source ",
                        "wavenumber must be nonzero for mode ",
                        alpha,
                        ".",
                    )
                end

                if wnm0 == 0.0
                    error(
                        "Error in initialize_rays!: Source vertical ",
                        "wavenumber is zero for mode ",
                        alpha,
                        ".",
                    )
                end

                # ------------------------------------------------------
                # Spectral packet widths.
                # ------------------------------------------------------

                if x_size == 1
                    dk_ini_nd = 0.0
                else
                    dk_ini_nd =
                        dkr_factor[alpha] * wnh0

                    if dk_ini_nd <= 0.0
                        error(
                            "Error in initialize_rays!: dk_ini_nd <= 0 ",
                            "for source mode ",
                            alpha,
                            ".",
                        )
                    end
                end

                if y_size == 1
                    dl_ini_nd = 0.0
                else
                    dl_ini_nd =
                        dlr_factor[alpha] * wnh0

                    if dl_ini_nd <= 0.0
                        error(
                            "Error in initialize_rays!: dl_ini_nd <= 0 ",
                            "for source mode ",
                            alpha,
                            ".",
                        )
                    end
                end

                dm_ini_nd =
                    dmr_factor[alpha] * abs(wnm0)

                if dm_ini_nd <= 0.0
                    error(
                        "Error in initialize_rays!: dm_ini_nd <= 0 ",
                        "for source mode ",
                        alpha,
                        ".",
                    )
                end

                # ------------------------------------------------------
                # Require the complete externally forced m interval to
                # remain on one side of m = 0.
                # ------------------------------------------------------

                m_lower = wnm0 - 0.5 * dm_ini_nd
                m_upper = wnm0 + 0.5 * dm_ini_nd

                if m_lower <= 0.0 <= m_upper
                    error(
                        "Error in initialize_rays!: ",
                        "ContinuousSpectralSource crosses m = 0 for mode ",
                        alpha,
                        ". Reduce dmr_factor.",
                    )
                end

                # ------------------------------------------------------
                # Spectral ray-volume centres.
                # ------------------------------------------------------

                rays.k[r, i, j, k] =
                    wnk0 -
                    0.5 * dk_ini_nd +
                    (ik - 0.5) * dk_ini_nd / nrk

                rays.l[r, i, j, k] =
                    wnl0 -
                    0.5 * dl_ini_nd +
                    (jl - 0.5) * dl_ini_nd / nrl

                rays.m[r, i, j, k] =
                    wnm0 -
                    0.5 * dm_ini_nd +
                    (km - 0.5) * dm_ini_nd / nrm

                # ------------------------------------------------------
                # Spectral ray-volume extents.
                # ------------------------------------------------------

                rays.dkray[r, i, j, k] =
                    dk_ini_nd / nrk

                rays.dlray[r, i, j, k] =
                    dl_ini_nd / nrl

                rays.dmray[r, i, j, k] =
                    dm_ini_nd / nrm

                # ------------------------------------------------------
                # Spectral phase-space volume represented by the packet.
                #
                # In the present 2-D single-column case this reduces to
                #
                #     pspvol = dm_ini_nd.
                # ------------------------------------------------------

                pspvol = dm_ini_nd

                if x_size > 1
                    pspvol *= dk_ini_nd
                end

                if y_size > 1
                    pspvol *= dl_ini_nd
                end

                # ------------------------------------------------------
                # Gaussian phase-space wave-action density.
                #
                # Identical normalization to GaussianDist.
                # ------------------------------------------------------

                rays.dens[r, i, j, k] =
                    wad0 / pspvol * m_weights[km]

                # ------------------------------------------------------
                # Check propagation direction and initialize WKB CFL
                # diagnostics using the same dispersion relation as
                # propagate_rays!.
                # ------------------------------------------------------

                n2r =
                    interpolate_stratification(zr, state, N2())

                if n2r < 0.0
                    error(
                        "Error in initialize_rays!: Negative ",
                        "stratification at source ray position.",
                    )
                end

                wnrk = rays.k[r, i, j, k]
                wnrl = rays.l[r, i, j, k]
                wnrm = rays.m[r, i, j, k]

                wnrh = sqrt(wnrk^2 + wnrl^2)

                omir =
                    branch *
                    sqrt(
                        n2r * wnrh^2 +
                        fc^2 * wnrm^2,
                    ) /
                    sqrt(wnrh^2 + wnrm^2)

                cgirz =
                    -wnrm *
                    (omir^2 - fc^2) /
                    (omir * (wnrh^2 + wnrm^2))

                # Every externally forced spectral component must enter
                # the physical domain through the lower boundary.
                if cgirz <= 0.0
                    error(
                        "Error in initialize_rays!: ",
                        "ContinuousSpectralSource contains a ",
                        "non-upward-propagating ray. Mode = ",
                        alpha,
                        ", m = ",
                        wnrm / lref,
                        " m^-1, cg_z = ",
                        cgirz * lref / tref,
                        " m s^-1.",
                    )
                end

                cgz_max[i, j, k] =
                    max(
                        cgz_max[i, j, k],
                        abs(cgirz),
                    )

                # Horizontal motion is not applied while a source ray
                # remains in k0 - 1. These group velocities are nevertheless
                # retained in the initial CFL diagnostics for general
                # multi-column configurations.
                if x_size > 1
                    cgirx =
                        wnrk *
                        (n2r - omir^2) /
                        (omir * (wnrh^2 + wnrm^2))

                    cgx_max[] =
                        max(cgx_max[], abs(cgirx))
                end

                if y_size > 1
                    cgiry =
                        wnrl *
                        (n2r - omir^2) /
                        (omir * (wnrh^2 + wnrm^2))

                    cgy_max[] =
                        max(cgy_max[], abs(cgiry))
                end
            end

            # Number of currently active source ray volumes in k0 - 1.
            nray[i, j, k] = r

            # Every possible source slot must have been visited, including
            # slots whose prescribed action is zero.
            if s != n_sfc
                error(
                    "Error in initialize_rays!: Number of source slots ",
                    s,
                    " != n_sfc = ",
                    n_sfc,
                    " at ",
                    (i, j, k),
                )
            end
        end
    end

    # ------------------------------------------------------------------
    # Global initial ray-volume count.
    #
    # On the lowest MPI subdomain this includes k0 - 1. All physical cells
    # are still empty. On other vertical subdomains it is simply zero.
    # ------------------------------------------------------------------

    kmin = ko == 0 ? k0 - 1 : k0

    @ivy local_sum =
        sum(nray[i0:i1, j0:j1, kmin:k1])

    global_sum =
        MPI.Allreduce(local_sum, +, comm)

    if master
        println("MS-GWaM:")
        println("Continuous spectral source initialized.")
        println("Global initial ray-volume count: ", global_sum)
        println("")
    end

    return
end

#=

function initialize_rays!(
    state::State,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)
    (; x_size, y_size) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (;
        nrx,
        nry,
        nrz,
        nrk,
        nrl,
        nrm,
        wave_modes,
        dkr_factor,
        dlr_factor,
        dmr_factor,
        wkb_mode,
        wave_modes,
        initial_wave_field,
    ) = state.namelists.wkb
    (; lref, tref, rhoref, uref) = state.constants
    (; comm, master, nxx, nyy, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, x, y, zc, jac) = state.grid
    (;
        nray_max,
        nray_wrk,
        n_sfc,
        nray,
        rays,
        surface_indices,
        cgx_max,
        cgy_max,
        cgz_max,
    ) = state.wkb

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref
    # Initialize local arrays.
    omi_ini = zeros(wave_modes, nxx, nyy, nzz)
    wnk_ini = zeros(wave_modes, nxx, nyy, nzz)
    wnl_ini = zeros(wave_modes, nxx, nyy, nzz)
    wnm_ini = zeros(wave_modes, nxx, nyy, nzz)
    wad_ini = zeros(wave_modes, nxx, nyy, nzz)

    # Compute initial wavenumbers, intrinsic frequencies and wave-action
    # densities with initial_wave_field.
    if wkb_mode != SteadyState()
        for k in k0:k1, j in j0:j1, i in i0:i1, alpha in 1:wave_modes
            (kdim, ldim, mdim, omegadim, adim) = initial_wave_field(
                alpha,
                x[i] * lref,
                y[j] * lref,
                zc[i, j, k] * lref,
            )
            wnk_ini[alpha, i, j, k] = kdim * lref
            wnl_ini[alpha, i, j, k] = ldim * lref
            wnm_ini[alpha, i, j, k] = mdim * lref
            omi_ini[alpha, i, j, k] = omegadim * tref
            wad_ini[alpha, i, j, k] = adim / rhoref / uref^2 / tref
        end
    else
        println(
            "Warning: MS-GWaM's steady-state mode currently ignores non-orographic initializations!",
        )
        println("")
    end

    # Add orographic wave modes.
    activate_orographic_source!(
        state,
        omi_ini,
        wnk_ini,
        wnl_ini,
        wnm_ini,
        wad_ini,
    )

    # Set initial spectral extents (these will be overwritten in the loop).
    dk_ini_nd = 0.0
    dl_ini_nd = 0.0
    dm_ini_nd = 0.0

    # Set vertical index bounds.
    kmin = ko == 0 ? k0 - 1 : k0
    kmax = k1

    # Loop over all grid cells with ray volumes.
    @ivy for k in kmin:kmax, j in j0:j1, i in i0:i1
        r = 0
        s = 0

        # Loop over all ray volumes within a spatial cell.
        for ix in 1:nrx,
            ik in 1:nrk,
            jy in 1:nry,
            jl in 1:nrl,
            kz in 1:nrz,
            km in 1:nrm,
            alpha in 1:wave_modes

            # Set ray-volume indices.
            if ko == 0 && k == k0 - 1
                s += 1

                # Set surface indices.
                surface_indices.ixs[s] = ix
                surface_indices.jys[s] = jy
                surface_indices.kzs[s] = kz
                surface_indices.iks[s] = ik
                surface_indices.jls[s] = jl
                surface_indices.kms[s] = km
                surface_indices.alphas[s] = alpha

                # Set surface ray-volume index.
                if wad_ini[alpha, i, j, k] == 0.0
                    surface_indices.rs[s, i, j] = -1
                    continue
                else
                    r += 1
                    surface_indices.rs[s, i, j] = r
                end
            else
                if wad_ini[alpha, i, j, k] == 0.0
                    continue
                end
                r += 1
            end

            # Set ray-volume positions.
            rays.x[r, i, j, k] = (x[i] - 0.5 * dx + (ix - 0.5) * dx / nrx)
            rays.y[r, i, j, k] = (y[j] - 0.5 * dy + (jy - 0.5) * dy / nry)
            rays.z[r, i, j, k] = (
                zc[i, j, k] - 0.5 * jac[i, j, k] * dz +
                (kz - 0.5) * jac[i, j, k] * dz / nrz
            )

            xr = rays.x[r, i, j, k]
            yr = rays.y[r, i, j, k]
            zr = rays.z[r, i, j, k]

            # Check if ray volume is too low.
            if zr < -dz
                error(
                    "Error in initialize_rays!: Ray volume",
                    r,
                    "at",
                    i,
                    j,
                    k,
                    "is too low!",
                )
            end

            # Compute local stratification.
            n2r = interpolate_stratification(zr, state, N2())

            # Set spatial extents.
            rays.dxray[r, i, j, k] = dx / nrx
            rays.dyray[r, i, j, k] = dy / nry
            rays.dzray[r, i, j, k] = jac[i, j, k] * dz / nrz

            wnk0 = wnk_ini[alpha, i, j, k]
            wnl0 = wnl_ini[alpha, i, j, k]
            wnm0 = wnm_ini[alpha, i, j, k]

            # Ensure correct wavenumber extents.
            if x_size > 1
                dk_ini_nd = dkr_factor[alpha] * sqrt(wnk0^2 + wnl0^2)
            end
            if y_size > 1
                dl_ini_nd = dlr_factor[alpha] * sqrt(wnk0^2 + wnl0^2)
            end
            if wnm0 == 0.0
                error("Error in WKB: wnm0 = 0!")
            else
                dm_ini_nd = dmr_factor[alpha] * abs(wnm0)
            end

            # Set ray-volume wavenumbers.
            rays.k[r, i, j, k] =
                (wnk0 - 0.5 * dk_ini_nd + (ik - 0.5) * dk_ini_nd / nrk)
            rays.l[r, i, j, k] =
                (wnl0 - 0.5 * dl_ini_nd + (jl - 0.5) * dl_ini_nd / nrl)
            rays.m[r, i, j, k] =
                (wnm0 - 0.5 * dm_ini_nd + (km - 0.5) * dm_ini_nd / nrm)

            # Set spectral extents.
            rays.dkray[r, i, j, k] = dk_ini_nd / nrk
            rays.dlray[r, i, j, k] = dl_ini_nd / nrl
            rays.dmray[r, i, j, k] = dm_ini_nd / nrm

            # Set spectral volume.
            pspvol = dm_ini_nd
            if x_size > 1
                pspvol = pspvol * dk_ini_nd
            end
            if y_size > 1
                pspvol = pspvol * dl_ini_nd
            end

            # Set phase-space wave-action density.
            rays.dens[r, i, j, k] = wad_ini[alpha, i, j, k] / pspvol

            # Interpolate winds to ray-volume position.
            uxr = interpolate_mean_flow(xr, yr, zr, state, U())
            vyr = interpolate_mean_flow(xr, yr, zr, state, V())
            wzr = interpolate_mean_flow(xr, yr, zr, state, W())

            wnrk = rays.k[r, i, j, k]
            wnrl = rays.l[r, i, j, k]
            wnrm = rays.m[r, i, j, k]
            wnrh = sqrt(wnrk^2 + wnrl^2)
            omir = omi_ini[alpha, i, j, k]

            # Compute maximum group velocities.
            cgirx = wnrk * (n2r - omir^2) / (omir * (wnrh^2 + wnrm^2))
            if abs(uxr + cgirx) > abs(cgx_max[])
                cgx_max[] = abs(uxr + cgirx)
            end
            cgiry = wnrl * (n2r - omir^2) / (omir * (wnrh^2 + wnrm^2))
            if abs(vyr + cgiry) > abs(cgy_max[])
                cgy_max[] = abs(vyr + cgiry)
            end
            cgirz = -wnrm * (omir^2 - fc^2) / (omir * (wnrh^2 + wnrm^2))
            if abs(wzr + cgirz) > abs(cgz_max[i, j, k])
                cgz_max[i, j, k] = max(cgz_max[i, j, k], abs(wzr + cgirz))
            end
        end

        # Set ray-volume count.
        nray[i, j, k] = r
        if r > nray_wrk
            error(
                "Error in initialize_rays!: nray",
                [i, j, k],
                " > nray_wrk =",
                nray_wrk,
            )
        end

        # Check if surface ray-volume count is correct.
        if ko == 0 && k == k0 - 1
            if s != n_sfc
                error(
                    "Error in initialize_rays!: s =",
                    s,
                    "/= n_sfc =",
                    n_sfc,
                    "at (i, j, k) = ",
                    (i, j, k),
                )
            end
        end
    end

    # Compute global ray-volume count.
    @ivy local_sum = sum(nray[i0:i1, j0:j1, kmin:kmax])
    global_sum = MPI.Allreduce(local_sum, +, comm)

    # Print information.
    if master
        println("MS-GWaM:")
        println("Global ray-volume count: ", global_sum)
        println("Maximum number of ray volumes per cell: ", nray_max)
        println("")
    end

    return
end
=#