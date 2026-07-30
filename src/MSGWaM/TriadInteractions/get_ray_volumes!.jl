function get_ray_volumes! end

function get_ray_volumes!(
    state::State,
    triad_mode::Triad2D;
    print_action_diagnostic::Bool = true,
)
    (; domain) = state
    (; master, comm, i0, i1, j0, j1, k0, k1) = domain

    (; x_size, y_size) = state.namelists.domain
    (; action_rel_tol, discarded_action_fraction_tol) = state.namelists.triad

    (; nray, rays, spec_tend) = state.wkb
    (; wavespectrum, was_pred, was_ray_signature, action_ref) = spec_tend
    (; kp, m, kpc, mc, delkp, delm) = spec_tend.spec_grid

    (; dx, dy, dz, x, y, zctilde, jac) = state.grid

    if !(0.0 < discarded_action_fraction_tol <= 1.0)
        throw(
            ArgumentError(
                "discarded_action_fraction_tol must satisfy " *
                "0 < discarded_action_fraction_tol <= 1.",
            ),
        )
    end

    if !isfinite(action_ref[]) || action_ref[] <= 0.0
        error("Invalid action_ref = ", action_ref[])
    end

    action_floor = action_rel_tol * action_ref[]

    local_total_action = 0.0
    local_discarded_action = 0.0

    for k in (k0 - 1):(k1 + 1),
        j in (j0 - 1):(j1 + 1),
        i in (i0 - 1):(i1 + 1)

        if all(iszero, @view wavespectrum[i, j, k, :, :])
            continue
        end

        #------------------------------------------------------
        # Redistribute the updated Eulerian wave spectrum among
        # the existing ray volumes.
        #------------------------------------------------------

        for r in 1:nray[i, j, k]
            xr = rays.x[r, i, j, k]
            yr = rays.y[r, i, j, k]
            zr = rays.z[r, i, j, k]

            dxr = rays.dxray[r, i, j, k]
            dyr = rays.dyray[r, i, j, k]
            dzr = rays.dzray[r, i, j, k]

            kr = rays.k[r, i, j, k]
            mr = rays.m[r, i, j, k]

            dkr = rays.dkray[r, i, j, k]
            dmr = rays.dmray[r, i, j, k]

            kpr = abs(kr)
            dkpr = dkr

            wadr = rays.dens[r, i, j, k]

            # The density is reconstructed from the updated
            # Eulerian wave spectrum.
            rays.dens[r, i, j, k] = 0.0

            imin, imax, jmin, jmax =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            kpmin, kpmax, mmin, mmax =
                compute_spectral_cell_indices(state, kpr, mr, dkpr, dmr)

            for iray in imin:imax
                if x_size > 1
                    dxi =
                        min(xr + dxr / 2, x[iray] + dx / 2) -
                        max(xr - dxr / 2, x[iray] - dx / 2)
                else
                    dxi = 1.0
                    dxr = 1.0
                end

                for jray in jmin:jmax
                    if y_size > 1
                        dyi =
                            min(yr + dyr / 2, y[jray] + dy / 2) -
                            max(yr - dyr / 2, y[jray] - dy / 2)
                    else
                        dyi = 1.0
                        dyr = 1.0
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
                            min(zr + dzr / 2, zctilde[iray, jray, kray]) -
                            max(zr - dzr / 2, zctilde[iray, jray, kray - 1])

                        for kpray in kpmin:kpmax
                            dkpi =
                                min(kpr + dkr / 2, kpc[kpray + 1]) -
                                max(kpr - dkr / 2, kpc[kpray])

                            for mray in mmin:mmax
                                if mr >= 0.0
                                    dmi =
                                        min(mr + dmr / 2, mc[mray + 2]) -
                                        max(mr - dmr / 2, mc[mray + 1])
                                else
                                    dmi =
                                        min(mr + dmr / 2, mc[mray + 1]) -
                                        max(mr - dmr / 2, mc[mray])
                                end

                                was0 = was_pred[
                                    iray,
                                    jray,
                                    kray,
                                    kpray,
                                    mray,
                                ]

                                was1 = wavespectrum[
                                    iray,
                                    jray,
                                    kray,
                                    kpray,
                                    mray,
                                ]

                                if was0 <= 0.0 || was1 == 0.0
                                    continue
                                end

                                # Total five-dimensional ray
                                # phase-space volume.
                                ray_vol = dxr * dyr * dzr * dkpr * dmr

                                # Fraction of that volume overlapping
                                # this physical-spectral grid cell.
                                occupied_vol =
                                    dxi * dyi * dzi * dkpi * dmi

                                rays.dens[r, i, j, k] +=
                                    wadr *
                                    occupied_vol *
                                    was1 /
                                    was0 /
                                    ray_vol
                            end
                        end
                    end
                end
            end
        end

        #------------------------------------------------------
        # Physical-cell volume for the global action diagnostic.
        #
        # Halo cells are needed for ray remapping and launching,
        # but they must not be included in the global total.
        #------------------------------------------------------

        owned_cell =
            i0 <= i <= i1 &&
            j0 <= j <= j1 &&
            k0 <= k <= k1

        if owned_cell
            dx_cell = x_size > 1 ? dx : 1.0
            dy_cell = y_size > 1 ? dy : 1.0

            physical_cell_volume =
                dx_cell * dy_cell * jac[i, j, k] * dz
        else
            physical_cell_volume = 0.0
        end

        #------------------------------------------------------
        # Launch new ray volumes for newly generated or
        # spectrally dispersed wave action that is not represented
        # by an existing ray volume.
        #------------------------------------------------------

        for mi in eachindex(m), kpi in eachindex(kp)
            was = wavespectrum[i, j, k, kpi, mi]
            was_sig = was_ray_signature[i, j, k, kpi, mi]

            dkps = delkp[kpi]
            dms = delm[mi]

            spectral_cell_action = was * dkps * dms

            launch_new_ray =
                !was_sig &&
                spectral_cell_action > action_floor

            discarded =
                !was_sig &&
                0.0 < spectral_cell_action <= action_floor

            if launch_new_ray
                kps = kp[kpi]
                ms = m[mi]

                launch_new_ray_vol!(
                    state,
                    i,
                    j,
                    k,
                    kps,
                    ms,
                    dkps,
                    dms,
                    was,
                    triad_mode,
                )
            end

            if owned_cell
                integrated_action =
                    spectral_cell_action * physical_cell_volume

                local_total_action += integrated_action

                if discarded
                    local_discarded_action += integrated_action
                end
            end
        end
    end

    #----------------------------------------------------------
    # Obtain global totals over all MPI ranks with one collective
    # communication operation.
    #----------------------------------------------------------

    local_action = [
        local_total_action,
        local_discarded_action,
    ]

    global_action = MPI.Allreduce(local_action, +, comm)

    total_action = global_action[1]
    discarded_action = global_action[2]

    if !isfinite(total_action) || !isfinite(discarded_action)
        error(
            "Nonfinite ray-volume action diagnostic: ",
            "total_action = ",
            total_action,
            ", discarded_action = ",
            discarded_action,
        )
    end

    if total_action < 0.0
        error("Negative total wave action detected: ", total_action)
    end

    discarded_action_fraction =
        total_action > 0.0 ?
        discarded_action / total_action :
        0.0

    if master && print_action_diagnostic
        println("")
        println("New ray-volume launch diagnostic:")
        println("  Total wave action              = ", total_action)
        println("  Discarded wave action          = ", discarded_action)
        println("  Discarded wave-action fraction = ", discarded_action_fraction)
        println("  Allowed discarded fraction     = ", discarded_action_fraction_tol)
        println("")
    end

    # Every rank throws the same error after the global reduction.
    # This is safer than terminating only the master rank.
    if discarded_action_fraction >= discarded_action_fraction_tol
        error(
            "Discarded wave-action fraction exceeded the allowed threshold: ",
            discarded_action_fraction,
            " >= ",
            discarded_action_fraction_tol,
        )
    end

    return nothing
end

function get_ray_volumes!(state::State, 
    wavespectrum_copy::Array{<: AbstractFloat, 5}, 
    triad_mode::Triad3DIso)
    (; domain, grid) = state
    (; branch) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = domain
    (; nray, rays, spec_tend) = state.wkb
    (; kp, m) = spec_tend.spec_grid
    (;lref ) = state.constants
    (; dx, dy, dz, x, y, zc ,zctilde, jac) = state.grid

    
    rays.dens .= 0
   
        

    (ukp, lkp) = half_logwidth(kp)
    (um, lm) = half_logwidth(m)

    @ivy for k in (k0 - 1):(k1 + 1),
        j in (j0 - 1):(j1 + 1),
        i in (i0 - 1):(i1 + 1)

        @ivy for mi in eachindex(m),
            kpi in eachindex(kp)
            
            was_old = wavespectrum_copy[i, j, k, kpi, mi]
            was = spec_tend.wavespectrum[i, j, k, kpi, mi]
            """
            if (i, j, k, kpi, mi) == (4, 4, 28, 10, 10)
                was = 1.0
            end
            """
            
            rv = spec_tend.ray_vol_signature[i, j, k, kpi, mi]

            if was == 0 && length(rv) == 0

                continue

            elseif was != 0 && length(rv) == 0
                #println("new ray volume loop called \n new ray volume launched at ", 
                #(x[i]*lref, y[j]*lref, zc[i, j, k]*lref, kp[kpi]/lref, m[mi]/lref))
                kpr = kp[kpi]
                mr = m[mi]
                dkpr = ukp[kpi] - lkp[kpi]
                dmr = um[mi] - lm[mi]
                launch_new_ray_vol!(state, i, j, k, kpr, mr, dkpr, dmr, was, triad_mode)

            else
                for tup in rv
                    #println(rays.dens[tup[1], tup[2], tup[3], tup[4]])
                    rays.dens[tup[1], tup[2], tup[3], tup[4]] += tup[5] * tup[6] * was_old / was
                    #println(rays.dens[tup[1], tup[2], tup[3], tup[4]])  
                end
                 
            end

        end

    end
 
end


