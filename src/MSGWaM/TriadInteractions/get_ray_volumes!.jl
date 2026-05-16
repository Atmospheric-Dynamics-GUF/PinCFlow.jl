function get_ray_volumes! end

#=
function get_ray_volumes!(state::State, 
    wavespectrum_copy::Array{<: AbstractFloat, 5}, 
    triad_mode::Triad2D)
    (; domain, grid) = state
    (; branch) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = domain
    (; nray, rays, spec_tend) = state.wkb
    (; kp, m, kpc, mc) = spec_tend.spec_grid
    (;lref ) = state.constants
    (; dx, dy, dz, x, y, zc ,zctilde, jac) = state.grid

    #setting up wave action density equal to zero for all rays, they will be re-written in loop
    rays.dens .= 0



    @ivy for k in (k0 - 1):(k1 + 1),
        j in (j0 - 1):(j1 + 1),
        i in (i0 - 1):(i1 + 1)

         for mi in eachindex(m),
            kpi in eachindex(kp)
            
            was_old = wavespectrum_copy[i, j, k, kpi, mi]
            was = spec_tend.wavespectrum[i, j, k, kpi, mi]
            
            
            rv = spec_tend.ray_vol_signature[i, j, k, kpi, mi]

            if was == 0 && length(rv) == 0

                continue

            elseif was != 0 && length(rv) == 0
                #println("new ray volume loop called \n new ray volume launched at ", 
                #(x[i]*lref, y[j]*lref, zc[i, j, k]*lref, kp[kpi]/lref, m[mi]/lref))
                kpr = kp[kpi]
                mr = m[mi]
                dkpr = kpc[kpi + 1] - kpc[kpi]
                if mr > 0 
                    dmr = mc[mi + 2] - mc[mi + 1]
                else
                    dmr = mc[mi + 1] - mc[mi]
                end
                launch_new_ray_vol!(state, i, j, k, kpr, mr, dkpr, dmr, was, triad_mode)
            else
                for tup in rv                    
                    rays.dens[tup[1], tup[2], tup[3], tup[4]] += (tup[5] * tup[6] * was / was_old / tup[7])
                end
                 
            end

        end

    end
 
end
=#

function get_ray_volumes!(state::State,  
    triad_mode::Triad2D)
    (; domain, grid) = state
    (; x_size, y_size) = state.namelists.domain
    (; branch) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = domain
    (; nray, rays, spec_tend) = state.wkb
    (; wavespectrum, was_pred, was_ray_signature) = spec_tend
    (; kp, m, kpc, mc) = spec_tend.spec_grid
    (;lref ) = state.constants
    (; dx, dy, dz, x, y, zc ,zctilde, jac) = state.grid

    #setting up wave action density equal to zero for all rays, they will be re-written in loop
    #rays.dens .= 0



   @ivy for k in (k0 - 1):(k1 + 1),
        j in (j0 - 1):(j1 + 1),
        i in (i0 - 1):(i1 + 1)

        if all(iszero, wavespectrum[i, j, k, :, :])
            continue
        end
    
        for r in 1:nray[i, j, k]   #mapping existing rays on Eulerian Grid

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

            kpr = abs(kr)
            dkpr = dkr

            wadr = rays.dens[r, i, j, k]
            rays.dens[r, i, j, k] = 0

            
            (imin, imax, jmin, jmax) =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            (kpmin, kpmax, mmin, mmax) = compute_spectral_cell_indices(state, kpr, mr, dkpr, dmr)
            
            for iray in imin:imax
                if x_size > 1
                    dxi = (
                        min(xr + dxr / 2, x[iray] + dx / 2) -
                        max(xr - dxr / 2, x[iray] - dx / 2)
                    )

                    fcpspx =  dxi / dx
                else
                    dxi = 1.0
                    fcpspx = 1.0
                    dxr = 1.0
                end

                for jray in jmin:jmax
                    if y_size > 1
                        dyi = (
                            min(yr + dyr / 2, y[jray] + dy / 2) -
                            max(yr - dyr / 2, y[jray] - dy / 2)
                        )

                        fcpspy =  dyi / dy
                    else
                        dyi = 1.0
                        fcpspy = 1.0
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
                            min((zr + dzr / 2), zctilde[iray, jray, kray]) -
                            max((zr - dzr / 2), zctilde[iray, jray, kray - 1])

                        fcpspz =  dzi / jac[iray, jray, kray] / dz
                        
                        for kpray in kpmin:kpmax
                             dkpi = 
                                min(kpr + dkr / 2, kpc[kpray + 1]) -  #kpi_i > 0, always lies between kpc_{i+1} to kpc_{i}
                                max(kpr - dkr / 2, kpc[kpray] )
                            dkp = kpc[kpray + 1] - kpc[kpray]
                                      
                            #fcpspkp =  dkpi / dkpr
                            fcpspkp =  dkpi / dkp

                             for mray in mmin:mmax
                                if mr >= 0       #becuase for mr > 0,  m_i > 0 always lies between mc_{i+2} to mc_{i+1}
                                    dmi = 
                                        min(mr + dmr / 2, mc[mray + 2]) -
                                        max(mr - dmr / 2, mc[mray + 1])
                                    dm = mc[mray + 2] - mc[mray + 1]
                                    #fcpspm = dmi / dmr
                                    fcpspm = dmi / dm
                                else
                                    dmi = 
                                        min(mr + dmr / 2, mc[mray + 1]) -
                                        max(mr - dmr / 2, mc[mray])
                                    dm = mc[mray + 1] - mc[mray]
                                    #fcpspm = dmi / dmr
                                    fcpspm = dmi / dm
                                end
                                was0 = was_pred[iray, jray, kray, kpray, mray]
                                was1 = wavespectrum[iray, jray, kray, kpray, mray]

                                if was0 <= 0 || was1 == 0
                                    continue
                                end
                                ray_vol = dxr * dyr * dzr * dkpr * dmr #total volume of the ray in 5D ray phase space
                                accupied_vol = dxi * dyi * dzi * dkpi * dmi #total volume of the ray volume (r, i, j, k) accupied by the grid cell
                                rays.dens[r, i, j, k] += wadr * accupied_vol * was1 / 
                                                        was0 / ray_vol
                                #println(spec_tend.wavespectrum[iray, jray, kray, kpray, mray])

                             end

                        end

                    end

                end

            end

        end

        #lanching new valumes for the newly generated wave modes
        for mi in eachindex(m),
            kpi in eachindex(kp)
        
            was = wavespectrum[i, j, k, kpi, mi]
            
            
            was_sig = was_ray_signature[i, j, k, kpi, mi]

            if was != 0 && was_sig == false
                #println("new ray volume loop called \n new ray volume launched at ", 
                #(x[i]*lref, y[j]*lref, zc[i, j, k]*lref, kp[kpi]/lref, m[mi]/lref))
                kps = kp[kpi]
                ms = m[mi]
                dkps = kpc[kpi + 1] - kpc[kpi]
                if ms > 0 
                    dms = mc[mi + 2] - mc[mi + 1]
                else
                    dms = mc[mi + 1] - mc[mi]
                end
                launch_new_ray_vol!(state, i, j, k, kps, ms, dkps, dms, was, triad_mode)

            end

        end

    end
 
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


