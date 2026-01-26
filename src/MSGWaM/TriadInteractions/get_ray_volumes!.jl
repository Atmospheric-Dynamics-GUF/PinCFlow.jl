function get_ray_volumes! end



function get_ray_volumes!(state::State, wavespectrum_copy::Array{<: AbstractFloat, 5}, triad_mode::Triad3DIso)
    (; domain, grid) = state
    (; branch) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = domain
    (; nray, rays, spec_tend) = state.wkb
    (; kp, m) = spec_tend.spec_grid
    (;lref ) = state.constants
    (; dx, dy, dz, x, y, zc ,zctilde, jac) = state.grid

    
    rays.dens .= 0

    println("calling the get ray volume")
    println("maximum wave density before getting it back", maximum(rays.dens) )
   
        

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


function get_ray_volumes!(state::State, wavespectrum_copy::Array{<: AbstractFloat, 5}, triad_mode::Triad2D)
    (; domain, grid) = state
    (; branch) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = domain
    (; nray, rays, spec_tend) = state.wkb
    (; kp, m, kpc, mc) = spec_tend.spec_grid
    (;lref ) = state.constants
    (; dx, dy, dz, x, y, zc ,zctilde, jac) = state.grid

    #setting up wave action density equal to zero for all rays, they will be re-written in loop
    rays.dens .= 0

   println("\n Getting the modified wave action density for the existing ray volumes, also lanching new ray volumes if required")
    #println("maximum wave density before getting it back", maximum(rays.dens) )
   
        

    (ukp, lkp) = half_logwidth(kp)
    (um, lm) = half_logwidth(m)

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
                println("new ray volume loop called \n new ray volume launched at ", 
                (x[i]*lref, y[j]*lref, zc[i, j, k]*lref, kp[kpi]/lref, m[mi]/lref))
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

                    rays.dens[tup[1], tup[2], tup[3], tup[4]] += (tup[5] * tup[6] * was / was_old / tup[7] )

                end
                 
            end

        end

    end
 
end