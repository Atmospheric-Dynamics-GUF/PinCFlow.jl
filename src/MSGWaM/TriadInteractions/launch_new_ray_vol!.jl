function launch_new_ray_vol! end
#=
function launch_new_ray_vol!(state::State, i::Integer, j::Integer, k::Integer, kpr::AbstractFloat, mr::AbstractFloat,
     dkpr::AbstractFloat, dmr::AbstractFloat, was::AbstractFloat, triad_mode::Triad2D)

     (; nray, rays, nray_wrk) = state.wkb
     (; dx, dy, dz, x, y, zc ,zctilde, jac) = state.grid
     (;lref ) = state.constants

     
     ray_index = nray[i, j, k] + 1
    
     @ivy if ray_index > nray_wrk
            error(
                "Error in lanching new ray volumes!: nray",
                [i, j, k],
                " > nray_wrk = ",
                nray_wrk,
            )
    end


     rays.x[ray_index, i, j, k] = x[i]
     rays.y[ray_index, i, j, k] = y[j]
     rays.z[ray_index, i, j, k] = zc[i, j, k]
     rays.k[ray_index, i, j, k] =  kpr ##right now it is launching the rays with positive wave numbers only.
     rays.l[ray_index, i, j, k] = 0.0 ##need to include the negative wave numbers also
     rays.m[ray_index, i, j, k] = mr
     rays.dxray[ray_index, i, j, k] = dx
     rays.dyray[ray_index, i, j, k] = dy
     rays.dzray[ray_index, i, j, k] = dz
     rays.dkray[ray_index, i, j, k] = dkpr
     rays.dlray[ray_index, i, j, k] = 0.0
     rays.dmray[ray_index, i, j, k] = dmr
     rays.dens[ray_index, i, j, k] = was
     nray[i, j, k] += 1

    #println("\n with the non dimensional ray volume properties ", 
    #            (kpr, mr, dkpr, dmr, was))
end
=#
function launch_new_ray_vol!(state::State, i::Integer, j::Integer, k::Integer, kpr::AbstractFloat, mr::AbstractFloat,
     dkpr::AbstractFloat, dmr::AbstractFloat, was::AbstractFloat, triad_mode::Triad2D)

     (; nray, rays, nray_wrk, cgx_max, cgz_max) = state.wkb
     (; dx, dy, dz, x, y, zc ,zctilde, jac) = state.grid
     (;lref, tref) = state.constants
     (; coriolis_frequency) = state.namelists.atmosphere
     (; nrx, nrz) = state.namelists.wkb

     fc = coriolis_frequency * tref
     
     ray_index = nray[i, j, k] + 1
     
     @ivy if ray_index > nray_wrk
            error(
                "Error in lanching new ray volumes!: nray",
                [i, j, k],
                " > nray_wrk = ",
                nray_wrk,
            )
    end
    for ix in 1:nrx,
        kz in 1:nrz,
        # Set ray-volume positions.
        rays.x[ray_index, i, j, k] = (x[i] - 0.5 * dx + (ix - 0.5) * dx / nrx)
        rays.y[ray_index, i, j, k] = y[j]
        rays.z[ray_index, i, j, k] = (
            zc[i, j, k] - 0.5 * jac[i, j, k] * dz +
            (kz - 0.5) * jac[i, j, k] * dz / nrz
            )

        rays.k[ray_index, i, j, k] =  kpr ##right now it is launching the rays with positive wave numbers only.
        rays.l[ray_index, i, j, k] = 0.0 ##need to include the negative wave numbers also
        rays.m[ray_index, i, j, k] = mr
        rays.dxray[ray_index, i, j, k] = dx / nrx
        rays.dyray[ray_index, i, j, k] = dy
        rays.dzray[ray_index, i, j, k] = jac[i, j, k] * dz / nrz
        rays.dkray[ray_index, i, j, k] = dkpr
        rays.dlray[ray_index, i, j, k] = 0.0
        rays.dmray[ray_index, i, j, k] = dmr
        rays.dens[ray_index, i, j, k] = was
        xr = rays.x[ray_index, i, j, k]
        yr = rays.y[ray_index, i, j, k]
        zr = rays.z[ray_index, i, j, k]

        #to compute the max group velocity in x, z direction
        n2r = interpolate_stratification(zr, state, N2())
        uxr = interpolate_mean_flow(xr, yr, zr, state, U())
        wzr = interpolate_mean_flow(xr, yr, zr, state, W())

        wnrk = rays.k[ray_index, i, j, k]
        wnrl = rays.l[ray_index, i, j, k]
        wnrm = rays.m[ray_index, i, j, k]
        wnrh = sqrt(wnrk^2 + wnrl^2)
        omir = sqrt(n2r) * compute_omega_hat(wnrk, wnrm)

        # Compute maximum group velocities.
        cgirx = wnrk * (n2r - omir^2) / (omir * (wnrh^2 + wnrm^2))
        if abs(uxr + cgirx) > abs(cgx_max[])
            cgx_max[] = abs(uxr + cgirx)
        end
        
        cgirz = -wnrm * (omir^2 - fc^2) / (omir * (wnrh^2 + wnrm^2))
        if abs(wzr + cgirz) > abs(cgz_max[i, j, k])
            cgz_max[i, j, k] = max(cgz_max[i, j, k], abs(wzr + cgirz))
        end
        ray_index += 1
        nray[i, j, k] += 1
    end
    
end


function launch_new_ray_vol!(state::State, i::Integer, j::Integer, k::Integer, kpr::AbstractFloat, mr::AbstractFloat,
     dkpr::AbstractFloat, dmr::AbstractFloat, was::AbstractFloat, triad_mode::Triad3DIso)

     (; nray, rays, nray_wrk) = state.wkb
     (; dx, dy, dz, x, y, zc ,zctilde, jac) = state.grid

     
     ray_index = nray[i, j, k] + 1
    
     @ivy if ray_index > nray_wrk
            error(
                "Error in lanching new ray volumes!: nray",
                [i, j, k],
                " > nray_wrk = ",
                nray_wrk,
            )
    end


     rays.x[ray_index, i, j, k] = x[i]
     rays.y[ray_index, i, j, k] = y[j]
     rays.z[ray_index, i, j, k] = zc[i, j, k]
     rays.k[ray_index, i, j, k] = sqrt(2) * kpr ##right now it is launching the rays with positive wave numbers only.
     rays.l[ray_index, i, j, k] = sqrt(2) * kpr ##need to include the negative wave numbers also
     rays.m[ray_index, i, j, k] = mr
     rays.dxray[ray_index, i, j, k] = dx
     rays.dyray[ray_index, i, j, k] = dy
     rays.dzray[ray_index, i, j, k] = dz
     rays.dkray[ray_index, i, j, k] = sqrt(2) * dkpr
     rays.dlray[ray_index, i, j, k] = sqrt(2) * dkpr
     rays.dmray[ray_index, i, j, k] = dmr
     rays.dens[ray_index, i, j, k] = was
     nray[i, j, k] += 1


end

